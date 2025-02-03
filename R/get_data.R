set_root_for_phyloseq <- function(physeq) {
  tree.unrooted <- phy_tree(physeq)
  treeDT <-
    cbind(data.table(tree.unrooted$edge),
          data.table(length = tree.unrooted$edge.length))[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  phy_tree(physeq) <-
    ape::root(tree.unrooted,
              outgroup = new.outgroup,
              resolve.root = TRUE)
  return(physeq)
}


get_possible_mitochondrial_otu <- function(physeq, seq_df){
  tax_df <- tax_table(physeq) %>% as.data.frame()
  tax_mito_df <- tax_df %>% filter(Domain %in% c("Archaea", "Bacteria")) %>% filter(Phylum == "?") %>%
    rownames_to_column(var="sequence_name") %>% left_join(seq_df, by="sequence_name")
  return(tax_mito_df)
}


remove_mitochondrial_otu_physeq <- function(physeq, otu_ids){
  physeq_filtered <- prune_taxa(!taxa_names(physeq) %in% otu_ids, physeq)
  return(physeq_filtered)
}


read_sig_metabolite <- function(file_path){
  df <- read_csv(file_path)
  df$comparison <- df$comparison %>%
    str_replace_all("A. fumigatus delta_pptA \\+ Corisone acetate", "PptA") %>%
    str_replace_all("A. fumigatus wt \\+ Cortisone acetate \\+ Voriconazol", "Tx") %>%
    str_replace_all("PBS \\+ Cortisone acetate \\+ Voriconazol", "VOR") %>%
    str_replace_all("A. fumigatus wt \\+ Cortisone acetate", "WT") %>%
    str_replace_all("PBS \\+ Cortisone acetate", "CA") %>%
    str_replace_all("PBS", "CTR")
  return(df)
}


# load lung phyloseq object and change OTU ID
get_lung_physeq <- function(phyloseq_path, meta_df) {
  library(microViz)
  load(phyloseq_path)
  sample_names(physeq) <- str_replace_all(sample_names(physeq),"@","-")
  physeq <- physeq %>%
    subset_taxa(Domain=="Bacteria") %>%
    subset_taxa(Genus != "Photobacterium") %>%
    subset_taxa(Genus != "Methanothermobacter") %>%
    subset_taxa(Genus != "Bradyrhizobium") %>%
    #subset_taxa(Genus != "Paracoccus") %>%
    subset_taxa(Genus != "Vulcaniibacterium") %>%
    subset_taxa(Genus != "Kineococcus") %>%
    #subset_taxa(Genus != "Klenkia") %>%
    subset_taxa(Genus != "Buchnera") %>%
    subset_taxa(Genus != "Streptophyta") %>%
    subset_taxa(Genus != "Sporolactobacillaceae_incertae_sedis") %>%
    subset_taxa(Genus != "Spirosoma") %>%
    subset_taxa(Genus != "Dyadobacter") %>%
    subset_taxa(Genus != "Candidatus Cloacamonas") %>%
    subset_taxa(Genus != "Candidatus Xiphinematobacter")
    
  new_sample_data_df <- sample_data(physeq) %>%
    data.frame() %>%
    rownames_to_column() %>%
    left_join(meta_df, by = c("rowname" = "sample")) %>%
    column_to_rownames() %>%
    filter(!sampleGroup %in% c(7,8)) %>%
    filter(!sampleType %in% c("PCR-POS"))
  sample_data(physeq) <- new_sample_data_df
  physeq <- physeq %>% tax_fix()
  return(physeq)
}


load_gut_baseline_physeq <- get_gut_physeq <- function(phyloseq_path, meta_df) {
  library(microViz)
  load(phyloseq_path)
  new_sample_data_df <- sample_data(physeq) %>%
    data.frame() %>%
    rownames_to_column() %>%
    filter(!str_starts(rowname,"F2")) %>%
    filter(!str_starts(rowname,"F3")) %>%
    filter(!str_starts(rowname,"Kot")) %>%
    filter(!str_starts(rowname,"pos")) %>%
    mutate(rowname_temp = ifelse(str_starts(rowname, "gm"), str_replace(rowname,"-","-gm"), rowname)) %>%
    mutate(id = ifelse(str_starts(rowname_temp, "F1"), str_remove(rowname_temp,"F1-"), rowname_temp)) %>%
    mutate(id = str_split(id,"-", simplify = T)[,2]) %>%
    left_join(meta_df, by = c("id" = "sampleName")) %>%
    mutate(sampleGroup = ifelse(str_starts(rowname_temp, "gm"), 7, sampleGroup)) %>%
    mutate(sampleGroup = ifelse(str_starts(rowname_temp, "SPF"), 8, sampleGroup)) %>%
    mutate(sampleGroup = ifelse(str_starts(rowname_temp, "neg"), "EXT-NEG", sampleGroup)) %>%
    mutate(sampleGroup2 = sampleGroup) %>%
    dplyr::rename(sampleName=id) %>%
    dplyr::select(-rowname_temp) %>%
    filter(!is.na(sampleGroup)) %>%
    column_to_rownames()%>%
    filter(!sampleGroup %in% c(7,8))
  sample_data(physeq) <- new_sample_data_df
  physeq <- physeq %>% tax_fix()
  return(physeq)
}
# load gut phyloseq object and change OTU ID
get_gut_physeq <- function(phyloseq_path, meta_df) {
  library(microViz)
  load(phyloseq_path)
  new_sample_data_df <- sample_data(physeq) %>%
    data.frame() %>%
    rownames_to_column() %>%
    filter(!str_starts(rowname,"F1")) %>%
    filter(!str_starts(rowname,"F2")) %>%
    filter(!str_starts(rowname,"Kot")) %>%
    filter(!str_starts(rowname,"pos")) %>%
    mutate(rowname_temp = ifelse(str_starts(rowname, "gm"), str_replace(rowname,"-","-gm"), rowname)) %>%
    mutate(id = ifelse(str_starts(rowname_temp, "F3"), str_remove(rowname_temp,"F3-"), rowname_temp)) %>%
    mutate(id = str_split(id,"-", simplify = T)[,2]) %>%
    left_join(meta_df, by = c("id" = "sampleName")) %>%
    mutate(sampleGroup = ifelse(str_starts(rowname_temp, "gm"), 7, sampleGroup)) %>%
    mutate(sampleGroup = ifelse(str_starts(rowname_temp, "SPF"), 8, sampleGroup)) %>%
    mutate(sampleGroup = ifelse(str_starts(rowname_temp, "neg"), "EXT-NEG", sampleGroup)) %>%
    mutate(sampleGroup2 = sampleGroup) %>%
    dplyr::rename(sampleName=id) %>%
    dplyr::select(-rowname_temp) %>%
    filter(!is.na(sampleGroup)) %>%
    column_to_rownames()%>%
    filter(!sampleGroup %in% c(7,8))
  sample_data(physeq) <- new_sample_data_df 
  physeq <- physeq %>% tax_fix()
  return(physeq)
}
# create metabolome phyloseq object and change OTU ID
get_metabolome_physeq <- function(otu_path, meta_df) {
  otu_df <- read_csv(otu_path)
  meta_df<- otu_df %>%
    dplyr::select(sample) %>%
    left_join(meta_df, by=c("sample"="sampleName")) %>%
    dplyr::select(-sample.y) %>%
    mutate(sampleName=sample) %>% 
    data.frame()
  rownames(meta_df) <- paste0("s",meta_df$sample)
  otu_df <- otu_df %>%
    mutate(sample=paste0("s",sample)) %>% 
    column_to_rownames(var="sample") %>%
    rotate_df()
  taxa_df <- data.frame("metabolite" = rownames(otu_df))
  rownames(taxa_df) <- taxa_df$metabolite
  META <- sample_data(meta_df)
  TAX <- tax_table(as.matrix(taxa_df))
  OTU <- otu_table(otu_df,taxa_are_rows=T)
  physeq <- phyloseq(META, TAX, OTU) %>% 
    transform_sample_counts(function(OTU)log2(OTU))
  return(physeq)
}
remove_spurious_taxa <- function(physeq, abundance_cutoff=0.0025,freq_cutoff=0){
  TMM_physeq <-
    transform_sample_counts(physeq, function(OTU)
      OTU / sum(OTU))
  kept_taxa <- filter_taxa(TMM_physeq, function(x) sum(x > abundance_cutoff) > freq_cutoff,F)
  filtered_physeq <- filtered_physeq <- prune_taxa(kept_taxa, physeq)
  return(filtered_physeq)
}
# remove contamination
get_decontaminated_lung_physeq <-
  function(physeq, output_dir_path = "") {
    bac_physeq <- physeq %>%
      subset_taxa(Domain %in% c("Archaea", "Bacteria"))
    if(!is.null(phy_tree(physeq, F))){
      bac_physeq <- bac_physeq %>%
        set_root_for_phyloseq()
    } 

    ## remove samples of pcr-pos,and negative control for SPF mice(sampleName is C1)
    bac_physeq %<>% subset_samples(!sampleType %in% c("PCR-POS")) %>%
      subset_samples(!sampleName %in% c("C1"))
    ## EXT
    bac_ext_physeq <-
      bac_physeq %>% subset_samples(!sampleType %in% c("PCR-NEG","ENV-Ctrl","7","8"))
    bac_ext_decom_otu <- otu_table(bac_ext_physeq) %>%
      as.matrix() %>%
      t()
    bac_ext_decom_taxa <-
      tax_table(bac_ext_physeq) %>% data.frame()
    bac_ext_neg_flag <- sample_data(bac_ext_physeq) %>%
      data.frame() %>%
      mutate(neg_flag = sampleType %in% c("ENV-Ctrl", "EXT-NEG", "PCR-NEG")) %>%
      pull(neg_flag)
    bac_ext_conc <-  sample_data(bac_ext_physeq) %>%
      data.frame() %>%
      pull(concLib)
    stopifnot(
      sample_data(bac_ext_physeq) %>% data.frame() %>% rownames() == rownames(bac_ext_decom_otu)
    )
    bac_ext_otu <- bac_ext_decom_otu[bac_ext_neg_flag,] %>% t
    bac_ext_contam_by_ext <- rowSums(bac_ext_otu>10) > 0
    bac_ext_contam.freq <-
      isContaminant(
        bac_ext_decom_otu,
        conc = bac_ext_conc,
        method = "either",
        threshold = 0.1,
        neg = bac_ext_neg_flag
      ) %>% merge(bac_ext_decom_taxa, by = 0)
    bac_ext_contam_by_ext <- bac_ext_contam_by_ext[bac_ext_contam.freq$Row.names]
    stopifnot(bac_ext_contam.freq$Row.names == names(bac_ext_contam_by_ext))
    bac_ext_contam.freq$contaminant <-   bac_ext_contam.freq$contaminant | bac_ext_contam_by_ext
    ## PCR
    bac_pcr_physeq <-
      bac_physeq %>% subset_samples(!sampleType %in% c("EXT-NEG","ENV-Ctrl","7","8"))
    bac_pcr_decom_otu <- otu_table(bac_pcr_physeq) %>%
      as.matrix() %>%
      t()
    bac_pcr_decom_taxa <-
      tax_table(bac_pcr_physeq) %>% data.frame()
    bac_pcr_neg_flag <- sample_data(bac_pcr_physeq) %>%
      data.frame() %>%
      mutate(neg_flag = sampleType %in% c("ENV-Ctrl", "EXT-NEG", "PCR-NEG")) %>%
      pull(neg_flag)
    bac_pcr_conc <-  sample_data(bac_pcr_physeq) %>%
      data.frame() %>%
      pull(concLib)
    stopifnot(
      sample_data(bac_pcr_physeq) %>% data.frame() %>% rownames() == rownames(bac_pcr_decom_otu)
    )
    bac_pcr_otu <- bac_pcr_decom_otu[bac_pcr_neg_flag,] %>% t
    bac_pcr_contam_by_pcr <- rowSums(bac_pcr_otu>10) > 0
    bac_pcr_contam.freq <-
      isContaminant(
        bac_pcr_decom_otu,
        conc = bac_pcr_conc,
        method = "either",
        threshold = 0.1,
        neg = bac_pcr_neg_flag
      ) %>% merge(bac_pcr_decom_taxa, by = 0)
    bac_pcr_contam_by_pcr <- bac_pcr_contam_by_pcr[bac_pcr_contam.freq$Row.names]
    stopifnot(bac_pcr_contam.freq$Row.names == names(bac_pcr_contam_by_pcr))
    bac_pcr_contam.freq$contaminant <-   bac_pcr_contam.freq$contaminant | bac_pcr_contam_by_pcr
    
    bac_selected_reg_contam <-
      bind_rows(bac_ext_contam.freq %>% filter(contaminant), bac_pcr_contam.freq %>% filter(contaminant))
    bac_env_physeq <-
      bac_physeq %>% subset_samples(!sampleType %in% c("EXT-NEG", "PCR-NEG", "7","8"))
    bac_env_decom_otu <- otu_table(bac_env_physeq) %>%
      as.matrix() %>%
      t()
    bac_env_decom_taxa <-
      tax_table(bac_env_physeq) %>% data.frame()
    bac_env_neg_flag <- sample_data(bac_env_physeq) %>%
      data.frame() %>%
      mutate(neg_flag = sampleType %in% c("ENV-Ctrl", "EXT-NEG", "PCR-NEG")) %>%
      pull(neg_flag)
    bac_env_conc <- sample_data(bac_env_physeq) %>%
      data.frame() %>%
      pull(concLib)
    stopifnot(
      sample_data(bac_env_physeq) %>% data.frame() %>% rownames() == rownames(bac_env_decom_otu)
    )
    bac_env_otu <- bac_env_decom_otu[bac_env_neg_flag,] %>% t
    bac_env_contam_by_env <- rowSums(bac_env_otu>10) > 0
    bac_env_contam.freq <-
      isContaminant(
        bac_env_decom_otu,
        conc=bac_env_conc,
        method = "either",
        threshold = 0.1,
        neg = bac_env_neg_flag
      ) %>% merge(bac_env_decom_taxa, by = 0)
    bac_env_contam_by_env <- bac_env_contam_by_env[bac_env_contam.freq$Row.names]
    stopifnot(bac_env_contam.freq$Row.names == names(bac_env_contam_by_env))
    bac_env_contam.freq$contaminant <-   bac_env_contam.freq$contaminant | bac_env_contam_by_env
    bac_selected_env_contam <-
      bac_env_contam.freq %>% filter(contaminant)
    
    bac_decomed_physeq <-
      prune_taxa(
        !taxa_names(bac_physeq) %in% c(
          bac_selected_env_contam$Row.names,
          bac_selected_reg_contam$Row.names
        ),
        bac_physeq
      ) %>%
      filter_taxa(function(x)
        sum(x > 0) > 0, TRUE)
    bac_decomed_physeq <- tax_glom(bac_decomed_physeq, taxrank = "Genus")
    ## remove unneeded taxa
    taxa_df <- tax_table(bac_decomed_physeq)
    tax_table(bac_decomed_physeq) <-
      taxa_df[, colSums(is.na(taxa_df)) < nrow(taxa_df)]
    ## change OTU ID
    taxa_names(bac_decomed_physeq) <-
      tax_table(bac_decomed_physeq) %>% as.data.frame %>% unite("OTU", rank_names(bac_decomed_physeq)) %>% pull(OTU)
    if (output_dir_path != "") {
      write_csv(bac_ext_contam.freq,
                file.path(output_dir_path, "lung_ext_contam.csv"))
      write_csv(bac_pcr_contam.freq,
                file.path(output_dir_path, "lung_pcr_contam.csv"))
      write_csv(bac_env_contam.freq,
                file.path(output_dir_path, "lung_env_contam.csv"))
      write.csv(
        as.data.frame(otu_table(bac_decomed_physeq)),
        file.path(output_dir_path, "lung_decontamed_genus.csv")
      )
    }
    
    return(bac_decomed_physeq)
  }


get_decontaminated_gut_physeq <-
  function(physeq, output_dir_path = "") {
    bac_physeq <- physeq %>%
      subset_taxa(Domain %in% c("Archaea", "Bacteria")) %>%
      set_root_for_phyloseq()
    
    bac_reg_physeq <-bac_physeq 
    bac_reg_decom_otu <- otu_table(bac_reg_physeq) %>%
      as.matrix() %>%
      t()
    bac_reg_decom_taxa <-
      tax_table(bac_reg_physeq) %>% data.frame()
    bac_reg_neg_flag <- sample_data(bac_reg_physeq) %>%
      data.frame() %>%
      mutate(neg_flag = sampleGroup %in% c("ENV-Ctrl", "EXT-NEG", "PCR-NEG")) %>%
      pull(neg_flag)
    stopifnot(
      sample_data(bac_reg_physeq) %>% data.frame() %>% rownames() == rownames(bac_reg_decom_otu)
    )
    bac_reg_contam.freq <-
      isContaminant(
        bac_reg_decom_otu,
        method = "prevalence",
        threshold = 0.1,
        neg = bac_reg_neg_flag
      ) %>% merge(bac_reg_decom_taxa, by = 0)
    bac_selected_reg_contam <-
      bac_reg_contam.freq %>% filter(contaminant)
    
    bac_decomed_physeq <-
      prune_taxa(
        !taxa_names(bac_physeq) %in% bac_selected_reg_contam$Row.names,
        bac_physeq
      ) %>%
      filter_taxa(function(x)
        sum(x > 0) > 0, TRUE)
    bac_decomed_physeq <- tax_glom(bac_decomed_physeq, taxrank = "Genus")
    ## remove unneeded taxa
    taxa_df <- tax_table(bac_decomed_physeq)
    tax_table(bac_decomed_physeq) <-
      taxa_df[, colSums(is.na(taxa_df)) < nrow(taxa_df)]
    ## change OTU ID
    taxa_names(bac_decomed_physeq) <-
      tax_table(bac_decomed_physeq) %>% as.data.frame %>% unite("OTU", rank_names(bac_decomed_physeq)) %>% pull(OTU)
    if (output_dir_path != "") {
      write_csv(bac_reg_contam.freq,
                file.path(output_dir_path, "gut_ext_contam.csv"))
      write.csv(
        as.data.frame(otu_table(bac_decomed_physeq)),
        file.path(output_dir_path, "gut_decontamed_genus.csv")
      )
    }
    
    return(bac_decomed_physeq)
  }


# preprocess metadata, change group id to group name
preprocess_lung_nanopore_metadata <- function(metadata_path, lung_short_read_meta_df) {
  group_name_map_df <-
    tibble(
      groupName_long = c(
        "PBS",
        "PBS_Cortisone acetate",
        "PBS_Cortisone acetate_Voriconazol",
        "A fumigatus wt_Cortisone acetate",
        "A fumigatus wt_Cortisone acetate_Voriconazol",
        "A fumigatus delta_pptA_Corisone acetate"
      ),
      groupName_short = c(
        "CTR",
        "CA",
        "VOR",
        "WT",
        "Tx",
        "PptA"
      )
    )
  nanopore_meta_df <- read_csv(metadata_path) %>%
    dplyr::mutate(nanopore_sample = paste0(LIBpool,"_",seqID)) %>% 
    left_join(lung_short_read_meta_df[c("sample","group","sampleTreatment","sampleName","sampleType","sampleGroup","sampleGroup2")], by="sampleName")
  nanopore_meta_df$sample[is.na(nanopore_meta_df$sample)] <- nanopore_meta_df$sampleName[is.na(nanopore_meta_df$sample)] 
  nanopore_meta_df$sampleType[str_detect(nanopore_meta_df$sample, "EXTNEG")] <- "EXT-NEG"
  nanopore_meta_df$sampleType[str_detect(nanopore_meta_df$sample, "EXTNEG")] <- "EXT-NEG"
  nanopore_meta_df$sampleGroup[str_detect(nanopore_meta_df$sample, "EXTNEG")] <-   nanopore_meta_df$sample[str_detect(nanopore_meta_df$sample, "EXTNEG")] 
  nanopore_meta_df$sampleGroup2[str_detect(nanopore_meta_df$sample, "EXTNEG")] <- "EXT-NEG"
  
  nanopore_meta_df$sample[is.na(nanopore_meta_df$sample)] <- nanopore_meta_df$sampleName[is.na(nanopore_meta_df$sample)] 
  nanopore_meta_df$sampleType[str_detect(nanopore_meta_df$sample, "PCRNEG")] <- "PCR-NEG"
  nanopore_meta_df$sampleType[str_detect(nanopore_meta_df$sample, "PCRNEG")] <- "PCR-NEG"
  nanopore_meta_df$sampleGroup[str_detect(nanopore_meta_df$sample, "PCRNEG")] <-   nanopore_meta_df$sample[str_detect(nanopore_meta_df$sample, "PCRNEG")] 
  nanopore_meta_df$sampleGroup2[str_detect(nanopore_meta_df$sample, "PCRNEG")] <- "PCR-NEG"
  
  return(nanopore_meta_df)
}


# load lung nanopore phyloseq object and change OTU ID
get_lung_nanopore_physeq <- function(otu_path, meta_df) {
  otu_df <- read_tsv(otu_path) %>% as.data.frame()
  rownames(otu_df) <- otu_df %>% mutate(OTU = paste(superkingdom, phylum, class, order, family, genus, species, sep = "_")) %>% pull(OTU)
  taxa_df <- otu_df %>% dplyr::select(superkingdom, phylum, class, order, family, genus, species)
  colnames(taxa_df) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  otu_df <- otu_df %>% dplyr::select(-superkingdom, -phylum, -class, -order, -family, -genus, -species)
  otu_df[is.na(otu_df)] <- 0
  colnames(otu_df) <- meta_df[match(colnames(otu_df), meta_df$nanopore_sample),]$sample
  meta_df <- data.frame(meta_df)
  rownames(meta_df) <- meta_df$sample
  
  META<- sample_data(meta_df)
  TAX <- tax_table(as.matrix(taxa_df))
  OTU <- otu_table(otu_df,taxa_are_rows=T)
  physeq <-phyloseq(META, TAX, OTU)
  
  physeq <- physeq %>%
    subset_taxa(Domain=="Bacteria") %>%
    subset_taxa(Genus != "Photobacterium") %>%
    subset_taxa(Genus != "Methanothermobacter") %>%
    subset_taxa(Genus != "Bradyrhizobium") %>%
    #subset_taxa(Genus != "Paracoccus") %>%
    subset_taxa(Genus != "Vulcaniibacterium") %>%
    subset_taxa(Genus != "Kineococcus") %>%
    #subset_taxa(Genus != "Klenkia") %>%
    subset_taxa(Genus != "Buchnera") %>%
    subset_taxa(Genus != "Streptophyta") %>%
    subset_taxa(Genus != "Sporolactobacillaceae_incertae_sedis") %>%
    subset_taxa(Genus != "Spirosoma") %>%
    subset_taxa(Genus != "Dyadobacter") %>%
    subset_taxa(Genus != "Candidatus Cloacamonas") %>%
    subset_taxa(Genus != "Candidatus Xiphinematobacter")
  
  new_sample_data_df <- sample_data(physeq) %>%
    data.frame() %>%
    filter(!sampleGroup %in% c(7,8)) %>%
    filter(!sampleType %in% c("PCR-POS"))
  sample_data(physeq) <- new_sample_data_df
  physeq <- physeq  %>%
    filter_taxa(function(x)
      sum(x > 0) > 0, TRUE) %>% 
    prune_samples(sample_sums(.)>0, .)
  return(physeq)
}


get_decontaminated_lung_nanopore_physeq <-
  function(physeq, output_dir_path = "") {
    bac_physeq <- physeq %>%
      subset_taxa(Domain %in% c("Archaea", "Bacteria"))
    if(!is.null(phy_tree(physeq, F))){
      bac_physeq <- bac_physeq %>%
        set_root_for_phyloseq()
    } 
    
    ## remove samples of pcr-pos,and negative control for SPF mice(sampleName is C1)
    bac_physeq %<>% subset_samples(!sampleType %in% c("PCR-POS")) %>%
      subset_samples(!sampleName %in% c("C1"))
    ## EXT
    bac_ext_physeq <-
      bac_physeq %>% subset_samples(!sampleType %in% c("PCR-NEG","ENV-Ctrl","7","8"))
    bac_ext_decom_otu <- otu_table(bac_ext_physeq) %>%
      as.matrix() %>%
      t()
    bac_ext_decom_taxa <-
      tax_table(bac_ext_physeq) %>% data.frame()
    bac_ext_neg_flag <- sample_data(bac_ext_physeq) %>%
      data.frame() %>%
      mutate(neg_flag = sampleType %in% c("ENV-Ctrl", "EXT-NEG", "PCR-NEG")) %>%
      pull(neg_flag)
    bac_ext_conc <-  sample_data(bac_ext_physeq) %>%
      data.frame() %>%
      pull(concLib)
    stopifnot(
      sample_data(bac_ext_physeq) %>% data.frame() %>% rownames() == rownames(bac_ext_decom_otu)
    )
    bac_ext_otu <- bac_ext_decom_otu[bac_ext_neg_flag,] %>% t
    bac_ext_contam_by_ext <- rowSums(bac_ext_otu>10) > 0
    bac_ext_contam.freq <-
      isContaminant(
        bac_ext_decom_otu,
        conc = bac_ext_conc,
        method = "either",
        threshold = 0.1,
        neg = bac_ext_neg_flag
      ) %>% merge(bac_ext_decom_taxa, by = 0)
    bac_ext_contam_by_ext <- bac_ext_contam_by_ext[bac_ext_contam.freq$Row.names]
    stopifnot(bac_ext_contam.freq$Row.names == names(bac_ext_contam_by_ext))
    bac_ext_contam.freq$contaminant <-   bac_ext_contam.freq$contaminant | bac_ext_contam_by_ext
    ## PCR
    bac_pcr_physeq <-
      bac_physeq %>% subset_samples(!sampleType %in% c("EXT-NEG","ENV-Ctrl","7","8"))
    bac_pcr_decom_otu <- otu_table(bac_pcr_physeq) %>%
      as.matrix() %>%
      t()
    bac_pcr_decom_taxa <-
      tax_table(bac_pcr_physeq) %>% data.frame()
    bac_pcr_neg_flag <- sample_data(bac_pcr_physeq) %>%
      data.frame() %>%
      mutate(neg_flag = sampleType %in% c("ENV-Ctrl", "EXT-NEG", "PCR-NEG")) %>%
      pull(neg_flag)
    if(sum(bac_pcr_neg_flag>0)){
      bac_pcr_conc <-  sample_data(bac_pcr_physeq) %>%
        data.frame() %>%
        pull(concLib)
      stopifnot(
        sample_data(bac_pcr_physeq) %>% data.frame() %>% rownames() == rownames(bac_pcr_decom_otu)
      )
      bac_pcr_otu <- bac_pcr_decom_otu[bac_pcr_neg_flag,] %>% t
      bac_pcr_contam_by_pcr <- rowSums(bac_pcr_otu>10) > 0
      bac_pcr_contam.freq <-
        isContaminant(
          bac_pcr_decom_otu,
          conc = bac_pcr_conc,
          method = "either",
          threshold = 0.1,
          neg = bac_pcr_neg_flag
        ) %>% merge(bac_pcr_decom_taxa, by = 0)
      bac_pcr_contam_by_pcr <- bac_pcr_contam_by_pcr[bac_pcr_contam.freq$Row.names]
      stopifnot(bac_pcr_contam.freq$Row.names == names(bac_pcr_contam_by_pcr))
      bac_pcr_contam.freq$contaminant <-   bac_pcr_contam.freq$contaminant | bac_pcr_contam_by_pcr
    } else{
      bac_pcr_contam.freq <- bac_ext_contam.freq %>% dplyr::select(Row.names,contaminant)
      bac_pcr_contam.freq$contaminant <- FALSE
    }
    bac_selected_reg_contam <-
      bind_rows(bac_ext_contam.freq %>% filter(contaminant), bac_pcr_contam.freq %>% filter(contaminant))
    bac_env_physeq <-
      bac_physeq %>% subset_samples(!sampleType %in% c("EXT-NEG", "PCR-NEG", "7","8"))
    bac_env_decom_otu <- otu_table(bac_env_physeq) %>%
      as.matrix() %>%
      t()
    bac_env_decom_taxa <-
      tax_table(bac_env_physeq) %>% data.frame()
    bac_env_neg_flag <- sample_data(bac_env_physeq) %>%
      data.frame() %>%
      mutate(neg_flag = sampleType %in% c("ENV-Ctrl", "EXT-NEG", "PCR-NEG")) %>%
      pull(neg_flag)
    bac_env_conc <- sample_data(bac_env_physeq) %>%
      data.frame() %>%
      pull(concLib)
    stopifnot(
      sample_data(bac_env_physeq) %>% data.frame() %>% rownames() == rownames(bac_env_decom_otu)
    )
    bac_env_otu <- bac_env_decom_otu[bac_env_neg_flag,] %>% t
    bac_env_contam_by_env <- rowSums(bac_env_otu>10) > 0
    bac_env_contam.freq <-
      isContaminant(
        bac_env_decom_otu,
        conc=bac_env_conc,
        method = "either",
        threshold = 0.1,
        neg = bac_env_neg_flag
      ) %>% merge(bac_env_decom_taxa, by = 0)
    bac_env_contam_by_env <- bac_env_contam_by_env[bac_env_contam.freq$Row.names]
    stopifnot(bac_env_contam.freq$Row.names == names(bac_env_contam_by_env))
    bac_env_contam.freq$contaminant <-   bac_env_contam.freq$contaminant | bac_env_contam_by_env
    bac_selected_env_contam <-
      bac_env_contam.freq %>% filter(contaminant)
    
    bac_decomed_physeq <-
      prune_taxa(
        !taxa_names(bac_physeq) %in% c(
          bac_selected_env_contam$Row.names,
          bac_selected_reg_contam$Row.names
        ),
        bac_physeq
      ) %>%
      filter_taxa(function(x)
        sum(x > 0) > 0, TRUE)
    ## remove unneeded taxa
    taxa_df <- tax_table(bac_decomed_physeq)
    tax_table(bac_decomed_physeq) <-
      taxa_df[, colSums(is.na(taxa_df)) < nrow(taxa_df)]
    ## change OTU ID
    taxa_names(bac_decomed_physeq) <-
      tax_table(bac_decomed_physeq) %>% as.data.frame %>% unite("OTU", rank_names(bac_decomed_physeq)) %>% pull(OTU)
    bac_decomed_physeq <- bac_decomed_physeq %>%  prune_samples(sample_sums(.)>0,.)
    bac_decomed_physeq <- bac_decomed_physeq %>% transform_sample_counts(function(x) x/sum(x))
      
    if (output_dir_path != "") {
      write_csv(bac_ext_contam.freq,
                file.path(output_dir_path, "lung_nanopore_ext_contam.csv"))
      write_csv(bac_pcr_contam.freq,
                file.path(output_dir_path, "lung_nanopore_pcr_contam.csv"))
      write_csv(bac_env_contam.freq,
                file.path(output_dir_path, "lung_nanopore_env_contam.csv"))
      write.csv(
        as.data.frame(otu_table(bac_decomed_physeq)),
        file.path(output_dir_path, "lung_nanopore_decontamed_genus.csv")
      )
    }
    
    return(bac_decomed_physeq)
  }


# get case samples data
get_case_physeq <- function(physeq) {
  case_physeq <- physeq %>%
    subset_samples(sampleGroup %in% 1:6) %>%
    filter_taxa(function(x)
      sum(x > 0) > 0, TRUE)
  new_sample_data_df <- sample_data(case_physeq) %>%
    data.frame() %>%
    mutate(group = factor(group, levels = c(
      "CTR",
      "CA",
      "VOR",
      "WT",
      "Tx",
      "PptA"
    )))
  sample_data(case_physeq) <- new_sample_data_df 
  return(case_physeq)
}

# get control samples data
get_control_physeq <- function(physeq) {
  control_physeq <- physeq %>%
    subset_samples(sampleGroup %in% 7:8) %>%
    filter_taxa(function(x)
      sum(x > 0) > 0, TRUE)
  return(control_physeq)
}

# get negative_control samples data
get_negative_control_physeq <- function(physeq) {
  control_physeq <- physeq %>%
    subset_samples(!sampleGroup %in% 1:8) %>%
    filter_taxa(function(x)
      sum(x > 0) > 0, TRUE)
  return(control_physeq)
}
# preprocess metadata, change group id to group name
preprocess_metadata <- function(metadata_path) {
  group_name_map_df <-
    tibble(
      groupName_long = c(
        "PBS",
        "PBS_Cortisone acetate",
        "PBS_Cortisone acetate_Voriconazol",
        "A fumigatus wt_Cortisone acetate",
        "A fumigatus wt_Cortisone acetate_Voriconazol",
        "A fumigatus delta_pptA_Corisone acetate"
      ),
      groupName_short = c(
        "CTR",
        "CA",
        "VOR",
        "WT",
        "Tx",
        "PptA"
      )
    )
  meta_df <- read_csv(metadata_path) %>%
    # filter(sampleGroup %in% 1:6) %>%
    mutate(seqID = paste0("S_", seqID)) %>%
    dplyr::rename(sample = seqID) %>%
    mutate(group = group_name_map_df[match(.$sampleTreatment, group_name_map_df$groupName_long),]$groupName_short)
  return(meta_df)
}


subset_phyloseq_by_comparison <- function(physeq, G1, G2) {
  new_bac_sample_data_df <- sample_data(physeq) %>%
    data.frame() %>%
    filter(group %in% c(G1, G2)) %>%
    mutate(group = factor(group, levels = c(G1, G2)))
  sample_data(physeq) <- new_bac_sample_data_df
  return(physeq)
}

get_omics_data <- function(lung_physeq, gut_physeq, lung_metabolome_physeq, plasma_metabolome_physeq, abundance_cutoff=0,prevalence_cutoff=0, transform = ""){
  source("R/utils.R")
  if(abundance_cutoff!=0|prevalence_cutoff!=0){
    lung_physeq <- lung_physeq %>% utils_filter_taxa(abundance_cutoff, prevalence_cutoff)
    gut_physeq <- gut_physeq %>% utils_filter_taxa(abundance_cutoff, prevalence_cutoff)
  }
  if(transform!=""){
    lung_physeq <- microbiome::transform(lung_physeq,transform=transform)
    gut_physeq <- microbiome::transform(gut_physeq,transform=transform)
  }
  lung_meta_df <- sample_data(lung_physeq) %>% data.frame()
  gut_meta_df <- sample_data(gut_physeq) %>% data.frame()
  lung_metabolome_meta_df <- sample_data(lung_metabolome_physeq) %>% data.frame()
  plasma_metabolome_meta_df <- sample_data(plasma_metabolome_physeq) %>% data.frame()
  
  kept_samples <- lung_meta_df %>% filter(sampleGroup %in% 1:6) %>% 
    pull(sampleName) %>%
    intersect(gut_meta_df$sampleName) %>%
    intersect(lung_metabolome_meta_df$sampleName) %>%
    intersect(plasma_metabolome_meta_df$sampleName)
  id_map <- lung_meta_df %>%
    dplyr::select(sampleName) %>%
    rownames_to_column(var="id1") %>%
    left_join(gut_meta_df %>% rownames_to_column(var="id2") , by="sampleName") %>%
    dplyr::select(id1, id2, sampleName,sampleGroup,group) %>%
    filter(sampleName %in% kept_samples)
  
  processed_lung_physeq <-  prune_samples(rownames(sample_data(lung_physeq)) %in% id_map$id1 ,lung_physeq) %>% set_sample_order("sampleName")
  processed_gut_physeq <-  prune_samples(rownames(sample_data(gut_physeq)) %in% id_map$id2 ,gut_physeq) %>% set_sample_order("sampleName")
  processed_lung_metabolome_physeq <-  prune_samples(get_variable(lung_metabolome_physeq,"sampleName") %in% id_map$sampleName ,lung_metabolome_physeq) %>% set_sample_order("sampleName")
  processed_plasma_metabolome_physeq <-  prune_samples(get_variable(plasma_metabolome_physeq,"sampleName") %in% id_map$sampleName ,plasma_metabolome_physeq) %>% set_sample_order("sampleName")
  stopifnot(get_variable(processed_lung_physeq,"sampleName") == get_variable(processed_gut_physeq,"sampleName"))
  stopifnot(get_variable(processed_lung_physeq,"sampleName")== get_variable(processed_lung_metabolome_physeq,"sampleName"))
  stopifnot(get_variable(processed_lung_physeq,"sampleName") == get_variable(processed_plasma_metabolome_physeq,"sampleName"))
  out_tbl <- tibble(name=c("gut_microbiome","lung_microbiome","lung_metabolome","plasma_metabolome"),
                    data=list(processed_gut_physeq,processed_lung_physeq,processed_lung_metabolome_physeq,processed_plasma_metabolome_physeq),
                    meta=list(id_map))
}


get_distribution_of_decontamination <- function(pre_physeq, pos_physeq) {
  n=10
  
  pre_physeq <- pre_physeq %>% 
    transform_sample_counts(function(OTU)OTU / sum(OTU)*100)
  pre_case_physeq <- pre_physeq %>% 
    subset_samples(sampleGroup %in% 1:6)
  pre_bs_physeq <- pre_physeq %>% 
    subset_samples(sampleType %in% c("ENV-Ctrl", "EXT-NEG", "PCR-NEG")) %>% 
    subset_samples(sampleName != "C1")
  pre_case_mean_abundance <- pre_case_physeq %>% 
    melt_phyloseq() %>% 
    replace(is.na(.), 0) %>% 
    dplyr::select(Abundance,Genus) %>% 
    group_by(Genus) %>% 
    summarise(case_mean_abundance=mean(Abundance))
  pre_bs_mean_abundance <- pre_bs_physeq %>% 
    melt_phyloseq() %>% 
    replace(is.na(.), 0) %>% 
    dplyr::select(Abundance,Genus) %>% 
    group_by(Genus) %>% 
    summarise(bs_mean_abundance=mean(Abundance))
  pre_all_mean_abundance <- full_join(pre_case_mean_abundance,pre_bs_mean_abundance) %>% 
    replace(is.na(.), 0)
  pre_ranked_case_barplot <- pre_all_mean_abundance %>% 
    slice_max(case_mean_abundance,n=n) %>%
    mutate(Genus = fct_reorder(Genus,case_mean_abundance)) %>% 
    pivot_longer(!Genus,values_to = "abundance",names_to = "group") %>% 
    mutate(group = ifelse(str_detect(group,"case"),"Lung samples","Blank control")) %>% 
    ggplot(aes(x=Genus, y=abundance))+
    geom_bar(stat="identity") +
    theme_pubr()+
    theme(axis.text.x = element_text(angle = 60, vjust = 0.5),
          plot.title = element_text(hjust = 0.5))+
    facet_wrap(vars(group),scales="free_y",ncol=1) +
    ylab("Abundance (%)")+
    xlab("")+
    labs(title="Ranked by tissue specimens before decontamination")
  pre_ranked_bs_barplot <- pre_all_mean_abundance %>% 
    slice_max(bs_mean_abundance,n=n) %>%
    mutate(Genus = fct_reorder(Genus,bs_mean_abundance)) %>% 
    pivot_longer(!Genus,values_to = "abundance",names_to = "group") %>% 
    mutate(group = ifelse(str_detect(group,"case"),"Lung samples","Blank control")) %>% 
    ggplot(aes(x=Genus, y=abundance))+
    geom_bar(stat="identity") +
    theme_pubr()+
    theme(axis.text.x = element_text(angle = 60, vjust = 0.5),
          plot.title = element_text(hjust = 0.5))+
    facet_wrap(vars(group),scales="free_y",ncol=1)+
    ylab("Abundance (%)")+
    xlab("")+
    labs(title="Ranked by blank control before decontamination")
  
  
  pos_physeq <- pos_physeq %>% 
    transform_sample_counts(function(OTU)OTU / sum(OTU)*100)
  pos_case_physeq <- pos_physeq %>% 
    subset_samples(sampleGroup %in% 1:6)
  pos_bs_physeq <- pos_physeq %>% 
    subset_samples(sampleType %in% c("ENV-Ctrl", "EXT-NEG", "PCR-NEG")) %>% 
    subset_samples(sampleName != "C1")
  pos_case_mean_abundance <- pos_case_physeq %>% 
    melt_phyloseq() %>% 
    replace(is.na(.), 0) %>% 
    dplyr::select(Abundance,Genus) %>% 
    group_by(Genus) %>% 
    summarise(case_mean_abundance=mean(Abundance))
  pos_bs_mean_abundance <- pos_bs_physeq %>% 
    melt_phyloseq() %>% 
    replace(is.na(.), 0) %>% 
    dplyr::select(Abundance,Genus) %>% 
    group_by(Genus) %>% 
    summarise(bs_mean_abundance=mean(Abundance))
  pos_all_mean_abundance <- full_join(pos_case_mean_abundance,pos_bs_mean_abundance) %>% 
    replace(is.na(.), 0)
  pos_ranked_case_barplot <- pre_all_mean_abundance %>% 
    filter(Genus %in% 
             (pos_all_mean_abundance %>% slice_max(case_mean_abundance,n=n) %>% pull(Genus))) %>%
    mutate(Genus = factor(Genus, levels =pos_all_mean_abundance %>% slice_max(case_mean_abundance,n=n) %>% pull(Genus) ),
           Genus = fct_rev(Genus)) %>%
    pivot_longer(!Genus,values_to = "abundance",names_to = "group") %>% 
    mutate(group = ifelse(str_detect(group,"case"),"Lung samples","Blank control")) %>% 
    ggplot(aes(x=Genus, y=abundance))+
    geom_bar(stat="identity") +
    theme_pubr()+
    ylab("Abundance (%)")+
    xlab("")+
    theme(axis.text.x = element_text(angle = 60, vjust = 0.5),
          plot.title = element_text(hjust = 0.5))+
    facet_wrap(vars(group),scales="free_y",ncol=1) +
    ylab("Abundance (%)")+
    labs(title="Ranked by tissue specimens after decontamination")
  pos_ranked_bs_barplot <-  pre_all_mean_abundance %>% 
    filter(Genus %in% 
             (pos_all_mean_abundance %>% slice_max(bs_mean_abundance,n=n) %>% pull(Genus))) %>% 
    mutate(Genus = factor(Genus, levels =pos_all_mean_abundance %>% slice_max(bs_mean_abundance,n=n) %>% pull(Genus) ),
           Genus = fct_rev(Genus)) %>%
    pivot_longer(!Genus,values_to = "abundance",names_to = "group") %>% 
    mutate(group = ifelse(str_detect(group,"case"),"Lung samples","Blank control")) %>% 
    ggplot(aes(x=Genus, y=abundance))+
    geom_bar(stat="identity") +
    theme_pubr()+
    theme(axis.text.x = element_text(angle = 60, vjust = 0.5),
          plot.title = element_text(hjust = 0.5))+
    facet_wrap(vars(group),scales="free_y",ncol=1)+
    ylab("Abundance (%)")+
    xlab("")+
    labs(title="Ranked by blank control after decontamination")
  
  
  wrap_plots(pre_ranked_case_barplot,pos_ranked_case_barplot,pre_ranked_bs_barplot,pos_ranked_bs_barplot)
}
