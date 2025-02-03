get_microbial_loading <- function(physeq, output_dir_path) {
  otu_df <- otu_table(physeq) %>% as.data.frame()
  meta_df <-
    sample_data(physeq) %>% data.frame() %>% rownames_to_column()

  ## read count by sample
  loading_by_sample <-
    colSums(otu_df) %>% as.data.frame %>% rownames_to_column()
  colnames(loading_by_sample) <- c("rowname", "read_count")
  loading_by_sample <-
    left_join(loading_by_sample, meta_df, by = "rowname") %>%
    dplyr::select(rowname, read_count, sampleName, sampleGroup2)

  return(loading_by_sample)
}

get_otu_before_decontam <- function(physeq){
  otu_long_df <- melt_phyloseq(physeq)
  otu_df <- otu_long_df %>%
    dplyr::select(Sample, OTU, Abundance, Domain,Phylum, Class,Order,Family,Genus,Species) %>%
    pivot_wider(names_from = Sample, values_from = Abundance)
  return(otu_df)
}

get_genus_overlap_of_case_control <-
  function(case_physeq,
           control_physeq,
           output_dir_path) {
    case_otu_df <- otu_table(case_physeq) %>% as.data.frame()
    case_meta_df <-
      sample_data(case_physeq) %>% data.frame() %>% rownames_to_column()

    control_otu_df <- otu_table(control_physeq) %>% as.data.frame()
    control_meta_df <-
      sample_data(control_physeq) %>% data.frame() %>% rownames_to_column()

    control_g7_otu_df <-
      otu_table(
        control_physeq %>% subset_samples(sampleGroup %in% c(7)) %>%
          filter_taxa(function(x)
            sum(x > 0) > 0, TRUE)
      ) %>%
      as.data.frame()

    control_g8_otu_df <-
      otu_table(
        control_physeq %>% subset_samples(sampleGroup %in% c(8)) %>%
          filter_taxa(function(x)
            sum(x > 0) > 0, TRUE)
      ) %>%
      as.data.frame()

    x <-
      list(
        G1_to_6 = rownames(case_otu_df),
        G7 = rownames(control_g7_otu_df),
        G8 = rownames(control_g8_otu_df)
      )
    out_vennplot <-
      ggVennDiagram(x, set_size = 5) + scale_fill_distiller(palette = "Reds", direction = 1)
    ggsave(
      file.path(output_dir_path, "genus_overlap_of_case_control.pdf"),
      out_vennplot,
      width = 8,
      height = 8
    )
  }

get_hetamap_detected_taxa <- function(physeq, count=T){
  ## remove unneeded taxa
  taxa_df <- tax_table(physeq)
  tax_table(physeq) <-
    taxa_df[, colSums(is.na(taxa_df)) < nrow(taxa_df)]
  physeq <- physeq %>% subset_samples(!sampleType %in% c("PCR-POS")) %>%
    subset_samples(!sampleName %in% c("C1"))
  meta_df <- data.frame(sample_data(physeq)) %>%  mutate(sampleGroup2=factor(sampleGroup2) ) %>%arrange(sampleGroup2)

  taxa_df <- data.frame(tax_table(physeq))
  otu_df <- as.data.frame(otu_table(physeq)) %>% merge(taxa_df, by=0) %>% dplyr::select(-`Row.names`,-Domain,-Phylum,-Class,-Order,-Family)
  otu_mat <- otu_df %>% dplyr::select(-Genus, -Species) %>% data.matrix()
  rownames(otu_mat) <- paste0(otu_df$Genus,"_", otu_df$Species)
  otu_mat <- otu_mat[,rownames(meta_df)]
  stopifnot(rownames(meta_df) == colnames(otu_mat))
  otu_mat <- otu_mat[order(rowSums(otu_mat>0),decreasing=T),][1:100,]
  if(count){
  col_fun <-
    colorRamp2(c(0, 100), c("#ffffff", "#FF0000"))
  }else{
    col_fun <-
      colorRamp2(c(0, 0.1), c("#ffffff", "#FF0000"))
  }
  split_vec <- meta_df$sampleGroup2
  ha <- HeatmapAnnotation(
    group = anno_block(gp = gpar(fill = 1:length(levels(split_vec))), labels = levels(split_vec))
  )
  heatmap_plot <-
    Heatmap(
      name="read count",
      otu_mat,
      col = col_fun,
      cluster_rows = F,
      cluster_columns = F,
      top_annotation = ha,column_split = split_vec,
      column_names_gp = gpar(fontsize=10),
      row_names_gp = gpar(fontsize=10)
      )
  return(heatmap_plot)
}


get_taxa_composition_plot <- function(physeq, n_taxa=40,tax_level="Genus",group_by="group"){
  plot_list <- physeq  %>%
    transform_sample_counts(function(x) x/sum(x)) %>%
    subset_taxa(!is.na(Domain)) %>%
    tax_fix(unknowns = c("?")) %>%
    comp_barplot(n_taxa = n_taxa, tax_level = tax_level, group_by = group_by)
  patch <- wrap_plots(plot_list, nrow = 2, guides = "collect")
  return(patch)
}