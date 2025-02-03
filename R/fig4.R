get_fig4_comparison_df <- function() {
  comparison_df <-data.frame(
    G1 = c(
      "WT",
      "VOR"

    ),
    G2 = c(
      "Tx",
      "Tx"
    )
  )
  comparison_df$comparison <-
    paste(comparison_df$G1, "vs", comparison_df$G2)
  rownames(comparison_df) <- comparison_df$comparison
  return(comparison_df)
}


# get all fig4 samples data
get_fig4_physeq <- function(physeq) {
  filtered_physeq <- physeq %>%
    subset_samples(sampleGroup %in% c(3,4,5)) %>%
    filter_taxa(function(x)
      sum(x > 0) > 0, TRUE)
  return(filtered_physeq)
}


# get extended fig4 samples data
get_extended_fig4_physeq <- function(physeq) {
  filtered_physeq <- physeq %>%
    subset_samples(sampleGroup %in% c(2,3,4,5)) %>%
    filter_taxa(function(x)
      sum(x > 0) > 0, TRUE)
  return(filtered_physeq)
}


plot_fig4_alpha <- function(physeq,plot_title=""){
  alpha_diversity_data <- calculate_alpha_diversity_data(physeq)
  color <- get_colormap_and_charmap() %>%
    filter(group %in% alpha_diversity_data$group) %>% pull(color)
  # chao1
  chao1_result <- alpha_diversity_data %>%
    dplyr::select(sample, group, chao1) %>%
    arrange(group)
  
  chao_alpha_boxplot <- ggboxplot(
    chao1_result,
    x = "group",
    y = "chao1",
    add = "boxplot",
    fill = "group",
  ) + stat_compare_means(method = "wilcox.test",comparisons=list(c("VOR","Tx"),c("WT","VOR"),c("WT","Tx")),label = "p.format",
                         label.x = 1.5,
                         size = 6.4) + theme_pubr() + ylab("Chao index") + xlab("") +
    theme(plot.title = element_text(size = 1, hjust = 0.5),
          text = element_text(22),
          legend.title = element_blank(),
          legend.key.size = unit(0.7, 'cm'),
          legend.spacing.x = unit(0.5, 'cm'),
          legend.text = element_text(size=22),
          axis.title = element_text(size = 27),
          axis.text.x=element_blank(),
          axis.ticks.length=unit(.1, "cm"),
          axis.text=element_text(size=22))+
    scale_fill_manual(values=color)
  
  # shannon
  shannon_result <- alpha_diversity_data %>%
    dplyr::select(sample, group, diversity_shannon) %>%
    arrange(group)
  
  shannon__alpha_boxplot <- ggboxplot(
    shannon_result,
    x = "group",
    y = "diversity_shannon",
    add = "boxplot",
    fill = "group",
  ) + stat_compare_means(method = "wilcox.test",comparisons=list(c("VOR","Tx"),c("WT","VOR"),c("WT","Tx")),label = "p.format",
                         label.x = 1.5,
                         size = 6.4) + theme_pubr() + ylab("Shannon index") + xlab("") +
    theme(plot.title = element_text(size = 1, hjust = 0.5),
          text = element_text(22), legend.title = element_blank() ,
          axis.title = element_text(size = 27),
          axis.text.x=element_blank(),
          legend.position="none",
          axis.ticks.length=unit(0.1, "cm"),
          axis.text=element_text(size=24),
          legend.text=element_text(size=22))+
    scale_fill_manual(values=color)
  
  
  out_plot <- ggarrange(chao_alpha_boxplot,shannon__alpha_boxplot, common.legend = T,align="hv",legend="top")
  out_plot <- annotate_figure(out_plot, top = text_grob(plot_title, size = 30))
  return(out_plot)
}


plot_fig4_microbiome_distance <- function(physeq, plot_title="") {
  distance = "unifrac"
  weighted = T
  meta_df <- data.frame(sample_data(physeq))
  ord_mat <-
    ordinate(physeq,
             method = "PCoA",
             distance = distance,
             weighted = weighted)
  dist_mat1 <-
    phyloseq::distance(physeq %>% subset_samples(group %in% c("WT","Tx")), method = distance, weighted = weighted)
  meta_df1 <-  data.frame(sample_data(physeq %>% subset_samples(group %in% c("WT","Tx"))))
  stopifnot(all(rownames(meta_df1)==names(dist_mat1)))
  dist_to_centroids_df1 <- dist_to_centroids(dist_mat1, meta_df1$group) %>%
    mutate(CentroidGroup ="WT vs Tx")

  dist_mat2 <-
    phyloseq::distance(physeq %>% subset_samples(group %in% c("VOR","Tx")), method = distance, weighted = weighted)
  meta_df2<-  data.frame(sample_data(physeq %>% subset_samples(group %in% c("VOR","Tx"))))
  stopifnot(all(rownames(meta_df2)==names(dist_mat2)))
  dist_to_centroids_df2 <- dist_to_centroids(dist_mat2, meta_df2$group) %>%
    mutate(CentroidGroup ="VOR vs Tx")

  dist_to_centroids_df <- bind_rows(dist_to_centroids_df1, dist_to_centroids_df2) %>%
    dplyr::rename("group"="CentroidGroup") %>%
    dplyr::rename("distance"="CentroidDistance") %>%
    mutate(group = factor(group))

  p <- wilcox_test(dist_to_centroids_df, distance ~ group)$p
    if(p<0.001){
    p <- format(p, scientific = T)
    }else {
       p = round(p, 3)
    }
  statistic <- wilcox_test(dist_to_centroids_df, distance ~ group)$statistic
  mu <- ddply(dist_to_centroids_df, "group", summarise, grp.median=median(distance))
  beta_boxplot <-ggplot(dist_to_centroids_df, aes(x=distance, color=group, fill=group)) +
    geom_histogram(aes(y=after_stat(density)), alpha=0.5, bins = 20,
                   position="identity")+
    geom_density(alpha=0.2)+
    geom_vline(data=mu, aes(xintercept=grp.median, color=group),
               linetype="dashed",size=1.5)+
    theme_pubr() +
    theme(
      plot.title = element_text(size = 30, hjust = 0.5),
      plot.subtitle = element_markdown(size = 22, hjust = 0.5),
      text = element_text(size = 24),
      axis.title=element_text(size=28),
      legend.text = element_text(size = 28),
      legend.title = element_blank(),
      legend.position = "top") +
    labs(subtitle = paste("p =", p), title = plot_title) +
    xlab("Distance") +
    ylab("Density") +
    scale_fill_aaas()+
    scale_color_aaas()
  return(beta_boxplot)
}


plot_fig4_metabolome_distance  <- function(physeq, plot_title="") {
  meta_df <- data.frame(sample_data(physeq))
  otu_df <- as.data.frame(otu_table(physeq)) %>% rotate_df()

  dist_mat1 <-
    phyloseq::distance(physeq %>% subset_samples(group %in% c("WT","Tx")), method = "euclidean")
  meta_df1 <-  data.frame(sample_data(physeq %>% subset_samples(group %in% c("WT","Tx"))))
  stopifnot(all(rownames(meta_df1)==names(dist_mat1)))
  dist_to_centroids_df1 <- dist_to_centroids(dist_mat1, meta_df1$group) %>%
    mutate(CentroidGroup ="WT vs Tx")

  dist_mat2 <-
    phyloseq::distance(physeq %>% subset_samples(group %in% c("VOR","Tx")), method = "euclidean")
  meta_df2<-  data.frame(sample_data(physeq %>% subset_samples(group %in% c("VOR","Tx"))))
  set.seed(42)
  stopifnot(all(rownames(meta_df2)==names(dist_mat2)))
  dist_to_centroids_df2 <- dist_to_centroids(dist_mat2, meta_df2$group) %>%
    mutate(CentroidGroup ="VOR vs Tx")

  dist_to_centroids_df <- bind_rows(dist_to_centroids_df1, dist_to_centroids_df2) %>%
    dplyr::rename("group"="CentroidGroup") %>%
    dplyr::rename("distance"="CentroidDistance") %>%
    mutate(group = factor(group))

  p <- wilcox_test(dist_to_centroids_df, distance ~ group)$p
  if(p<0.001){
    p <- format(p, scientific = T)
    }else {
       p = round(p, 3)
    }
  statistic <- wilcox_test(dist_to_centroids_df, distance ~ group)$statistic
  mu <- ddply(dist_to_centroids_df, "group", summarise, grp.median=median(distance))
  beta_boxplot <- ggplot(dist_to_centroids_df, aes(x=distance, color=group, fill=group)) +
    geom_histogram(aes(y=after_stat(density)), alpha=0.5, bins = 20,
                   position="identity")+
    geom_density(alpha=0.2)+
    geom_vline(data=mu, aes(xintercept=grp.median, color=group),
               linetype="dashed",size=1.5)+
    theme_pubr() +
    theme(
      plot.title = element_text(size = 30, hjust = 0.5),
      plot.subtitle = element_markdown(size = 22, hjust = 0.5),
      text = element_text(size = 24),
      axis.title=element_text(size=28),
      legend.text = element_text(size = 28),
      legend.title = element_blank(),
      legend.position = "top") +
    labs(subtitle = paste("p =", p), title = plot_title) +
    xlab("Distance") +
    ylab("Density") +
    scale_fill_aaas()+
    scale_color_aaas()
  return(beta_boxplot)
}


plot_fig4_microbiome_volcano <-
  function(lung_meta_sig_df, gut_meta_sig_df) {
    lung_df <-  lung_meta_sig_df %>%
      filter(comparison == "WT vs Tx") %>%
      mutate(label = ifelse(padj<0.1&abs(log2FC)>log2(1.5), Genus, NA)) %>%
      mutate(type = log2FC <0)
    lung_volcano_plot <- ggplot(data = lung_df, aes(x = log2FC, y = -log10(padj), col = type,label = label)) +
      geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.1), col = "gray", linetype = 'dashed') +
      geom_point(size = 2) +
      theme_pubclean() +
      geom_text_repel(max.overlaps = Inf) +
      scale_color_manual(values = c("#00AFBB", "#bb0c00")) +
      labs(
           x = expression("log"[2]*"FC"), y = expression("-log"[10]*"(q-value)"),
           title = "Lung microbiome" )  +
      theme(legend.position = "none")

    gut_df <-  gut_meta_sig_df %>%
      filter(comparison == "WT vs Tx") %>%
      mutate(label = ifelse(padj<0.1&abs(log2FC)>log2(1.5), Genus, NA)) %>%
      mutate(type = log2FC <0)
    gut_volcano_plot <- ggplot(data = gut_df, aes(x = log2FC, y = -log10(padj), col = type,label = label)) +
      geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.1), col = "gray", linetype = 'dashed') +
      geom_point(size = 2) +
      theme_pubclean() +
      geom_text_repel(max.overlaps = Inf) +
      scale_color_manual(values = c("#00AFBB", "#bb0c00")) +
      labs(
        x = expression("log"[2]*"FC"), y = expression("-log"[10]*"(q-value)"),
        title = "Gut microbiome" )  +
      theme(legend.position = "none")
    return(list(lung_volcano_plot, gut_volcano_plot))

  }


plot_fig4_metabolome_volcano <-
  function(lung_meta_sig_df, gut_meta_sig_df) {
    lung_df <-  lung_meta_sig_df %>%
      filter(comparison == "WT vs Tx") %>%
      mutate(label = ifelse(padj<0.1&abs(log2FC)>log2(1.5), metabolite, NA)) %>%
      mutate(type = log2FC <0)
    lung_volcano_plot <- ggplot(data = lung_df, aes(x = log2FC, y = -log10(padj), col = type,label = label)) +
      geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.1), col = "gray", linetype = 'dashed') +
      geom_point(size = 2) +
      theme_pubclean() +
      geom_text_repel(max.overlaps = Inf) +
      scale_color_manual(values = c("#00AFBB", "#bb0c00")) +
      labs(
        x = expression("log"[2]*"FC"), y = expression("-log"[10]*"(q-value)"),
        title = "Lung metabolome" )  +
      theme(legend.position = "none")

    gut_df <-  gut_meta_sig_df %>%
      filter(comparison == "WT vs Tx") %>%
      mutate(label = ifelse(padj<0.1&abs(log2FC)>log2(1.5), metabolite, NA)) %>%
      mutate(type = log2FC <0)
    gut_volcano_plot <- ggplot(data = gut_df, aes(x = log2FC, y = -log10(padj), col = type,label = label)) +
      geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.1), col = "gray", linetype = 'dashed') +
      geom_point(size = 2) +
      theme_pubclean() +
      geom_text_repel(max.overlaps = Inf) +
      scale_color_manual(values = c("#00AFBB", "#bb0c00")) +
      labs(
        x = expression("log"[2]*"FC"), y = expression("-log"[10]*"(q-value)"),
        title = "Plasma metabolome" )  +
      theme(legend.position = "none")
    return(list(lung_volcano_plot, gut_volcano_plot))

  }

plot_fig4_microbiome_venn <-
  function(lung_meta_sig_df, gut_meta_sig_df) {
    library(ggvenn)
    lung_infection_taxa <-  lung_meta_sig_df %>%
      filter(comparison == "CA vs WT") %>%
      filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5)) %>% pull(OTU)
    lung_treatment_taxa <-  lung_meta_sig_df %>%
      filter(comparison == "WT vs Tx") %>%
      filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5)) %>% pull(OTU)
    lung_vor_taxa <-  lung_meta_sig_df %>%
      filter(comparison == "CA vs VOR") %>%
      filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5)) %>% pull(OTU)

    lung_treatment_taxa <- lung_treatment_taxa[!lung_treatment_taxa %in% lung_vor_taxa]
    lung_venn_data <- list_to_data_frame(list("CA vs WT"=lung_infection_taxa, "WT vs Tx"=lung_treatment_taxa))
    lung_venn_plot <- ggplot(lung_venn_data, aes(A=`CA vs WT`,B= `WT vs Tx`))+
    geom_venn(fill_color=c("#7FC07E","#FFD17E"),set_name_size=7,text_size=6)+
      labs(title = "Lung microbiome") +
      theme_void() +
      coord_fixed(clip="off") +
      theme(
      plot.title = element_text(size = 25, hjust = 0.5))

    gut_infection_taxa <-  gut_meta_sig_df %>%
      filter(comparison == "CA vs WT") %>%
      filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5)) %>% pull(OTU)
    gut_treatment_taxa <-  gut_meta_sig_df %>%
      filter(comparison == "WT vs Tx") %>%
      filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5)) %>% pull(OTU)
    gut_vor_taxa <-  gut_meta_sig_df %>%
      filter(comparison == "CA vs VOR") %>%
      filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5)) %>% pull(OTU)
    gut_treatment_taxa <- gut_treatment_taxa[!gut_treatment_taxa %in% gut_vor_taxa]
    gut_venn_data <- list_to_data_frame(list("CA vs WT"=gut_infection_taxa, "WT vs Tx"=gut_treatment_taxa))
    gut_venn_plot <- ggplot(gut_venn_data, aes(A=`CA vs WT`,B= `WT vs Tx`))+
    geom_venn(fill_color=c("#7FC07E","#FFD17E"),set_name_size=7,text_size=6)+
      labs(title = "Gut microbiome") +
      theme_void() +
      coord_fixed(clip="off") +
      theme(
      plot.title = element_text(size = 25, hjust = 0.5))

    return(list(lung_venn_plot, gut_venn_plot))

  }


plot_fig4_metabolome_venn <-
function(lung_meta_sig_df, plasma_meta_sig_df) {
  library(ggvenn)
  lung_infection_taxa <-  lung_meta_sig_df %>%
    filter(comparison == "CA vs WT") %>%
    filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5)) %>% pull(metabolite)
  lung_treatment_taxa <-  lung_meta_sig_df %>%
    filter(comparison == "WT vs Tx") %>%
    filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5)) %>% pull(metabolite)
  lung_vor_taxa <-  lung_meta_sig_df %>%
    filter(comparison == "CA vs VOR") %>%
    filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5)) %>% pull(metabolite)
  lung_treatment_taxa <- lung_treatment_taxa[!lung_treatment_taxa %in% lung_vor_taxa]
  lung_venn_data <- list_to_data_frame(list("CA vs WT"=lung_infection_taxa, "WT vs Tx"=lung_treatment_taxa))
  lung_venn_plot <- ggplot(lung_venn_data, aes(A=`CA vs WT`,B= `WT vs Tx`))+
    geom_venn(fill_color=c("#7FC07E","#FFD17E"),set_name_size=7,text_size=6)+
    labs(title = "Lung metabolome")  +
      theme_void() + 
      coord_fixed(clip="off") +
      theme(
      plot.title = element_text(size = 25, hjust = 0.5))
  plasma_infection_taxa <-  plasma_meta_sig_df %>%
    filter(comparison == "CA vs WT") %>%
    filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5)) %>% pull(metabolite)
  plasma_treatment_taxa <-  plasma_meta_sig_df %>%
    filter(comparison == "WT vs Tx") %>%
    filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5)) %>% pull(metabolite)
  plasma_vor_taxa <-  plasma_meta_sig_df %>%
    filter(comparison == "CA vs VOR") %>%
    filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5)) %>% pull(metabolite)
  plasma_treatment_taxa <- plasma_treatment_taxa[!plasma_treatment_taxa %in% plasma_vor_taxa]
  plasma_venn_data <- list_to_data_frame(list("CA vs WT"=plasma_infection_taxa, "WT vs Tx"=plasma_treatment_taxa))
  plasma_venn_plot <- ggplot(plasma_venn_data, aes(A=`CA vs WT`,B= `WT vs Tx`))+
    geom_venn(fill_color=c("#7FC07E","#FFD17E"),set_name_size=7,text_size=6)+
    labs(title = "Plasma metabolome")  +
      theme_void() +
      coord_fixed(clip="off") +
      theme(
      plot.title = element_text(size = 25, hjust = 0.5))
  return(list(lung_venn_plot, plasma_venn_plot))

}


get_inter_group_beta_dist <-
  function(physeq, index, weighted=F, diversity_name = "") {
    meta_df <- data.frame(sample_data(physeq))
    if (!is.factor(meta_df$group)) {
      stop("the group column must be factor")
    }
    if (diversity_name == "") {
      diversity_name <- index
    }
    dist_mat <-
      phyloseq::distance(physeq, method = index, weighted = weighted)
    dist_df <- dist_mat %>% as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column("sample1") %>%
      pivot_longer(!sample1, names_to = "sample2") %>%
      mutate(G1 = meta_df[match(sample1,rownames(meta_df)),]$group, G2 = meta_df[match(sample2,rownames(meta_df)),]$group) %>%
      mutate(comparison = paste0(G1, " vs ", G2)) %>%
      filter(as.numeric(G1)!=as.numeric(G2)) %>%
      filter(startsWith(comparison,"Tx")) %>%
      mutate(index = diversity_name)
    return(dist_df)
  }


plot_microbiome_feast <- function(feast_df){
  feast_long_df <- feast_df %>%
    pivot_longer(!sample, names_to = "sample_group", values_to ="proportion")
  if(all(str_starts(feast_df$sample, "F3"))){
    plot_lab = "Gut microbiome"
    feast_long_df <- feast_long_df %>%
      mutate(group = gsub("^[^_]+_", "\\1",sample_group ))
  } else{
    plot_lab = "Lung microbiome"
    feast_long_df <- feast_long_df %>%
      mutate(group = gsub("^[^_]+_[^_]+_", "\\1",sample_group ))
  }
  plot_df <- feast_long_df %>%
    dplyr::select(sample, group, proportion) %>%
    group_by(sample, group) %>%
    summarise("Source proportion" = sum(proportion)) %>%
    ungroup() %>%
    mutate(group = factor(group, level=c(as.character(get_colormap_and_charmap()$group),"Unknown")))
  color <- get_colormap_and_charmap() %>%
    filter(group %in% plot_df$group) %>% pull(color)
  feast_plot <-  ggplot(data=plot_df, mapping = aes(x=group,y=`Source proportion`,fill=group))+
    geom_violin()+
    geom_boxplot(width=0.1)+
    stat_compare_means(method = "wilcox.test",comparisons=list(c("WT","VOR")),label = "p.format",
                         label.x = 1.5,
                         size = 4)+
    theme_pubr()+
    scale_fill_manual(values=c(color, "grey")) +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          text = element_text(16),
          legend.title = element_blank(),
          legend.key.size = unit(.5, 'cm'),
          legend.spacing.x = unit(.5, 'cm'),
          legend.text = element_text(size=12),
          axis.title = element_text(size = 16),
          axis.text.x=element_blank(),
          axis.ticks.length=unit(.1, "cm"),
          axis.text=element_text(size=16))+
    labs(title = plot_lab) + xlab("")
  return(feast_plot)
}

plot_metabolome_intergroup_dist <- function(dist_df, plot_lab=""){
  plot_df <- dist_df %>%
    filter(G2!="PptA")
  color <- get_colormap_and_charmap() %>%
    filter(group %in% plot_df$G2) %>% pull(color)
  feast_plot <-  ggplot(data=plot_df, mapping = aes(x=G2,y=1/value,fill=G2))+
    geom_violin()+
    geom_boxplot(width=0.1)+
    stat_compare_means(method = "wilcox.test",comparisons=list(c("WT","VOR")),label = "p.format",
                       label.x = 1.5,
                       size = 4)+
    theme_pubr()+
    scale_fill_manual(values=c(color, "grey")) +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          text = element_text(16),
          legend.title = element_blank(),
          legend.key.size = unit(.5, 'cm'),
          legend.spacing.x = unit(.5, 'cm'),
          legend.text = element_text(size=12),
          axis.title = element_text(size = 16),
          axis.text.x=element_blank(),
          axis.ticks.length=unit(.1, "cm"),
          axis.text=element_text(size=16))+
    labs(title = plot_lab) + xlab("") + ylab("1/Euclidean distance")
  return(feast_plot)
}

get_constant_taxa <- function(sig_tbl) {
  if("metabolite" %in% colnames(sig_tbl)){
    sig_tbl <- sig_tbl %>% dplyr::rename(OTU=metabolite)
  }
  filtered_taxa <- sig_tbl %>%
    filter(comparison %in% c("CA vs WT", "WT vs Tx")) %>%
    filter(pval < 0.05) %>%
    filter(abs(log2FC) > log2(1.5)) %>%
    dplyr::count(OTU, sort = T) %>%
    filter(n==2) %>%
    pull(OTU)
  filtered_tbl = tibble()
  tryCatch({
  filtered_tbl <- sig_tbl %>%
    filter(OTU  %in% filtered_taxa) %>%
    filter(comparison %in% c("CA vs WT", "WT vs Tx")) %>%
    dplyr::select(OTU, comparison, log2FC) %>%
    pivot_wider(names_from = comparison, values_from = log2FC) %>%
    filter(`CA vs WT`*`WT vs Tx`<0)},
  error = function(e){
    print(e)
    })
  return(filtered_tbl)
}

plot_constant_taxa <- function(physeq, constant_taxa_tbl,top_n=999,plot_lab="") {
  if(nrow(constant_taxa_tbl)<1){
    return(ggplot())
  }
  constant_taxa_tbl <- constant_taxa_tbl %>%
    arrange(desc(abs(`CA vs WT`)+abs(`WT vs Tx`))) %>%
    head(top_n)
  filtered_physeq <- prune_taxa(constant_taxa_tbl$OTU, physeq)
  otu_long_df <-  melt_phyloseq(filtered_physeq) %>%
    filter(group %in% c("CA", "WT", "Tx"))
  if("Genus" %in% colnames(otu_long_df)){
    otu_long_df <- otu_long_df %>%
      filter(Genus!="?") %>%
      mutate(OTU=paste0("*",Genus,"*"),
             OTU=ifelse(str_detect(OTU,"Family"),str_replace(OTU," Family\\*","* Family"),OTU),
             OTU=ifelse(str_detect(OTU,"Order"),str_replace(OTU," Order\\*","* Order"),OTU),
             OTU=ifelse(str_detect(OTU,"Class"),str_replace(OTU," Class\\*","* Class"),OTU),
             OTU=ifelse(str_detect(OTU,"Phylum"),str_replace(OTU," Phylum\\*","* Phylum"),OTU),
             OTU=ifelse(str_detect(OTU,"Domain"),str_replace(OTU," Domain\\*","* Domain"),OTU)
             )
  }
  otu_long_df <- otu_long_df %>%
    mutate(OTU = str_replace(OTU,"3-\\(propan-2-yl\\)-octahydropyrrolo\\[1,2-a\\]pyrazine-1,4-dione","CID:6992260"),
           OTU = str_replace(OTU,"3,4-Dihydroxybenzenesulfonic acid","CID:189003"),
           OTU = str_replace(OTU,"3-\\(2-Hydroxyethyl\\)indole","CID:10685"),
           OTU = str_replace(OTU," ","<br>"))
  means <- otu_long_df %>% group_by(OTU, group) %>%
    dplyr::summarize(med = median(Abundance))
  color <- get_colormap_and_charmap() %>%
    filter(group %in% means$group) %>% pull(color)
  ggplot(otu_long_df, aes(x=group, y=Abundance,fill=group)) +
    geom_boxplot() +
    facet_wrap(~OTU, scales = "free_y")+
    geom_line(data = means,
              mapping = aes(x=group, y=med,group=OTU), color="grey50") +
    theme_pubr()+
    theme(plot.title = element_text(size = 30, hjust = 0.5),
          text = element_text(20),
          strip.text = element_markdown(size = 24),
          strip.background = element_blank(),
          legend.title = element_blank(),
          legend.key.size = unit(1, 'cm'),
          legend.spacing.x = unit(.5, 'cm'),
          legend.text = element_text(size=30),
          axis.title = element_text(size = 24),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text=element_text(size=20))+
    labs(title = plot_lab) + xlab("") +
    scale_fill_manual(values=color)
}

get_constant_taxa_cor <- function(physeq_taxa_tbl){
  library(psych)
  filtered_physeq_taxa_tbl = tibble()
  tryCatch({
  filtered_physeq_taxa_tbl <- physeq_taxa_tbl %>%
    mutate(filtered_physeq = map2(physeq,taxa,~prune_taxa(.y$OTU,.x)))
  },error = function(e){
    print(e)
  })
  if(nrow(filtered_physeq_taxa_tbl)<1) return(tibble())
  lung_micro_long_tbl <- filtered_physeq_taxa_tbl %>%
    filter(omics=="lung_microbiome") %>% pull(filtered_physeq) %>% `[[`(1) %>%
    melt_phyloseq %>%
    dplyr::select(Genus, Abundance,sampleName) %>%
    dplyr::rename(OTU=Genus) %>%
    filter(OTU!="?") %>%
    mutate(cate="lung microbiome") %>%
    mutate(id=paste0(cate,'@',OTU))
  gut_micro_long_tbl <- filtered_physeq_taxa_tbl %>%
    filter(omics=="gut_microbiome") %>% pull(filtered_physeq) %>% `[[`(1) %>%
    transform_sample_counts(function(OTU) OTU / sum(OTU)) %>%
    melt_phyloseq %>%
    dplyr::select(Genus, Abundance,sampleName) %>%
    dplyr::rename(OTU=Genus) %>%
    filter(OTU!="?") %>%
    mutate(cate="gut microbiome") %>%
    mutate(id=paste0(cate,'@',OTU))
  lung_metabo_long_tbl <- filtered_physeq_taxa_tbl %>%
    filter(omics=="lung_metabolome") %>% pull(filtered_physeq) %>% `[[`(1) %>%
    melt_phyloseq %>%
    dplyr::select(metabolite, Abundance,sampleName) %>%
    dplyr::rename(OTU=metabolite) %>%
    mutate(cate="lung metabolome") %>%
    mutate(id=paste0(cate,'@',OTU))
  plasma_metabo_long_tbl <- filtered_physeq_taxa_tbl %>%
    filter(omics=="plasma_metabolome") %>% pull(filtered_physeq) %>% `[[`(1) %>%
    melt_phyloseq %>%
    dplyr::select(metabolite, Abundance,sampleName) %>%
    dplyr::rename(OTU=metabolite) %>%
    mutate(cate="plasma metabolome") %>%
    mutate(id=paste0(cate,'@',OTU))

  merged_otu_table <- bind_rows(lung_micro_long_tbl,gut_micro_long_tbl,lung_metabo_long_tbl,plasma_metabo_long_tbl) %>%
    dplyr::select(id,Abundance,sampleName) %>%
    pivot_wider(names_from = "id", values_from = "Abundance") %>%
    mutate(sampleName=paste0("s",sampleName)) %>%
    column_to_rownames(var="sampleName")
  merged_taxa_table <-  bind_rows(lung_micro_long_tbl,gut_micro_long_tbl,lung_metabo_long_tbl,plasma_metabo_long_tbl) %>%
    dplyr::select(id,cate,OTU) %>%
    unique() %>%
    column_to_rownames(var="id")
  merged_meta_table <-  filtered_physeq_taxa_tbl %>%
    filter(omics=="lung_microbiome") %>% pull(filtered_physeq) %>% `[[`(1) %>%
    sample_data() %>%
    as.data.frame()
  rownames(merged_meta_table) <- paste0("s",merged_meta_table$sampleName)
  META <- sample_data(merged_meta_table)
  TAX <- tax_table(as.matrix(merged_taxa_table))
  OTU <- otu_table(merged_otu_table,taxa_are_rows=F)
  physeq <- phyloseq(META, TAX, OTU)

  cor_obj <- corr.test(merged_otu_table,method="spearman", adjust="none", ci=F)
  cor_df <- cor_obj$r %>% as.data.frame %>% rownames_to_column(var="x") %>% pivot_longer(!x, names_to = "y", values_to = "r") %>%
    filter(x>y)
  pval_df <- cor_obj$p %>% as.data.frame %>% rownames_to_column(var="x") %>% pivot_longer(!x, names_to = "y", values_to = "pval") %>%
    filter(x>y)
  out_df <-cor_df %>% full_join(pval_df, by=c("x","y")) %>%  mutate(padj=p.adjust(pval,method="BH"))
  return(out_df)
}

plot_constant_taxa_cor_network <- function(cor_tbl, plot_title=""){
  library(igraph)
  library(ggnetwork)
  if(nrow(cor_tbl)<1) return(ggplot())
  filtered_cor_tbl <- cor_tbl %>%
    #filter(str_split_i(x,"@",1)!=str_split_i(y,"@",1)) %>%
    filter(pval <= 0.05)
  links <- filtered_cor_tbl %>%
    dplyr::select(x,y, r,pval,padj) %>%
    mutate(direction=ifelse(r<0, "negative","positive"))

  nodes <- unique(c(filtered_cor_tbl$x, filtered_cor_tbl$y)) %>%
    tibble() %>%
    dplyr::rename("name"=".") %>%
    mutate(temp=name) %>%
    separate(temp, c("omics","taxa"), sep = "@") %>% 
    arrange(name)
  network <- graph_from_data_frame(d=links, vertices=nodes, directed=F) %>% ggnetwork(layout=layout_in_circle(.))
  #network$p <- factor(network$p,level=c("p <= 0.05","p > 0.05"))
  #network$direction <- factor(network$direction,level=c("negative","positive"))
  ggplot(network, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(color=direction)) +
    geom_nodes(aes(color = omics),size=6) +
    #geom_nodelabel_repel(aes(color = omics, label = taxa),
    #                     fontface = "bold", box.padding = unit(1, "lines"))+
    theme_blank() +
    guides(colour = guide_legend(nrow = 1))+
    scale_color_brewer(palette = "Paired")+
    labs(title = plot_title) +
    theme(plot.title = element_text(size = 24,hjust = 0.5))
}
