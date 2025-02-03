get_fig3_comparison_df <- function() {
  comparison_df <-data.frame(
    G1 = c(
      "CA",
      "CA"
    ),
    G2 = c(
      "WT",
      "PptA"
    )
  )
  comparison_df$comparison <-
    paste(comparison_df$G1, "vs", comparison_df$G2)
  rownames(comparison_df) <- comparison_df$comparison
  return(comparison_df)
}


get_fig3_physeq <- function(physeq) {
  filtered_physeq <- physeq %>%
    subset_samples(sampleGroup %in% c(2,4,6)) %>%
    filter_taxa(function(x)
      sum(x > 0) > 0, TRUE)
  return(filtered_physeq)
}


get_pbs_ca_physeq <- function(physeq) {
  filtered_physeq <- physeq %>%
    subset_samples(sampleGroup %in% c(2)) %>%
    filter_taxa(function(x)
      sum(x > 0) > 0, TRUE)
  return(filtered_physeq)
}


get_af_wt_physeq <- function(physeq) {
  filtered_physeq <- physeq %>%
    subset_samples(sampleGroup %in% c(4)) %>%
    filter_taxa(function(x)
      sum(x > 0) > 0, TRUE)
  return(filtered_physeq)
}


get_af_ppta_physeq <- function(physeq) {
  filtered_physeq <- physeq %>%
    subset_samples(sampleGroup %in% c(6)) %>%
    filter_taxa(function(x)
      sum(x > 0) > 0, TRUE)
  return(filtered_physeq)
}


plot_fig3_alpha <- function(physeq,plot_title=""){
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
  ) + stat_compare_means(method = "wilcox.test",comparisons=list(c("CA","WT"),c("CA","PptA"),c("WT","PptA")),label = "p.format",
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
  ) + stat_compare_means(method = "wilcox.test",comparisons=list(c("CA","WT"),c("CA","PptA"),c("WT","PptA")),label = "p.format",
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

plot_fig3_beta <- function(physeq, distance="unifrac",weighted = T, plot_title="") {z
  meta_df <- data.frame(sample_data(physeq))
  color <- get_colormap_and_charmap() %>%
    filter(group %in% meta_df$group) %>% pull(color)
  ord_mat <-
    ordinate(physeq,
             method = "PCoA",
             distance = distance,
             weighted = weighted)
  dist_mat1 <-
    phyloseq::distance(physeq %>% subset_samples(group %in% c("CA","WT")), method = distance, weighted = weighted)
  meta_df1 <-  data.frame(sample_data(physeq %>% subset_samples(group %in% c("CA","WT"))))
  set.seed(42)
  adonis_result1 <-
    adonis2(dist_mat1 ~ group, data = meta_df1, permutations = 999)
  p1 <- adonis_result1$`Pr(>F)`[1]
  r1 <- round(adonis_result1$`R2`[1],3)

  dist_mat2 <-
    phyloseq::distance(physeq %>% subset_samples(group %in% c("CA","PptA")), method = distance, weighted = weighted)
  meta_df2<-  data.frame(sample_data(physeq %>% subset_samples(group %in% c("CA","PptA"))))
  set.seed(42)
  adonis_result2 <-
    adonis2(dist_mat2 ~ group, data = meta_df2, permutations = 999)
  p2 <- adonis_result2$`Pr(>F)`[1]
  r2 <- round(adonis_result2$`R2`[1],3)

  dist_mat3 <-
    phyloseq::distance(physeq %>% subset_samples(group %in% c("WT","PptA")), method = distance, weighted = weighted)
  meta_df3<-  data.frame(sample_data(physeq %>% subset_samples(group %in% c("WT","PptA"))))
  set.seed(42)
  adonis_result3 <-
    adonis2(dist_mat3 ~ group, data = meta_df3, permutations = 999)
  p3 <- adonis_result3$`Pr(>F)`[1]
  r3 <- round(adonis_result3$`R2`[1],3)

  beta_boxplot <- plot_ordination(physeq, ord_mat, color = "group",) +
    geom_point(size=5)+
    stat_ellipse(type = "norm",
                 linetype = 1,
                 lwd = 0.5) +
    stat_mean(aes(color=group),size=8,shape=13)+
    theme_pubr() +
    theme(
      plot.title = element_text(size = 30, hjust = 0.5),
      plot.subtitle = element_markdown(size = 24, hjust = 0.5),
      text = element_text(size = 24),
      legend.text = element_text(size = 24),
      legend.title = element_blank(),
      legend.position = "top",
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 1
      )
    ) +
    labs(subtitle = paste("CA vs WT: r<sup>2</sup> =", r1,", p = ", p1,
                          "<br />CA vs PptA: r<sup>2</sup> =", r2,", p = ", p2,
                          "<br />WT vs PptA: r<sup>2</sup> =", r3,", p = ", p3), title = plot_title) +  scale_color_manual(values=color)
  return(beta_boxplot)
}

plot_fig3_plsda <- function(physeq, plot_title="") {
  meta_df <- data.frame(sample_data(physeq))
  color <- get_colormap_and_charmap() %>%
    filter(group %in% meta_df$group) %>% pull(color)
  otu_df <- as.data.frame(otu_table(physeq)) %>% rotate_df()

  plsda = opls(otu_df,meta_df$group,predI=2,scaleC="standard",fig.pdfC="none")
  q2 =  plsda@summaryDF$`Q2(cum)`
  
  adonis_result <- adonis2(otu_df~meta_df$group, method="euclidean")

  dist_mat1 <-
    phyloseq::distance(physeq %>% subset_samples(group %in% c("CA","WT")), method = "euclidean")
  meta_df1 <-  data.frame(sample_data(physeq %>% subset_samples(group %in% c("CA","WT"))))
  set.seed(42)
  adonis_result1 <-
    adonis2(dist_mat1 ~ group, data = meta_df1, permutations = 999)
  p1 <- adonis_result1$`Pr(>F)`[1]
  r1 <- round(adonis_result1$`R2`[1],3)

  dist_mat2 <-
    phyloseq::distance(physeq %>% subset_samples(group %in% c("CA","PptA")), method = "euclidean")
  meta_df2<-  data.frame(sample_data(physeq %>% subset_samples(group %in% c("CA","PptA"))))
  set.seed(42)
  adonis_result2 <-
    adonis2(dist_mat2 ~ group, data = meta_df2, permutations = 999)
  p2 <- adonis_result2$`Pr(>F)`[1]
  r2 <- round(adonis_result2$`R2`[1],3)

  dist_mat3 <-
    phyloseq::distance(physeq %>% subset_samples(group %in% c("WT","PptA")), method = "euclidean", weighted = F)
  meta_df3<-  data.frame(sample_data(physeq %>% subset_samples(group %in% c("WT","PptA"))))
  set.seed(42)
  adonis_result3 <-
    adonis2(dist_mat3 ~ group, data = meta_df3, permutations = 999)
  p3 <- adonis_result3$`Pr(>F)`[1]
  r3 <- round(adonis_result3$`R2`[1],3)

  sample.score = plsda@scoreMN %>%
    as.data.frame() %>%
    mutate(group = meta_df$group)

  beta_boxplot = ggplot(sample.score, aes(p1, p2, color = group)) +
    geom_point(size=5) +
    stat_ellipse(type = "norm",
                 linetype = 1,
                 lwd = 0.5) +
    stat_mean(aes(color=group),size=8,shape=13)+
    annotate("text", x = 0.8*ceiling(min(sample.score$p1)), y = 1.1*floor(min(sample.score$p2)), label = paste("Q2(cum)=",q2), size = 7) +
    theme_pubr() +
    theme(
      plot.title = element_text(size = 30, hjust = 0.5),
      plot.subtitle = element_markdown(size = 24, hjust = 0.5),
      text = element_text(size = 24),
      legend.text = element_text(size = 30),
      legend.title = element_blank(),
      legend.position = "top",
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 1
      )
    ) +
    labs(subtitle = paste("CA vs WT: r<sup>2</sup> =", r1,", p = ", p1,
                          "<br />CA vs PptA: r<sup>2</sup> =", r2,", p = ", p2,
                          "<br />WT vs PptA: r<sup>2</sup> =", r3,", p = ", p3), title = plot_title) +
    labs(x = paste0("P1=",plsda@modelDF$R2X[1] *100,"%") ,y = paste0("P1=",plsda@modelDF$R2X[2] *100,"%"))+
    scale_color_manual(values=color)
  return(beta_boxplot)
}

plot_fig3_microbiome_heatmap <-
  function(lung_micro_physeq, lung_micro_deseq2_df, gut_micro_physeq, gut_micro_deseq2_df) {
    lung_micro_TMM_physeq<-
      transform_sample_counts(lung_micro_physeq, function(OTU) OTU / sum(OTU)) %>%
      subset_samples(sampleGroup %in% c(2,4,6)) %>%
      set_sample_order(c("group","sampleName"))
    gut_micro_TMM_physeq<-
      transform_sample_counts(gut_micro_physeq, function(OTU) OTU / sum(OTU)) %>%
      subset_samples(sampleGroup %in% c(2,4,6)) %>%
      set_sample_order(c("group","sampleName"))
    lung_micro_deseq2_filtered_taxa <- lung_micro_deseq2_df %>%
      filter(Genus != "?") %>%
      filter(pval < 0.05) %>%
      filter(abs(log2FC) > log2(1.5)) %>%
      pull(Genus)
    lung_micro_deseq2_filtered_df <- lung_micro_deseq2_df %>%
      filter(comparison %in% c("CA vs WT", "CA vs PptA" )) %>%
      filter(Genus %in% lung_micro_deseq2_filtered_taxa) %>%
      mutate(significance = pval < 0.05) %>% 
      mutate(FDR = padj < 0.1) %>%
      arrange(desc(comparison), Genus)
    gut_micro_deseq2_filtered_taxa <- gut_micro_deseq2_df %>%
      filter(Genus != "?") %>%
      filter(pval < 0.05) %>%
      filter(abs(log2FC) > log2(1.5)) %>%
      pull(Genus)
    gut_micro_deseq2_filtered_df <- gut_micro_deseq2_df %>%
      filter(comparison %in% c("CA vs WT", "CA vs PptA" )) %>%
      filter(Genus %in% gut_micro_deseq2_filtered_taxa) %>% 
      mutate(significance = pval < 0.05) %>%
      mutate(FDR = padj < 0.1) %>%
      arrange(desc(comparison), Genus)

    ## for lung microbiome
    lung_filtered_genus <- lung_micro_deseq2_filtered_df$Genus
    lung_micro_TMM_filtered_physeq <- lung_micro_TMM_physeq %>%
      prune_taxa(data.frame(tax_table(.))$Genus %in% lung_filtered_genus,.)
    lung_micro_TMM_otu_df <- melt_phyloseq(lung_micro_TMM_filtered_physeq) %>%
      dplyr::select(Genus,group,Abundance) %>%
      group_by(Genus) %>%
      mutate(Abundance = scale(Abundance)) %>%
      group_by(Genus,group) %>%
      summarise(Abundance = mean(Abundance)) %>%
      pivot_wider(names_from = group, values_from = Abundance)
    lung_mat <- lung_micro_TMM_otu_df %>%
      as.data.frame() %>%
      arrange(Genus) %>%
      column_to_rownames(var="Genus") %>%
      as.matrix()
    lung_micro_deseq2_filtered_df_subset1 <- lung_micro_deseq2_filtered_df[lung_micro_deseq2_filtered_df$comparison=="CA vs WT",] %>% filter(significance)
    lung_micro_deseq2_filtered_df_subset2 <-lung_micro_deseq2_filtered_df[lung_micro_deseq2_filtered_df$comparison=="CA vs PptA",] %>% filter(significance)
    lung_mat_info <- data.frame(Genus=rownames(lung_mat)) %>%
      mutate(comparison1 = Genus %in% lung_micro_deseq2_filtered_df_subset1$Genus) %>%
      mutate(comparison1_FDR = lung_micro_deseq2_filtered_df_subset1[match(Genus, lung_micro_deseq2_filtered_df_subset1$Genus),]$FDR) %>%
      mutate(comparison1_log2FC = lung_micro_deseq2_filtered_df_subset1[match(Genus, lung_micro_deseq2_filtered_df_subset1$Genus),]$log2FC) %>%
      mutate(comparison2 = Genus %in% lung_micro_deseq2_filtered_df_subset2$Genus) %>%
      mutate(comparison2_FDR = lung_micro_deseq2_filtered_df_subset2[match(Genus, lung_micro_deseq2_filtered_df_subset2$Genus),]$FDR) %>%
      mutate(comparison2_log2FC = lung_micro_deseq2_filtered_df_subset2[match(Genus, lung_micro_deseq2_filtered_df_subset2$Genus),]$log2FC) %>%
      arrange(Genus)
    lung_mat_info[is.na(lung_mat_info)] <- FALSE
    lung_mat_info <- lung_mat_info %>%
      mutate(comparison1_mark = ifelse(comparison1_log2FC >0,ifelse(comparison1,ifelse(comparison1_FDR,"⬆","⇧"),""),ifelse(comparison1,ifelse(comparison1_FDR,"⬇","⇩"),"")),
             comparison2_mark = ifelse(comparison2_log2FC >0,ifelse(comparison2,ifelse(comparison2_FDR,"⬆","⇧"),""),ifelse(comparison2,ifelse(comparison2_FDR,"⬇","⇩"),"")))

    lung_group <- sample_data(lung_micro_TMM_filtered_physeq) %>% as.data.frame() %>% pull(group) %>% rev %>% factor(levels = c("CA","WT","PptA"))
    lung_color <- get_colormap_and_charmap() %>%
      filter(group %in% levels(lung_group)) %>% pull(color)
    names(lung_color) <- levels(lung_group)
    lung_col_fun <-colorRamp2(c(-1, 0, 1), c("#0b5394","#ffffff","#d9b611"))

    lung_mat <- lung_mat %>% t
    lung_mat_info <- lung_mat_info %>%
      dplyr::select(comparison1_mark,comparison2_mark) %>%
      t

    lung_cell_fun = function(j, i, x, y, width, height, fill) {
      mark=""
      if(i==2){
        mark= lung_mat_info["comparison1_mark",j]
      }
      if(i==3){
        mark= lung_mat_info["comparison2_mark",j]
      }
        grid.text(mark,x = x, y = y,
          gp = gpar(col = "black", fill = NA, fontfamily = "Noto Sans JP"))
    }
    lung_heatmap_plot <-
      Heatmap(
        lung_mat,
        name = "lung",
        col = lung_col_fun,
        cell_fun = lung_cell_fun,
        border_gp = gpar(col = "black", lty = 2),
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_columns = F,
        column_title= "Lung",
        column_title_gp = gpar(fontsize=14),
        column_names_side = "top",
        column_names_gp = gpar(fontsize=12,fontface="italic"),
        column_names_rot = 60,
        cluster_rows = F,
        row_title_gp = gpar(fontsize=12),
        row_names_side = "left",
        row_names_gp = gpar(fontsize=14),
        width = ncol(lung_mat)*unit(5, "mm"),
        height = nrow(lung_mat)*unit(5, "mm"),
        heatmap_legend_param = list(title = "Z-score"),
        show_heatmap_legend = T
      )
    ## for gut microbiome
    gut_filtered_genus <- gut_micro_deseq2_filtered_df$Genus
    gut_micro_TMM_filtered_physeq <- gut_micro_TMM_physeq %>%
      prune_taxa(data.frame(tax_table(.))$Genus %in% gut_filtered_genus,.)
    gut_micro_TMM_otu_df <- melt_phyloseq(gut_micro_TMM_filtered_physeq) %>%
      dplyr::select(Genus,group,Abundance) %>%
      group_by(Genus) %>%
      mutate(Abundance = scale(Abundance)) %>%
      group_by(Genus,group) %>%
      summarise(Abundance = mean(Abundance)) %>%
      pivot_wider(names_from = group, values_from = Abundance)
    gut_mat <- gut_micro_TMM_otu_df %>% 
      as.data.frame() %>%
      arrange(Genus) %>%
      column_to_rownames(var="Genus") %>%
      as.matrix()
    gut_micro_deseq2_filtered_df_subset1 <- gut_micro_deseq2_filtered_df[gut_micro_deseq2_filtered_df$comparison=="CA vs WT",] %>% filter(significance)
    gut_micro_deseq2_filtered_df_subset2 <-gut_micro_deseq2_filtered_df[gut_micro_deseq2_filtered_df$comparison=="CA vs PptA",] %>% filter(significance)
    gut_mat_info <- data.frame(Genus=rownames(gut_mat)) %>%
      mutate(comparison1 = Genus %in% gut_micro_deseq2_filtered_df_subset1$Genus) %>%
      mutate(comparison1_FDR = gut_micro_deseq2_filtered_df_subset1[match(Genus, gut_micro_deseq2_filtered_df_subset1$Genus),]$FDR) %>%
      mutate(comparison1_log2FC = gut_micro_deseq2_filtered_df_subset1[match(Genus, gut_micro_deseq2_filtered_df_subset1$Genus),]$log2FC) %>%
      mutate(comparison2 = Genus %in% gut_micro_deseq2_filtered_df_subset2$Genus) %>%
      mutate(comparison2_FDR = gut_micro_deseq2_filtered_df_subset2[match(Genus, gut_micro_deseq2_filtered_df_subset2$Genus),]$FDR) %>%
      mutate(comparison2_log2FC = gut_micro_deseq2_filtered_df_subset2[match(Genus, gut_micro_deseq2_filtered_df_subset2$Genus),]$log2FC) %>%
      arrange(Genus)
    gut_mat_info[is.na(gut_mat_info)] <- FALSE
    gut_mat_info <- gut_mat_info %>%
      mutate(comparison1_mark = ifelse(comparison1_log2FC >0,ifelse(comparison1,ifelse(comparison1_FDR,"⬆","⇧"),""),ifelse(comparison1,ifelse(comparison1_FDR,"⬇","⇩"),"")),
             comparison2_mark = ifelse(comparison2_log2FC >0,ifelse(comparison2,ifelse(comparison2_FDR,"⬆","⇧"),""),ifelse(comparison2,ifelse(comparison2_FDR,"⬇","⇩"),"")))

    gut_group <- sample_data(gut_micro_TMM_filtered_physeq) %>% as.data.frame() %>% pull(group) %>% rev %>% factor(levels = c("CA","WT","PptA"))
    gut_color <- get_colormap_and_charmap() %>%
      filter(group %in% levels(gut_group)) %>% pull(color)
    names(gut_color) <- levels(gut_group)
    gut_col_fun <-colorRamp2(c(-1, 0, 1), c("#0b5394","#ffffff","#d9b611"))
    
    gut_mat <- gut_mat %>% t
    gut_mat_info <- gut_mat_info %>%
      dplyr::select(comparison1_mark,comparison2_mark) %>%
      t

    gut_cell_fun = function(j, i, x, y, width, height, fill) {
      mark=""
      if(i==2){
        mark= gut_mat_info["comparison1_mark",j]
      }
      if(i==3){
        mark= gut_mat_info["comparison2_mark",j]
      }
        grid.text(mark,x = x, y = y,
          gp = gpar(col = "black", fill = NA, fontfamily = "Noto Sans JP"))
    }
    gut_heatmap_plot <-
      Heatmap(
        gut_mat,
        name = "gut",
        col = gut_col_fun,
        cell_fun = gut_cell_fun,
        border_gp = gpar(col = "black", lty = 2),
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_columns = F,
        column_title= "Gut",
        column_title_gp = gpar(fontsize=14),
        column_names_side = "top",
        column_names_gp = gpar(fontsize=12,fontface="italic"),
        column_names_rot = 60,
        cluster_rows = F,
        row_title_gp = gpar(fontsize=12),
        row_names_side = "left",
        row_names_gp = gpar(fontsize=14),
        width = ncol(gut_mat)*unit(5, "mm"),
        height = nrow(gut_mat)*unit(5, "mm"),
        heatmap_legend_param = list(title = "Z-score"),
        show_heatmap_legend = F
      )
    ht_list <- lung_heatmap_plot + gut_heatmap_plot
    return(ht_list)
  }

plot_fig3_metabolome_heatmap <-
  function(lung_meta_physeq, lung_meta_sig_df, plasma_meta_physeq, plasma_meta_sig_df, n=10) {
    lung_meta_physeq<-lung_meta_physeq %>%
      subset_samples(sampleGroup %in% c(2,4,6)) %>%
      set_sample_order(c("group","sampleName"))
    plasma_meta_physeq<-plasma_meta_physeq %>%
      subset_samples(sampleGroup %in% c(2,4,6)) %>%
      set_sample_order(c("group","sampleName"))
    lung_meta_deseq2_filtered_taxa <- lung_meta_sig_df %>%
      filter(pval < 0.05) %>%
      filter(abs(log2FC) > log2(1.5)) %>%
      filter(comparison %in% c("CA vs WT", "CA vs PptA" )) %>%
      group_by(comparison) %>%
      arrange(pval) %>%
      slice_head(n=n) %>%
      pull(metabolite)
    lung_meta_deseq2_filtered_df <- lung_meta_sig_df %>%
      filter(comparison %in% c("CA vs WT", "CA vs PptA" )) %>%
      filter(metabolite %in% lung_meta_deseq2_filtered_taxa) %>% 
      mutate(significance = pval < 0.05) %>% 
      mutate(FDR = padj < 0.1) %>%
      arrange(desc(comparison), metabolite)

    plasma_meta_deseq2_filtered_taxa <- plasma_meta_sig_df %>%
      filter(pval < 0.05) %>%
      filter(abs(log2FC) > log2(1.5)) %>%
      filter(comparison %in% c("CA vs WT", "CA vs PptA" )) %>%
      group_by(comparison) %>%
      arrange(pval) %>%
      slice_head(n=n) %>%
      pull(metabolite)
    plasma_meta_deseq2_filtered_df <- plasma_meta_sig_df %>%
      filter(comparison %in% c("CA vs WT", "CA vs PptA" )) %>%
      filter(metabolite %in% plasma_meta_deseq2_filtered_taxa) %>% 
      mutate(significance = pval < 0.05) %>% 
      mutate(FDR = padj < 0.1) %>%
      arrange(desc(comparison), metabolite)

    ## for lung metabolome
    lung_filtered_metabolite <- lung_meta_deseq2_filtered_df$metabolite
    lung_meta_log2_filtered_physeq <- lung_meta_physeq %>%
      prune_taxa(lung_filtered_metabolite,.)
    lung_meta_log2_otu_df <- melt_phyloseq(lung_meta_log2_filtered_physeq) %>%
      dplyr::select(metabolite,group,Abundance) %>%
      group_by(metabolite) %>%
      mutate(Abundance = scale(Abundance)) %>%
      group_by(metabolite,group) %>%
      summarise(Abundance = mean(Abundance)) %>%
      pivot_wider(names_from = group, values_from = Abundance)

    lung_mat <- lung_meta_log2_otu_df %>%
      as.data.frame() %>%
      arrange(metabolite) %>%
      column_to_rownames(var="metabolite") %>%
      as.matrix()
    lung_meta_deseq2_filtered_df_subset1 <- lung_meta_deseq2_filtered_df[lung_meta_deseq2_filtered_df$comparison=="CA vs WT",] %>% filter(significance)
    lung_meta_deseq2_filtered_df_subset2 <-lung_meta_deseq2_filtered_df[lung_meta_deseq2_filtered_df$comparison=="CA vs PptA",] %>% filter(significance)
    lung_mat_info <- data.frame(metabolite=rownames(lung_mat)) %>%
      mutate(comparison1 = metabolite %in% lung_meta_deseq2_filtered_df_subset1$metabolite) %>%
      mutate(comparison1_FDR = lung_meta_deseq2_filtered_df_subset1[match(metabolite, lung_meta_deseq2_filtered_df_subset1$metabolite),]$FDR) %>%
      mutate(comparison1_log2FC = lung_meta_deseq2_filtered_df_subset1[match(metabolite, lung_meta_deseq2_filtered_df_subset1$metabolite),]$log2FC) %>%
      mutate(comparison2 = metabolite %in% lung_meta_deseq2_filtered_df_subset2$metabolite) %>%
      mutate(comparison2_FDR = lung_meta_deseq2_filtered_df_subset2[match(metabolite, lung_meta_deseq2_filtered_df_subset2$metabolite),]$FDR) %>%
      mutate(comparison2_log2FC = lung_meta_deseq2_filtered_df_subset2[match(metabolite, lung_meta_deseq2_filtered_df_subset2$metabolite),]$log2FC) %>%
      arrange(metabolite)
    lung_mat_info[is.na(lung_mat_info)] <- FALSE
        lung_mat_info <- lung_mat_info %>%
      mutate(comparison1_mark = ifelse(comparison1_log2FC >0,ifelse(comparison1,ifelse(comparison1_FDR,"⬆","⇧"),""),ifelse(comparison1,ifelse(comparison1_FDR,"⬇","⇩"),"")),
             comparison2_mark = ifelse(comparison2_log2FC >0,ifelse(comparison2,ifelse(comparison2_FDR,"⬆","⇧"),""),ifelse(comparison2,ifelse(comparison2_FDR,"⬇","⇩"),"")))

    lung_group <- sample_data(lung_meta_log2_filtered_physeq) %>% as.data.frame() %>% pull(group) %>% rev %>% factor(levels = c("CA","WT","PptA"))
    lung_color <- get_colormap_and_charmap() %>%
      filter(group %in% levels(lung_group)) %>% pull(color)
    names(lung_color) <- levels(lung_group)
    lung_col_ha <- HeatmapAnnotation(Group =anno_block(labels =levels(lung_group),gp=gpar(fill=lung_color),labels_gp= gpar(fontsize = 28)))
    lung_col_fun <-colorRamp2(c(-1, 0, 1), c( "#0b5394","#ffffff","#d9b611"))
    rownames(lung_mat) <- rownames(lung_mat) %>%
      str_replace("1,4a-Dimethyl-6-methylene-5-\\[2-\\(2-oxo-2,5-dihydro-3-furanyl\\)ethyl\\]decahydro-1-naphthalenecarboxylic acid", "CID:14806212") %>%
      str_replace("5,6-dimethyl-4-oxo-4H-pyran-2-carboxylic acid", "CID:2786724") %>%
      str_replace("\\(8aR,12S,12aR\\)-12-Hydroxy-4-methyl-4,5,6,7,8,8a,12,12a-octahydro-2H-3-benzoxecine-2,9\\(1H\\)-dione", "CID:637324") %>%
      str_replace("3-\\(propan-2-yl\\)-octahydropyrrolo\\[1,2-a\\]pyrazine-1,4-dione", "CID:98951") %>%
      str_replace("2-\\[\\(1S\\)-1-Hydroxyethyl\\]-4\\(1H\\)-quinazolinone", "CID:11116825") %>%
      str_replace("10-\\(Hydroxymethyl\\)-4-methyl-8,11-dioxa-2,6-diazatricyclo\\[7.2.1.02,7\\]dodeca-3,6-dien-5-one", "CID:457982")
    
    lung_mat <- lung_mat %>% t
    lung_mat_info <- lung_mat_info %>%
      dplyr::select(comparison1_mark,comparison2_mark) %>%
      t

    lung_cell_fun <- function(j, i, x, y, width, height, fill) {
      mark=""
      if(i==2){
        mark= lung_mat_info["comparison1_mark",j]
      }
      if(i==3){
        mark= lung_mat_info["comparison2_mark",j]
      }
        grid.text(mark,x = x, y = y,
          gp = gpar(col = "black", fill = NA, fontfamily = "Noto Sans JP"))
    }
    lung_heatmap_plot <-
      Heatmap(
        lung_mat,
        name = "lung",
        col = lung_col_fun,
        cell_fun = lung_cell_fun,
        border_gp = gpar(col = "black", lty = 2),
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_columns = F,
        column_title= "Lung",
        column_title_gp = gpar(fontsize=14),
        column_names_side = "top",
        column_names_gp = gpar(fontsize=12),
        column_names_rot = 60,
        cluster_rows = F,
        row_title_gp = gpar(fontsize=14),
        row_names_side = "left",
        row_names_gp = gpar(fontsize=12),
        width = ncol(lung_mat)*unit(5, "mm"),
        height = nrow(lung_mat)*unit(5, "mm"),
        heatmap_legend_param = list(title = "Normalized abundance"),
        show_heatmap_legend = T
      )

    ## for plasma metabolome
    plasma_filtered_metabolite <- plasma_meta_deseq2_filtered_df$metabolite
    plasma_meta_log2_filtered_physeq <- plasma_meta_physeq %>%
      prune_taxa(plasma_filtered_metabolite,.)
    plasma_meta_log2_otu_df <- melt_phyloseq(plasma_meta_log2_filtered_physeq) %>%
      dplyr::select(metabolite,group,Abundance) %>%
      group_by(metabolite) %>%
      mutate(Abundance = scale(Abundance)) %>%
      group_by(metabolite,group) %>%
      summarise(Abundance = mean(Abundance)) %>%
      pivot_wider(names_from = group, values_from = Abundance)

    plasma_mat <- plasma_meta_log2_otu_df %>%
      as.data.frame() %>%
      arrange(metabolite) %>%
      column_to_rownames(var="metabolite") %>%
      as.matrix()
    plasma_meta_deseq2_filtered_df_subset1 <- plasma_meta_deseq2_filtered_df[plasma_meta_deseq2_filtered_df$comparison=="CA vs WT",] %>% filter(significance)
    plasma_meta_deseq2_filtered_df_subset2 <-plasma_meta_deseq2_filtered_df[plasma_meta_deseq2_filtered_df$comparison=="CA vs PptA",] %>% filter(significance)
    plasma_mat_info <- data.frame(metabolite=rownames(plasma_mat)) %>%
      mutate(comparison1 = metabolite %in% plasma_meta_deseq2_filtered_df_subset1$metabolite) %>%
      mutate(comparison1_FDR = plasma_meta_deseq2_filtered_df_subset1[match(metabolite, plasma_meta_deseq2_filtered_df_subset1$metabolite),]$FDR) %>%
      mutate(comparison1_log2FC = plasma_meta_deseq2_filtered_df_subset1[match(metabolite, plasma_meta_deseq2_filtered_df_subset1$metabolite),]$log2FC) %>%
      mutate(comparison2 = metabolite %in% plasma_meta_deseq2_filtered_df_subset2$metabolite) %>%
      mutate(comparison2_FDR = plasma_meta_deseq2_filtered_df_subset2[match(metabolite, plasma_meta_deseq2_filtered_df_subset2$metabolite),]$FDR) %>%
      mutate(comparison2_log2FC = plasma_meta_deseq2_filtered_df_subset2[match(metabolite, plasma_meta_deseq2_filtered_df_subset2$metabolite),]$log2FC) %>%
      arrange(metabolite)
    plasma_mat_info[is.na(plasma_mat_info)] <- FALSE
        plasma_mat_info <- plasma_mat_info %>%
      mutate(comparison1_mark = ifelse(comparison1_log2FC >0,ifelse(comparison1,ifelse(comparison1_FDR,"⬆","⇧"),""),ifelse(comparison1,ifelse(comparison1_FDR,"⬇","⇩"),"")),
             comparison2_mark = ifelse(comparison2_log2FC >0,ifelse(comparison2,ifelse(comparison2_FDR,"⬆","⇧"),""),ifelse(comparison2,ifelse(comparison2_FDR,"⬇","⇩"),"")))

    plasma_group <- sample_data(plasma_meta_log2_filtered_physeq) %>% as.data.frame() %>% pull(group) %>% rev %>% factor(levels = c("CA","WT","PptA"))
    plasma_color <- get_colormap_and_charmap() %>%
      filter(group %in% levels(plasma_group)) %>% pull(color)
    names(plasma_color) <- levels(plasma_group)
    plasma_col_ha <- HeatmapAnnotation(Group =anno_block(labels =levels(plasma_group),gp=gpar(fill=plasma_color),labels_gp= gpar(fontsize = 28)))
    plasma_col_fun <-colorRamp2(c(-1, 0, 1), c( "#0b5394","#ffffff","#d9b611"))
    rownames(plasma_mat) <- rownames(plasma_mat) %>%
      str_replace("1,4a-Dimethyl-6-methylene-5-\\[2-\\(2-oxo-2,5-dihydro-3-furanyl\\)ethyl\\]decahydro-1-naphthalenecarboxylic acid", "CID:14806212") %>%
      str_replace("5,6-dimethyl-4-oxo-4H-pyran-2-carboxylic acid", "CID:2786724") %>%
      str_replace("\\(8aR,12S,12aR\\)-12-Hydroxy-4-methyl-4,5,6,7,8,8a,12,12a-octahydro-2H-3-benzoxecine-2,9\\(1H\\)-dione", "CID:637324") %>%
      str_replace("3-\\(propan-2-yl\\)-octahydropyrrolo\\[1,2-a\\]pyrazine-1,4-dione", "CID:98951") %>%
      str_replace("2-\\[\\(1S\\)-1-Hydroxyethyl\\]-4\\(1H\\)-quinazolinone", "CID:11116825") %>%
      str_replace("3-\\(1-hydroxyethyl\\)-2,3,6,7,8,8a-hexahydropyrrolo\\[1,2-a\\]pyrazine-1,4-dione", "CID:73146465") %>%
      str_replace("α-aminobutyric acid\\/γ-Aminobutyric acid", "α\\/γ-Aminobutyric acid")

    
    plasma_mat <- plasma_mat %>% t
    plasma_mat_info <- plasma_mat_info %>%
      dplyr::select(comparison1_mark,comparison2_mark) %>%
      t

    plasma_cell_fun <- function(j, i, x, y, width, height, fill) {
      mark=""
      if(i==2){
        mark= plasma_mat_info["comparison1_mark",j]
      }
      if(i==3){
        mark= plasma_mat_info["comparison2_mark",j]
      }
        grid.text(mark,x = x, y = y,
          gp = gpar(col = "black", fill = NA, fontfamily = "Noto Sans JP"))
    }
    plasma_heatmap_plot <-
      Heatmap(
        plasma_mat,
        name = "plasma",
        col = plasma_col_fun,
        cell_fun = plasma_cell_fun,
        border_gp = gpar(col = "black", lty = 2),
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_columns = F,
        column_title= "Plasma",
        column_title_gp = gpar(fontsize=14),
        column_names_side = "top",
        column_names_gp = gpar(fontsize=12),
        column_names_rot = 60,
        cluster_rows = F,
        row_title_gp = gpar(fontsize=14),
        row_names_side = "left",
        row_names_gp = gpar(fontsize=12),
        width = ncol(plasma_mat)*unit(5, "mm"),
        height = nrow(plasma_mat)*unit(5, "mm"),
        heatmap_legend_param = list(title = "Normalized abundance"),
        show_heatmap_legend = F
      )
    ht_list <- lung_heatmap_plot + plasma_heatmap_plot
    return(ht_list)
  }


get_correlation_between_omics <- function(omics_tbl) {
  lung_microbiome_otu_df <- omics_tbl %>% filter(name=="lung_microbiome") %>% pull(data) %>% `[[`(1) %>% otu_table() %>% as.data.frame() %>% t
  gut_microbiome_otu_df <- omics_tbl %>% filter(name=="gut_microbiome") %>%pull(data) %>% `[[`(1) %>% otu_table() %>% as.data.frame() %>% t
  lung_metabolome_otu_df <- omics_tbl %>% filter(name=="lung_metabolome") %>%pull(data) %>% `[[`(1) %>% otu_table() %>% as.data.frame() %>% t
  plasma_metabolome_otu_df <- omics_tbl %>% filter(name=="plasma_metabolome") %>%pull(data) %>% `[[`(1) %>% otu_table() %>% as.data.frame() %>% t

  out_tbl <- tibble(X=c("lung_microbiome_otu_df","gut_microbiome_otu_df","lung_metabolome_otu_df","plasma_metabolome_otu_df"),
                    Y=c("lung_microbiome_otu_df","gut_microbiome_otu_df","lung_metabolome_otu_df","plasma_metabolome_otu_df")) %>%
    tidyr::expand(X,Y) %>%
    filter(X>Y) %>%
    arrange(X,Y)

  out_tbl <- out_tbl %>% rowwise() %>%
    dplyr::mutate(result = list(corr.test(get(X), get(Y),  method = "spearman",adjust="none"))) %>% ungroup()
  return(out_tbl)
}


compare_correlation_between_omics <- function(ctl_cor_tbl,exp_cor_tbl) {
  exp_cor_r_tbl <- exp_cor_tbl %>%
    mutate(result = map(result, ~.x$r%>% as.data.frame %>%  rownames_to_column(var="x") %>% pivot_longer(!x, names_to = "y", values_to ="r"))) %>%
    unnest(cols = c(result)) %>%
    mutate(temp = paste0(X,Y,x,y,sep="@"))
  ctl_cor_r_tbl <- ctl_cor_tbl %>%
    mutate(result = map(result, ~.x$r%>% as.data.frame %>%  rownames_to_column(var="x") %>% pivot_longer(!x, names_to = "y", values_to ="r"))) %>%
    unnest(cols = c(result)) %>%
    mutate(temp = paste0(X,Y,x,y,sep="@")) %>%
    filter(temp %in% exp_cor_r_tbl$temp)
  exp_cor_r_tbl <- exp_cor_r_tbl %>%
    filter(temp %in% ctl_cor_r_tbl$temp)
  merged_tbl <- bind_rows(exp_cor_r_tbl, ctl_cor_r_tbl) %>%
    mutate(r = abs(r)) %>%
    arrange(group, temp) %>%
    dplyr::select(-temp)
  G1= unique(ctl_cor_tbl$group)
  G2=unique(exp_cor_tbl$group)
  log2FC_df <- merged_tbl %>%
    mutate(r = abs(r)) %>%
    group_by(X, Y, group) %>%
    dplyr::summarize(median_r = sum(r)) %>%
    pivot_wider(names_from = group, values_from = median_r) %>%
    mutate(log2FC = log2(get(G2) / get(G1))) %>%
    ungroup()
  res <- merged_tbl %>%
    group_by(X, Y) %>%
    group_modify(~tidy(ks.test(r ~ group, data = .x)))
  sigtab <- res %>%
    dplyr::rename(pval=p.value) %>%
    mutate(G1=G1, G2=G2) %>%
    mutate(comparison = paste0(G1, " vs ", G2)) %>%
    left_join(log2FC_df, by=c("X","Y"))
  return(sigtab)
}


plot_correlation_between_omics_chord <- function(cor_tbl) {
  library(ggnetwork)
  library(sna)
  links <- cor_tbl %>%
    dplyr::select(X, Y, log2FC, pval) %>%
    mutate(X=str_to_title(X), Y=str_to_title(Y)) %>%
    mutate(X=str_remove_all(X, "_otu_df"), Y=str_remove_all(Y, "_otu_df")) %>%
    mutate(X=str_replace_all(X, "_", " "), Y=str_replace_all(Y, "_", " "))  %>%
    mutate(X=str_replace_all(X, "microbiome", "microbiome"), Y=str_replace_all(Y, "microbiome", "microbiome")) %>%
    mutate(`abs (log2FC)`=abs(log2FC)) %>%
    mutate(direction = factor(ifelse(log2FC > 0, "increase", "decrease"), levels=c("increase", "decrease"))) %>%
    mutate(significance = factor(ifelse(pval <= 0.05, "p ≤ 0.05", "p > 0.05"), levels=c("p ≤ 0.05", "p > 0.05")))
  network <- links %>% ggnetwork(layout = "circle")
  network_plot <- ggplot(network, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(linetype=significance, linewidth=`abs (log2FC)`, color=direction)) +
    geom_nodes(color = "grey", size=15) +
    geom_nodelabel_repel(aes(label = vertex.names),
                   fontface = "bold", size=5) +
    theme_blank()
  return(network_plot)
}


get_modima_result <- function(lung_physeq, gut_physeq, lung_metabolome_physeq, plasma_metabolome_physeq){
  lung_physeq <- lung_physeq  %>% set_sample_order("sampleName")
  gut_physeq <- gut_physeq  %>% set_sample_order("sampleName")
  lung_metabolome_physeq <- lung_metabolome_physeq %>% set_sample_order("sampleName")
  plasma_metabolome_physeq <- plasma_metabolome_physeq %>% set_sample_order("sampleName")
  stopifnot(get_variable(lung_physeq,"sampleName") == get_variable(gut_physeq,"sampleName"))
  stopifnot(get_variable(lung_physeq,"sampleName") == get_variable(lung_metabolome_physeq,"sampleName"))
  stopifnot(get_variable(lung_physeq,"sampleName") == get_variable(plasma_metabolome_physeq,"sampleName"))

  lung_micro_dist <- phyloseq::distance(lung_physeq,"wunifrac")
  gut_micro_dist <- phyloseq::distance(gut_physeq,"wunifrac")
  lung_meta_dist <- phyloseq::distance(lung_metabolome_physeq,"euclidean")
  plasma_meta_dist <- phyloseq::distance(plasma_metabolome_physeq,"euclidean")

  out_tbl <- tibble(X=c("lung_micro_dist","gut_micro_dist","lung_meta_dist","plasma_meta_dist"),
                    M=c("lung_micro_dist","gut_micro_dist","lung_meta_dist","plasma_meta_dist"),
                    Y=c("lung_micro_dist","gut_micro_dist","lung_meta_dist","plasma_meta_dist")) %>%
    tidyr::expand(X,M,Y) %>%
    filter(!X==Y) %>%
    filter(!X==M) %>%
    filter(!M==Y) %>%
    arrange(X,M,Y)
  set.seed(42)
  out_tbl <- out_tbl %>% rowwise() %>%
    dplyr::mutate(result = list(modima(get(X), get(M), get(Y), nrep=9999))) %>%
    dplyr::mutate(X_M = result$estimates["exposure–mediator bcdcor"]) %>%
    dplyr::mutate(M_Y = result$estimates["mediator–response bcdcor"]) %>%
    dplyr::mutate(TE = result$estimates["exposure–response bcdcor"]) %>%
    dplyr::mutate(ACME = result$statistic) %>% #statistic = exposure–mediator bcdcor * mediator–response–exposure pdcor
    dplyr::mutate(M_Prop = round(ACME/TE,4)) %>%
    dplyr::mutate(pvalue = result$p.value)
  return(out_tbl)
}


plot__modima_result <- function(pbs_modiam_tbl,wt_modima_tbl, ppta_modima_tbl) {
  pbs_modiam_tbl <- pbs_modiam_tbl %>%
    mutate(group="CA")
  wt_modima_tbl <- wt_modima_tbl %>%
    mutate(group="WT")
  ppta_modima_tbl <- ppta_modima_tbl %>%
    mutate(group="PptA")
  merge_tbl <- bind_rows(pbs_modiam_tbl, wt_modima_tbl, ppta_modima_tbl) %>%
    mutate(path = paste(str_to_title(X),str_to_title(M),str_to_title(Y), sep = " -> ")) %>%
    mutate(path= str_remove_all(path, "_dist")) %>%
    mutate(path = str_replace_all(path, "_", " ")) %>%
    mutate(path = str_replace_all(path, "micro", "microbiome")) %>%
    mutate(path = str_replace_all(path, "meta", "metabolome")) %>%
    filter(path %in%  c("Gut microbiome -> Plasma metabolome -> Lung microbiome",
                        "Gut microbiome -> Plasma metabolome -> Lung metabolome",
                        "Lung microbiome -> Plasma metabolome -> Gut microbiome",
                        "Lung metabolome -> Plasma metabolome -> Gut microbiome")) %>%
    mutate(path = factor(path, levels= c(
      "Lung metabolome -> Plasma metabolome -> Gut microbiome",
      "Lung microbiome -> Plasma metabolome -> Gut microbiome",
      "Gut microbiome -> Plasma metabolome -> Lung microbiome",
      "Gut microbiome -> Plasma metabolome -> Lung metabolome"))) %>%
    mutate(group = factor(group, levels=c("CA","WT","PptA")))
  color <- get_colormap_and_charmap() %>%
    filter(group %in% merge_tbl$group) %>% pull(color)
  bar_plot <- ggplot(merge_tbl, aes(fill=group, y=M_Prop, x=path)) +
    geom_col(position = "dodge") +
    geom_text(
      aes(label = paste0("p=",round(pvalue,3))),
      colour = "black", size = 5,
      vjust = -0.2, position = position_dodge(0.9)) +
    ylab("Mediation proportion") +
    xlab("") +
    theme_pubr() +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5),
      text = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      legend.position = "top"
    )+
    scale_fill_manual(values=color) +
    coord_flip()
}


get_hima_result <- function(omics_tbl,sig_object_tbl=NULL){
  library(future)
  library(furrr)
  library(hdmed)

  lung_microbiome_otu_df <- omics_tbl %>% filter(name=="lung_microbiome") %>% pull(data) %>% `[[`(1) %>% otu_table() %>% as.data.frame()
  gut_microbiome_otu_df <- omics_tbl %>% filter(name=="gut_microbiome") %>%pull(data) %>% `[[`(1) %>% otu_table() %>% as.data.frame()
  lung_metabolome_otu_df <- omics_tbl %>% filter(name=="lung_metabolome") %>%pull(data) %>% `[[`(1) %>% otu_table() %>% as.data.frame()
  plasma_metabolome_otu_df <- omics_tbl %>% filter(name=="plasma_metabolome") %>%pull(data) %>% `[[`(1) %>% otu_table() %>% as.data.frame()

  if(!is.null(sig_object_tbl)){
    lung_microbiome_sig_df <- sig_object_tbl %>% filter(name=="lung_microbiome") %>% pull(data) %>% `[[`(1) %>%
      filter(comparison %in% c("CA vs WT","CA vs PptA") ) %>%
      filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5))
    gut_microbiome_sig_df <- sig_object_tbl %>% filter(name=="gut_microbiome") %>%pull(data) %>% `[[`(1) %>%
      filter(comparison %in% c("CA vs WT","CA vs PptA") ) %>%
      filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5))
    lung_metabolome_sig_df <- sig_object_tbl %>% filter(name=="lung_metabolome") %>%pull(data) %>% `[[`(1) %>% dplyr::rename(OTU=metabolite) %>%
      filter(comparison %in% c("CA vs WT","CA vs PptA") ) %>%
       filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5))
    plasma_metabolome_sig_df <- sig_object_tbl %>% filter(name=="plasma_metabolome") %>%pull(data) %>% `[[`(1) %>%  dplyr::rename(OTU=metabolite) %>%
      filter(comparison %in% c("CA vs WT","CA vs PptA") ) %>%
      filter(pval < 0.05) %>% filter(abs(log2FC) > log2(1.5))

    lung_microbiome_otu_df <- lung_microbiome_otu_df %>%
      rownames_to_column(var="OTU") %>% filter(OTU %in% lung_microbiome_sig_df$OTU) %>% column_to_rownames(var='OTU')
    gut_microbiome_otu_df <- gut_microbiome_otu_df %>%
      rownames_to_column(var="OTU") %>% filter(OTU %in% gut_microbiome_sig_df$OTU) %>% column_to_rownames(var='OTU')
    lung_metabolome_otu_df <- lung_metabolome_otu_df %>%
      rownames_to_column(var="OTU") %>% filter(OTU %in% lung_metabolome_sig_df$OTU) %>% column_to_rownames(var='OTU')
    plasma_metabolome_otu_df <- plasma_metabolome_otu_df %>%
      rownames_to_column(var="OTU") %>% filter(OTU %in% plasma_metabolome_sig_df$OTU)%>% column_to_rownames(var='OTU')
  }

  lung_microbiome_otu_df <- lung_microbiome_otu_df %>% t
  gut_microbiome_otu_df <- gut_microbiome_otu_df %>% t
  lung_metabolome_otu_df <- lung_metabolome_otu_df %>% t
  plasma_metabolome_otu_df <- plasma_metabolome_otu_df %>% t

  tidy_mediation_result <- function(mediation_result){
    re <- mediation_result$effects
    ACME <- as.numeric(re[1,"estimate"])
    ADE <- as.numeric(re[2,"estimate"])
    TE <- as.numeric(re[3,"estimate"])
    med_prop <-ACME/TE
    clean_re <- tibble(ACME=ACME, ADE=ADE, TE=TE, med_prop=med_prop)
  }
  set.seed(42)
  plan(multisession, workers = 10)
  out_tbl1 <-crossing(colnames(gut_microbiome_otu_df),colnames(lung_microbiome_otu_df)) %>%
    dplyr::rename(x=`colnames(gut_microbiome_otu_df)`, y=`colnames(lung_microbiome_otu_df)`) %>%
    rowwise() %>%
    mutate(x_data=list(gut_microbiome_otu_df[,x]),y_data=list(lung_microbiome_otu_df[,y]),M="Plasma metabolome",m_data=list(plasma_metabolome_otu_df)) %>%
    ungroup() %>%
    mutate(results = future_pmap(list(x_data,m_data,y_data), mediate_hima, .options=furrr_options(seed = 42))) %>%
    mutate(global = future_map(results, tidy_mediation_result, .options=furrr_options(seed = 42))) %>%
    mutate(selected_features = future_map(results, ~.$contributions, .options=furrr_options(seed = 42))) %>%
    dplyr::select(-x_data,-m_data,-y_data, -results) %>%
    unnest(cols = c(global)) %>%
    mutate(X="Gut microbiome", Y="Lung microbiome")
  set.seed(42)
  plan(multisession, workers = 10)
  out_tbl2 <-crossing(colnames(lung_microbiome_otu_df),colnames(gut_microbiome_otu_df)) %>%
    dplyr::rename(x=`colnames(lung_microbiome_otu_df)`, y=`colnames(gut_microbiome_otu_df)`) %>%
    rowwise() %>%
    mutate(x_data=list(lung_microbiome_otu_df[,x]),y_data=list(gut_microbiome_otu_df[,y]),M="Plasma metabolome",m_data=list(plasma_metabolome_otu_df)) %>%
    ungroup() %>%
    mutate(results = future_pmap(list(x_data,m_data,y_data), mediate_hima, .options=furrr_options(seed = 42))) %>%
    mutate(global = future_map(results, tidy_mediation_result, .options=furrr_options(seed = 42))) %>%
    mutate(selected_features = future_map(results, ~.$contributions, .options=furrr_options(seed = 42))) %>%
    dplyr::select(-x_data,-m_data,-y_data, -results) %>%
    unnest(cols = c(global)) %>%
    mutate(X="Lung microbiome", Y="Gut microbiome")
  set.seed(42)
  plan(multisession, workers = 10)
  out_tbl3 <-crossing(colnames(gut_microbiome_otu_df),colnames(lung_metabolome_otu_df)) %>%
    dplyr::rename(x=`colnames(gut_microbiome_otu_df)`, y=`colnames(lung_metabolome_otu_df)`) %>%
    rowwise() %>%
    mutate(x_data=list(gut_microbiome_otu_df[,x]),y_data=list(lung_metabolome_otu_df[,y]),M="Plasma metabolome",m_data=list(plasma_metabolome_otu_df)) %>%
    ungroup() %>%
    mutate(results = future_pmap(list(x_data,m_data,y_data), mediate_hima, .options=furrr_options(seed = 42))) %>%
    mutate(global = future_map(results, tidy_mediation_result, .options=furrr_options(seed = 42))) %>%
    mutate(selected_features = future_map(results, ~.$contributions, .options=furrr_options(seed = 42))) %>%
    dplyr::select(-x_data,-m_data,-y_data, -results) %>%
    unnest(cols = c(global)) %>%
    mutate(X="Gut microbiome", Y="Lung metabolome")
  set.seed(42)
  plan(multisession, workers = 10)
  out_tbl4 <-crossing(colnames(lung_metabolome_otu_df),colnames(gut_microbiome_otu_df)) %>%
    dplyr::rename(x=`colnames(lung_metabolome_otu_df)`, y=`colnames(gut_microbiome_otu_df)`) %>%
    rowwise() %>%
    mutate(x_data=list(lung_metabolome_otu_df[,x]),y_data=list(gut_microbiome_otu_df[,y]),M="Plasma metabolome",m_data=list(plasma_metabolome_otu_df)) %>%
    ungroup() %>%
    mutate(results = future_pmap(list(x_data,m_data,y_data), mediate_hima, .options=furrr_options(seed = 42))) %>%
    mutate(global = future_map(results, tidy_mediation_result, .options=furrr_options(seed = 42))) %>%
    mutate(selected_features = future_map(results, ~.$contributions, .options=furrr_options(seed = 42))) %>%
    dplyr::select(-x_data,-m_data,-y_data, -results) %>%
    unnest(cols = c(global)) %>%
    mutate(X="Lung metabolome", Y="Gut microbiome")
  out_tbl <- bind_rows(out_tbl1,out_tbl2,out_tbl3,out_tbl4)
  return(out_tbl)
}


extract_hima_mediation <- function(wt_mediation_tbl, ppta_mediation_tbl){
  calculate_adjusted_pvalue <- . %>% mutate(padj_group = p.adjust(pval, method = "BH"))
  wt_mediation_tbl <- wt_mediation_tbl %>%
    mutate(path=paste(str_to_title(X),str_to_title(M),str_to_title(Y),sep=" -> ")) %>%
    unnest(cols=c(selected_features)) %>%
    dplyr::rename(pval=ab_pv) %>%
    group_by(path) %>%
    nest() %>%
    rowwise() %>%
    mutate(data = list(calculate_adjusted_pvalue(data))) %>%
    unnest()%>%
    ungroup() %>%
    mutate(padj_all = p.adjust(pval, method = "BH")) %>%
    mutate(group="WT")
  ppta_mediation_tbl <- ppta_mediation_tbl %>%
    mutate(path=paste(str_to_title(X),str_to_title(M),str_to_title(Y),sep=" -> ")) %>%
    unnest(cols=c(selected_features)) %>%
    dplyr::rename(pval=ab_pv) %>%
    group_by(path) %>%
    nest() %>%
    rowwise() %>%
    mutate(data = list(calculate_adjusted_pvalue(data))) %>%
    unnest() %>%
    ungroup() %>%
    mutate(padj_all = p.adjust(pval, method = "BH")) %>%
    mutate(group="PptA")
  merged_top_tbl <- bind_rows(wt_mediation_tbl,ppta_mediation_tbl)
  return(merged_top_tbl)
}

plot_hima_box <- function(pbs_mediation_tbl, wt_mediation_tbl, ppta_mediation_tbl, scale_med_prop=F){
  pbs_mediation_long_tbl <- pbs_mediation_tbl %>%
    mutate(path=paste(str_to_title(X),str_to_title(M),str_to_title(Y),sep=" -> ")) %>%
    mutate(id=paste0(x,y,M,X,Y)) %>%
    mutate(group="CA")
  wt_mediation_long_tbl <- wt_mediation_tbl %>%
    mutate(path=paste(str_to_title(X),str_to_title(M),str_to_title(Y),sep=" -> ")) %>%
    mutate(id=paste0(x,y,M,X,Y)) %>%
    mutate(group="WT")
  ppta_mediation_long_tbl <- ppta_mediation_tbl %>%
    mutate(path=paste(str_to_title(X),str_to_title(M),str_to_title(Y),sep=" -> ")) %>%
    mutate(id=paste0(x,y,M,X,Y)) %>%
    mutate(group="PptA")

  pbs_mediation_long_tbl <- pbs_mediation_long_tbl %>%
    filter(med_prop>=0)%>%
    filter(med_prop <=1)
  wt_mediation_long_tbl <- wt_mediation_long_tbl %>%
    filter(med_prop>=0)%>%
    filter(med_prop <=1)
  ppta_mediation_long_tbl <- ppta_mediation_long_tbl %>%
    filter(med_prop>=0) %>%
    filter(med_prop <=1)

  merged_long_tbl <- bind_rows(pbs_mediation_long_tbl, wt_mediation_long_tbl,ppta_mediation_long_tbl)
  merged_long_tbl <- merged_long_tbl %>%
    mutate(path=str_replace_all(path,"Gut Microbiome -> Plasma Metabolome -> Lung Metabolome", "gMicro -> pMeta -> lMeta")) %>%
    mutate(path=str_replace_all(path,"Gut Microbiome -> Plasma Metabolome -> Lung Microbiome", "gMicro -> pMeta -> lMicro")) %>%
    mutate(path=str_replace_all(path,"Lung Metabolome -> Plasma Metabolome -> Gut Microbiome", "lMeta -> pMeta -> gMicro")) %>%
    mutate(path=str_replace_all(path,"Lung Microbiome -> Plasma Metabolome -> Gut Microbiome", "lMicro -> pMeta -> gMicro")) %>%
    mutate(group = factor(group, levels=c("CA","WT","PptA")))

  if(scale_med_prop){
    merged_long_tbl <-   merged_long_tbl %>%
      mutate(med_prop = ifelse(med_prop <0, 0, med_prop)) %>%
      mutate(med_prop = ifelse(med_prop>1, 1, med_prop))
  }

  color <- get_colormap_and_charmap() %>%
    filter(group %in% unique(merged_long_tbl$group)) %>% pull(color)

  mediation_plot <-ggplot(merged_long_tbl, aes(x=group, y=med_prop, fill=group)) +
  geom_boxplot()+
  facet_wrap(~path, scales="free")  +
    stat_compare_means(comparisons=list(c("CA","WT"),c("CA","PptA"),c("WT","PptA")), method = "wilcox.test",label = "p.format", size = 4)+
    ylab("Mediation proportion") +
    xlab("") +
    theme_pubr() +
    theme(
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      plot.title = element_text(size = 24, hjust = 0.5),
      strip.text.x = element_text(size = 18, color = "black"),
      strip.background = element_blank(),
      text = element_text(size = 18),
      legend.text = element_text(size = 18),
      legend.title = element_blank(),
      legend.position = "top"
    )+
    ylim(0,1.5)+
    coord_cartesian(clip="off")+
    scale_fill_manual(values=color)
}

get_pcoa_reduced_mediation_combination <- function(omics_tbl,k=5){
  lung_microbiome_otu_df <- omics_tbl %>% filter(name=="lung_microbiome") %>% pull(data) %>% `[[`(1) %>% phyloseq::distance(method="unifrac", weighted=T) %>%
    cmdscale(k = k) %>% as.data.frame
  gut_microbiome_otu_df <- omics_tbl %>% filter(name=="gut_microbiome") %>%pull(data) %>% `[[`(1) %>% phyloseq::distance(method="unifrac", weighted=T) %>%
    cmdscale(k = k) %>% as.data.frame
  lung_metabolome_otu_df <- omics_tbl %>% filter(name=="lung_metabolome") %>%pull(data) %>% `[[`(1) %>% phyloseq::distance(method="euclidean",) %>%
    cmdscale(k = k) %>% as.data.frame
  plasma_metabolome_otu_df <- omics_tbl %>% filter(name=="plasma_metabolome") %>%pull(data) %>% `[[`(1) %>% phyloseq::distance(method="euclidean") %>%
    cmdscale(k = k) %>% as.data.frame
  # Cumulative corrected relative eigenvalues
  lung_microbiome_cum_eig <- omics_tbl %>% filter(name=="lung_microbiome") %>%pull(data) %>% `[[`(1) %>%
    phyloseq::distance(method="unifrac", weighted=T) %>% pcoa() %>% `[[`("values") %>%  `[[`("Cum_corr_eig") %>% `[[`(k)
  gut_microbiome_cum_eig <- omics_tbl %>% filter(name=="gut_microbiome") %>%pull(data) %>% `[[`(1) %>%
    phyloseq::distance(method="unifrac", weighted=T) %>% pcoa() %>% `[[`("values") %>%  `[[`("Cum_corr_eig") %>% `[[`(k)
  lung_metabolome_cum_eig <- omics_tbl %>% filter(name=="lung_metabolome") %>%pull(data) %>% `[[`(1) %>%
    phyloseq::distance(method="euclidean") %>% pcoa() %>% `[[`("values") %>%  `[[`("Cumul_eig") %>% `[[`(k)
  plasma_metabolome_cum_eig <- omics_tbl %>% filter(name=="plasma_metabolome") %>%pull(data) %>% `[[`(1) %>%
    phyloseq::distance(method="euclidean") %>% pcoa() %>% `[[`("values") %>%  `[[`("Cumul_eig") %>% `[[`(k)

  out_tbl <- tibble(X=c("lung_microbiome_otu_df","gut_microbiome_otu_df","lung_metabolome_otu_df","plasma_metabolome_otu_df"),
                    M=c("lung_microbiome_otu_df","gut_microbiome_otu_df","lung_metabolome_otu_df","plasma_metabolome_otu_df"),
                    Y=c("lung_microbiome_otu_df","gut_microbiome_otu_df","lung_metabolome_otu_df","plasma_metabolome_otu_df")) %>%
    tidyr::expand(X,M,Y) %>%
    filter(!X==Y) %>%
    filter(!X==M) %>%
    filter(!M==Y) %>%
    arrange(X,M,Y) %>%
    rowwise() %>%
    mutate(X_data=list(get(X)),M_data=list(get(M)),Y_data=list(get(Y))) %>%
    mutate(X=str_remove(X,"_otu_df"), M=str_remove(M,"_otu_df"),Y=str_remove(Y,"_otu_df"))

  get_all_possible_combination <- function(X_data, M_data, Y_data){
    combined_df <- expand_grid(x = colnames(X_data), m = colnames(M_data), y = colnames(Y_data)) %>%
      mutate(
        X_vec = map2(x, list(X_data), ~ .y[[.x]]),
        M_vec = map2(m, list(M_data), ~ .y[[.x]]),
        Y_vec = map2(y, list(Y_data), ~ .y[[.x]])
      ) %>%
      rowwise() %>%
      mutate(data=list(tibble(x=X_vec,m=M_vec,y=Y_vec))) %>%
      dplyr::select(-X_vec,-M_vec,-Y_vec) %>%
      ungroup()
  }
  out_tbl <- out_tbl %>% rowwise() %>% mutate(data=list(get_all_possible_combination(X_data, M_data, Y_data))) %>% ungroup() %>%
    mutate(cum_eig=list(list(lung_microbiome=lung_microbiome_cum_eig,
                         gut_microbiome=gut_microbiome_cum_eig,
                         lung_metabolome=lung_metabolome_cum_eig,
                         plasma_metabolome=plasma_metabolome_cum_eig)))
}


extract_pcoa_reduced_mediation <- function(re_tbl){
  re_tbl_extraction <- re_tbl %>%
    dplyr::select(!contains("data"),-cum_eig) %>%
    unnest(results) %>%
    dplyr::select(-x,-m,-y)
  return(re_tbl_extraction)
}


plot_pcoa_reduced_mediation_bar <- function(pbs_re_tbl, wt_re_tbl, ppta_re_tbl){
  merged_tbl <- wt_re_tbl %>%
    mutate(group="WT") %>%
    bind_rows(ppta_re_tbl %>% mutate(group="PptA")) %>%
    bind_rows(pbs_re_tbl %>% mutate(group="CA")) %>%
    #filter(str_detect(M,"metabolome")) %>%
    #mutate(med_prop=ifelse(med_prop <0, 0, med_prop)) %>%
    mutate(path=paste(X,M,Y, sep="->"))
  ggplot(data=merged_tbl,  aes(x=group, y=med_prop, fill=group)) +
    geom_boxplot() +
    stat_compare_means(comparisons=list(c("CA","WT"),c("CA","PptA"),c("WT","PptA")))+
    facet_wrap(~path, scale="free")
}


get_ligilactobacillus_mediation_combination <- function(omics_tbl,sig_object_tbl){
  gut_microbiome_sig_df <- sig_object_tbl %>% filter(name=="gut_microbiome") %>%pull(data) %>% `[[`(1) %>%
    filter(comparison %in% c("CA vs WT","CA vs PptA") ) %>% filter(pval < 0.05)
  lung_metabolome_sig_df <- sig_object_tbl %>% filter(name=="lung_metabolome") %>%pull(data) %>% `[[`(1) %>% dplyr::rename(OTU=metabolite) %>% 
    filter(comparison %in% c("CA vs WT","CA vs PptA") ) %>% filter(pval < 0.05)
  plasma_metabolome_sig_df <- sig_object_tbl %>% filter(name=="plasma_metabolome") %>%pull(data) %>% `[[`(1) %>%  dplyr::rename(OTU=metabolite) %>% 
    filter(comparison %in% c("CA vs WT","CA vs PptA") ) %>% filter(pval < 0.05)

  lung_microbiome_otu_df <- omics_tbl %>% filter(name=="lung_microbiome") %>% pull(data) %>% `[[`(1) %>%
    otu_table() %>% as.data.frame() %>% rownames_to_column(var="OTU") %>%
    filter(OTU=="Bacteria_Firmicutes_Bacilli_Lactobacillales_Lactobacillaceae_Ligilactobacillus") %>%
    column_to_rownames(var='OTU') %>% rotate_df()
  gut_microbiome_otu_df <- omics_tbl %>% filter(name=="gut_microbiome") %>%pull(data) %>% `[[`(1) %>% otu_table() %>% as.data.frame() %>% 
    rownames_to_column(var="OTU") %>% filter(OTU %in% gut_microbiome_sig_df$OTU) %>% column_to_rownames(var='OTU') %>% rotate_df()
  lung_metabolome_otu_df <- omics_tbl %>% filter(name=="lung_metabolome") %>%pull(data) %>% `[[`(1) %>% otu_table() %>% as.data.frame() %>% 
    rownames_to_column(var="OTU") %>% filter(OTU %in% lung_metabolome_sig_df$OTU) %>% column_to_rownames(var='OTU') %>% rotate_df()
  plasma_metabolome_otu_df <- omics_tbl %>% filter(name=="plasma_metabolome") %>%pull(data) %>% `[[`(1) %>% otu_table() %>% as.data.frame() %>% 
    rownames_to_column(var="OTU") %>% filter(OTU %in% plasma_metabolome_sig_df$OTU)%>% column_to_rownames(var='OTU') %>% rotate_df()
  out_tbl <- tibble(X=c("lung_microbiome_otu_df","gut_microbiome_otu_df"),
                    M=c("plasma_metabolome_otu_df","plasma_metabolome_otu_df"),
                    Y=c("gut_microbiome_otu_df","lung_microbiome_otu_df")) %>%
    rowwise() %>%
    mutate(X_data=list(get(X)),M_data=list(get(M)),Y_data=list(get(Y))) %>% 
    mutate(X=str_remove(X,"_otu_df"), M=str_remove(M,"_otu_df"),Y=str_remove(Y,"_otu_df")) %>%
    mutate(path=paste(str_to_title(X),str_to_title(M),str_to_title(Y),sep=" -> "))

  get_all_possible_combination <- function(X_data, M_data, Y_data){
    combined_df <- expand_grid(x = colnames(X_data), m = colnames(M_data), y = colnames(Y_data)) %>%
      mutate(
        X_vec = map2(x, list(X_data), ~ .y[[.x]]),
        M_vec = map2(m, list(M_data), ~ .y[[.x]]),
        Y_vec = map2(y, list(Y_data), ~ .y[[.x]])
      ) %>%
      rowwise() %>%
      mutate(data=list(tibble(x=X_vec,m=M_vec,y=Y_vec))) %>%
      dplyr::select(-X_vec,-M_vec,-Y_vec) %>%
      ungroup()
  }
  out_tbl <- out_tbl %>% rowwise() %>% mutate(data=list(get_all_possible_combination(X_data, M_data, Y_data))) %>% ungroup()
  return(out_tbl)
}


get_mediation <- function(data_tbl){
  library(future)
  library(furrr)
  get_mediation_one <- function(dat){
    x2y_beta=NA
    x2y_p=NA
    x2m_beta=NA
    x2m_p=NA
    xm2y_beta=NA
    xm2y_p=NA
    med_total = NA
    med_ade=NA
    med_acmd=NA
    med_p=NA
    med_prop=NA
    tryCatch({
      x2y.fit <- lm(y ~ x, data = dat)
      x2y_beta <- summary(x2y.fit)$coefficients[2,1]
      x2y_p <- summary(x2y.fit)$coefficients[2,4]
      x2m.fit <- lm(m ~ x, data = dat)
      x2m_beta <- summary(x2m.fit)$coefficients[2,1]
      x2m_p <- summary(x2m.fit)$coefficients[2,4]
      xm2y.fit <- lm(y ~ x + m, data = dat)
      xm2y_beta <- summary(xm2y.fit)$coefficients[3,1]
      xm2y_p <- summary(xm2y.fit)$coefficients[3,4]
      med.out <- mediate(x2m.fit, xm2y.fit, treat = "x", mediator = "m")
      med_prop <- med.out$n.avg
      med_total <-med.out$tau.coef
      med_ade <- med.out$z.avg
      med_acmd <- med.out$d.avg
      med_p <- med.out$d.avg.p
    },error = function(e){
      print(e)
    })
    out_tbl <- tibble(x2y_beta=x2y_beta,x2y_p=x2y_p,
                      x2m_beta=x2m_beta,x2m_p=x2m_p,
                      xm2y_beta=xm2y_beta,xm2y_p=xm2y_p,
                      med_total = med_total,med_ade=med_ade,med_acmd=med_acmd,
                      med_p=med_p,med_prop=med_prop)
    return(out_tbl)
  }
  set.seed(42)
  plan(multisession, workers = 10)
  out_tbl <- data_tbl %>% ungroup() %>% mutate(results = future_map(data, get_mediation_one, .options=furrr_options(seed = 42))) %>%
    dplyr::select(-data) %>%
    unnest(results)
  return(out_tbl)
}

extract_mediation <- function(wt_mediation_tbl, ppta_mediation_tbl){
  calculate_adjusted_pvalue <- . %>% mutate(padj_group = p.adjust(pval, method = "BH"))
  wt_mediation_tbl <- wt_mediation_tbl %>%
    mutate(path=paste(str_to_title(X),str_to_title(M),str_to_title(Y),sep=" -> ")) %>%
    unnest(cols=c(results)) %>%
    mutate(pval=med_p) %>%
    group_by(path) %>%
    nest() %>%
    rowwise() %>%
    mutate(data = list(calculate_adjusted_pvalue(data))) %>%
    unnest()%>%
    ungroup() %>%
    mutate(padj_all = p.adjust(pval, method = "BH")) %>%
    dplyr::select(-X_data,-M_data,-Y_data,-data) %>%
    mutate(group="WT") 
  ppta_mediation_tbl <- ppta_mediation_tbl %>%
    mutate(path=paste(str_to_title(X),str_to_title(M),str_to_title(Y),sep=" -> ")) %>%
    unnest(cols=c(results)) %>%
    mutate(pval=med_p) %>%
    group_by(path) %>%
    nest() %>%
    rowwise() %>%
    mutate(data = list(calculate_adjusted_pvalue(data))) %>%
    unnest() %>%
    ungroup() %>%
    mutate(padj_all = p.adjust(pval, method = "BH")) %>%
    dplyr::select(-X_data,-M_data,-Y_data,-data) %>%
    mutate(group="PptA")
  merged_top_tbl <- bind_rows(wt_mediation_tbl,ppta_mediation_tbl) %>%
    filter(med_prop>=0)%>%
    filter(med_prop <=1)
  return(merged_top_tbl)
}