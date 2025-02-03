get_fig2_comparison_df <- function() {
  comparison_df <-data.frame(
    G1 = c(
      "CTR",
      "CA"
    ),
    G2 = c(
      "CA",
      "VOR"
    )
  )
  comparison_df$comparison <-
    paste(comparison_df$G1, "vs", comparison_df$G2)
  rownames(comparison_df) <- comparison_df$comparison
  return(comparison_df)
}

# get all fig2 samples data
get_fig2_physeq <- function(physeq) {
  filtered_physeq <- physeq %>%
    subset_samples(sampleGroup %in% 1:3) %>%
    filter_taxa(function(x)
      sum(x > 0) > 0, TRUE)
  return(filtered_physeq)
}

plot_fig2_alpha <- function(physeq,plot_title=""){
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
  ) + stat_compare_means(method = "wilcox.test",comparisons=list(c("CTR","CA"),c("CA","VOR")),label = "p.format",
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
  ) + stat_compare_means(method = "wilcox.test",comparisons=list(c("CTR","CA"),c("CA","VOR")),label = "p.format",
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

plot_fig2_beta <- function(physeq, distance="unifrac", weighted=T, plot_title="") {
  meta_df <- data.frame(sample_data(physeq))
  color <- get_colormap_and_charmap() %>%
    filter(group %in% meta_df$group) %>% pull(color)
  ord_mat <-
    ordinate(physeq,
             method = "PCoA",
             distance = distance,
             weighted = weighted)
  dist_mat1 <-
    phyloseq::distance(physeq %>% subset_samples(group %in% c("CTR","CA")), method = distance, weighted = weighted)
  meta_df1 <-  data.frame(sample_data(physeq %>% subset_samples(group %in% c("CTR","CA"))))
  set.seed(42)
  adonis_result1 <-
    adonis2(dist_mat1 ~ group, data = meta_df1, permutations = 999)
  p1 <- adonis_result1$`Pr(>F)`[1]
  r1 <- round(adonis_result1$`R2`[1],3)

  dist_mat2 <-
    phyloseq::distance(physeq %>% subset_samples(group %in% c("CA","VOR")), method = distance, weighted = weighted)
  meta_df2<-  data.frame(sample_data(physeq %>% subset_samples(group %in% c("CA","VOR"))))
  set.seed(42)
  adonis_result2 <-
    adonis2(dist_mat2 ~ group, data = meta_df2, permutations = 999)
  p2 <- adonis_result2$`Pr(>F)`[1]
  r2 <- round(adonis_result2$`R2`[1],3)

  beta_boxplot <- plot_ordination(physeq, ord_mat, color = "group") +
    geom_point(size=5)+
    stat_ellipse(type = "norm",
                 linetype = 1,
                 lwd = 0.5) +
    stat_mean(aes(color=group),size=8,shape=13)+
    theme_pubr() +
    theme(
      plot.title = element_text(size = 1, hjust = 0.5),
      plot.subtitle = element_markdown(size = 22, hjust = 0.5),
      text = element_text(size = 22),
      axis.title.x=element_text(colour='black', size=25),
      axis.title.y=element_text(colour='black', size=25),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1
      ),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      legend.position = "top"
    ) +
    labs(subtitle = paste("CTR vs CA: r<sup>2</sup> =", r1,", p = ", p1,
                     "<br />CA vs VOR: r<sup>2</sup> =", r2,", p = ", p2), title = plot_title) +  scale_color_manual(values=color)
  return(beta_boxplot)
}

plot_fig2_plsda <- function(physeq, plot_title="") {
  meta_df <- data.frame(sample_data(physeq))
  color <- get_colormap_and_charmap() %>%
    filter(group %in% meta_df$group) %>% pull(color)
  otu_df <- as.data.frame(otu_table(physeq)) %>% rotate_df()

  plsda = opls(otu_df,meta_df$group,predI=2,scaleC="standard",fig.pdfC="none")
  q2 =  plsda@summaryDF$`Q2(cum)`

  adonis_result <- adonis2(otu_df~meta_df$group, method="euclidean")

  dist_mat1 <-
    phyloseq::distance(physeq %>% subset_samples(group %in% c("CTR","CA")), method = "euclidean")
  meta_df1 <-  data.frame(sample_data(physeq %>% subset_samples(group %in% c("CTR","CA"))))
  set.seed(42)
  adonis_result1 <-
    adonis2(dist_mat1 ~ group, data = meta_df1, permutations = 999)
  p1 <- adonis_result1$`Pr(>F)`[1]
  r1 <- round(adonis_result1$`R2`[1],3)

  dist_mat2 <-
    phyloseq::distance(physeq %>% subset_samples(group %in% c("CA","VOR")), method = "euclidean")
  meta_df2<-  data.frame(sample_data(physeq %>% subset_samples(group %in% c("CA","VOR"))))
  set.seed(42)
  adonis_result2 <-
    adonis2(dist_mat2 ~ group, data = meta_df2, permutations = 999)
  p2 <- adonis_result2$`Pr(>F)`[1]
  r2 <- round(adonis_result2$`R2`[1],3)

  sample.score = plsda@scoreMN %>%
    as.data.frame() %>%
    mutate(group = meta_df$group)

  beta_boxplot = ggplot(sample.score, aes(p1, p2, color = group)) +
    geom_point(size=5) +
    geom_point(aes(-10,-10), color = 'white') +
    stat_ellipse(type = "norm",
                 linetype = 1,
                 lwd = 0.5) +
    stat_mean(aes(color=group),size=8,shape=13)+
    annotate("text", x = 0.8*ceiling(min(sample.score$p1)), y = -10, label = paste("Q2(cum)=",q2), size = 7) +
    theme_pubr() +
    theme(
      plot.title = element_text(size = 1, hjust = 0.5),
      plot.subtitle = element_markdown(size = 22, hjust = 0.5),
      text = element_text(size = 22),
      axis.title.x=element_text(colour='black', size=25),
      axis.title.y=element_text(colour='black', size=25),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      legend.position = "top",
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 1
      )
    ) +
    labs(subtitle = paste("CTR vs CA: r<sup>2</sup> =", r1,", p = ", p1,
                     "<br />CA vs VOR: r<sup>2</sup> =", r2,", p = ", p2), title = plot_title) +
    labs(x = paste0("P1=",plsda@modelDF$R2X[1] *100,"%") ,y = paste0("P1=",plsda@modelDF$R2X[2] *100,"%"))+
    scale_color_manual(values=color)
  return(beta_boxplot)
}

plot_stacked_sig_taxa_barplot <- function(deseq2_df,n=10){
  deseq2_filter_df <- deseq2_df %>%
    filter(Genus != "?") %>%
    filter(pval < 0.05) %>%
    filter(abs(log2FC) > log2(1.5)) %>%
    group_by(comparison) %>%
    arrange(pval) %>%
    slice_head(n=n)  %>%
    ungroup() %>%
    mutate(color = log2FC < 0)%>%
    mutate(mark = ifelse(padj < 0.1, "★", ""))

  deseq2_filter_df_subset1 <- deseq2_filter_df %>%
    filter(comparison == "CTR vs CA") %>%
    arrange(log2FC, Genus) %>%
    mutate(Genus=paste0("*",Genus,"*"),
           Genus=ifelse(str_detect(Genus,"Family"),str_replace(Genus," Family\\*","* Family"),Genus),
           Genus=ifelse(str_detect(Genus,"Order"),str_replace(Genus," Order\\*","* Order"),Genus),
           Genus=ifelse(str_detect(Genus,"Class"),str_replace(Genus," Class\\*","* Class"),Genus),
           Genus=ifelse(str_detect(Genus,"Phylum"),str_replace(Genus," Phylum\\*","* Phylum"),Genus),
           Genus=ifelse(str_detect(Genus,"Domain"),str_replace(Genus," Domain\\*","* Domain"),Genus)
    ) %>% 
    mutate(Genus= factor(Genus, levels = Genus))

  plot_subset1 <- ggplot(deseq2_filter_df_subset1, aes(x = Genus, y = log2FC, fill=color, label=mark)) +
    geom_bar(stat = "identity",
             show.legend = FALSE) +
    geom_text(aes(hjust = 0.5 - sign(log2FC)/1.5),size=7)+
    coord_flip(ylim=c(round(1.6*min(deseq2_filter_df_subset1$log2FC)),round(1.6*max(deseq2_filter_df_subset1$log2FC)))) +
    theme_pubr() +
    theme(plot.title = element_text(size = 26, hjust = 0.5),
          text = element_text(size = 22),
          axis.text.y=element_markdown(),
          axis.title.x=element_text(colour='black', size=25),
          axis.title.y=element_text(colour='black', size=25))+
    scale_fill_manual(values=c('#C466A0','#66C48A'))+
    labs(title="CTR vs CA") + xlab(paste0("Genus"))

  deseq2_filter_df_subset2 <- deseq2_filter_df %>%
    filter(comparison == "CA vs VOR") %>%
    arrange(log2FC, Genus) %>%
    mutate(Genus=paste0("*",Genus,"*"),
           Genus=ifelse(str_detect(Genus,"Family"),str_replace(Genus," Family\\*","* Family"),Genus),
           Genus=ifelse(str_detect(Genus,"Order"),str_replace(Genus," Order\\*","* Order"),Genus),
           Genus=ifelse(str_detect(Genus,"Class"),str_replace(Genus," Class\\*","* Class"),Genus),
           Genus=ifelse(str_detect(Genus,"Phylum"),str_replace(Genus," Phylum\\*","* Phylum"),Genus),
           Genus=ifelse(str_detect(Genus,"Domain"),str_replace(Genus," Domain\\*","* Domain"),Genus)
    ) %>% 
    mutate(Genus= factor(Genus, levels = Genus))
  plot_subset2 <- ggplot(deseq2_filter_df_subset2, aes(x = Genus, y = log2FC, fill=color,label=mark)) +
    geom_bar(stat = "identity",
             show.legend = FALSE) +
    geom_text(aes(hjust = 0.5 - sign(log2FC)/1.5),size=7)+
    coord_flip(ylim=c(round(1.6*min(deseq2_filter_df_subset2$log2FC)),round(1.6*max(deseq2_filter_df_subset2$log2FC)))) +
    theme_pubr() +
    theme(plot.title = element_text(size = 26, hjust = 0.5),
          text = element_text(size = 22),
          axis.text.y=element_markdown(),
          axis.title.x=element_text(colour='black', size=25),
          axis.title.y=element_text(colour='black', size=25))+
    scale_fill_manual(values=c('#C466A0','#66C48A'))+
    labs(title="CA vs VOR")+ xlab(paste0("Genus"))
  return(list(plot_subset1,plot_subset2))
}

plot_stacked_sig_metabolite_barplot <- function(sig_metabolites_df,n=10){
  filter_df <- sig_metabolites_df %>%
    filter(pval < 0.05) %>%
    filter(abs(log2FC) > log2(1.5)) %>%
    group_by(comparison) %>%
    arrange(pval) %>%
    slice_head(n=n) %>%
    ungroup() %>%
    mutate(color = log2FC < 0)%>%
    mutate(mark = ifelse(padj < 0.1, "★", ""))

  filter_df_subset1 <- filter_df %>%
    filter(comparison == "CTR vs CA") %>%
    arrange(log2FC, metabolite) %>%
    mutate(metabolite = str_replace(metabolite, "1,4a-Dimethyl-6-methylene-5-\\[2-\\(2-oxo-2,5-dihydro-3-furanyl\\)ethyl\\]decahydro-1-naphthalenecarboxylic acid", "CID:14806212")) %>%
    mutate(metabolite = str_replace(metabolite, "5,6-dimethyl-4-oxo-4H-pyran-2-carboxylic acid", "CID:2786724")) %>%
    mutate(metabolite= factor(metabolite, levels = metabolite))

  plot_subset1 <- ggplot(filter_df_subset1, aes(x = metabolite, y = log2FC, fill=color, label=mark)) +
    geom_bar(stat = "identity",
             show.legend = FALSE) +
    geom_text(aes(hjust = 0.5 - sign(log2FC)/1.5),size=7)+
    coord_flip(ylim=c(round(1.6*min(filter_df_subset1$log2FC)),round(1.6*max(filter_df_subset1$log2FC)))) +
    theme_pubr() +
    theme(plot.title = element_text(size = 26, hjust = 0.5),
          text = element_text(size = 22),
          axis.title.x=element_text(colour='black', size=25),
          axis.title.y=element_text(colour='black', size=25))+
    scale_fill_manual(values=c('#C466A0','#66C48A'))+
    labs(title="CTR vs CA")+ xlab(paste0("Metabolite"))

  filter_df_subset2 <- filter_df %>%
    filter(comparison == "CA vs VOR") %>%
    arrange(log2FC, metabolite) %>%
    mutate(metabolite = str_replace(metabolite, "\\(8aR,12S,12aR\\)-12-Hydroxy-4-methyl-4,5,6,7,8,8a,12,12a-octahydro-2H-3-benzoxecine-2,9\\(1H\\)-dione", "CID:637324")) %>%
    mutate(metabolite = str_replace(metabolite, "3-\\(propan-2-yl\\)-octahydropyrrolo\\[1,2-a\\]pyrazine-1,4-dione", "CID:98951")) %>%
    mutate(metabolite = str_replace(metabolite, "5,6-dimethyl-4-oxo-4H-pyran-2-carboxylic acid", "CID:2786724")) %>%
    mutate(metabolite = str_replace(metabolite, "2-\\[\\(1S\\)-1-Hydroxyethyl\\]-4\\(1H\\)-quinazolinone", "CID:11116825")) %>%
    mutate(metabolite= factor(metabolite, levels = metabolite))

  plot_subset2 <- ggplot(filter_df_subset2, aes(x = metabolite, y = log2FC, fill=color,label=mark)) +
    geom_bar(stat = "identity",
             show.legend = FALSE) +
    geom_text(aes(hjust = 0.5 - sign(log2FC)/1.5),size=7)+
    coord_flip(ylim=c(round(1.6*min(filter_df_subset2$log2FC)),round(1.6*max(filter_df_subset2$log2FC)))) +
    theme_pubr() +
    theme(plot.title = element_text(size = 26, hjust = 0.5),
          text = element_text(size = 22),
          axis.title.x=element_text(colour='black', size=25),
          axis.title.y=element_text(colour='black', size=25))+
    scale_fill_manual(values=c('#C466A0','#66C48A'))+
    labs(title="CA vs VOR")+ xlab(paste0("Metabolite"))
  return(list(plot_subset1,plot_subset2))
}

get_omics_similarity <- function(omics_tbl, omics_clr_tbl){
  lung_physeq <- omics_tbl %>% filter(name=="lung_microbiome") %>%pull(data) %>% `[[`(1)
  gut_physeq <- omics_tbl %>% filter(name=="gut_microbiome") %>%pull(data) %>% `[[`(1)
  lung_metabolome_physeq <- omics_tbl %>% filter(name=="lung_metabolome") %>%pull(data) %>% `[[`(1)
  plasma_metabolome_physeq <- omics_tbl %>% filter(name=="plasma_metabolome") %>%pull(data) %>% `[[`(1)

  lung_clr_otu_df <- omics_clr_tbl %>% filter(name=="lung_microbiome") %>%pull(data) %>% `[[`(1) %>% otu_table() %>% data.frame() %>% rotate_df()
  gut_clr_otu_df <- omics_clr_tbl %>% filter(name=="gut_microbiome") %>%pull(data) %>% `[[`(1) %>% otu_table() %>% data.frame() %>% rotate_df()
  lung_clr_metabolome_otu_df <- omics_clr_tbl %>% filter(name=="lung_metabolome") %>%pull(data) %>% `[[`(1) %>% otu_table() %>% data.frame() %>% rotate_df()
  plasma_clr_metabolome_otu_df <- omics_clr_tbl %>% filter(name=="plasma_metabolome") %>%pull(data) %>% `[[`(1) %>% otu_table() %>% data.frame() %>% rotate_df()

  lung_dist <- phyloseq::distance(lung_physeq,method = "wunifrac")
  gut_dist <- phyloseq::distance(gut_physeq, method = "wunifrac")
  lung_metabolome_dist <- phyloseq::distance(lung_metabolome_physeq,method = "euclidean")
  plasma_metabolome_dist <- phyloseq::distance(plasma_metabolome_physeq, method = "euclidean")

  comparison_dist_tbl <- list(lung_dist,gut_dist,lung_metabolome_dist,plasma_metabolome_dist) %>%
    combn(2) %>% t %>% as_tibble()
  comparison_clr_otu_tbl <- list(lung_clr_otu_df,gut_clr_otu_df,lung_clr_metabolome_otu_df,plasma_clr_metabolome_otu_df) %>%
    combn(2) %>% t %>% as_tibble()
  comparison_name_tbl <- c("lung_microbiome","gut_microbiome","lung_metabolome","plasma_metabolome") %>%
    combn(2) %>% t %>% as_tibble()

  comparison_tbl <- bind_cols(comparison_name_tbl, comparison_dist_tbl,comparison_clr_otu_tbl)
  colnames(comparison_tbl) <- c("G1_name", "G2_name", "G1_dist", "G2_dist", "G1_clr_otu", "G2_clr_otu")

  out_tbl <- comparison_tbl %>% rowwise() %>% mutate(mantel_result= list(mantel(G1_dist,G2_dist,method="spearman")),
                                                      RV_result=list(RV2(G1_clr_otu,G2_clr_otu)),
                                                      RV_p=RV2.p(G1_clr_otu,G2_clr_otu),
                                                      plots = list(plot_omics_similarity(G1_dist,G2_dist,
                                                      G1_clr_otu,G2_clr_otu,
                                                      G1_name, G2_name)))
  out_tbl$RV_cor <- sapply(out_tbl$RV_result, `[`,1)
  out_tbl$mantel_p <- sapply(out_tbl$mantel_result, `[[`, "signif")
  out_tbl$mantel_cor <- sapply(out_tbl$mantel_result, `[[`, "statistic")
  return(out_tbl)
}


plot_omics_similarity <- function(G1_dist,G2_dist, G1_clr_otu,G2_clr_otu, G1_name, G2_name){
  cor <-  round(RV2(G1_clr_otu,G2_clr_otu),3)
  p  <- round(RV2.p(G1_clr_otu,G2_clr_otu),3)
  mds.s <- monoMDS(G1_dist)
  mds.r <- monoMDS(G2_dist)

  pro.s.r <- procrustes(mds.s,mds.r)

  X <- data.frame(pro.s.r$X) %>% mutate(ID =1:nrow(.))
  colnames(X) <- c("X1","X2","ID")
  Y <- data.frame(pro.s.r$Yrot) %>% mutate(ID =1:nrow(.))
  colnames(Y) <- c("X1","X2","ID")
  plot_df <- bind_rows(X,Y)

  color1 = "#C6551F"
  color2 = "#8494FF"
  if (G1_name=="lung_microbiome"){
    G1_name <- "Microbiome (lung)"
    color1 <- "#C6551F"
    }
  if (G1_name=="gut_microbiome") {
    G1_name <- "Microbiome (gut)"
    color1 <- "#C6551F"}
  if (G2_name=="lung_metabolome"){
    G2_name <- "Metabolome (lung)"
    color2 <- "#8494FF"
    }
  if (G2_name=="plasma_metabolome"){
    G2_name <- "Metabolome (plasma)"
    color2 <- "#8494FF"
    }
  plot_df <- plot_df %>%
    mutate(omics=c(rep(G1_name,nrow(X)),rep(G2_name,nrow(Y))),
      color=c(rep(color1,nrow(X)),rep(color2,nrow(Y))))
  p <- ggplot(plot_df,aes(X1, X2,fill=omics)) +
      geom_point(size = 4, shape = 21) +
      geom_line(aes(group = ID), alpha = 0.6, colour = "grey50") +
      theme_pubr()+
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 0, hjust = 0.5),
          plot.subtitle = element_text(size = 22, hjust = 0.5),
          text = element_text(size = 22),
          #legend.position="right",
          legend.spacing.x = unit(0.1, 'cm'),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=25),
          axis.title.y=element_text(colour='black', size=25),
          axis.text=element_text(colour='black',size=22)) +
    labs(x = 'PC 1', y = 'PC 2') +
    labs(subtitle = paste("RV: r =", cor,", p = ", p)) +
    geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
    geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
    geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
    geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
    scale_fill_manual(values=c(color1,color2))
}

get_bacteria_metabolite_correlation <- function(omics_data) {
  gut_micro_tbl <- omics_data %>%
    filter(name == "gut_microbiome") %>%
    pull(data) %>% .[[1]] %>%
    transform_sample_counts(function(x) x/sum(x)) %>%
    melt_phyloseq() %>% 
    dplyr::select(sampleName, Genus, Abundance)
  plasma_metab_tbl <- omics_data %>% 
    filter(name == "plasma_metabolome") %>%
    pull(data) %>% .[[1]] %>%
    melt_phyloseq() %>% 
    dplyr::select(sampleName, metabolite, Abundance)
  correlation_tbl <- gut_micro_tbl %>%
    inner_join(plasma_metab_tbl, by = "sampleName",relationship = "many-to-many") %>%
    group_by(Genus, metabolite) %>%
    mutate(cor= cor(Abundance.x, Abundance.y,method = "spearman"),
           pval = cor.test(Abundance.x, Abundance.y,method = "spearman")$p.value) %>%
    ungroup() %>%
    dplyr::select(Genus, metabolite, cor, pval) %>%
    distinct() %>%
    mutate(cor = round(cor,3),
           group = "gut_microbiome and plasma_metabolome")
  
  lung_micro_tbl <- omics_data %>%
    filter(name == "lung_microbiome") %>%
    pull(data) %>% .[[1]] %>%
    transform_sample_counts(function(x) x/sum(x)) %>%
    melt_phyloseq() %>% 
    dplyr::select(sampleName, Genus, Abundance)
  lung_metab_tbl <- omics_data %>%
    filter(name == "lung_metabolome") %>%
    pull(data) %>% .[[1]] %>%
    melt_phyloseq() %>% 
    dplyr::select(sampleName, metabolite, Abundance)
  correlation_tbl2 <- lung_micro_tbl %>%
    inner_join(lung_metab_tbl, by = "sampleName",relationship = "many-to-many") %>%
    group_by(Genus, metabolite) %>%
    mutate(cor= cor(Abundance.x, Abundance.y,method = "spearman"),
           pval = cor.test(Abundance.x, Abundance.y,method = "spearman")$p.value) %>%
    ungroup() %>%
    dplyr::select(Genus, metabolite, cor, pval) %>%
    distinct() %>%
    mutate(cor = round(cor,3),
           group = "lung_microbiome and lung_metabolome")
  return(bind_rows(correlation_tbl,correlation_tbl2))
}