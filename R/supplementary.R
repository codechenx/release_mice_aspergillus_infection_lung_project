get_gut_baseline_beta_plot <- function(physeq, plot_title="") {
  distance = "unifrac"
  weighted = T
  physeq <- physeq %>% subset_samples(sampleGroup %in% 1:6)
  meta_df <- data.frame(sample_data(physeq))
  color <- get_colormap_and_charmap() %>%  pull(color)
  ord_mat <-
    ordinate(physeq,
             method = "PCoA",
             distance = distance,
             weighted = weighted)
  dist_mat<-
    phyloseq::distance(physeq, method = distance, weighted = weighted)
  set.seed(42)
  adonis_result <-
    adonis2(dist_mat ~ group, data = meta_df, permutations = 999)
  p <- adonis_result$`Pr(>F)`[1]
  r <- round(adonis_result$`R2`[1],3)

  beta_boxplot <- plot_ordination(physeq, ord_mat, color = "group") +
    geom_point(size=5) +
    stat_ellipse(type = "norm",
                 linetype = 1,
                 lwd = 0.5) +
    stat_mean(aes(color=group),size=8,shape=13)+
    theme_pubr() +
    theme(
      plot.title = element_text(size = 40, hjust = 0.5),
      plot.subtitle = element_markdown(size = 28, hjust = 0.5),
      text = element_text(size = 36),
      legend.text = element_text(size = 32),
      legend.title = element_blank(),
      legend.position = "top",
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 1
      )
    ) +
    labs(subtitle = paste("all: r<sup>2</sup> =", r,", p = ", p), title = plot_title)+scale_color_manual(values=color)
  return(beta_boxplot)
}


plot_taxa_composition <- function(physeq,taxrank="Genus",n_taxa=15) {
  physeq <- physeq %>% transform_sample_counts(function(x) x/sum(x)*100)
  if(taxrank=="Genus"){
    tax_table(physeq) <- tax_table(physeq) %>% as.data.frame %>%
      mutate(Genus=paste0("*",Genus,"*"),
             Genus=ifelse(str_detect(Genus,"Family"),str_replace(Genus," Family\\*","* Family"),Genus),
             Genus=ifelse(str_detect(Genus,"Order"),str_replace(Genus," Order\\*","* Order"),Genus),
             Genus=ifelse(str_detect(Genus,"Class"),str_replace(Genus," Class\\*","* Class"),Genus),
             Genus=ifelse(str_detect(Genus,"Phylum"),str_replace(Genus," Phylum\\*","* Phylum"),Genus),
             Genus=ifelse(str_detect(Genus,"Domain"),str_replace(Genus," Domain\\*","* Domain"),Genus)
      ) %>% 
      as.matrix() %>%
      tax_table
  } else if(taxrank=="Species"){
    tax_table(physeq) <- tax_table(physeq) %>% as.data.frame %>%
      mutate(Species=paste0("*",Species,"*")) %>%
      as.matrix() %>%
      tax_table
  }
  customPal <- tax_palette(
    data = physeq, rank = taxrank, pal = "brewerPlus", n = n_taxa, add = c(Other = "white")
  )
  plot_all <- physeq  %>%
    subset_taxa(!is.na(Domain)) %>%
    set_sample_order("group") %>%
    comp_barplot(n_taxa = n_taxa, tax_level = taxrank, group_by = "group",merge_other = F,palette = customPal)
  for(i in 1:length(plot_all)){
    plot_all[[i]] <- plot_all[[i]] + theme(plot.margin = margin(0, .8, 0, 0, "pt"))
  }
  composition_plot <- wrap_plots(plot_all, nrow = 1, guides = "collect",axes = "collect",axis_titles = "collect")&
    theme(plot.title=element_text(hjust = 0.5),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.key.size = unit(.3,"cm"),
          legend.title = element_text(size=10),
          legend.text = element_markdown(size=8),
    )&ylab("Abundance (%)")
  return(composition_plot)
}

get_power_analysis <- function(physeq){
  library(powmic)
  library(SpiecEasi)
  library(lqmm)
  #source("utils.R")
  otu_tbl <- physeq %>% 
    utils_filter_taxa(abundance_cutoff = 0.0001,prevalence_cutoff = 0.1) %>% 
    #transform_sample_counts(function(x) x/sum(x)*100) %>%
    otu_table %>%
    as.data.frame %>% 
    t
  distrib='NB'
  params=estParams(otu_tbl,Sigma=NULL,method='CCLasso',distrib=distrib)
  lmu0=log(params$mu)
  lphi0=log(params$phi)
  #lp0=log(params$p0)
  Sigma=params$Sigma %>% make.positive.definite(tol=1e-3)
  params.sim=setParams.NB(nTaxa=1000,p.DA=0.05,Sigma=Sigma,lmu0=lmu0,lphi0=lphi0,lfc.mu=1.5,sim.seed = 42)
  powmic.out=powmic(n1s=c(150), n2s=c(150),
                    params=params.sim,distrib=distrib,DAmethod='wilcox',nsims = 20)
  assess.out = assess(powmic.out, alpha.type='pval',alpha.level=0.05,stratify.type='prevalence')
  summaryAssess(assess.out,assess.type='overall')
}


plot_alpha <- function(physeq,plot_title=""){
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
  ) + stat_compare_means(method = "wilcox.test",comparisons=list(c("CTR","CA"),c("CA","VOR"),c("CA","WT"),c("CA","PptA"),c("WT","PptA"),c("VOR","Tx"),c("WT","VOR"),c("WT","Tx")),label = "p.format",
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
  ) + stat_compare_means(method = "wilcox.test",comparisons=list(c("CTR","CA"),c("CA","VOR"),c("CA","WT"),c("CA","PptA"),c("WT","PptA"),c("VOR","Tx"),c("WT","VOR"),c("WT","Tx")),label = "p.format",
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


plot_cage_effect <- function(physeq,plot_title=""){
  plot_beta_plot <- function(subset_physeq,plot_title="") {
    distance = "unifrac"
    weighted = T
    meta_df <- data.frame(sample_data(subset_physeq))
    ord_mat <-
      ordinate(subset_physeq,
               method = "PCoA",
               distance = distance,
               weighted = weighted)
    dist_mat<-
      phyloseq::distance(subset_physeq, method = distance, weighted = weighted)
    set.seed(42)
    adonis_result <-
      adonis2(dist_mat ~ cage, data = meta_df, permutations = 999)
    p <- adonis_result$`Pr(>F)`[1]
    r <- round(adonis_result$`R2`[1],3)
    
    beta_boxplot <- plot_ordination(subset_physeq, ord_mat, color = "cage") +
      geom_point(size=2) +
      stat_ellipse(type = "norm",
                   linetype = 1,
                   lwd = 0.5) +
      stat_mean(aes(color=cage),size=3,shape=13)+
      theme_pubr() +
      theme(
        plot.title = element_text(size = 10, hjust = 0.5),
        plot.subtitle = element_markdown(size = 7, hjust = 0.5),
        text = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.position = "top",
        panel.border = element_rect(
          colour = "black",
          fill = NA,
          size = 1
        )
      ) +
      labs(subtitle = paste("r<sup>2</sup> =", r,", p = ", p), title = plot_title)
    return(beta_boxplot)
  }
  
  sample_data(physeq) <- sample_data(physeq) %>% 
    data.frame() %>%
    mutate(cage =(as.numeric(sampleName)-1)%/%5,
           cage = as.factor(cage))%>% 
    sample_data()

  
  combined_tbl <- tibble(subset_group = c("CTR","CA","VOR","WT","Tx","PptA")) %>% 
    mutate(subset_physeq = map(subset_group,~subset_physeq_by_group(physeq, .x))) %>% 
    mutate(beta_plot = map2(subset_physeq,subset_group, ~plot_beta_plot(.x,.y)))
    
  wrap_plots(combined_tbl$beta_plot, nrow = 2)&
    theme(plot.title=element_text(hjust = 0.5),
          legend.key.size = unit(.1,"cm"),
          legend.title = element_text(size=5),
          legend.text = element_markdown(4),
    )
}


commbine_all_omics_similarity_plot <- function(omics_similarity_data){
  CTR_plots <-   omics_similarity_data %>% 
    filter(subset_group=="CTR") %>%
    pull(similarity_data) %>%
    `[[`(1) %>% 
    pull(plots)
  CTR_plot <- wrap_plots(list(CTR_plots[[2]],CTR_plots[[5]]))+plot_annotation(title = "CTR")&
    theme(plot.title=element_text(hjust = 0.5,size = 30))
  
  CA_plots <-   omics_similarity_data %>%
    filter(subset_group=="CA") %>%
    pull(similarity_data) %>%
    `[[`(1) %>% 
    pull(plots)
  CA_plot <- wrap_plots(list(CA_plots[[2]],CA_plots[[5]]))+plot_annotation(title = "CA")&
    theme(plot.title=element_text(hjust = 0.5,size = 30))
  
  VOR_plots <-   omics_similarity_data %>%
    filter(subset_group=="VOR") %>%
    pull(similarity_data) %>%
    `[[`(1) %>% 
    pull(plots)
  VOR_plot <- wrap_plots(list(VOR_plots[[2]],VOR_plots[[5]]))+plot_annotation(title = "VOR")&
    theme(plot.title=element_text(hjust = 0.5,size = 30))
  
  WT_plots <-   omics_similarity_data %>%
    filter(subset_group=="WT") %>%
    pull(similarity_data) %>%
    `[[`(1) %>% 
    pull(plots)
  WT_plot <- wrap_plots(list(WT_plots[[2]],WT_plots[[5]]))+plot_annotation(title = "WT")&
    theme(plot.title=element_text(hjust = 0.5,size = 30))
  
  Tx_plots <-   omics_similarity_data %>%
    filter(subset_group=="Tx") %>%
    pull(similarity_data) %>%
    `[[`(1) %>% 
    pull(plots)
  Tx_plot <- wrap_plots(list(Tx_plots[[2]],Tx_plots[[5]]))+plot_annotation(title = "Tx")&
    theme(plot.title=element_text(hjust = 0.5,size = 30))
  
  PptA_plots <-   omics_similarity_data %>%
    filter(subset_group=="PptA") %>%
    pull(similarity_data) %>%
    `[[`(1) %>% 
    pull(plots)
  PptA_plot <- wrap_plots(list(PptA_plots[[2]],PptA_plots[[5]]))+plot_annotation(title = "PptA")&
    theme(plot.title=element_text(hjust = 0.5,size = 30))
  
  
  ggarrange(plotlist = list(CTR_plot,CA_plot,VOR_plot,WT_plot,Tx_plot,PptA_plot),ncol = 2,nrow = 3)
}


commbine_all_omics_similarity_plot2 <- function(omics_similarity_data){
  CTR_CA_plots <-   omics_similarity_data %>% 
    .[["similarity_data"]] %>% 
    .[[1]] %>% 
    pull(plots)
  CTR_CA_plot <- wrap_plots(list(CTR_CA_plots[[2]],CTR_CA_plots[[5]]))+plot_annotation(title = "CTR + CA")&
    theme(plot.title=element_text(hjust = 0.5,size = 30))
  
  CA_VOR_plots <-   omics_similarity_data %>%
    .[["similarity_data"]] %>% 
    .[[2]] %>% 
    pull(plots)
  CA_VOR_plot <- wrap_plots(list(CA_VOR_plots[[2]],CA_VOR_plots[[5]]))+plot_annotation(title = "CA + VOR")&
    theme(plot.title=element_text(hjust = 0.5,size = 30))
  
  CA_WT_plots <-   omics_similarity_data %>%
    .[["similarity_data"]] %>% 
    .[[3]] %>% 
    pull(plots)
  CA_WT_plot <- wrap_plots(list(CA_WT_plots[[2]],CA_WT_plots[[5]]))+plot_annotation(title = "CA + WT")&
    theme(plot.title=element_text(hjust = 0.5,size = 30))
  
  CA_PptA_plots <-   omics_similarity_data %>%
    .[["similarity_data"]] %>% 
    .[[4]] %>% 
    pull(plots)
  CA_PptA_plot <- wrap_plots(list(CA_PptA_plots[[2]],CA_PptA_plots[[5]]))+plot_annotation(title = "CA + PptA")&
    theme(plot.title=element_text(hjust = 0.5,size = 30))
  
  WT_Tx_plots <-   omics_similarity_data %>%
    .[["similarity_data"]] %>% 
    .[[5]] %>% 
    pull(plots)
  WT_Tx_plot <- wrap_plots(list(WT_Tx_plots[[2]],WT_Tx_plots[[5]]))+plot_annotation(title = "WT + Tx")&
    theme(plot.title=element_text(hjust = 0.5,size = 30))
  
  ggarrange(plotlist = list(CTR_CA_plot,CA_VOR_plot,CA_WT_plot,CA_PptA_plot,WT_Tx_plot),ncol = 2,nrow = 3)
}