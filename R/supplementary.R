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