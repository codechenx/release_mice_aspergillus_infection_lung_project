plot_taxa_composition <- function(physeq) {
n_taxa=15
physeq <- physeq %>% transform_sample_counts(function(x) x/sum(x))
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
customPal <- tax_palette(
  data = physeq, rank = "Genus", pal = "kelly", n = n_taxa, add = c(Other = "white")
)
customPal[1] <- "royalblue"
customPal[2] <- "#f3c300"
customPal[3] <- "aquamarine4"
customPal["*Bacteria* Domain"] <- "seagreen2"
plot_all <- physeq  %>%
  subset_taxa(!is.na(Domain)) %>%
  set_sample_order("group") %>%
  comp_barplot(n_taxa = n_taxa, tax_level = "Genus", group_by = "group",merge_other = F,palette = customPal)
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
  )&ylab("Relative abundance")
return(composition_plot)
}