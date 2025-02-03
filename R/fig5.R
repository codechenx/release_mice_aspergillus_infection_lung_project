
get_obligate_aerobe_anaerobe <- function() {
  tbl <- tibble(
    genus = c("Brevibacterium","Pseudoclavibacter","Acinetobacter","Sphingomonas",
    "Phascolarctobacterium","Romboutsia","Anaerococcus","Olsenella"),
    oxygen = c("aerobe","aerobe","aerobe","aerobe",
    "anaerobe","anaerobe","anaerobe","anaerobe")
  )
  return(tbl)
}


get_fig5_physeq <- function(physeq) {
  filtered_physeq <- physeq %>%
    subset_samples(sampleGroup %in% c(2,4,5)) %>%
    filter_taxa(function(x)
      sum(x > 0) > 0, TRUE)
  return(filtered_physeq)
}


get_qpcr_tbl <- function(file_path,meta_tbl) {
  qpcr_tbl <- read_csv(file_path) %>%
    left_join(meta_tbl %>% mutate(sampleName=as.numeric(sampleName)), by = c("id" = "sampleName"))
  return(qpcr_tbl)
}


plot_qpcr_box <- function(qpcr_tbl) {
  plot_tbl <- qpcr_tbl %>%
    filter(group %in% c("CA","WT","Tx","PptA")) %>%
    dplyr::select(group,L_murinus,A_fumigatus) %>%
    gather(key="species",value="abundance",L_murinus,A_fumigatus) %>%
    mutate(species=str_replace(species,"_",". "))
  color <- get_colormap_and_charmap() %>%
    filter(group %in% plot_tbl$group) %>% pull(color)
  qpcr_box <- ggplot(plot_tbl, aes(x=group, y=abundance,fill=group)) +
    geom_boxplot() +
    facet_wrap(~species)+
    stat_compare_means(comparisons=list(c("CA","WT"),c("CA","Tx"),c("WT","Tx"),c("CA","PptA")),method = "wilcox.test",label = "p.format", size = 5)+
    theme_pubr()+
    theme(text = element_text(size=18),
          strip.text = element_text(size=20,face = "italic"),
          strip.background = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          legend.key.size = unit(.5, 'cm'),
          legend.spacing.x = unit(.4, 'cm'),
          legend.text = element_text(size=18),
          axis.title = element_text(size = 18),
          axis.text=element_text(size=16))+
    ylab("Abundance (pseudo-log)") +
    xlab("")+
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks =c(1,10,100,500,1000,2000,3000,4000))+
    scale_fill_manual(values=color)
  return(qpcr_box)
}


plot_qpcr_cor <- function(qpcr_tbl) {
  plot_tbl <- qpcr_tbl %>%
    filter(group %in% c("WT","Tx")) %>%
    filter(L_murinus!=0) %>%
    mutate(group=fct_drop(group)) %>%
    dplyr::select(group,L_murinus,A_fumigatus) %>%
    mutate(L_murinus=rank(L_murinus),A_fumigatus=rank(A_fumigatus))
  color <- get_colormap_and_charmap() %>%
    filter(group %in% plot_tbl$group) %>% pull(color)
  qpcr_cor_plot <- ggscatter(plot_tbl, x = "A_fumigatus", y = "L_murinus",scales = "free",
    add = "reg.line", conf.int = F) +
    stat_cor(method = "spearman",size=5,p.accuracy = 0.001, r.accuracy = 0.001,label.y.npc=0.9) +
    xlab(expression("Rank ("~italic("A. fumigatus")~")")) +
    ylab(expression("Rank ("~italic("L. murinus")~")")) +
    theme(text = element_text(size=18),
          strip.text = element_text(size=20),
          strip.background = element_rect(color="black", fill="white"),
          legend.title = element_blank(),
          legend.position = "top",
          legend.key.size = unit(.5, 'cm'),
          legend.spacing.x = unit(.4, 'cm'),
          legend.text = element_text(size=18),
          axis.title = element_text(size = 18),
          axis.text=element_text(size=16))+
    scale_color_manual(values=color)
  return(qpcr_cor_plot)
}


plot_lung_aerobe_anaerobe <- function(physeq, tbl,plot_lab="") {
  if(nrow(tbl)<1){
    return(ggplot())
  }
  physeq<- transform_sample_counts(physeq, function(OTU) OTU / sum(OTU)*100)
  otu_long_tbl <-  melt_phyloseq(physeq) %>%
    filter(Genus %in% tbl$genus) %>%
    mutate(oxygen = tbl$oxygen[match(Genus, tbl$genus)]) %>%
    dplyr::select(Abundance, oxygen, group, Genus, Sample)
  
  aerobe_anaerobe_tbl <- otu_long_tbl %>%
    group_by(group, oxygen,Sample) %>%
    summarize(Abundance = sum(Abundance)) 

  mean_ratio_tbl <- otu_long_tbl %>%
    dplyr::select(Abundance, oxygen, group, Genus) %>%
    group_by(group, oxygen) %>%
    summarize(Abundance = sum(Abundance)) %>%
    ungroup() %>%
    spread(oxygen, Abundance) %>%
    mutate(ratio = anaerobe/aerobe)

  ratio_tbl <- otu_long_tbl %>%
    group_by(group, oxygen,Sample) %>%
    summarize(Abundance = sum(Abundance)) %>%
    ungroup() %>%
    spread(oxygen, Abundance) %>%
    mutate(ratio = anaerobe/aerobe) %>%
    filter(!is.nan(ratio)) %>%
    mutate(ratio=ifelse(is.infinite(ratio),max(ratio[is.finite(ratio)]),ratio))

  color <- get_colormap_and_charmap() %>%
    filter(group %in% otu_long_tbl$group) %>% pull(color)

  ratio_barplot <- ggplot(ratio_tbl, aes(x=group, y=ratio, fill=group)) +
    geom_boxplot() +
    theme_pubr()+
    stat_compare_means(method = "wilcox.test",comparisons=list(c("CA","WT"),c("CA","VOR"),c("WT","VOR")),label = "p.format", size = 4)+
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          text = element_text(size=18),
          strip.text = element_text(size=20),
          legend.title = element_blank(),
          legend.position = "top",
          legend.key.size = unit(.5, 'cm'),
          legend.spacing.x = unit(.4, 'cm'),
          legend.text = element_text(size=18),
          axis.title = element_text(size = 18),
          axis.title.x=element_blank(),
          axis.line.x = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text=element_text(size=18))+
    labs(title = plot_lab) + ylab("Anaerobe/aerobe ratio") +
    scale_fill_manual(values=color)

  mean_ratio_barplot <- ggplot(mean_ratio_tbl, aes(x=group, y=ratio, fill=group)) +
    geom_bar(stat="identity") +
    theme_pubr()+
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          text = element_text(size=18),
          strip.text = element_text(size=20),
          legend.title = element_blank(),
          legend.position = "top",
          legend.key.size = unit(.5, 'cm'),
          legend.spacing.x = unit(.4, 'cm'),
          legend.text = element_text(size=18),
          axis.title = element_text(size = 18),
          axis.title.x=element_blank(),
          axis.line.x = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text=element_text(size=18))+
    labs(title = plot_lab) + ylab("Anaerobe/aerobe ratio") +
    scale_fill_manual(values=color)
abundance_jitter <- ggplot(aerobe_anaerobe_tbl, aes(x=group, y=Abundance, color=oxygen)) +
  geom_jitter(size=2)  +
  theme_pubr()+
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        text = element_text(size=18),
        strip.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(.8, 'cm'),
        legend.spacing.x = unit(.4, 'cm'),
        legend.text = element_text(size=18),
        axis.title = element_text(size = 18),
        axis.text.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.x=element_blank(),
        axis.line.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text=element_text(size=18))+
  labs(title = plot_lab) + ylab("Abundance (%)")+  xlab("") +
  scale_y_continuous(
    trans = trans_reverser("pseudo_log"),
    breaks = c(0, 1, 10, 80)
  )+
  scale_color_manual(values=c('#1b7c3d','#f16c23'))
mean_ratio_barplot+abundance_jitter + plot_layout(ncol=1, heights=c(1,1.5))
}