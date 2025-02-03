library(targets)
library(tarchetypes)
library(tidyverse)
library(crew)
source("R/utils.R")
source("R/get_data.R")
source("R/get_colormap.R")
source("R/diversity.R")
source("R/significant_taxa.R")
source("R/run_feast.R")
source('R/fig4.R')
tar_option_set(
  packages = c(
    "ggpubr","ggplot2","ggsci","sjmisc", "ggtext","ggrepel","scales","rstatix","circlize","usedist",
    "ComplexHeatmap","data.table","showtext","ape","phylosmith","phyloseq","vegan","ropls","plyr","tidyverse"),
    controller = crew_controller_local(workers = tar_config_get("workers"))

)

list(
  # output path
  tar_target(output_dir, "results/fig4"),
  tar_target(final_output_dir, "results/final_figures"),

  # data path
  tar_target(lung_decontam_genus_physeq_file, "processed_data/lung_decontam_genus_physeq.rds", format = "file"),
  tar_target(gut_decontam_genus_physeq_file, "processed_data/gut_decontam_genus_physeq.rds", format = "file"),
  tar_target(lung_metabolome_physeq_file, "processed_data/lung_metabolome_physeq.rds", format = "file"),
  tar_target(plasma_metabolome_physeq_file, "processed_data/plasma_metabolome_physeq.rds", format = "file"),
  tar_target(lung_metabolome_sig_file, "raw_data/lung_dif_0.05_fdr_log.csv", format = "file"),
  tar_target(plasma_metabolome_sig_file, "raw_data/plasma_dif_0.05_fdr_log.csv", format = "file"),


  # load data
  tar_target(lung_decontam_genus_physeq, readRDS(lung_decontam_genus_physeq_file)),
  tar_target(lung_decontam_genus_rl_physeq, transform_sample_counts(lung_decontam_genus_physeq, function(OTU)OTU / sum(OTU))),
  tar_target(gut_decontam_genus_physeq, readRDS(gut_decontam_genus_physeq_file)),
  tar_target(gut_decontam_genus_rl_physeq, transform_sample_counts(gut_decontam_genus_physeq, function(OTU)OTU / sum(OTU))),
  tar_target(lung_metabolome_physeq, readRDS(lung_metabolome_physeq_file)),
  tar_target(plasma_metabolome_physeq, readRDS(plasma_metabolome_physeq_file)),
  tar_target(lung_metabolome_sig_df, read_sig_metabolite(lung_metabolome_sig_file)),
  tar_target(plasma_metabolome_sig_df, read_sig_metabolite(plasma_metabolome_sig_file)),
  ## load deseq2 result from fig3(infection)
  tar_target(fig3_lung_combined_different_genus_deseq2, targets::tar_read(lung_combined_different_genus_deseq2, store = "cache/3.fig3")),
  tar_target(fig3_gut_combined_different_genus_deseq2, targets::tar_read(gut_combined_different_genus_deseq2, store = "cache/3.fig3")),
  ## get phyloseq objects for deseq2 of fig4
  tar_target(fig4_comparison_df, get_fig4_comparison_df()),
  tar_target(lung_subset_decontam_genus_physeq, subset_phyloseq_by_comparison(lung_decontam_genus_physeq, G1 = fig4_comparison_df$G1, G2 = fig4_comparison_df$G2),
             pattern = fig4_comparison_df,
             iteration = "list"
  ),
  tar_target(gut_subset_decontam_genus_physeq, subset_phyloseq_by_comparison(gut_decontam_genus_physeq, G1 = fig4_comparison_df$G1, G2 = fig4_comparison_df$G2),
             pattern = fig4_comparison_df,
             iteration = "list"
  ),
  tar_target(lung_subset_metabolome_physeq, subset_phyloseq_by_comparison(lung_metabolome_physeq, G1 = fig4_comparison_df$G1, G2 = fig4_comparison_df$G2),
             pattern = fig4_comparison_df,
             iteration = "list"
  ),
  tar_target(plasma_subset_metabolome_physeq, subset_phyloseq_by_comparison(plasma_metabolome_physeq, G1 = fig4_comparison_df$G1, G2 = fig4_comparison_df$G2),
             pattern = fig4_comparison_df,
             iteration = "list"
  ),
  ## get fig4 phyloseq objects
  tar_target(fig4_lung_decontam_genus_physeq, get_fig4_physeq(lung_decontam_genus_physeq)),
  tar_target(fig4_gut_decontam_genus_physeq, get_fig4_physeq(gut_decontam_genus_physeq)),
  tar_target(fig4_lung_metabolome_physeq, get_fig4_physeq(lung_metabolome_physeq)),
  tar_target(fig4_plasma_metabolome_physeq, get_fig4_physeq(plasma_metabolome_physeq)),
  ## get extended fig4 phyloseq objects including PBS_CA(CA), AF_wt_CA(WT), and AF_wt_CA_Vor(Tx)
  tar_target(extended_fig4_lung_decontam_genus_physeq, get_extended_fig4_physeq(lung_decontam_genus_physeq)),
  tar_target(extended_fig4_gut_decontam_genus_physeq, get_extended_fig4_physeq(gut_decontam_genus_physeq)),
  tar_target(extended_fig4_lung_metabolome_physeq, get_extended_fig4_physeq(lung_metabolome_physeq)),
  tar_target(extended_fig4_plasma_metabolome_physeq, get_extended_fig4_physeq(plasma_metabolome_physeq)),

  
  # alpha diversity
  tar_target(lung_alpha_plot, plot_fig4_alpha(fig4_lung_decontam_genus_physeq,plot_title="")),
  tar_target(gut_alpha_plot, plot_fig4_alpha(fig4_gut_decontam_genus_physeq,plot_title="")),
  
  
  # the distance between Vor vs Tx and WT vs Tx diversity
  tar_target(lung_beta_plot, plot_fig4_microbiome_distance(fig4_lung_decontam_genus_physeq,plot_title="Lung microbiome")),
  tar_target(gut_beta_plot, plot_fig4_microbiome_distance(fig4_gut_decontam_genus_physeq,plot_title="Gut microbiome")),
  tar_target(lung_plsda_plot, plot_fig4_metabolome_distance(fig4_lung_metabolome_physeq,plot_title="Lung metabolome")),
  tar_target(plasma_plsda_plot, plot_fig4_metabolome_distance(fig4_plasma_metabolome_physeq,plot_title="Plasma metabolome")),


  # significant taxa
  ## filter genus by abundance and prevalence
  tar_target(fig4_lung_filtered_decontam_genus_physeq,utils_filter_taxa(fig4_lung_decontam_genus_physeq,abundance_cutoff = 0.0001,prevalence_cutoff = 0.1)),
  tar_target(fig4_gut_filtered_decontam_genus_physeq,utils_filter_taxa(fig4_gut_decontam_genus_physeq,abundance_cutoff = 0.0001,prevalence_cutoff = 0.1)),
  tar_target(lung_subset_filtered_decontam_genus_physeq, subset_phyloseq_by_comparison(fig4_lung_filtered_decontam_genus_physeq, G1 = fig4_comparison_df$G1, G2 = fig4_comparison_df$G2),
             pattern = fig4_comparison_df,
             iteration = "list"
  ),
  tar_target(gut_subset_filtered_decontam_genus_physeq, subset_phyloseq_by_comparison(fig4_gut_filtered_decontam_genus_physeq, G1 = fig4_comparison_df$G1, G2 = fig4_comparison_df$G2),
             pattern = fig4_comparison_df,
             iteration = "list"
  ),
  ## DESeq2
  tar_target(lung_different_genus_deseq2,get_different_taxa_deseq2(lung_subset_filtered_decontam_genus_physeq),pattern = lung_subset_filtered_decontam_genus_physeq,iteration = "list"),
  tar_target(lung_combined_different_genus_deseq2,bind_rows(lung_different_genus_deseq2)),
  tar_target(export_lung_deseq2_results, tar_write_csv(lung_combined_different_genus_deseq2,file.path(output_dir, "lung_different_genus_deseq2.csv"))),
  tar_target(gut_different_genus_deseq2,get_different_taxa_deseq2(gut_subset_filtered_decontam_genus_physeq),pattern = gut_subset_filtered_decontam_genus_physeq,iteration = "list"),
  tar_target(gut_combined_different_genus_deseq2,bind_rows(gut_different_genus_deseq2)),
  tar_target(export_gut_deseq2_results, tar_write_csv(gut_combined_different_genus_deseq2,file.path(output_dir, "gut_different_genus_deseq2.csv"))),
  ## volcano plot
  tar_target(microbiome_volcano, plot_fig4_microbiome_volcano(lung_combined_different_genus_deseq2, gut_combined_different_genus_deseq2)),
  tar_target(metabolome_volcano, plot_fig4_metabolome_volcano(lung_metabolome_sig_df, plasma_metabolome_sig_df)),
  ## venn plot
  tar_target(microbiome_venn, plot_fig4_microbiome_venn(bind_rows(lung_combined_different_genus_deseq2, fig3_lung_combined_different_genus_deseq2),bind_rows(gut_combined_different_genus_deseq2, fig3_gut_combined_different_genus_deseq2))),
  tar_target(metabolome_venn, plot_fig4_metabolome_venn(lung_metabolome_sig_df, plasma_metabolome_sig_df)),


  # constant changing taxa(fig4c)
  tar_target(lung_all_different_microbiome_deseq2,bind_rows(fig3_lung_combined_different_genus_deseq2,lung_combined_different_genus_deseq2)),
  tar_target(gut_all_different_microbiome_deseq2,bind_rows(fig3_gut_combined_different_genus_deseq2,gut_combined_different_genus_deseq2)),
  tar_target(lung_all_different_metabolome_deseq2,lung_metabolome_sig_df),
  tar_target(plasma_all_different_metabolome_deseq2,plasma_metabolome_sig_df),

  tar_target(lung_microbiome_constant_changing_taxa, get_constant_taxa(lung_all_different_microbiome_deseq2)),
  tar_target(export_lung_microbiome_constant_changing_taxa, tar_write_csv(lung_microbiome_constant_changing_taxa,file.path(output_dir, "lung_microbiome_constant_changing_taxa.csv"))),
  tar_target(gut_microbiome_constant_changing_taxa, get_constant_taxa(gut_all_different_microbiome_deseq2)),
  tar_target(export_gut_microbiome_constant_changing_taxa, tar_write_csv(gut_microbiome_constant_changing_taxa,file.path(output_dir, "gut_microbiome_constant_changing_taxa.csv"))),
  tar_target(lung_metabolome_constant_changing_taxa, get_constant_taxa(lung_all_different_metabolome_deseq2)),
  tar_target(export_lung_metabolome_constant_changing_taxa, tar_write_csv(lung_metabolome_constant_changing_taxa,file.path(output_dir, "lung_metabolome_constant_changing_taxa.csv"))),
  tar_target(plasma_metabolome_constant_changing_taxa, get_constant_taxa(plasma_all_different_metabolome_deseq2)),
  tar_target(export_plasma_all_different_metabolome_deseq2, tar_write_csv(plasma_metabolome_constant_changing_taxa,file.path(output_dir, "plasma_metabolome_constant_changing_taxa.csv"))),


  tar_target(lung_microbiome_constant_changing_plot, plot_constant_taxa(lung_decontam_genus_rl_physeq, lung_microbiome_constant_changing_taxa, plot_lab="Lung microbiome")),
  tar_target(gut_microbiome_constant_changing_plot, plot_constant_taxa(gut_decontam_genus_rl_physeq, gut_microbiome_constant_changing_taxa, plot_lab="Gut microbiome")),
  tar_target(lung_metabolome_constant_changing_plot, plot_constant_taxa(lung_metabolome_physeq, lung_metabolome_constant_changing_taxa, plot_lab="Lung metabolome")),
  tar_target(plasma_metabolome_constant_changing_plot, plot_constant_taxa(plasma_metabolome_physeq, plasma_metabolome_constant_changing_taxa,top_n=6, plot_lab="Plasma metabolome")),


  # create a integrated network for constant changing genus and metabolites
  ## all
  tar_target(constant_changing_taxa_cor, get_constant_taxa_cor(tibble(omics=c("lung_microbiome","gut_microbiome","lung_metabolome","plasma_metabolome"),
                                                                                 physeq=list(microbiome::transform(extended_fig4_lung_decontam_genus_physeq,"rclr"),
                                                                                             microbiome::transform(extended_fig4_gut_decontam_genus_physeq,"rclr"),
                                                                                             extended_fig4_lung_metabolome_physeq,
                                                                                             extended_fig4_plasma_metabolome_physeq),
                                                                                 taxa=list(lung_microbiome_constant_changing_taxa,
                                                                                           gut_microbiome_constant_changing_taxa,
                                                                                           lung_metabolome_constant_changing_taxa,
                                                                                           plasma_metabolome_constant_changing_taxa)))),
  tar_target(constant_changing_taxa_network, plot_constant_taxa_cor_network(constant_changing_taxa_cor)),
  ## PBS_CA(CA)
  tar_target(pbs_ca_constant_changing_taxa_cor, get_constant_taxa_cor(tibble(omics=c("lung_microbiome","gut_microbiome","lung_metabolome","plasma_metabolome"),
                                                                          physeq=list(microbiome::transform(extended_fig4_lung_decontam_genus_physeq %>% subset_samples(group=="CA"),"rclr"),
                                                                                      microbiome::transform(extended_fig4_gut_decontam_genus_physeq %>% subset_samples(group=="CA"),"rclr"),
                                                                                      extended_fig4_lung_metabolome_physeq %>% subset_samples(group=="CA"),
                                                                                      extended_fig4_plasma_metabolome_physeq %>% subset_samples(group=="CA")),
                                                                          taxa=list(lung_microbiome_constant_changing_taxa,
                                                                                    gut_microbiome_constant_changing_taxa,
                                                                                    lung_metabolome_constant_changing_taxa,
                                                                                    plasma_metabolome_constant_changing_taxa)))),
  tar_target(pbs_ca_constant_changing_taxa_network, plot_constant_taxa_cor_network(pbs_ca_constant_changing_taxa_cor,"CA")),
  ## pbs_ca_vor(VOR)
  tar_target(pbs_ca_vor_constant_changing_taxa_cor, get_constant_taxa_cor(tibble(omics=c("lung_microbiome","gut_microbiome","lung_metabolome","plasma_metabolome"),
                                                                             physeq=list(microbiome::transform(extended_fig4_lung_decontam_genus_physeq %>% subset_samples(group=="VOR"),"rclr"),
                                                                                         microbiome::transform(extended_fig4_gut_decontam_genus_physeq %>% subset_samples(group=="VOR"),"rclr"),
                                                                                         extended_fig4_lung_metabolome_physeq %>% subset_samples(group=="VOR"),
                                                                                         extended_fig4_plasma_metabolome_physeq %>% subset_samples(group=="VOR")),
                                                                             taxa=list(lung_microbiome_constant_changing_taxa,
                                                                                       gut_microbiome_constant_changing_taxa,
                                                                                       lung_metabolome_constant_changing_taxa,
                                                                                       plasma_metabolome_constant_changing_taxa)))),
  tar_target(pbs_ca_vor_constant_changing_taxa_network, plot_constant_taxa_cor_network(af_wt_ca_constant_changing_taxa_cor,"VOR")),
  ## AF_wt_CA(WT)
  tar_target(af_wt_ca_constant_changing_taxa_cor, get_constant_taxa_cor(tibble(omics=c("lung_microbiome","gut_microbiome","lung_metabolome","plasma_metabolome"),
                                                                             physeq=list(microbiome::transform(extended_fig4_lung_decontam_genus_physeq %>% subset_samples(group=="WT"),"rclr"),
                                                                                         microbiome::transform(extended_fig4_gut_decontam_genus_physeq %>% subset_samples(group=="WT"),"rclr"),
                                                                                         extended_fig4_lung_metabolome_physeq %>% subset_samples(group=="WT"),
                                                                                         extended_fig4_plasma_metabolome_physeq %>% subset_samples(group=="WT")),
                                                                             taxa=list(lung_microbiome_constant_changing_taxa,
                                                                                       gut_microbiome_constant_changing_taxa,
                                                                                       lung_metabolome_constant_changing_taxa,
                                                                                       plasma_metabolome_constant_changing_taxa)))),
  tar_target(af_wt_ca_constant_changing_taxa_network, plot_constant_taxa_cor_network(af_wt_ca_constant_changing_taxa_cor,"WT")),
  ## AF_wt_CA_Vor(Tx)
  tar_target(af_wt_ca_vor_constant_changing_taxa_cor, get_constant_taxa_cor(tibble(omics=c("lung_microbiome","gut_microbiome","lung_metabolome","plasma_metabolome"),
                                                                               physeq=list(microbiome::transform(extended_fig4_lung_decontam_genus_physeq %>% subset_samples(group=="Tx"),"rclr"),
                                                                                           microbiome::transform(extended_fig4_gut_decontam_genus_physeq %>% subset_samples(group=="Tx"),"rclr"),
                                                                                           extended_fig4_lung_metabolome_physeq %>% subset_samples(group=="Tx"),
                                                                                           extended_fig4_plasma_metabolome_physeq %>% subset_samples(group=="Tx")),
                                                                               taxa=list(lung_microbiome_constant_changing_taxa,
                                                                                         gut_microbiome_constant_changing_taxa,
                                                                                         lung_metabolome_constant_changing_taxa,
                                                                                         plasma_metabolome_constant_changing_taxa)))),
  tar_target(af_wt_ca_vor_constant_changing_taxa_network, plot_constant_taxa_cor_network(af_wt_ca_vor_constant_changing_taxa_cor,"Tx")),

  # feast(supplementary)
  ## microbiome
  tar_target(lung_genus_feat_result, get_feast(lung_decontam_genus_physeq)),
  tar_target(gut_genus_feat_result, get_feast(gut_decontam_genus_physeq)),
  tar_target(lung_genus_feat_plot, plot_microbiome_feast(lung_genus_feat_result)),
  tar_target(gut_genus_feat_plot, plot_microbiome_feast(gut_genus_feat_result)),
  ## metabolome
  tar_target(lung_metabolome_dist_result, get_inter_group_beta_dist(lung_metabolome_physeq, index = "euclidean")),
  tar_target(plasma_metabolome_dist_result, get_inter_group_beta_dist(plasma_metabolome_physeq, index = "euclidean")),
  tar_target(lung_metabolome_intergroup_plot, plot_metabolome_intergroup_dist(lung_metabolome_dist_result, "Lung metabolome")),
  tar_target(plasma_metabolome_intergroup_plot, plot_metabolome_intergroup_dist(plasma_metabolome_dist_result, "Plasma metabolome")),

  # combine all figures
  tar_target(export_lung_alpha_plot, tar_ggsave(lung_alpha_plot, file.path(output_dir, "lung_alpha_diversity.pdf"))),
  tar_target(export_gut_alpha_plot, tar_ggsave(gut_alpha_plot, file.path(output_dir, "gut_alpha_diversity.pdf"))),
  ## fig4a
  tar_target(fig4a, ggarrange(plotlist = list(lung_beta_plot, gut_beta_plot,lung_plsda_plot,plasma_plsda_plot),nrow = 1,align = "hv", common.legend = T, legend="top")),
  tar_target(export_fig4a, tar_ggsave(fig4a,file.path(output_dir, "fig4a.pdf"),width=25,height=6), format = "file"),
  ## fig4b
  tar_target(fig4b, ggarrange(plotlist = c(microbiome_venn,metabolome_venn),ncol = 2,nrow = 2,align = "hv")),
  tar_target(export_fig4b, tar_ggsave(fig4b,file.path(output_dir, "fig4b.pdf"),width=13,height=13), format = "file"),
  ## fig4c
  tar_target(fig4c, ggarrange(plotlist = list(lung_microbiome_constant_changing_plot, gut_microbiome_constant_changing_plot,lung_metabolome_constant_changing_plot,plasma_metabolome_constant_changing_plot),ncol = 2,nrow = 2,align = "hv", common.legend = T, legend="top")),
  tar_target(export_fig4c, tar_ggsave(fig4c,file.path(output_dir, "fig4c.pdf"),width=25,height=14), format = "file"),
  ## fig4d
  tar_target(fig4d, ggarrange(plotlist = list(pbs_ca_constant_changing_taxa_network,af_wt_ca_constant_changing_taxa_network,af_wt_ca_vor_constant_changing_taxa_network),ncol = 3,nrow = 1,align = "hv", common.legend = T, legend = "top")),
  tar_target(export_fig4d, tar_ggsave(fig4d,file.path(output_dir, "fig4d.pdf"), width=10,height=4), format = "file"),
  ## fig4_s1
  tar_target(fig4_s1, ggarrange(plotlist = list(lung_genus_feat_plot, gut_genus_feat_plot,lung_metabolome_intergroup_plot,plasma_metabolome_intergroup_plot),ncol = 2,nrow = 2,align = "hv", common.legend = T)),
  tar_target(export_fig4_s1, tar_ggsave(fig4_s1,file.path(output_dir, "fig4_s1.pdf"),width=12,height=10), format = "file"),
  ## fig4_s2
  tar_target(fig4_s2, ggarrange(plotlist = c(microbiome_volcano,metabolome_volcano),nrow = 1,align = "hv", common.legend = T)),
  tar_target(export_fig4_s2, tar_ggsave(fig4_s2,file.path(output_dir, "fig4_s2.pdf"),width=40,height=6), format = "file"),
  ## final figures
  tar_target(fig4, ggarrange(plotlist = list(fig4a,fig4b,fig4c,fig4d),nrow = 4,ncol=1,align = "hv")),
  tar_target(export_fig4, tar_ggsave(fig4,file.path(final_output_dir, "fig4.pdf"),width=30,height=40), format = "file")
)