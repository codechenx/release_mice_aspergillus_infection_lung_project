library(targets)
library(tarchetypes)
library(tidyverse)
library(crew)
source("R/utils.R")
source("R/get_data.R")
source("R/get_colormap.R")
source("R/diversity.R")
source("R/significant_taxa.R")
source("R/fig2.R")
tar_option_set(
  packages = c(
    "ggpubr","ggplot2","patchwork","microViz","circlize","ComplexHeatmap","ggtext","showtext","ropls","sjmisc","MatrixCorrelation",
    "rstatix","picante","DESeq2","data.table","microbiome","decontam","ape","phyloseq","phylosmith","vegan","tidyverse"),
  controller = crew_controller_local(workers = tar_config_get("workers"))
)

list(
  # output path
  tar_target(output_dir, "results/fig2"),
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
  tar_target(gut_decontam_genus_physeq, readRDS(gut_decontam_genus_physeq_file)),
  tar_target(lung_metabolome_physeq, readRDS(lung_metabolome_physeq_file)),
  tar_target(plasma_metabolome_physeq, readRDS(plasma_metabolome_physeq_file)),
  tar_target(lung_metabolome_sig_df, read_sig_metabolite(lung_metabolome_sig_file)),
  tar_target(plasma_metabolome_sig_df, read_sig_metabolite(plasma_metabolome_sig_file)),

  tar_target(fig2_comparison_df, get_fig2_comparison_df()),
  tar_target(lung_subset_decontam_genus_physeq, subset_phyloseq_by_comparison(lung_decontam_genus_physeq, G1 = fig2_comparison_df$G1, G2 = fig2_comparison_df$G2),
             pattern = fig2_comparison_df,
             iteration = "list"
  ),
  tar_target(gut_subset_decontam_genus_physeq, subset_phyloseq_by_comparison(gut_decontam_genus_physeq, G1 = fig2_comparison_df$G1, G2 = fig2_comparison_df$G2),
             pattern = fig2_comparison_df,
             iteration = "list"
  ),
  tar_target(lung_subset_metabolome_physeq, subset_phyloseq_by_comparison(lung_metabolome_physeq, G1 = fig2_comparison_df$G1, G2 = fig2_comparison_df$G2),
             pattern = fig2_comparison_df,
             iteration = "list"
  ),
  tar_target(plasma_subset_metabolome_physeq, subset_phyloseq_by_comparison(plasma_metabolome_physeq, G1 = fig2_comparison_df$G1, G2 = fig2_comparison_df$G2),
             pattern = fig2_comparison_df,
             iteration = "list"
  ),

  tar_target(fig2_lung_decontam_genus_physeq, get_fig2_physeq(lung_decontam_genus_physeq)),
  tar_target(fig2_gut_decontam_genus_physeq, get_fig2_physeq(gut_decontam_genus_physeq)),
  tar_target(fig2_lung_metabolome_physeq, get_fig2_physeq(lung_metabolome_physeq)),
  tar_target(fig2_plasma_metabolome_physeq, get_fig2_physeq(plasma_metabolome_physeq)),


  # alpha diversity
  tar_target(lung_alpha_plot, plot_fig2_alpha(fig2_lung_decontam_genus_physeq,plot_title="")),
  tar_target(gut_alpha_plot, plot_fig2_alpha(fig2_gut_decontam_genus_physeq,plot_title="")),


  # beta diversity
  tar_target(lung_beta_plot, plot_fig2_beta(fig2_lung_decontam_genus_physeq,plot_title="")),
  tar_target(gut_beta_plot, plot_fig2_beta(fig2_gut_decontam_genus_physeq,plot_title="")),


  # plsda
  tar_target(lung_plsda_plot, plot_fig2_plsda(fig2_lung_metabolome_physeq,plot_title="")),
  tar_target(plasma_plsda_plot, plot_fig2_plsda(fig2_plasma_metabolome_physeq,plot_title="")),
  tar_target(lung_metabolome_beta_plot,plot_fig2_beta(fig2_lung_metabolome_physeq,distance="euclidean",weighted = F,plot_title="Lung metabolome")),
  tar_target(plasma_metabolome_beta_plot,plot_fig2_beta(fig2_plasma_metabolome_physeq,distance="euclidean",weighted = F,plot_title="Plasma metabolome")),


  # significant taxa
  ### filter genus by abundance and prevalence
  tar_target(lung_subset_filtered_decontam_genus_physeq,utils_filter_taxa(lung_subset_decontam_genus_physeq,abundance_cutoff = 0.0001,prevalence_cutoff = 0.1),pattern = lung_subset_decontam_genus_physeq,iteration = "list"),
  tar_target(gut_subset_filtered_decontam_genus_physeq,utils_filter_taxa(gut_subset_decontam_genus_physeq,abundance_cutoff = 0.0001,prevalence_cutoff = 0.1),pattern = gut_subset_decontam_genus_physeq,iteration = "list"),
  ### DESeq2
  tar_target(lung_different_genus_deseq2,get_different_taxa_deseq2(lung_subset_filtered_decontam_genus_physeq),pattern = lung_subset_filtered_decontam_genus_physeq,iteration = "list"),
  tar_target(lung_combined_different_genus_deseq2,bind_rows(lung_different_genus_deseq2)),
  tar_target(export_lung_deseq2_results, tar_write_csv(lung_combined_different_genus_deseq2,file.path(output_dir, "lung_different_genus_deseq2.csv"))),
  tar_target(gut_different_genus_deseq2,get_different_taxa_deseq2(gut_subset_filtered_decontam_genus_physeq),pattern = gut_subset_filtered_decontam_genus_physeq,iteration = "list"),
  tar_target(gut_combined_different_genus_deseq2,bind_rows(gut_different_genus_deseq2)),
  tar_target(export_gut_deseq2_results, tar_write_csv(gut_combined_different_genus_deseq2,file.path(output_dir, "gut_different_genus_deseq2.csv"))),
  tar_target(lung_combined_different_genus_deseq2_plot_list, plot_stacked_sig_taxa_barplot(lung_combined_different_genus_deseq2)),
  tar_target(gut_combined_different_genus_deseq2_plot_list, plot_stacked_sig_taxa_barplot(gut_combined_different_genus_deseq2)),


  # significant metabolites
  tar_target(plasma_metabolome_sig_plot_list, plot_stacked_sig_metabolite_barplot(plasma_metabolome_sig_df)),
  tar_target(lung_metabolome_sig_plot_list, plot_stacked_sig_metabolite_barplot(lung_metabolome_sig_df)),


  # omics similarity
  tar_target(omics_data, get_omics_data(fig2_lung_decontam_genus_physeq, fig2_gut_decontam_genus_physeq, fig2_lung_metabolome_physeq,fig2_plasma_metabolome_physeq, transform="")),
  tar_target(omics_clr_data, get_omics_data(fig2_lung_decontam_genus_physeq, fig2_gut_decontam_genus_physeq, fig2_lung_metabolome_physeq,fig2_plasma_metabolome_physeq, transform="rclr")),
  tar_target(omics_similarity_data, get_omics_similarity(omics_data,omics_clr_data)),


  # correlation between microbiome and metabolome
  ## CTR and CA
  tar_target(CTR_CA_omics_data, get_omics_data(lung_decontam_genus_physeq %>% subset_samples(group %in% c("CTR","CA")),
                                               gut_decontam_genus_physeq %>% subset_samples(group %in% c("CTR","CA")),
                                               lung_metabolome_physeq %>% subset_samples(group %in% c("CTR","CA")),
                                               fig2_plasma_metabolome_physeq %>% subset_samples(group %in% c("CTR","CA")))),
  tar_target(metabolome_microbiome_correlation, get_bacteria_metabolite_correlation(CTR_CA_omics_data)),
  tar_target(exp_metabolome_microbiome_correlation, tar_write_csv(metabolome_microbiome_correlation,file.path(output_dir, "metabolome_microbiome_correlation.csv"))),
  
  # combine and export all figures
  tar_target(export_lung_metabolome_beta_plot, tar_ggsave(lung_metabolome_beta_plot, file.path(output_dir, "lung_metabolome_beta_diversity.pdf"))),
  tar_target(export_plasma_metabolome_beta_plot, tar_ggsave(plasma_metabolome_beta_plot, file.path(output_dir, "plasma_metabolome_beta_diversity.pdf"))),
  
  tar_target(fig2a, ggarrange(plotlist = list(lung_alpha_plot, lung_beta_plot,lung_combined_different_genus_deseq2_plot_list[[1]],lung_combined_different_genus_deseq2_plot_list[[2]]),nrow = 1,heights=1.5,align = "h")),
  tar_target(export_fig2a, tar_ggsave(fig2a,file.path(output_dir, "fig2a.pdf"),width=27,height=6), format = "file"),

  tar_target(fig2b, ggarrange(plotlist = list(gut_alpha_plot, gut_beta_plot,gut_combined_different_genus_deseq2_plot_list[[1]],gut_combined_different_genus_deseq2_plot_list[[2]]),nrow = 1,align = "h")),
  tar_target(export_fig2b, tar_ggsave(fig2b,file.path(output_dir, "fig2b.pdf"),width=27,height=6), format = "file"),

  tar_target(fig2c, ggarrange(plotlist = list(omics_similarity_data$plots[[2]], lung_plsda_plot,lung_metabolome_sig_plot_list[[1]],lung_metabolome_sig_plot_list[[2]]),nrow = 1,align = "h")),
  tar_target(export_fig2c, tar_ggsave(fig2c,file.path(output_dir, "fig2c.pdf"),width=27,height=6), format = "file"),

  tar_target(fig2d, ggarrange(plotlist = list(omics_similarity_data$plots[[5]], plasma_plsda_plot,plasma_metabolome_sig_plot_list[[1]],plasma_metabolome_sig_plot_list[[2]]),nrow = 1,align = "h")),
  tar_target(export_fig2d, tar_ggsave(fig2d,file.path(output_dir, "fig2d.pdf"),width=27,height=6), format = "file"),
  ## final figures
  tar_target(fig2, ggarrange(plotlist = list(
    ggarrange(lung_alpha_plot,gut_alpha_plot,omics_similarity_data$plots[[2]],omics_similarity_data$plots[[5]],nrow = 4,ncol=1),
    ggarrange(lung_beta_plot,gut_beta_plot,lung_plsda_plot,plasma_plsda_plot,nrow = 4,ncol=1,align = "v"),
    ggarrange(lung_combined_different_genus_deseq2_plot_list[[1]],gut_combined_different_genus_deseq2_plot_list[[1]],lung_metabolome_sig_plot_list[[1]],plasma_metabolome_sig_plot_list[[1]],nrow = 4,ncol=1,align = "v"),
    ggarrange(lung_combined_different_genus_deseq2_plot_list[[2]],gut_combined_different_genus_deseq2_plot_list[[2]],lung_metabolome_sig_plot_list[[2]],plasma_metabolome_sig_plot_list[[2]],nrow = 4,ncol=1,align = "v")
    ),
    nrow = 1,ncol=4,align = "h")),
  tar_target(export_fig2, tar_ggsave(fig2,file.path(final_output_dir, "fig2.pdf"),width=27,height=25), format = "file")
)
