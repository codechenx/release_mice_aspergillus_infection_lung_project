library(targets)
library(tarchetypes)
library(tidyverse)
library(crew)
source("R/utils.R")
source("R/get_data.R")
source("R/get_colormap.R")
source("R/supplementary.R")
tar_option_set(
  packages = c(
    "ggpubr","ggplot2","sjmisc","rstatix","data.table","microViz","patchwork","phylosmith","ggtext","showtext","phyloseq","vegan","tidyverse"),
  controller = crew_controller_local(workers = tar_config_get("workers"))
)

list(
  # output path
  tar_target(output_dir, "results/supplementary"),
  tar_target(final_output_dir, "results/final_figures"),


  # data path
  tar_target(gut_decontam_genus_physeq_file, "processed_data/gut_decontam_genus_physeq.rds", format = "file"),
  tar_target(gut_baseline_decontam_genus_physeq_file, "processed_data/gut_baseline_decontam_genus_physeq.rds", format = "file"),
  tar_target(lung_nanopore_decontam_genus_physeq_file, "processed_data/lung_nanopore_decontam_genus_physeq.rds", format = "file"),
  
  
  # load data
  tar_target(gut_baseline_decontam_genus_physeq, readRDS(gut_baseline_decontam_genus_physeq_file)),
  tar_target(gut_decontam_genus_physeq, readRDS(gut_decontam_genus_physeq_file)),
  tar_target(lung_nanopore_physeq, readRDS(lung_nanopore_decontam_genus_physeq_file)),
  
  
  # beta diversity
  tar_target(gut_baseline_beta_plot, get_gut_baseline_beta_plot(gut_baseline_decontam_genus_physeq,plot_title="gut microbiome")),
  
  
  # plot composition of bacteria in gut
  tar_target(gut_taxa_composition_plot, plot_taxa_composition(gut_decontam_genus_physeq)),
  
  # plot composition of bacteria in lung nanopore
  tar_target(lung_nanopore_taxa_composition_plot, plot_taxa_composition(lung_nanopore_physeq,taxrank="Species",n_taxa=40)),
  
  # combine all figures
  ## supplementary figures
  tar_target(export_gut_baseline_beta_plot, tar_ggsave(gut_baseline_beta_plot,file.path(output_dir, "gut_baseline_beta.pdf"),width=10,height=10), format = "file"),
  ## final figures
  tar_target(export_fig5S1, tar_ggsave(gut_baseline_beta_plot,file.path(final_output_dir, "s1.pdf"),width=10,height=10), format = "file"),
  tar_target(export_fig5S2, tar_ggsave(gut_taxa_composition_plot,file.path(final_output_dir, "s2.pdf"),width=10,height=3), format = "file"),
  tar_target(export_fig5S3, tar_ggsave(lung_nanopore_taxa_composition_plot,file.path(final_output_dir, "s3.pdf"),width=15,height=3), format = "file")
)
