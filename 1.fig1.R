library(targets)
library(tarchetypes)
library(tidyverse)
library(crew)
source("R/utils.R")
source("R/get_colormap.R")
source("R/fig1.R")
tar_option_set(
  packages = c(
    "ggpubr","ggplot2","patchwork","ggtext","showtext","microViz","phyloseq","phylosmith","vegan","tidyverse"),
  controller = crew_controller_local(workers = tar_config_get("workers"))
)

list(
  # output path
  tar_target(output_dir, "results/fig1"),
  tar_target(final_output_dir, "results/final_figures"),


  # data path
  tar_target(lung_decontam_genus_physeq_file, "processed_data/lung_decontam_genus_physeq.rds", format = "file"),


  # load data
  tar_target(lung_decontam_genus_physeq, readRDS(lung_decontam_genus_physeq_file)),


  # plot composition of bacteria in lung
  tar_target(lung_taxa_composition_plot, plot_taxa_composition(lung_decontam_genus_physeq)),


  # combine and export all figures
  tar_target(export_fig1a, tar_ggsave(lung_taxa_composition_plot,file.path(output_dir, "fig1a.pdf"),width=10,height=3), format = "file"),
  ## final figures
  tar_target(export_fig1, tar_ggsave(lung_taxa_composition_plot,file.path(final_output_dir, "fig1.pdf"),width=10,height=3), format = "file")


)
