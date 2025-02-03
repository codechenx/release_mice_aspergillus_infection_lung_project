library(targets)
library(tarchetypes)
library(tidyverse)
library(crew)
source("R/utils.R")
source("R/get_data.R")
source("R/get_colormap.R")
source("R/fig5.R")
tar_option_set(
  packages = c(
    "ggpubr","ggplot2","ggforce","patchwork","ggsci","showtext","sjmisc","rstatix","ggbreak",
    "psych","phylosmith","phyloseq","vegan","tidyverse"),
  controller = crew_controller_local(workers = tar_config_get("workers"))
)


list(
  # output path
  tar_target(output_dir, "results/fig5"),
  tar_target(final_output_dir, "results/final_figures"),


  # data path
  tar_target(lung_decontam_genus_physeq_file, "processed_data/lung_decontam_genus_physeq.rds", format = "file"),
  tar_target(qpcr_file, "raw_data/qpcr.csv", format = "file"),

  # load data
  tar_target(lung_decontam_genus_physeq, readRDS(lung_decontam_genus_physeq_file)),
  tar_target(qpcr_tbl, get_qpcr_tbl(qpcr_file,data.frame(sample_data(lung_decontam_genus_physeq)))),

  # get data for fig5
  tar_target(obligate_aerobe_anaerobe_tbl, get_obligate_aerobe_anaerobe()),
  tar_target(fig5_lung_decontam_genus_physeq, get_fig5_physeq(lung_decontam_genus_physeq)),

  # plot qPCR of L. murinus and A. fumigatus
  tar_target(fig5a1, plot_qpcr_box(qpcr_tbl)),
  tar_target(fig5a2, plot_qpcr_cor(qpcr_tbl)),
  # plot obligate aerobic and anaerobic bacteria
  tar_target(fig5b, plot_lung_aerobe_anaerobe(fig5_lung_decontam_genus_physeq,obligate_aerobe_anaerobe_tbl)),

  # combine all figures
  ## fig5a
  tar_target(fig5a, fig5a1 + fig5a2),
  tar_target(export_fig5a, tar_ggsave(fig5a,file.path(output_dir, "fig5a.pdf"),width=10,height=7), format = "file"),
  ## fig5b
  tar_target(export_fig5b, tar_ggsave(fig5b,file.path(output_dir, "fig5b.pdf"),width=4,height=9), format = "file")
  ## final figures
)