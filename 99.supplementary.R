library(targets)
library(tarchetypes)
library(tidyverse)
library(crew)
source("R/utils.R")
source("R/get_data.R")
source("R/get_colormap.R")
source("R/diversity.R")
source("R/fig2.R")
source("R/supplementary.R")
tar_option_set(
  packages = c(
    "ggpubr","ggplot2","sjmisc","rstatix","MatrixCorrelation","data.table","microViz","patchwork","phylosmith","ggtext","showtext","phyloseq","vegan","tidyverse"),
  controller = crew_controller_local(workers = tar_config_get("workers"))
)

list(
  # output path
  tar_target(output_dir, "results/supplementary"),
  tar_target(final_output_dir, "results/final_figures"),


  # data path
  tar_target(lung_decontam_genus_physeq_file, "processed_data/lung_decontam_genus_physeq.rds", format = "file"),
  tar_target(gut_decontam_genus_physeq_file, "processed_data/gut_decontam_genus_physeq.rds", format = "file"),
  tar_target(gut_baseline_decontam_genus_physeq_file, "processed_data/gut_baseline_decontam_genus_physeq.rds", format = "file"),
  tar_target(lung_nanopore_decontam_genus_physeq_file, "processed_data/lung_nanopore_decontam_genus_physeq.rds", format = "file"),
  tar_target(lung_metabolome_physeq_file, "processed_data/lung_metabolome_physeq.rds", format = "file"),
  tar_target(plasma_metabolome_physeq_file, "processed_data/plasma_metabolome_physeq.rds", format = "file"),
  
  # load data
  tar_target(lung_decontam_genus_physeq, readRDS(lung_decontam_genus_physeq_file)),
  tar_target(gut_baseline_decontam_genus_physeq, readRDS(gut_baseline_decontam_genus_physeq_file)),
  tar_target(gut_decontam_genus_physeq, readRDS(gut_decontam_genus_physeq_file)),
  tar_target(lung_nanopore_physeq, readRDS(lung_nanopore_decontam_genus_physeq_file)),
  tar_target(lung_metabolome_physeq, readRDS(lung_metabolome_physeq_file)),
  tar_target(plasma_metabolome_physeq, readRDS(plasma_metabolome_physeq_file)),
  
  # alpha diversity of each group
  tar_target(lung_alpha_plot, plot_alpha(lung_decontam_genus_physeq,plot_title="lung microbiome")),
  tar_target(gut_alpha_plot, plot_alpha(gut_decontam_genus_physeq,plot_title="gut microbiome")),
  
  
  # omics similarity of each group
  tar_target(omics_data,tibble(subset_group = c("CTR","CA","VOR","WT","Tx","PptA")) %>%
    mutate(subset_lung_genus_physeq = map(subset_group,~subset_physeq_by_group(lung_decontam_genus_physeq, .x)),
           subset_gut_genus_physeq = map(subset_group,~subset_physeq_by_group(gut_decontam_genus_physeq, .x)),
           subset_lung_metabolome_physeq = map(subset_group,~subset_physeq_by_group(lung_metabolome_physeq, .x)),
           subset_plasma_metabolome_physeq = map(subset_group,~subset_physeq_by_group(plasma_metabolome_physeq, .x)),
           data = pmap(list(subset_lung_genus_physeq,subset_gut_genus_physeq,subset_lung_metabolome_physeq,subset_plasma_metabolome_physeq),
                             ~get_omics_data(..1,..2,..3,..4, transform="")),
           clr_data = pmap(list(subset_lung_genus_physeq,subset_gut_genus_physeq,subset_lung_metabolome_physeq,subset_plasma_metabolome_physeq),
                                 ~get_omics_data(..1,..2,..3,..4, transform="rclr"))) %>% 
      dplyr::select(subset_group,data,clr_data)),
  #tar_target(omics_similarity_data, get_omics_similarity(omics_data,omics_clr_data)),
  tar_target(omics_similarity_data,omics_data %>% mutate(similarity_data= map2(data,clr_data,~get_omics_similarity(.x,.y)))),
  tar_target(all_omics_similarity_plot,commbine_all_omics_similarity_plot(omics_similarity_data)),
  
  
  # omics similarity of each group
  tar_target(omics_data2,tibble(subset_group = list(c("CTR","CA"),c("CA","VOR"),c("CA","WT"),c("CA","PptA"),c("WT","Tx"))) %>%
               mutate(subset_lung_genus_physeq = map(subset_group,~subset_physeq_by_group(lung_decontam_genus_physeq, .x)),
                      subset_gut_genus_physeq = map(subset_group,~subset_physeq_by_group(gut_decontam_genus_physeq, .x)),
                      subset_lung_metabolome_physeq = map(subset_group,~subset_physeq_by_group(lung_metabolome_physeq, .x)),
                      subset_plasma_metabolome_physeq = map(subset_group,~subset_physeq_by_group(plasma_metabolome_physeq, .x)),
                      data = pmap(list(subset_lung_genus_physeq,subset_gut_genus_physeq,subset_lung_metabolome_physeq,subset_plasma_metabolome_physeq),
                                  ~get_omics_data(..1,..2,..3,..4, transform="")),
                      clr_data = pmap(list(subset_lung_genus_physeq,subset_gut_genus_physeq,subset_lung_metabolome_physeq,subset_plasma_metabolome_physeq),
                                      ~get_omics_data(..1,..2,..3,..4, transform="rclr"))) %>%
               dplyr::select(subset_group,data,clr_data)),
  tar_target(omics_similarity_data2,omics_data2 %>% mutate(similarity_data= map2(data,clr_data,~get_omics_similarity(.x,.y)))),
  tar_target(all_omics_similarity_plot2,commbine_all_omics_similarity_plot2(omics_similarity_data2)),
  
  
  # beta diversity of gut baseline
  tar_target(gut_baseline_beta_plot, get_gut_baseline_beta_plot(gut_baseline_decontam_genus_physeq,plot_title="gut microbiome")),
  
  
  # plot composition of bacteria in gut
  tar_target(gut_taxa_composition_plot, plot_taxa_composition(gut_decontam_genus_physeq)),
  
  # plot composition of bacteria in lung nanopore
  tar_target(lung_nanopore_taxa_composition_plot, plot_taxa_composition(lung_nanopore_physeq,taxrank="Species",n_taxa=40)),
  
  # power analysis
  #tar_target(lung_power_analysis, get_power_analysis(lung_decontam_genus_physeq)),
  
  # cage effect
  tar_target(lung_cage_effect, plot_cage_effect(lung_decontam_genus_physeq)),
  
  # combine all figures
  tar_target(export_all_omics_similarity_plot, tar_ggsave(all_omics_similarity_plot,file.path(output_dir, "all_omics_similarity.pdf"),width=30,height=24), format = "file"),
  tar_target(export_all_omics_similarity_plot2, tar_ggsave(all_omics_similarity_plot2,file.path(output_dir, "all_omics_similarity2.pdf"),width=30,height=24), format = "file"),
  ## supplementary figures
  tar_target(export_lung_alpha_plot, tar_ggsave(lung_alpha_plot,file.path(output_dir, "lung_alpha.pdf"),width=10,height=10), format = "file"),
  tar_target(export_gut_alpha_plot, tar_ggsave(gut_alpha_plot,file.path(output_dir, "gut_alpha.pdf"),width=10,height=10), format = "file"),
  tar_target(export_gut_baseline_beta_plot, tar_ggsave(gut_baseline_beta_plot,file.path(output_dir, "gut_baseline_beta.pdf"),width=10,height=10), format = "file"),
  tar_target(export_lung_cage_effect, tar_ggsave(lung_cage_effect,file.path(output_dir, "lung_cage_effect.pdf"),width=8,height=6), format = "file"),
  ## final figures
  tar_target(export_fig5S1, tar_ggsave(gut_baseline_beta_plot,file.path(final_output_dir, "s1.pdf"),width=10,height=10), format = "file"),
  tar_target(export_fig5S2, tar_ggsave(gut_taxa_composition_plot,file.path(final_output_dir, "s2.pdf"),width=10,height=3), format = "file"),
  tar_target(export_fig5S3, tar_ggsave(lung_nanopore_taxa_composition_plot,file.path(final_output_dir, "s3.pdf"),width=15,height=3), format = "file")
)
