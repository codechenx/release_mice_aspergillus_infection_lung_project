library(targets)
library(tarchetypes)
library(tidyverse)
library(crew)
source("R/utils.R")
source("R/basic_stat.R")
source("R/get_data.R")

tar_option_set(
  packages = c(
    "ggpubr","ggplot2","patchwork","microViz","showtext","ggVennDiagram","circlize","ComplexHeatmap","Biostrings",
    "data.table","sjmisc","decontam","ape","phyloseq", "phylosmith","vegan","tidyverse"),
    controller = crew_controller_local(workers = tar_config_get("workers"))
)
list(
  # output path
  tar_target(data_output_dir, "processed_data"),
  tar_target(basic_stata_output_dir, "results/basic_stat"),
  
  # load meta data
  tar_target(metadata_file, "raw_data/metadata.csv", format = "file"),
  tar_target(meta_df, preprocess_metadata(metadata_file)),
  tar_target(lung_nanopore_metadata_file, "raw_data/lung_nanopore_metadata.csv", format = "file"),
  tar_target(lung_nanopore_meta_df, preprocess_lung_nanopore_metadata(lung_nanopore_metadata_file, meta_df)),
  
  # load lung metabolome
  tar_target(lung_metabolome_physeq_file, "raw_data/lung_mebo_pre.csv", format = "file"),
  tar_target(lung_metabolome_all_physeq, get_metabolome_physeq(lung_metabolome_physeq_file, meta_df)),
  tar_target(lung_metabolome_physeq, get_case_physeq(lung_metabolome_all_physeq)),
  tar_target(save_lung_metabolome_physeq, tar_saveRDS(lung_metabolome_physeq,file.path(data_output_dir, "lung_metabolome_physeq.rds")), format = "file"),
  
  # load plasma metabolome
  tar_target(plasma_metabolome_physeq_file, "raw_data/plasma_mebo_pre.csv", format = "file"),
  tar_target(plasma_metabolome_all_physeq, get_metabolome_physeq(plasma_metabolome_physeq_file, meta_df)),
  tar_target(plasma_metabolome_physeq, get_case_physeq(plasma_metabolome_all_physeq)),
  tar_target(save_plasma_metabolome_physeq, tar_saveRDS(plasma_metabolome_physeq,file.path(data_output_dir, "plasma_metabolome_physeq.rds")), format = "file"),
  
  # load lung data
  tar_target(lung_physeq_file, "raw_data/lung_phyloseq.Rdata", format = "file"),
  tar_target(lung_OTU_fna_file, "raw_data/lung_OTU.fna", format = "file"),
  tar_target(lung_genus_physeq, get_lung_physeq(lung_physeq_file, meta_df)),
  tar_target(lung_mito,get_possible_mitochondrial_otu(lung_genus_physeq, fasta2dataframe(lung_OTU_fna_file))),
  tar_target(save_lung_mito, tar_write_csv(lung_mito, file.path(data_output_dir, "lung_mito.csv")), format = "file"),
  tar_target(export_lung_mito, dataframe2fasta(lung_mito, file.path(data_output_dir, "lung_mito.fna")), format = "file"),
  tar_target(lung_genus_no_mito_physeq, remove_mitochondrial_otu_physeq(lung_genus_physeq, c("OTU1","OTU6","OTU131","OTU108", "OTU89","OTU138"))),
  tar_target(lung_clean_otu_physeq, remove_spurious_taxa(lung_genus_no_mito_physeq, abundance_cutoff=0.0025,freq_cutoff=1)),
  tar_target(save_lung_otu_table, tar_write_csv(get_otu_before_decontam(lung_clean_otu_physeq),file.path(basic_stata_output_dir,"lung_otu_table_before_decontam.csv"))),
  tar_target(lung_heatmap_detected_otu, get_hetamap_detected_taxa(lung_genus_no_mito_physeq)),
  tar_target(export_lung_heatmap_detected_otu, export_heatmap(lung_heatmap_detected_otu,file.path(basic_stata_output_dir, "lung_heatmap_detected_otu.pdf"), width = 50, height = 30)),
  tar_target(lung_decontam_all_genus_physeq, get_decontaminated_lung_physeq(lung_clean_otu_physeq, output_dir_path = basic_stata_output_dir)),
  tar_target(lung_decontam_genus_physeq, get_case_physeq(lung_decontam_all_genus_physeq)),
  tar_target(save_lung_decontam_genus_physeq, tar_saveRDS(lung_decontam_genus_physeq,file.path(data_output_dir, "lung_decontam_genus_physeq.rds")), format = "file"),
  tar_target(lung_distribution_of_decontamination_barplot, get_distribution_of_decontamination(lung_clean_otu_physeq %>% tax_glom("Genus"), lung_decontam_all_genus_physeq)),
  tar_target(exprot_lung_distribution_of_decontamination_barplot, tar_ggsave(lung_distribution_of_decontamination_barplot,file.path(basic_stata_output_dir, "lung_distribution_of_decontamination_barplot.pdf"),width=12,height=10), format = "file"),
  # basic statistics of lung data
  tar_target(lung_read_count_data, get_microbial_loading(lung_clean_otu_physeq)),
  tar_target(export_lung_read_count_data, tar_write_csv(lung_read_count_data, file.path(basic_stata_output_dir, "lung_read_count_by_sample.csv")), format = "file"),
  tar_target(lung_taxa_composition_plot, get_taxa_composition_plot(lung_decontam_genus_physeq)),
  tar_target(export_lung_taxa_composition_plot, tar_ggsave(lung_taxa_composition_plot,file.path(basic_stata_output_dir, "lung_composition.pdf"),width=20,height=10), format = "file"),
  
  tar_target(lung_negative_control_genus_physeq,tax_glom(get_negative_control_physeq(lung_clean_otu_physeq),taxrank = "Genus")),
  tar_target(lung_negative_control_genus_composition_plot, get_taxa_composition_plot(lung_negative_control_genus_physeq,group_by ="sampleGroup2")),
  tar_target(export_lung_negative_control_genus_composition_plot, tar_ggsave(lung_negative_control_genus_composition_plot,file.path(basic_stata_output_dir, "lung_negative_control_genus_composition.pdf"),width=20,height=20), format = "file"),

  
  # load gut data
  tar_target(gut_physeq_file, "raw_data/gut_phyloseq.Rdata", format = "file"),
  tar_target(gut_OTU_fna_file, "raw_data/gut_OTU.fna", format = "file"),
  tar_target(gut_genus_physeq, get_gut_physeq(gut_physeq_file, meta_df)),
  tar_target(gut_mito,get_possible_mitochondrial_otu(gut_genus_physeq, fasta2dataframe(gut_OTU_fna_file))),
  tar_target(save_gut_mito, tar_write_csv(gut_mito, file.path(data_output_dir, "gut_mito.csv")), format = "file"),
  tar_target(export_gut_mito, dataframe2fasta(gut_mito, file.path(data_output_dir, "gut_mito.fna")), format = "file"),
  tar_target(gut_genus_no_mito_physeq, remove_mitochondrial_otu_physeq(gut_genus_physeq, c("OTU415"))),
  tar_target(gut_clean_otu_physeq, remove_spurious_taxa(gut_genus_no_mito_physeq, abundance_cutoff=0.0025,freq_cutoff=0)),
  tar_target(save_gut_otu_table, tar_write_csv(get_otu_before_decontam(gut_clean_otu_physeq),file.path(basic_stata_output_dir,"gut_otu_table_before_decontam.csv"))),
  tar_target(gut_heatmap_detected_otu, get_hetamap_detected_taxa(gut_genus_no_mito_physeq)),
  tar_target(export_gut_heatmap_detected_otu, export_heatmap(gut_heatmap_detected_otu,file.path(basic_stata_output_dir, "gut_heatmap_detected_otu.pdf"), width = 50, height = 30)),  
  tar_target(gut_decontam_all_genus_physeq, get_decontaminated_gut_physeq(gut_clean_otu_physeq, output_dir_path = basic_stata_output_dir)),
  tar_target(gut_decontam_genus_physeq, get_case_physeq(gut_decontam_all_genus_physeq)),
  tar_target(save_gut_decontam_genus_physeq, tar_saveRDS(gut_decontam_genus_physeq,file.path(data_output_dir, "gut_decontam_genus_physeq.rds")), format = "file"),
  ## load gut baseline
  tar_target(gut_baseline_genus_physeq, load_gut_baseline_physeq(gut_physeq_file, meta_df)),
  tar_target(gut_baseline_genus_no_mito_physeq, remove_mitochondrial_otu_physeq(gut_baseline_genus_physeq, c("OTU415"))),
  tar_target(gut_baseline_clean_genus_physeq, remove_spurious_taxa(gut_baseline_genus_no_mito_physeq, abundance_cutoff=0.0025,freq_cutoff=0)),
  tar_target(gut_baseline_decontam_all_genus_physeq, get_decontaminated_gut_physeq(gut_baseline_clean_genus_physeq, output_dir_path = basic_stata_output_dir)),
  tar_target(save_gut_baseline_decontam_genus_physeq, tar_saveRDS(gut_baseline_decontam_all_genus_physeq,file.path(data_output_dir, "gut_baseline_decontam_genus_physeq.rds")), format = "file"),
  ## basic statistics of gut data
  tar_target(gut_read_count_data, get_microbial_loading(gut_clean_otu_physeq)),
  tar_target(export_gut_read_count_data, tar_write_csv(gut_read_count_data, file.path(basic_stata_output_dir, "gut_read_count_by_sample.csv")), format = "file"),
  tar_target(gut_taxa_composition_plot, get_taxa_composition_plot(gut_decontam_genus_physeq)),
  tar_target(export_gut_taxa_composition_plot, tar_ggsave(gut_taxa_composition_plot,file.path(basic_stata_output_dir, "gut_composition.pdf"),width=20,height=10), format = "file"),
  
  tar_target(gut_negative_control_genus_physeq,tax_glom(get_negative_control_physeq(gut_clean_otu_physeq),taxrank = "Genus")),
  tar_target(gut_negative_control_taxa_composition_plot, get_taxa_composition_plot(gut_negative_control_genus_physeq,group_by ="sampleGroup2")),
  tar_target(export_gut_negative_control_taxa_composition_plot, tar_ggsave(gut_negative_control_taxa_composition_plot,file.path(basic_stata_output_dir, "gut_negative_control_taxa_composition.pdf"),width=20,height=10), format = "file"),
  
  
  # load lung nanopore data
  tar_target(lung_nanopore_physeq_file, "raw_data/lung_nanopore_OTU.tsv", format = "file"),
  tar_target(lung_nanopore_otu_physeq, get_lung_nanopore_physeq(lung_nanopore_physeq_file, lung_nanopore_meta_df)),
  tar_target(lung_nanopore_clean_otu_physeq, remove_spurious_taxa(lung_nanopore_otu_physeq, abundance_cutoff=0.0001,freq_cutoff=1)),
  tar_target(save_lung_nanopore_otu_table, tar_write_csv(get_otu_before_decontam(lung_nanopore_clean_otu_physeq),file.path(basic_stata_output_dir,"lung_nanopore_otu_table_before_decontam.csv"))),
  tar_target(lung_nanopore_heatmap_detected_otu, get_hetamap_detected_taxa(lung_nanopore_clean_otu_physeq, count=F)),
  tar_target(export_lung_nanopore_heatmap_detected_otu, export_heatmap(lung_nanopore_heatmap_detected_otu,file.path(basic_stata_output_dir, "lung_nanopore_heatmap_detected_otu.pdf"), width = 50, height = 30)),
  tar_target(lung_nanopore_decontam_all_genus_physeq, get_decontaminated_lung_nanopore_physeq(lung_nanopore_clean_otu_physeq, output_dir_path = basic_stata_output_dir)),
  tar_target(lung_nanopore_decontam_genus_physeq, get_case_physeq(lung_nanopore_decontam_all_genus_physeq)),
  tar_target(save_lung_nanopore_decontam_genus_physeq, tar_saveRDS(lung_nanopore_decontam_genus_physeq,file.path(data_output_dir, "lung_nanopore_decontam_genus_physeq.rds")), format = "file"),
  tar_target(save_lung_nanopore_decontam_otu_table, lung_nanopore_decontam_genus_physeq %>% melt_phyloseq() %>% dplyr::select(Sample, OTU, Domain:Species, Abundance) %>% pivot_wider(names_from = Sample, values_from = Abundance) %>% tar_write_csv(file.path(basic_stata_output_dir,"lung_nanopore_otu_table_after_decontam.csv")))
)
