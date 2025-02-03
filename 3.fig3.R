library(targets)
library(tarchetypes)
library(tidyverse)
library(crew)
source("R/utils.R")
source("R/get_data.R")
source("R/get_colormap.R")
source("R/diversity.R")
source("R/significant_taxa.R")
source("R/modima.R")
source("R/fig3.R")
tar_option_set(
  packages = c(
    "ggpubr","ggplot2","ggsci","microViz","circlize","ComplexHeatmap","ggtext",
    "showtext","ropls","sjmisc","MatrixCorrelation","rstatix","data.table",
    "psych","ggfittext","ggalluvial","energy","mediation",
    "microbiome","decontam","ape","phylosmith","phyloseq","vegan","tidyverse"),
  controller = crew_controller_local(workers = tar_config_get("workers"))
)


list(
  # output path
  tar_target(output_dir, "results/fig3"),
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

  tar_target(fig3_comparison_df, get_fig3_comparison_df()),
  tar_target(lung_subset_decontam_genus_physeq, subset_phyloseq_by_comparison(lung_decontam_genus_physeq, G1 = fig3_comparison_df$G1, G2 = fig3_comparison_df$G2),
             pattern = fig3_comparison_df,
             iteration = "list"
  ),
  tar_target(gut_subset_decontam_genus_physeq, subset_phyloseq_by_comparison(gut_decontam_genus_physeq, G1 = fig3_comparison_df$G1, G2 = fig3_comparison_df$G2),
             pattern = fig3_comparison_df,
             iteration = "list"
  ),
  tar_target(lung_subset_metabolome_physeq, subset_phyloseq_by_comparison(lung_metabolome_physeq, G1 = fig3_comparison_df$G1, G2 = fig3_comparison_df$G2),
             pattern = fig3_comparison_df,
             iteration = "list"
  ),
  tar_target(plasma_subset_metabolome_physeq, subset_phyloseq_by_comparison(plasma_metabolome_physeq, G1 = fig3_comparison_df$G1, G2 = fig3_comparison_df$G2),
             pattern = fig3_comparison_df,
             iteration = "list"
  ),

  # split data by group
  ## Fig3 data
  tar_target(fig3_lung_decontam_genus_physeq, get_fig3_physeq(lung_decontam_genus_physeq)),
  tar_target(fig3_gut_decontam_genus_physeq, get_fig3_physeq(gut_decontam_genus_physeq)),
  tar_target(fig3_lung_metabolome_physeq, get_fig3_physeq(lung_metabolome_physeq)),
  tar_target(fig3_plasma_metabolome_physeq, get_fig3_physeq(plasma_metabolome_physeq)),
  ## filter genus by abundance and prevalence
  tar_target(fig3_lung_filtered_decontam_genus_physeq,utils_filter_taxa(fig3_lung_decontam_genus_physeq,abundance_cutoff = 0.0001,prevalence_cutoff = 0.1)),
  tar_target(fig3_gut_filtered_decontam_genus_physeq,utils_filter_taxa(fig3_gut_decontam_genus_physeq,abundance_cutoff = 0.0001,prevalence_cutoff = 0.1)),
  ## PBS_CA(CA)
  tar_target(pbs_lung_decontam_genus_physeq, get_pbs_ca_physeq(fig3_lung_filtered_decontam_genus_physeq)),
  tar_target(pbs_gut_decontam_genus_physeq, get_pbs_ca_physeq(fig3_gut_filtered_decontam_genus_physeq)),
  tar_target(pbs_lung_metabolome_physeq, get_pbs_ca_physeq(lung_metabolome_physeq)),
  tar_target(pbs_plasma_metabolome_physeq, get_pbs_ca_physeq(plasma_metabolome_physeq)),
  ## transform PBS_CA(CA)
  tar_target(pbs_omics_data, get_omics_data(pbs_lung_decontam_genus_physeq, pbs_gut_decontam_genus_physeq, pbs_lung_metabolome_physeq, pbs_plasma_metabolome_physeq)),
  tar_target(pbs_omics_filtered_data, get_omics_data(pbs_lung_decontam_genus_physeq, pbs_gut_decontam_genus_physeq, pbs_lung_metabolome_physeq, pbs_plasma_metabolome_physeq)),
  tar_target(pbs_omics_clr_data, get_omics_data(pbs_lung_decontam_genus_physeq, pbs_gut_decontam_genus_physeq, pbs_lung_metabolome_physeq, pbs_plasma_metabolome_physeq, transform="rclr")),
  tar_target(pbs_omics_filtered_clr_data, get_omics_data(pbs_lung_decontam_genus_physeq, pbs_gut_decontam_genus_physeq, pbs_lung_metabolome_physeq, pbs_plasma_metabolome_physeq,transform="rclr")),
  
  ## wild type AF(WT)
  tar_target(wt_lung_decontam_genus_physeq, get_af_wt_physeq(fig3_lung_filtered_decontam_genus_physeq)),
  tar_target(wt_gut_decontam_genus_physeq, get_af_wt_physeq(fig3_gut_filtered_decontam_genus_physeq)),
  tar_target(wt_lung_metabolome_physeq, get_af_wt_physeq(lung_metabolome_physeq)),
  tar_target(wt_plasma_metabolome_physeq, get_af_wt_physeq(plasma_metabolome_physeq)),
  ## transform wild type AF(WT)
  tar_target(wt_omics_data, get_omics_data(wt_lung_decontam_genus_physeq, wt_gut_decontam_genus_physeq, wt_lung_metabolome_physeq, wt_plasma_metabolome_physeq)),
  tar_target(wt_omics_filtered_data, get_omics_data(wt_lung_decontam_genus_physeq, wt_gut_decontam_genus_physeq, wt_lung_metabolome_physeq, wt_plasma_metabolome_physeq)),
  tar_target(wt_omics_clr_data, get_omics_data(wt_lung_decontam_genus_physeq, wt_gut_decontam_genus_physeq, wt_lung_metabolome_physeq, wt_plasma_metabolome_physeq, transform="rclr")),
  tar_target(wt_omics_filtered_clr_data, get_omics_data(wt_lung_decontam_genus_physeq, wt_gut_decontam_genus_physeq, wt_lung_metabolome_physeq, wt_plasma_metabolome_physeq,transform="rclr")),
  
  ## mutation type AF(PptA)
  tar_target(ppta_lung_decontam_genus_physeq, get_af_ppta_physeq(fig3_lung_filtered_decontam_genus_physeq)),
  tar_target(ppta_gut_decontam_genus_physeq, get_af_ppta_physeq(fig3_gut_filtered_decontam_genus_physeq)),
  tar_target(ppta_lung_metabolome_physeq, get_af_ppta_physeq(lung_metabolome_physeq)),
  tar_target(ppta_plasma_metabolome_physeq, get_af_ppta_physeq(plasma_metabolome_physeq)),
  ## transform mutation type AF(PptA)
  tar_target(ppta_omics_data, get_omics_data(ppta_lung_decontam_genus_physeq, ppta_gut_decontam_genus_physeq, ppta_lung_metabolome_physeq, ppta_plasma_metabolome_physeq)),
  tar_target(ppta_omics_filtered_data, get_omics_data(ppta_lung_decontam_genus_physeq, ppta_gut_decontam_genus_physeq, ppta_lung_metabolome_physeq, ppta_plasma_metabolome_physeq)),
  tar_target(ppta_omics_clr_data, get_omics_data(ppta_lung_decontam_genus_physeq, ppta_gut_decontam_genus_physeq, ppta_lung_metabolome_physeq, ppta_plasma_metabolome_physeq, transform="rclr")),
  tar_target(ppta_omics_filtered_clr_data, get_omics_data(ppta_lung_decontam_genus_physeq, ppta_gut_decontam_genus_physeq, ppta_lung_metabolome_physeq, ppta_plasma_metabolome_physeq,transform="rclr")),
  
  

  # alpha diversity
  tar_target(lung_alpha_plot, plot_fig3_alpha(fig3_lung_decontam_genus_physeq,plot_title="")),
  tar_target(gut_alpha_plot, plot_fig3_alpha(fig3_gut_decontam_genus_physeq,plot_title="")),


  # beta diversity
  tar_target(lung_beta_plot, plot_fig3_beta(fig3_lung_decontam_genus_physeq,plot_title="Lung microbiome")),
  tar_target(gut_beta_plot, plot_fig3_beta(fig3_gut_decontam_genus_physeq,plot_title="Gut microbiome")),
  # plsda
  tar_target(lung_plsda_plot, plot_fig3_plsda(fig3_lung_metabolome_physeq,plot_title="Lung metabolome")),
  tar_target(plasma_plsda_plot, plot_fig3_plsda(fig3_plasma_metabolome_physeq,plot_title="Plasma metabolome")),
  tar_target(lung_metabolome_beta_plot,plot_fig3_beta(fig3_lung_metabolome_physeq,distance="euclidean",weighted = F,plot_title="Lung metabolome")),
  tar_target(plasma_metabolome_beta_plot,plot_fig3_beta(fig3_plasma_metabolome_physeq,distance="euclidean",weighted = F,plot_title="Plasma metabolome")),
  
  # significant taxa(fig3b)
  ## filter genus by abundance and prevalence
  tar_target(lung_subset_filtered_decontam_genus_physeq, subset_phyloseq_by_comparison(fig3_lung_filtered_decontam_genus_physeq, G1 = fig3_comparison_df$G1, G2 = fig3_comparison_df$G2),
             pattern = fig3_comparison_df,
             iteration = "list"
  ),
  tar_target(gut_subset_filtered_decontam_genus_physeq, subset_phyloseq_by_comparison(fig3_gut_filtered_decontam_genus_physeq, G1 = fig3_comparison_df$G1, G2 = fig3_comparison_df$G2),
             pattern = fig3_comparison_df,
             iteration = "list"
  ),
  ## DESeq2
  tar_target(lung_different_genus_deseq2,get_different_taxa_deseq2(lung_subset_filtered_decontam_genus_physeq),pattern = lung_subset_filtered_decontam_genus_physeq,iteration = "list"),
  tar_target(lung_combined_different_genus_deseq2,bind_rows(lung_different_genus_deseq2)),
  tar_target(export_lung_deseq2_results, tar_write_csv(lung_combined_different_genus_deseq2,file.path(output_dir, "lung_different_genus_deseq2.csv"))),
  tar_target(gut_different_genus_deseq2,get_different_taxa_deseq2(gut_subset_filtered_decontam_genus_physeq),pattern = gut_subset_filtered_decontam_genus_physeq,iteration = "list"),
  tar_target(gut_combined_different_genus_deseq2,bind_rows(gut_different_genus_deseq2)),
  tar_target(export_gut_deseq2_results, tar_write_csv(gut_combined_different_genus_deseq2,file.path(output_dir, "gut_different_genus_deseq2.csv"))),
  ## heatmap
  tar_target(microbiome_heatmap, plot_fig3_microbiome_heatmap(lung_decontam_genus_physeq,lung_combined_different_genus_deseq2, gut_decontam_genus_physeq, gut_combined_different_genus_deseq2)),
  tar_target(metabolome_heatmap, plot_fig3_metabolome_heatmap(lung_metabolome_physeq,lung_metabolome_sig_df, plasma_metabolome_physeq, plasma_metabolome_sig_df,n=10)),


  # compare correlation change
  tar_target(pbs_correlation_results, get_correlation_between_omics(pbs_omics_filtered_clr_data)),
  tar_target(wt_correlation_results, get_correlation_between_omics(wt_omics_filtered_clr_data)),
  tar_target(ppta_correlation_results, get_correlation_between_omics(ppta_omics_filtered_clr_data)),
  ## wt and control
  tar_target(compare_correlation_between_wt_pbs_results, compare_correlation_between_omics(pbs_correlation_results %>% mutate(group="pbs"),
                                                                                           wt_correlation_results%>% mutate(group="wt"))),
  tar_target(export_compare_correlation_between_wt_pbs_results,tar_write_csv(compare_correlation_between_wt_pbs_results, file.path(output_dir,"compare_correlation_between_wt_pbs_results.csv")),format = "file"),
  tar_target(correlation_between_wt_pbs_chord, plot_correlation_between_omics_chord(compare_correlation_between_wt_pbs_results)),
  ## wt and mutation
  tar_target(compare_correlation_between_wt_ppa_results, compare_correlation_between_omics(ppta_correlation_results %>% mutate(group="ppta"),
                                                                                           wt_correlation_results%>% mutate(group="wt"))),
  tar_target(export_compare_correlation_between_wt_ppa_results,tar_write_csv(compare_correlation_between_wt_ppa_results, file.path(output_dir,"compare_correlation_between_wt_ppa_results.csv")),format = "file"),
  tar_target(correlation_between_wt_ppa_chord, plot_correlation_between_omics_chord(compare_correlation_between_wt_ppa_results)),


  # modiam analysis
  ## PBS_CA
  tar_target(pbs_modima_result, get_modima_result(pbs_lung_decontam_genus_physeq, pbs_gut_decontam_genus_physeq, pbs_lung_metabolome_physeq, pbs_plasma_metabolome_physeq)),
  tar_target(pbs_export_modima_results,tar_write_csv(pbs_modima_result, file.path(output_dir,"pbs_modima.csv")),format = "file"),
  ## wild type AF
  tar_target(wt_modima_result, get_modima_result(wt_lung_decontam_genus_physeq, wt_gut_decontam_genus_physeq, wt_lung_metabolome_physeq, wt_plasma_metabolome_physeq)),
  tar_target(export_wt_modima_results,tar_write_csv(wt_modima_result, file.path(output_dir,"wt_modima.csv")),format = "file"),
  ## mutation type AF
  tar_target(ppta_modima_result, get_modima_result(ppta_lung_decontam_genus_physeq, ppta_gut_decontam_genus_physeq, ppta_lung_metabolome_physeq, ppta_plasma_metabolome_physeq)),
  tar_target(export_ppta_modima_results,tar_write_csv(ppta_modima_result, file.path(output_dir,"ppta_modima.csv")),format = "file"),
  ## compare all
  tar_target(modima_barplot, plot__modima_result(pbs_modima_result,wt_modima_result,ppta_modima_result)),


  # # mediation analysis by PCOA reduced dimension data
  # ## PBS_CA
  # tar_target(pbs_pcoa_mediation_combination, get_pcoa_reduced_mediation_combination(pbs_omics_filtered_data,k=5)),
  # tar_target(pbs_pcoa_mediation_results,pbs_pcoa_mediation_combination %>% rowwise %>% mutate(results = list(get_mediation(data)))),
  # tar_target(pbs_pcoa_mediation_summary_tbl, extract_pcoa_reduced_mediation(pbs_pcoa_mediation_results)),
  # ## wild type AF
  # tar_target(wt_pcoa_mediation_combination, get_pcoa_reduced_mediation_combination(wt_omics_filtered_data,k=5)),
  # tar_target(wt_pcoa_mediation_results,wt_pcoa_mediation_combination %>% rowwise %>% mutate(results = list(get_mediation(data)))),
  # tar_target(wt_pcoa_mediation_summary_tbl, extract_pcoa_reduced_mediation(wt_pcoa_mediation_results)),
  # ## mutation type AF
  # tar_target(ppta_pcoa_mediation_combination, get_pcoa_reduced_mediation_combination(ppta_omics_filtered_data,k=5)),
  # tar_target(ppta_pcoa_mediation_results,ppta_pcoa_mediation_combination %>% rowwise %>% mutate(results = list(get_mediation(data)))),
  # tar_target(ppta_pcoa_mediation_summary_tbl, extract_pcoa_reduced_mediation(ppta_pcoa_mediation_results)),
  # ## plot and compare mediation
  # tar_target(ppta_pcoa_mediation_summary_barplot, plot_pcoa_reduced_mediation_bar(pbs_pcoa_mediation_summary_tbl,wt_pcoa_mediation_summary_tbl,ppta_pcoa_mediation_summary_tbl)),


  # mediation analysis by hima
  tar_target(sig_data, tibble(name=c("gut_microbiome","lung_microbiome","lung_metabolome","plasma_metabolome"),data=list(gut_combined_different_genus_deseq2,lung_combined_different_genus_deseq2,lung_metabolome_sig_df,plasma_metabolome_sig_df))),
  ## PBS_CA
  tar_target(pbs_mediation_result, get_hima_result(pbs_omics_filtered_clr_data)),
  tar_target(pbs_sig_taxa_mediation_result, get_hima_result(pbs_omics_filtered_clr_data, sig_data)),
  ## wild type AF
  tar_target(wt_mediation_result, get_hima_result(wt_omics_filtered_clr_data)),
  tar_target(wt_sig_taxa_mediation_result, get_hima_result(wt_omics_filtered_clr_data, sig_data)),
  ## mutation type AF
  tar_target(ppta_mediation_result, get_hima_result(ppta_omics_filtered_clr_data)),
  tar_target(ppta_sig_taxa_mediation_result, get_hima_result(ppta_omics_filtered_clr_data, sig_data)),
  ## compare mediation proportion between wild type and mutation type AF
  tar_target(overall_mediation_boxplot, plot_hima_box(pbs_mediation_result, wt_mediation_result, ppta_mediation_result,scale_med_prop=T)),
  tar_target(overall_sig_taxa_mediation_boxplot, plot_hima_box(pbs_sig_taxa_mediation_result, wt_sig_taxa_mediation_result, ppta_sig_taxa_mediation_result,scale_med_prop=T)),
  ## extract mediation
  tar_target(extracted_mediation, extract_hima_mediation(wt_mediation_result,ppta_mediation_result)),
  tar_target(export_extracted_mediation,tar_write_csv(extracted_mediation, file.path(output_dir,"hima_extracted_mediation.csv")),format = "file"),
  tar_target(extracted_sig_taxa_mediation, extract_hima_mediation(wt_sig_taxa_mediation_result,ppta_sig_taxa_mediation_result)),
  tar_target(export_extracted_sig_taxa_mediation,tar_write_csv(extracted_sig_taxa_mediation, file.path(output_dir,"hima_extracted_sig_taxa_mediation.csv")),format = "file"),


  # mediation analysis for Ligilactobacillus
  ## wild type AF
  tar_target(wt_ligilactobacillus_mediation_combination,get_ligilactobacillus_mediation_combination(wt_omics_filtered_clr_data, sig_data)),
  tar_target(wt_ligilactobacillus_mediation_results,wt_ligilactobacillus_mediation_combination %>% rowwise %>% mutate(results = list(get_mediation(data)))),
  ## mutation type AF
  tar_target(ppta_ligilactobacillus_mediation_combination,get_ligilactobacillus_mediation_combination(ppta_omics_filtered_clr_data, sig_data)),
  tar_target(ppta_ligilactobacillus_mediation_results,ppta_ligilactobacillus_mediation_combination %>% rowwise %>% mutate(results = list(get_mediation(data)))),
  ## extract mediation results
  tar_target(ligilactobacillus_mediation_summary_tbl, extract_mediation(wt_ligilactobacillus_mediation_results,ppta_ligilactobacillus_mediation_results)),
  tar_target(export_ligilactobacillus_mediation_summary_tbl,tar_write_csv(ligilactobacillus_mediation_summary_tbl, file.path(output_dir,"ligilactobacillus_mediation_summary.csv")),format = "file"),
  
  
  # combine all figures
  tar_target(export_lung_alpha_plot, tar_ggsave(lung_alpha_plot, file.path(output_dir, "lung_alpha_diversity.pdf"))),
  tar_target(export_gut_alpha_plot, tar_ggsave(gut_alpha_plot, file.path(output_dir, "gut_alpha_diversity.pdf"))),
  tar_target(export_lung_metabolome_beta_plot, tar_ggsave(lung_metabolome_beta_plot, file.path(output_dir, "lung_metabolome_beta_diversity.pdf"))),
  tar_target(export_plasma_metabolome_beta_plot, tar_ggsave(plasma_metabolome_beta_plot, file.path(output_dir, "plasma_metabolome_beta_diversity.pdf"))),
  ## fig3a
  tar_target(fig3a, ggarrange(plotlist = list(lung_beta_plot, gut_beta_plot,lung_plsda_plot,plasma_plsda_plot),nrow = 1,align = "hv", common.legend = TRUE, legend="top")),
  tar_target(export_fig3a, tar_ggsave(fig3a,file.path(output_dir, "fig3a.pdf"),width=25,height=7), format = "file"),
  ## fig3b
  tar_target(export_microbiome_heatmap, tar_export_complexheatmap(microbiome_heatmap,file.path(output_dir, "fig3b_microbiome_heatmap.pdf"),heatmap_legend_side = "bottom", width=12,height=4), format = "file"),
  tar_target(export_metabolome_heatmap, tar_export_complexheatmap(metabolome_heatmap,file.path(output_dir, "fig3b_metabolome_heatmap.pdf"),heatmap_legend_side = "bottom", width=12,height=4), format = "file"),
  ## fig3c
  tar_target(export_correlation_between_wt_ppta_chord,tar_ggsave(correlation_between_wt_ppa_chord, file.path(output_dir,"fig3_supplementary.pdf"), width=10,height=8),format = "file"),
  tar_target(export_overall_sig_taxa_mediation_boxplot,tar_ggsave(overall_sig_taxa_mediation_boxplot, file.path(output_dir,"fig3c.pdf"), width=8,height=8),format = "file"),
  tar_target(export_overall_mediation_boxplot,tar_ggsave(overall_mediation_boxplot, file.path(output_dir,"fig3c2.pdf"), width=8,height=8),format = "file"),
  tar_target(export_modima_barplot,tar_ggsave(modima_barplot, file.path(output_dir,"modiam.pdf"), width=18,height=10),format = "file")
)