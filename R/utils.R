tar_write_csv <- function(x, file, ...) {
  write_csv(x, file, ...)
  return(file)
}
tar_ggsave <- function(plot,filename,...){
  font_add_google("Noto Sans JP")
  showtext_auto()
  ggsave(filename=filename,plot=plot,bg = "transparent",...)
  showtext_auto(FALSE)
  return(filename)
}

tar_export_complexheatmap <- function(plot,filename,heatmap_legend_side="right",...){
  library(showtext)
  font_add_google("Noto Sans JP")
  showtext_auto()
  #if(!str_ends(filename, "pdf")|!str_ends(filename, "svg"))stop("wrong file format")
  pdf(file=filename,bg = "transparent", ...)
  draw(plot, background = "transparent",heatmap_legend_side=heatmap_legend_side)
  dev.off()
  showtext_auto(FALSE)
  return(filename)
}

tar_saveRDS<- function(object, file = "", ...) {
  saveRDS(object, file, ...)
  return(file)
}

utils_taxa_glom <- function(...){
  physeq <- tax_glom(...)
  ## remove unneeded taxa rank
  taxa_df <- tax_table(physeq)
  tax_table(physeq) <- taxa_df[,colSums(is.na(taxa_df))<nrow(taxa_df)]
  ## change OTU ID
  taxa_names(physeq) <- tax_table(physeq) %>% as.data.frame %>% unite("OTU",rank_names(physeq)) %>% pull(OTU)
  return(physeq)
}

utils_filter_taxa <- function(physeq, abundance_cutoff, prevalence_cutoff){
  TMM_physeq <-
    transform_sample_counts(physeq, function(OTU)
      OTU / sum(OTU))
  kept_taxa <- filter_taxa(TMM_physeq, function(x) sum(x > abundance_cutoff) > prevalence_cutoff*length(x),F)
  filtered_physeq <- prune_taxa(kept_taxa, physeq)
  return(filtered_physeq)
}

subset_physeq_by_group <- function(physeq, groups) {
  new_bac_sample_data_df <- sample_data(physeq) %>%
    data.frame() %>%
    filter(group %in% groups) %>%
    mutate(group = factor(group, levels = groups))
  sample_data(physeq) <- new_bac_sample_data_df
  return(physeq)
}


export_heatmap <- function(complexheatmap_objects, output_file_path, ...) {
  pdf(output_file_path, ...)
  draw(complexheatmap_objects)
  dev.off()
  return(output_file_path)
}

export_plots <- function(plots_list, plot_name="plot", output_dir, combine=F,ncol=NULL,nrow=NULL, ...) {
  if(is.ggplot(plots_list)){
    plot_name <- paste0(plot_name, ".pdf")
    output_path <- file.path(output_dir, plot_name)
    ggsave(output_path, plots_list, ...)
    return(output_path)
  }
  if(combine){
    plot_name <- paste0(plot_name, ".pdf")
    output_path <- file.path(output_dir, plot_name)
    combined_plot <- ggarrange(plotlist=plots_list,ncol=ncol,nrow=nrow)
    ggsave(output_path, combined_plot, ...)
    return(output_path)
  }
  if(length(plot_name)==1){
    plot_name <- paste0(plot_name, 1:length(plots_list))
    plot_name <- paste0(plot_name, ".pdf")
  }
  if(length(plots_list)!=length(plot_name))
    stop("Error: The plot_name and plots_list do not have the same length!")

  for(i in 1:length(plot_name)){
    output_path <- file.path(output_dir, plot_name[i])
    ggsave(output_path,
           plots_list[[i]],...)
  }
  return(output_path)
}

RV2.p <- function(mat1, mat2, permuation_n=999){
  stopifnot(nrow(mat1)==nrow(mat2))
  origianl_cor <- RV2(mat1, mat2)
  n <- nrow(mat1)
  set.seed(42)
  x <- replicate(permuation_n, {
    new_order <- sample(1:n)
    RV2(mat1[new_order,],mat2)
  })
  p<- (sum(abs(x)>abs(origianl_cor))+1)/(permuation_n+1)
  return(p)
}

fasta2dataframe <- function(fasta_path){
  fasta<- readDNAStringSet(fasta_path)
  df <- data.frame(sequence_name = names(fasta), sequence = as.character(fasta))
  return(df)
}

dataframe2fasta <- function(df, out_file_path){
  seqset <- DNAStringSet(as.character(df$sequence))
  names(seqset) <- df$sequence_name
  writeXStringSet(seqset, format="fasta", file=out_file_path)
}
