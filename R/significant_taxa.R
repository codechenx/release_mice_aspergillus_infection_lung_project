export_significant_taxa <- function(diff_taxa_df, pval_cutoff=0.05, padj_cutoff=1, output_dir){
  filtered_taxa_df <- diff_taxa_df %>%
    filter(pval < pval_cutoff & padj < padj_cutoff)
  file_name = "sig_diff_taxa.csv"
  output_path <- file.path(output_dir, file_name)
  write_csv(filtered_taxa_df,output_path)
  return(output_path)
}

# DESeq2
get_different_taxa_deseq2 <-
  function(physeq, output_dir_path = "") {
    library(DESeq2)
    taxrank <- rev(colnames(tax_table(physeq)))[1]
    meta_df <- data.frame(sample_data(physeq))
    if (!is.factor(meta_df$group)) {
      stop("the group column must be factor")
    }

    G1 <- levels(meta_df$group)[1]
    G2 <- levels(meta_df$group)[2]

    diagdds <- phyloseq_to_deseq2(physeq, ~ group)
    diagdds <- DESeq(diagdds, fitType = "local",sfType="poscounts")
    res <- results(diagdds, cooksCutoff = FALSE)
    sigtab <- res %>%
      as.data.frame %>%
      dplyr::rename(pval=pvalue) %>%
      dplyr::rename(log2FC=log2FoldChange) %>%
      mutate(FC=2**log2FC) %>%
      mutate(G1 = G1, G2 = G2) %>%
      mutate(comparison = paste0(G1, " vs ", G2)) %>%
      mutate(method = "DESeq2")
    sigtab <-
      cbind(sigtab, as(tax_table(physeq)[rownames(sigtab), ], "matrix")) %>%
      rownames_to_column(var = "OTU") %>%
      arrange(pval)
    if (output_dir_path != "") {
      file_name <-
        paste0(str_to_lower(taxrank), "_deseq2_", G1, " vs ", G2, ".csv")
      output_path <- file.path(output_dir_path, file_name)
      write_csv(sigtab, output_path)
    }
    return(sigtab)
  }
