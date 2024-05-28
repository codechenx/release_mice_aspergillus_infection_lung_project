get_feast <- function(physeq){
  library(FEAST)
  library(stringi)
  taxrank <- rev(colnames(tax_table(physeq)))[1]
  physeq <- transform_sample_counts(physeq, function(OTU) round(OTU))
  meta_df <- sample_data(physeq) %>%
    data.frame() %>%
    filter(group %in% c("CTR", "CA", "VOR", "WT", "Tx")) %>%
    mutate(SourceSink = ifelse(group == "Tx", "Sink", "Source")) %>%
    mutate(Env=as.character(group)) %>%
    mutate(id=1:length(SourceSink)) %>%
    mutate(id= ifelse(SourceSink=="Sink", id, NA)) %>%
    dplyr::select(Env, SourceSink, id)
  otu_df <- as.data.frame(otu_table(physeq))
  rownames(otu_df) <- paste0("taxa",1:nrow(otu_df))
  random_file_name = stri_rand_strings(1,20)
  FEAST(C = t(otu_df), metadata = meta_df,different_sources_flag=0,dir_path=tempdir(), outfile=random_file_name)
  FEAST_result <- read.table(file.path(tempdir(),paste0(random_file_name,"_source_contributions_matrix.txt")), sep = "\t") %>%
    rownames_to_column(var="sample") %>%
    as_tibble()
  colnames(FEAST_result) <- str_replace_all(colnames(FEAST_result),"\\.","-")
  return(FEAST_result)
}