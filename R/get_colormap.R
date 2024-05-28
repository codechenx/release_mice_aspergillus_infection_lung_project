# Get colormap for groups
get_colormap_and_charmap <- function(){
  library(RColorBrewer)
  col_char_map <- data.frame(
    group = factor(c(
      "CTR",
      "CA",
      "VOR",
      "WT",
      "Tx",
      "PptA"
    ), levels =c(
      "CTR",
      "CA",
      "VOR",
      "WT",
      "Tx",
      "PptA"
    )),
    color = c("#047F97","#F16522","#71C467","#94D2DC","#DC8676","#4876B9"),
    char = c(
      "#",
      "@",
      "$",
      "%",
      "&",
      "+"
    ))
  return(col_char_map)
}