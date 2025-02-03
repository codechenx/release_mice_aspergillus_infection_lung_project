# utils
extract_plot_diversity_for_comparision <-
  function(combined_diversity_df,G1,G2) {
    G1n = G1
    G2n = G2
    selected_df <-
      combined_diversity_df %>%
      mutate(comparison = paste0(G1, "_vs_", G2)) %>%
      group_by(index, comparison) %>%
      filter(G1==G1n&G2==G2n)
    alpha_plot <-
      ggarrange(plotlist = selected_df$plot, common.legend=T)
    return(alpha_plot)
  }

# Alpha diversity
calculate_alpha_diversity_data <- function(physeq) {
  meta_df <- sample_data(physeq) %>% as.data.frame()
  microbiome::alpha(physeq,index = c("Observed","Chao1","Shannon","Simpson","ACE","InvSimpson")) %>%
    rownames_to_column(var = "sample") %>%
    mutate(group = meta_df[.$sample] %>% pull(group))
}

get_all_alpha_diversity <- function(physeq) {
  alpha_diversity_data <- calculate_alpha_diversity_data(physeq)
  # chao1
  chao1_result <- alpha_diversity_data %>%
    dplyr::select(sample, group, chao1) %>%
    boxplot_alpha_diversity("chao1")
  # shannon
  shannon_result <- alpha_diversity_data %>%
    dplyr::select(sample, group, diversity_shannon) %>%
    boxplot_alpha_diversity("shannon")
  # combine all
  bind_rows(chao1_result, shannon_result)
}


boxplot_alpha_diversity <-
  function(alpha_diversity_df,
           diversity_name = "",
           paired = FALSE) {
    if (!is.factor(alpha_diversity_df$group)) {
      stop("the group column must be factor")
    }
    if (diversity_name == "") {
      diversity_name <- colnames(alpha_diversity_df)[3]
    }
    colnames(alpha_diversity_df) <- c("sample", "group", "diversity")
    G1 <- levels(alpha_diversity_df$group)[1]
    G2 <- levels(alpha_diversity_df$group)[2]
    color <- get_colormap_and_charmap() %>%
      filter(group %in% levels(alpha_diversity_df$group)) %>% pull(color)

    p <- alpha_diversity_df %>%
      wilcox_test(diversity ~ group, paired = paired) %>%
      pull(p)
    fc <-
      log2(
        mean(alpha_diversity_df %>% filter(group == G2) %>% pull(diversity)) / mean(alpha_diversity_df %>% filter(group == G1) %>% pull(diversity))
      )
    alpha_boxplot <- ggboxplot(
      alpha_diversity_df,
      x = "group",
      y = "diversity",
      add = "boxplot",
      fill = "group",
    ) + stat_compare_means(method = "wilcox.test",
                           label = "p.format",
                           label.x = 1.5,
                           size = 6) + xlab("") + theme_pubr() + ylab(diversity_name) +
      scale_fill_manual(values=color)
    return(tibble(
      index = diversity_name,
      G1 = G1,
      G2 = G2,
      fc = fc,
      p = p,
      plot = list(alpha_boxplot)
    ))
  }


# Beta diversity
get_all_beta_diversity <- function(physeq) {
  TMM_physeq <-
    transform_sample_counts(physeq, function(OTU)
      OTU / sum(OTU))
    # weighted UniFrac
    weighted_unifrac_result <-
      boxplot_beta_diversity(TMM_physeq,
                             "unifrac",
                             weighted = T,
                             diversity_name = "Weighted Unifrac")
  return(weighted_unifrac_result)
}


boxplot_beta_diversity <-
  function(physeq, index, weighted, diversity_name = "") {
    meta_df <- data.frame(sample_data(physeq))
    if (!is.factor(meta_df$group)) {
      stop("the group column must be factor")
    }
    G1 <- levels(meta_df$group)[1]
    G2 <- levels(meta_df$group)[2]
    color <- get_colormap_and_charmap() %>%
      filter(group %in% levels(meta_df$group)) %>% pull(color)
    if (diversity_name == "") {
      diversity_name <- index
    }
    ord_mat <-
      ordinate(physeq,
               method = "PCoA",
               distance = index,
               weighted = weighted)
    dist_mat <-
      phyloseq::distance(physeq, method = index, weighted = weighted)
    set.seed(42)
    adonis_result <-
      adonis2(dist_mat ~ group, data = meta_df, permutations = 999)
    p <- adonis_result$`Pr(>F)`[1]
    r2 <- round(adonis_result$`R2`[1],3)
    beta_boxplot <- plot_ordination(physeq, ord_mat, color = "group") +
      stat_ellipse(type = "norm",
                   linetype = 1,
                   lwd = 0.5) +
      theme_pubr() +
      theme(
        plot.title = element_text(size = 12, hjust = 0.5),
        text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "top",
        panel.border = element_rect(
          colour = "black",
          fill = NA,
          size = 1
        )
      ) +
      labs(tag = paste("r<sup>2</sup> =", r2,", p = ", p), title = diversity_name) +
      theme(plot.tag.position = c(0.5, 0.8), plot.tag = element_markdown()) +  scale_color_manual(values=color)
    return(tibble(
      index = diversity_name,
      G1 = G1,
      G2 = G2,
      r2 = r2,
      p = p,
      plot = list(beta_boxplot)
    ))
  }

# plsda
get_all_plsda <- function(physeq) {
  meta_df <- data.frame(sample_data(physeq))
  if (!is.factor(meta_df$group)) {
    stop("the group column must be factor")
  }
  G1 <- levels(meta_df$group)[1]
  G2 <- levels(meta_df$group)[2]
  color <- get_colormap_and_charmap() %>%
    filter(group %in% levels(meta_df$group)) %>% pull(color)
  otu_df <- as.data.frame(otu_table(physeq)) %>% rotate_df()

  plsda = opls(otu_df,meta_df$group,predI=2,scaleC="standard",fig.pdfC="none")

  adonis_result <- adonis2(otu_df~meta_df$group, method="euclidean")
  p <- adonis_result$`Pr(>F)`[1]
  r2 <- round(adonis_result$`R2`[1],3)

  sample.score = plsda@scoreMN %>%
    as.data.frame() %>%
    mutate(group = meta_df$group)

  boxplot = ggplot(sample.score, aes(p1, p2, color = group)) +
    geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) +
    geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
    geom_point() +
    geom_point(aes(-10,-10), color = 'white') +
    stat_ellipse(type = "norm",
                 linetype = 1,
                 lwd = 0.5) +
    theme_pubr() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      text = element_text(size = 10),
      legend.text = element_text(size = 12),
      legend.title = element_blank(),
      legend.position = "top",
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 1
      )
    ) +
    labs(tag = paste("r<sup>2</sup> =", r2,", p = ", p)) +
    labs(x = paste0("P1=",plsda@modelDF$R2X[1] *100,"%") ,y = paste0("P1=",plsda@modelDF$R2X[2] *100,"%"))+
    theme(plot.tag.position = c(0.5, 0.8), plot.tag = element_markdown()) +  scale_color_manual(values=color)
  return(tibble(
    index = "PLS-DA",
    G1 = G1,
    G2 = G2,
    r2=r2,
    p = p,
    plot = list(boxplot)
  ))
}
