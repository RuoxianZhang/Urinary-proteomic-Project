source(file.path(Sys.getenv("UR_FILE0624_PROJECT_DIR", unset = Sys.getenv("UR_PROJECT_DIR", unset = Sys.getenv("UR_FILE0623_PROJECT_DIR", unset = getwd()))), "code", "00_common_functions.R"))

ensure_output_dirs()

read_result <- function(...) read.csv(file.path(result_dir, ...), check.names = FALSE)
panel_data <- function(name, x) {
  write_csv_safe(x, file.path(dir_map$plotted, paste0(name, ".csv")))
  x
}

old_panels <- list.files(dir_map$plotted, pattern = "[.]csv$", full.names = TRUE)
if (length(old_panels) > 0) unlink(old_panels)

page_in <- function(width_pt, height_pt) c(width = width_pt / 72, height = height_pt / 72)
gold_pages <- list(
  Figure1 = page_in(595, 841),
  Figure2 = page_in(806, 1116),
  Figure3 = page_in(595, 1094),
  Figure4 = page_in(595, 1065),
  Figure5 = page_in(595, 1008),
  Figure6 = page_in(595.275591, 841.889764),
  Supplementary_Figure_1 = page_in(662, 255),
  Supplementary_Figure_2 = page_in(619, 871),
  Supplementary_Figure_3 = page_in(619, 295),
  Supplementary_Figure_4 = page_in(619, 457),
  Supplementary_Figure_5 = page_in(1238.4, 368.16)
)

gold_pal <- c(
  Healthy = "#67B7C4",
  LungCancer = "#C75A56",
  `Lung cancer` = "#C75A56",
  Lung_cancer = "#C75A56",
  Control = "#67B7C4",
  Other = "grey80",
  NSCLC = "#5A8FB0",
  SCLC = "#6977A5",
  LUAD = "#5BAA87",
  LUSC = "#A77E58",
  N0 = "#AEBBC5",
  Nplus = "#C75A56",
  Female = "#D8A4A1",
  Male = "#7EA9C2",
  Smoker = "#8F7A63",
  `Non-smoker` = "#B8C7D1",
  Drinker = "#A98962",
  `Non-drinker` = "#C8D1D8",
  No_or_not_recorded = "#B8C7D1",
  `Age <=60` = "#BFD7DC",
  `Age >60` = "#C9A4A1",
  Discovery = "#738B9C",
  `Discovery CV` = "#738B9C",
  `Discovery 5-fold out-of-fold` = "#738B9C",
  `Discovery 5-fold OOF` = "#738B9C",
  `Discovery OOB` = "#738B9C",
  Validation = "#C75A56",
  `Validation holdout` = "#C75A56",
  `Internal validation` = "#C75A56",
  `Zhang validation` = "#E08A7D",
  positive_up = "#C75A56",
  negative_up = "#67B7C4",
  not_selected = "grey82"
)

theme_gold <- function(base_size = 6.5) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 0.8),
      axis.text = element_text(color = "black", size = base_size - 0.4),
      axis.title = element_text(size = base_size),
      axis.line = element_line(linewidth = 0.22),
      axis.ticks = element_line(linewidth = 0.18),
      panel.border = element_rect(fill = NA, color = "grey72", linewidth = 0.22),
      strip.background = element_rect(fill = "grey97", color = "grey72", linewidth = 0.22),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      legend.key.size = unit(0.24, "cm"),
      legend.title = element_text(face = "bold", size = base_size - 0.4),
      legend.text = element_text(size = base_size - 0.7),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

save_gold_pdf <- function(plot, name) {
  dims <- gold_pages[[sub("[.]pdf$", "", name)]]
  save_pdf(plot, file.path(dir_map$figures, name), width = dims[["width"]], height = dims[["height"]])
}

empty_plot <- function(title = "No computed data") {
  ggplot() + annotate("text", 0, 0, label = title, size = 2.5) + theme_void()
}

format_dataset <- function(x) gsub("_", " ", x)

ring_data <- function(items) {
  bind_rows(lapply(seq_along(items), function(i) {
    df <- items[[i]]
    df$ring <- names(items)[i]
    df$ring_id <- i
    df
  })) %>%
    group_by(ring, ring_id) %>%
    arrange(category, .by_group = TRUE) %>%
    mutate(
      frac = n / sum(n),
      xmax = cumsum(frac),
      xmin = lag(xmax, default = 0),
      ymin = ring_id - 0.35,
      ymax = ring_id + 0.28
    ) %>%
    ungroup()
}

plot_ring <- function(df, title, center_label = NULL) {
  if (nrow(df) == 0) return(empty_plot(title))
  p <- ggplot(df, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, fill = category)) +
    geom_rect(color = "white", linewidth = 0.35) +
    coord_polar(theta = "x") +
    xlim(0, 1) +
    scale_y_continuous(limits = c(0.35, max(df$ymax) + 0.2)) +
    scale_fill_manual(values = gold_pal, na.value = "grey80") +
    labs(title = title, fill = NULL) +
    theme_void(base_size = 6.5) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 7.2),
      legend.position = "right",
      legend.text = element_text(size = 5.7),
      legend.key.size = unit(0.22, "cm"),
      plot.margin = ggplot2::margin(1, 1, 1, 1)
    )
  if (!is.null(center_label)) {
    p <- p + annotate("text", x = 0, y = 0.72, label = center_label, size = 2.4, fontface = "bold")
  }
  p
}

ring_panel <- function(df, view, title, center_label = NULL, show_legend = FALSE) {
  p <- plot_ring(df %>% filter(.data$view == view), title, center_label = center_label)
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

contiguous_track_segments <- function(track_df, track_levels) {
  track_df %>%
    arrange(track, x) %>%
    group_by(track) %>%
    mutate(run_id = cumsum(category != lag(category, default = dplyr::first(category)))) %>%
    group_by(track, run_id, category) %>%
    summarise(
      xmin = min(x) - 0.5,
      xmax = max(x) + 0.5,
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      track = factor(as.character(track), levels = track_levels),
      y = as.numeric(track),
      ymin = y - 0.36,
      ymax = y + 0.36,
      xmid = (xmin + xmax) / 2,
      label = ifelse(n >= 8, as.character(n), "")
    )
}

volcano_plot <- function(df, title, label_n = 8, label_genes = NULL) {
  plot_df <- df %>%
    mutate(
      minus_log10_p = -log10(pmax(P.Value, .Machine$double.xmin)),
      status = case_when(
        P.Value < 0.05 & logFC > 0.3 ~ "positive_up",
        P.Value < 0.05 & logFC < -0.3 ~ "negative_up",
        TRUE ~ "not_selected"
      )
    )
  if (is.null(label_genes)) {
    label_df <- plot_df %>% filter(status != "not_selected") %>% arrange(P.Value) %>% head(label_n)
  } else {
    label_df <- plot_df %>%
      filter(gene %in% label_genes) %>%
      mutate(label_order = match(gene, label_genes)) %>%
      arrange(label_order)
  }
  ggplot(plot_df, aes(logFC, minus_log10_p, color = status)) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey75", linewidth = 0.22) +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = 2, color = "grey75", linewidth = 0.22) +
    geom_point(size = 0.55, alpha = 0.72) +
    ggrepel::geom_text_repel(data = label_df, aes(label = gene), size = 1.8, max.overlaps = Inf, min.segment.length = 0) +
    scale_color_manual(values = gold_pal, guide = guide_legend(title = NULL, override.aes = list(size = 1.8))) +
    labs(title = title, x = "log2 fold change", y = "-log10(P)") +
    theme_gold() +
    theme(legend.position = "bottom")
}

forest_auc <- function(df, title, method_col = "method", dataset_filter = NULL) {
  plot_df <- df
  if (!is.null(dataset_filter)) plot_df <- plot_df %>% filter(dataset %in% dataset_filter)
  if (!method_col %in% names(plot_df)) method_col <- "model"
  plot_df <- plot_df %>%
    mutate(
      method_label = .data[[method_col]],
      method_label = stringr::str_replace_all(method_label, "_", " "),
      method_label = factor(method_label, levels = rev(unique(method_label)))
    )
  ggplot(plot_df, aes(x = auc, y = method_label, color = dataset)) +
    geom_segment(aes(x = ci_low, xend = ci_high, yend = method_label), color = "grey70", linewidth = 0.38) +
    geom_point(size = 1.45) +
    geom_vline(xintercept = 0.5, linetype = 2, color = "grey75", linewidth = 0.22) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0.5, 0.75, 1)) +
    scale_color_manual(values = gold_pal, na.value = "grey40") +
    labs(title = title, x = "AUC (95% CI)", y = NULL, color = NULL) +
    theme_gold() +
    theme(legend.position = "bottom")
}

roc_plot <- function(roc_df, perf_df, title) {
  lab <- perf_df %>%
    mutate(label = paste0(dataset, " AUC=", fmt_num(auc, 3))) %>%
    group_by(dataset) %>%
    summarise(label = dplyr::first(label), .groups = "drop") %>%
    mutate(fpr = 0.52, tpr = seq(0.25, 0.12, length.out = n()))
  ggplot(roc_df, aes(fpr, tpr, color = dataset)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey78", linewidth = 0.22) +
    geom_path(linewidth = 0.55) +
    geom_text(data = lab, aes(x = fpr, y = tpr, label = label), inherit.aes = FALSE, hjust = 0, size = 2.0) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_color_manual(values = gold_pal, na.value = "grey35") +
    labs(title = title, x = "1 - specificity", y = "Sensitivity", color = NULL) +
    theme_gold() +
    theme(legend.position = "bottom")
}

confusion_plot <- function(df, title) {
  if (nrow(df) == 0) return(empty_plot(title))
  if (all(c("Healthy", "LungCancer") %in% c(df$actual, df$predicted))) {
    actual_levels <- c("Healthy", "LungCancer")
    predicted_levels <- c("LungCancer", "Healthy")
  } else if (all(c("N0", "Nplus") %in% c(df$actual, df$predicted))) {
    actual_levels <- c("Nplus", "N0")
    predicted_levels <- c("Nplus", "N0")
  } else {
    actual_levels <- unique(df$actual)
    predicted_levels <- unique(df$predicted)
  }
  df <- df %>%
    mutate(
      actual = factor(actual, levels = actual_levels),
      predicted = factor(predicted, levels = predicted_levels)
    )
  ggplot(df, aes(actual, predicted, fill = n)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = n), size = 2.4) +
    facet_wrap(~dataset) +
    scale_fill_gradient(low = "#F3F4F5", high = "#C75A56") +
    labs(title = title, x = "Actual class", y = "Predicted class", fill = "n") +
    theme_gold()
}

make_expr_long <- function(expr_df, clinical, genes, samples, name, group_cols = c("group", "major_type", "histology", "n_status")) {
  genes <- intersect(genes, expr_df$genesymbol)
  samples <- intersect(samples, names(expr_df))
  if (length(genes) == 0 || length(samples) == 0) return(panel_data(name, data.frame()))
  mat <- as.matrix(expr_df[match(genes, expr_df$genesymbol), samples, drop = FALSE])
  rownames(mat) <- genes
  mode(mat) <- "numeric"
  z <- cap_value(row_z(mat), 2.5)
  data.frame(genesymbol = rownames(z), z, check.names = FALSE) %>%
    pivot_longer(cols = all_of(samples), names_to = "sample_id", values_to = "z_abundance") %>%
    left_join(clinical %>% dplyr::select(sample_id, any_of(group_cols)), by = "sample_id") %>%
    panel_data(name, .)
}

heatmap_plot <- function(df, title, sample_order = NULL, feature_order = NULL) {
  if (nrow(df) == 0) return(empty_plot(title))
  if (is.null(sample_order)) sample_order <- unique(df$sample_id)
  if (is.null(feature_order)) feature_order <- unique(df$genesymbol)
  df <- df %>%
    mutate(sample_id = factor(sample_id, levels = sample_order), genesymbol = factor(genesymbol, levels = rev(feature_order)))
  hm <- ggplot(df, aes(sample_id, genesymbol, fill = z_abundance)) +
    geom_tile() +
    scale_fill_gradient2(low = "#67B7C4", mid = "white", high = "#C75A56", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish) +
    labs(title = title, x = NULL, y = NULL, fill = "z") +
    theme_gold(base_size = 5.8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "right")
  ann_cols <- intersect(c("model_set", "dataset", "progression3", "group", "major_type", "histology", "n_status", "truth"), names(df))
  ann_cols <- ann_cols[vapply(ann_cols, function(col) dplyr::n_distinct(df[[col]][!is.na(df[[col]])]) > 1, logical(1))]
  ann_cols <- ann_cols[seq_len(min(length(ann_cols), 3))]
  if (length(ann_cols) == 0) return(hm)
  ann_df <- df %>%
    distinct(sample_id, across(all_of(ann_cols))) %>%
    pivot_longer(cols = all_of(ann_cols), names_to = "track", values_to = "category") %>%
    mutate(
      sample_id = factor(sample_id, levels = sample_order),
      category = ifelse(is.na(category) | category == "", "Control", as.character(category)),
      track = factor(track, levels = rev(ann_cols))
    )
  ann <- ggplot(ann_df, aes(sample_id, track, fill = category)) +
    geom_tile(height = 0.86) +
    scale_fill_manual(values = gold_pal, na.value = "grey80") +
    theme_void(base_size = 5.2) +
    theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(0, 2, 0, 2),
      axis.text.y = element_text(color = "black", size = 4.7),
      axis.text.x = element_blank()
    )
  ann / hm + plot_layout(heights = c(0.12 + 0.055 * length(ann_cols), 1))
}

module_boxplot <- function(df, xvar, title, fill_var = xvar, ncol = 3) {
  ggplot(df, aes(x = .data[[xvar]], y = module_score, fill = .data[[fill_var]])) +
    geom_boxplot(outlier.size = 0.28, linewidth = 0.18) +
    facet_wrap(~module, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = gold_pal, na.value = "grey75", guide = "none") +
    labs(title = title, x = NULL, y = "Module score") +
    theme_gold(base_size = 5.8) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

clinical <- read_result("processed_data", "clinical_processed.csv")
clinical$sample_id <- as.character(clinical$sample_id)
expr_df <- read_result("processed_data", "expression_processed.csv")
qc_sample <- read_result("processed_data", "QC_sample_metrics.csv")
qc_sample$sample_id <- as.character(qc_sample$sample_id)
qc_protein <- read_result("processed_data", "QC_protein_metrics.csv")
pca <- read_result("processed_data", "QC_pca_scores.csv")
pca$sample_id <- as.character(pca$sample_id)
disc_curve <- read_result("processed_data", "protein_discovery_curve_summary.csv")

gold_genes <- list(
  Figure3B_labels = c("ANXA11", "LYNX1", "HRG", "IGFBP2", "PDCD6", "ANXA7"),
  Figure3C_labels = c("FTH1", "ORM2", "LYVE1", "APOC2", "CFD", "APOA2"),
  Figure3D_heatmap = c(
    "SLC13A2", "ACE2", "ANXA5", "ATP1A1", "ATP1A3", "SLC4A4", "AQP1", "STOM",
    "S100A8", "TIMP1", "ANXA1", "LBP", "ORM2", "ORM1", "HP", "SERPINA3", "SAA1",
    "FTH1", "A2M", "APCS", "SERPINA1", "SAA2", "CRP", "FGB", "SERPINC1", "C4B",
    "FGG", "C4A", "CFH", "CFB", "C9", "FGA", "C3", "PLG", "CFI", "C1QB", "C1QC",
    "C2", "C5", "C6", "C7", "C8A", "C8B", "C8G"
  ),
  Figure3E_representative = c("SLC13A2", "ACE2", "ANXA5", "S100A8", "TIMP1", "ANXA1"),
  Figure4B_labels = c("FTH1", "LBP", "FTL", "LYVE1", "CTSZ", "FGB", "NDUFA9", "RTN4", "FGG"),
  Figure4D_LUAD_heatmap = c("NDUFA9", "LBP", "ORM1", "SLC25A5", "FTH1", "FTL", "DNPH1", "CTSZ", "CTSC", "LRG1", "LAMP1", "LYVE1", "FGB", "FGG", "RTN4"),
  Figure4E_LUAD_representative = c("FTH1", "LBP", "FTL", "FGB", "CTSZ", "FGG"),
  Figure4G_labels = c("HP", "GPC1", "LAMTOR1", "SLC35F6", "RAB32", "ORM2", "S100A7", "RAB17", "SH3BGRL3"),
  Figure4I_LUSC_heatmap = c("RAB17", "S100A7", "PTER", "HP", "ORM2", "CRP", "GPC1", "SPAG9", "LAMTOR1", "LAMTOR3", "SH3BGRL3", "AKR7A2", "CANX", "GANAB", "CETP"),
  Figure4J_LUSC_representative = c("HP", "LAMTOR1", "GPC1", "SLC35F6", "RAB32", "RAB17"),
  Figure5H_ModelB_features = c("CETP", "ORM1", "LBP", "ORM2", "C4B_2", "FGG", "SMIM22", "GLB1", "ATP5L", "FTH1", "HP", "APOC2", "HIST1H4K", "RBP5", "CP", "SLC5A8", "RTN4", "RAB7A"),
  Figure5J_representative = c("CETP", "ORM1", "LBP", "ORM2", "C4B_2", "FGG", "FTH1", "HP", "APOC2", "CP", "RTN4", "RAB7A"),
  Figure6DE_ModelC_features = c("HIST1H4A", "HIST1H4B", "HIST1H4C", "VDAC2", "VDAC3", "TMEM63A", "RAB7A", "NDRG1", "RPSA"),
  Supplementary2A_labels = c("LYNX1", "APOA2", "DCD", "IGFBP2", "CDH13", "KNG1", "PSAP", "SERPINF2", "PDCD6", "ANXA11"),
  Supplementary2B_shifts = c("SPP1", "ANXA11", "PTMA", "CAPN7", "CYSTM1", "PDCD6", "SAA1", "AQP1", "GPRC5C", "PROM1", "A1BG", "IGFBP2", "S100A8", "TMEM256", "AHSG", "B2M", "SPRR3", "DCD", "APOA2", "LYNX1"),
  Supplementary2C_heatmap = c(
    "ANXA11", "PDCD6", "PROM1", "CAPN7", "SLC5A12", "ANXA5", "CPNE8", "ATP1A3", "ACE2",
    "SLC6A19", "ANXA4", "ANXA7", "MAPK14", "GPRC5C", "UMOD", "ATP1A2", "MUC1", "PEF1",
    "PDCD6IP", "ROBO4", "DCD", "CDH13", "KNG1", "PSAP", "SULF2", "SERPINF2", "CD248",
    "B2M", "IGJ", "VASN", "MXRA8", "CLEC3B", "PCDH12", "NAPSA", "AHSG", "DSC2"
  )
)
keep_genes <- function(x) intersect(x, expr_df$genesymbol)
gold_gene_index <- bind_rows(lapply(names(gold_genes), function(panel) {
  data.frame(panel = panel, genesymbol = gold_genes[[panel]], present_in_expression = gold_genes[[panel]] %in% expr_df$genesymbol, check.names = FALSE)
}))
write_csv_safe(gold_gene_index, file.path(dir_map$metadata, "publication_display_gene_sets_used.csv"))

# Figure 1 --------------------------------------------------------------------
sample_order <- clinical %>%
  mutate(order_key = case_when(
    group == "Healthy" ~ 1,
    major_type == "NSCLC" & histology == "LUAD" & n_status == "N0" ~ 2,
    major_type == "NSCLC" & histology == "LUAD" & n_status == "Nplus" ~ 3,
    major_type == "NSCLC" & histology == "LUSC" & n_status == "N0" ~ 4,
    major_type == "NSCLC" & histology == "LUSC" & n_status == "Nplus" ~ 5,
    major_type == "SCLC" & n_status == "N0" ~ 6,
    major_type == "SCLC" & n_status == "Nplus" ~ 7,
    TRUE ~ 99
  )) %>%
  arrange(order_key, sample_id) %>%
  mutate(x = row_number())
track_df <- bind_rows(
  sample_order %>% transmute(sample_id, x, track = "Disease status", category = ifelse(group == "Healthy", "Healthy", "LungCancer")),
  sample_order %>% transmute(sample_id, x, track = "Histology", category = ifelse(group == "Healthy", "Control", as.character(histology))),
  sample_order %>% transmute(sample_id, x, track = "Nodal status", category = ifelse(is.na(n_status), "Control", as.character(n_status))),
  sample_order %>% transmute(sample_id, x, track = "Disease", category = ifelse(group == "Healthy", "Control", "LungCancer")),
  sample_order %>% transmute(sample_id, x, track = "Major type", category = as.character(major_type)),
  sample_order %>% transmute(sample_id, x, track = "NSCLC nodes", category = case_when(
    major_type == "NSCLC" & n_status == "N0" ~ "N0",
    major_type == "NSCLC" & n_status == "Nplus" ~ "Nplus",
    TRUE ~ "Control"
  ))
)
track_levels_f1 <- rev(c("Disease status", "Histology", "Nodal status", "Disease", "Major type", "NSCLC nodes"))
track_df$track <- factor(track_df$track, levels = track_levels_f1)
track_segments <- contiguous_track_segments(track_df, track_levels_f1)
sample_breaks <- sample_order %>%
  group_by(order_key) %>%
  summarise(xend = max(x) + 0.5, .groups = "drop") %>%
  filter(order_key < max(order_key))
f1a <- panel_data("Figure1A_corrected_study_cohort_tracks", track_segments)
p1a <- ggplot(f1a, aes(fill = category)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = "white", linewidth = 0.26) +
  geom_vline(data = sample_breaks, aes(xintercept = xend), inherit.aes = FALSE, color = "black", linewidth = 0.15) +
  geom_text(aes(x = xmid, y = y, label = label), size = 1.7, color = "black") +
  scale_y_continuous(breaks = seq_along(track_levels_f1), labels = track_levels_f1, expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = gold_pal, na.value = "grey80") +
  labs(title = "Corrected study cohort", x = NULL, y = NULL, fill = NULL) +
  theme_gold(base_size = 6.2) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  )

f1b <- panel_data("Figure1B_protein_detection_frequency_distribution", qc_protein)
p1b <- ggplot(f1b, aes(detection_fraction * 100)) +
  geom_histogram(bins = 44, fill = "#D9EEF3", color = "white", linewidth = 0.15) +
  geom_vline(xintercept = c(33, 75, 95), color = "#E6A09B", linewidth = 0.32) +
  labs(title = "Distribution of protein detection frequency", x = "Detection fraction across samples (%)", y = "Proteins") +
  theme_gold()

f1c <- panel_data("Figure1C_sample_level_protein_coverage", qc_sample)
p1c <- ggplot(f1c, aes(major_type, detected_protein_n, fill = major_type)) +
  geom_violin(width = 0.82, color = "grey35", linewidth = 0.18, alpha = 0.78) +
  geom_boxplot(width = 0.20, outlier.size = 0.25, linewidth = 0.18) +
  geom_jitter(width = 0.08, size = 0.35, alpha = 0.28) +
  scale_fill_manual(values = gold_pal, guide = "none") +
  labs(title = "Sample-level protein coverage", x = NULL, y = "Quantified proteins") +
  theme_gold()

f1d <- panel_data("Figure1D_pca_top1000_variable_proteins", pca %>% mutate(status = ifelse(is.na(n_status), "Healthy", as.character(n_status))))
p1d <- ggplot(f1d, aes(PC1, PC2, color = major_type)) +
  stat_ellipse(aes(fill = major_type), geom = "polygon", alpha = 0.08, color = NA) +
  stat_ellipse(linewidth = 0.3, alpha = 0.8) +
  geom_point(size = 0.9, alpha = 0.80) +
  scale_color_manual(values = gold_pal) +
  scale_fill_manual(values = gold_pal, guide = "none") +
  labs(title = "PCA of urinary proteome profiles", x = paste0("PC1 (", fmt_num(unique(f1d$PC1_var)[1], 1), "%)"), y = paste0("PC2 (", fmt_num(unique(f1d$PC2_var)[1], 1), "%)")) +
  theme_gold()

f1e <- panel_data("Figure1E_cumulative_protein_discovery_curve", disc_curve)
p1e <- ggplot(f1e, aes(sample_n, mean_cumulative, color = group, fill = group)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.55) +
  scale_color_manual(values = gold_pal) +
  scale_fill_manual(values = gold_pal) +
  labs(title = "Cumulative protein discovery curve", x = "Number of sampled individuals", y = "Cumulative proteins") +
  theme_gold() +
  theme(legend.position = "bottom")

fig1 <- (p1a / (p1b | p1c) / (p1d | p1e)) +
  plot_layout(heights = c(1.0, 1.45, 1.45)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 7))
save_gold_pdf(fig1, "Figure1.pdf")

# Figure 2 --------------------------------------------------------------------
ma_split <- read_result("modelA", "modelA_sample_split.csv")
ma_perf <- read_result("modelA", "modelA_final_performance.csv")
ma_roc <- read_result("modelA", "modelA_final_roc_coordinates.csv")
ma_conf <- read_result("modelA", "modelA_final_confusion_matrices.csv")
ma_scores <- read_result("modelA", "modelA_final_prediction_scores.csv")
ma_cand <- read_result("modelA", "modelA_candidate_model_performance.csv")
ma_features <- read_result("modelA", "modelA_feature_table_direction_weights.csv")
ma_cv <- read_result("modelA", "modelA_rf_feature_set_cross_validation.csv")
ma_heat <- read_result("modelA", "modelA_direction_adjusted_heatmap_data.csv")
zhang_clin <- read_result("processed_data", "zhang_clinical_LC_CTL.csv")
zhang_clin$sample_id <- as.character(zhang_clin$sample_id)

modela_ring <- function(df, view_label) {
  ring_data(list(
    Status = df %>% count(category = class, name = "n"),
    Histology = df %>% mutate(category = ifelse(class == "Healthy", "Control", as.character(histology))) %>% count(category, name = "n"),
    Nodal = df %>% mutate(category = ifelse(is.na(n_status), "Control", as.character(n_status))) %>% count(category, name = "n")
  )) %>% mutate(view = view_label)
}
f2a <- panel_data("Figure2A_modelA_cohort_overview_ring", bind_rows(
  modela_ring(ma_split %>% filter(model_set == "Discovery"), "Discovery"),
  modela_ring(ma_split %>% filter(model_set == "Internal validation"), "Internal validation"),
  ring_data(list(
    Status = zhang_clin %>% count(category = group, name = "n"),
    Dataset = zhang_clin %>% count(category = dataset, name = "n")
  )) %>% mutate(view = "Zhang validation")
))
p2a <- patchwork::wrap_plots(
  ring_panel(f2a, "Discovery", "Discovery set", center_label = paste0("n=", sum(ma_split$model_set == "Discovery")), show_legend = TRUE),
  ring_panel(f2a, "Internal validation", "Internal validation", center_label = paste0("n=", sum(ma_split$model_set == "Internal validation"))) /
    ring_panel(f2a, "Zhang validation", "Zhang validation", center_label = paste0("n=", nrow(zhang_clin))),
  ncol = 2,
  widths = c(1.45, 0.78)
) + plot_annotation(title = "Model A cohort overview") &
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 7.2))
f2b <- panel_data("Figure2B_discovery_set_algorithm_comparison", ma_cand %>% filter(dataset == "Discovery"))
p2b <- forest_auc(f2b, "Discovery-set algorithm comparison", dataset_filter = "Discovery")
f2c <- panel_data("Figure2C_discovery_random_forest_feature_set_comparison", ma_cv %>% mutate(dataset = "Discovery CV", method = feature_set, auc = mean_cv_auc, ci_low = low_2.5, ci_high = high_97.5))
p2c <- forest_auc(f2c, "Discovery random-forest feature-set comparison")
f2d <- panel_data("Figure2D_discovery_signed_five_protein_panel", ma_features %>% mutate(dataset = "Discovery", signed_contribution = discovery_mean_lung_cancer - discovery_mean_healthy))
p2d <- ggplot(f2d, aes(reorder(genesymbol, signed_contribution), signed_contribution, fill = direction_label)) +
  geom_col(color = "grey35", linewidth = 0.15) +
  coord_flip() +
  scale_fill_manual(values = c(`Higher in lung cancer` = gold_pal[["LungCancer"]], `Lower in lung cancer` = gold_pal[["Healthy"]]), guide = guide_legend(title = NULL)) +
  labs(title = "Discovery signed five-protein panel", x = NULL, y = "Discovery mean difference") +
  theme_gold()
f2e <- panel_data("Figure2E_five_protein_score_distribution", ma_scores)
p2e <- ggplot(f2e, aes(truth, modelA_score, fill = truth)) +
  geom_hline(yintercept = unique(ma_perf$threshold[ma_perf$dataset == "Discovery"])[1], linetype = 2, linewidth = 0.22) +
  geom_violin(color = "grey35", linewidth = 0.18, alpha = 0.78) +
  geom_boxplot(width = 0.18, outlier.size = 0.25, linewidth = 0.18) +
  facet_wrap(~dataset, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = gold_pal, guide = "none") +
  labs(title = "Five-protein score distributions", x = NULL, y = "Score") +
  theme_gold(base_size = 5.8)
f2f <- panel_data("Figure2F_final_modelA_auc_by_cohort", ma_perf %>% mutate(method = dataset))
p2f <- forest_auc(f2f, "AUC by cohort")
f2g <- panel_data("Figure2G_internal_validation_roc", ma_roc %>% filter(dataset == "Internal validation"))
p2g <- roc_plot(f2g, ma_perf %>% filter(dataset == "Internal validation"), "Internal ROC")
f2h <- panel_data("Figure2H_internal_validation_confusion", ma_conf %>% filter(dataset == "Internal validation"))
p2h <- confusion_plot(f2h, "Internal confusion")
f2i <- panel_data("Figure2I_zhang_validation_roc", ma_roc %>% filter(dataset == "Zhang validation"))
p2i <- roc_plot(f2i, ma_perf %>% filter(dataset == "Zhang validation"), "Zhang ROC")
f2j <- panel_data("Figure2J_zhang_validation_confusion", ma_conf %>% filter(dataset == "Zhang validation"))
p2j <- confusion_plot(f2j, "Zhang confusion")
f2k <- panel_data("Figure2K_modelA_direction_adjusted_heatmap", ma_heat %>% dplyr::rename(z_abundance = direction_adjusted_z))
p2k <- (
  patchwork::wrap_elements(full = heatmap_plot(f2k %>% filter(dataset == "Discovery"), "Discovery score-adjusted marker heatmap", sample_order = unique((f2k %>% filter(dataset == "Discovery"))$sample_id), feature_order = unique(f2k$genesymbol))) /
    patchwork::wrap_elements(full = heatmap_plot(f2k %>% filter(dataset == "Internal validation"), "Internal validation", sample_order = unique((f2k %>% filter(dataset == "Internal validation"))$sample_id), feature_order = unique(f2k$genesymbol))) /
    patchwork::wrap_elements(full = heatmap_plot(f2k %>% filter(dataset == "Zhang validation"), "Zhang independent", sample_order = unique((f2k %>% filter(dataset == "Zhang validation"))$sample_id), feature_order = unique(f2k$genesymbol)))
) + plot_layout(heights = c(1, 0.85, 0.85))

fig2 <- (p2a | p2b) / (p2c | p2d) / p2e / (p2f | p2g | p2h) / (p2i | p2j | plot_spacer()) / p2k +
  plot_layout(heights = c(0.95, 0.88, 0.78, 0.9, 0.9, 1.65)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 7))
save_gold_pdf(fig2, "Figure2.pdf")

# Figure 3 --------------------------------------------------------------------
nsclc_de <- read_result("differential_analysis", "NSCLC_Nplus_vs_N0_limma_all_proteins.csv")
n0_de <- read_result("differential_analysis", "NSCLC_N0_vs_Healthy_limma_all_proteins.csv")
prog_scores <- read_result("differential_analysis", "progression_curated_module_scores.csv")
nsclc_mod_stats <- read_result("differential_analysis", "NSCLC_curated_module_statistics.csv")
prog_clin <- clinical %>%
  mutate(progression3 = case_when(
    group == "Healthy" ~ "Healthy",
    major_type == "NSCLC" & n_status == "N0" ~ "NSCLC_N0",
    major_type == "NSCLC" & n_status == "Nplus" ~ "NSCLC_Nplus",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(progression3))
f3a <- panel_data("Figure3A_NSCLC_progression_cohort_ring", ring_data(list(
  Status = prog_clin %>% count(category = progression3, name = "n"),
  Histology = prog_clin %>% mutate(category = ifelse(group == "Healthy", "Control", as.character(histology))) %>% count(category, name = "n"),
  Nodal = prog_clin %>% mutate(category = ifelse(is.na(n_status), "Control", as.character(n_status))) %>% count(category, name = "n")
)))
p3a <- plot_ring(f3a, "NSCLC progression cohort overview", center_label = paste0("n=", nrow(prog_clin)))
f3b <- panel_data("Figure3B_NSCLC_N0_vs_Healthy_volcano", n0_de)
p3b <- volcano_plot(f3b, "NSCLC N0 vs Healthy", label_genes = gold_genes$Figure3B_labels)
f3c <- panel_data("Figure3C_NSCLC_Nplus_vs_N0_volcano", nsclc_de)
p3c <- volcano_plot(f3c, "NSCLC Nplus vs N0", label_genes = gold_genes$Figure3C_labels)
top_nsclc_genes <- keep_genes(gold_genes$Figure3D_heatmap)
f3d <- make_expr_long(expr_df, prog_clin, top_nsclc_genes, prog_clin$sample_id, "Figure3D_selected_protein_groups_across_NSCLC_progression", group_cols = c("progression3"))
p3d <- patchwork::wrap_elements(full = heatmap_plot(f3d, "Selected protein groups across NSCLC progression", sample_order = unique(f3d$sample_id[order(f3d$progression3)]), feature_order = top_nsclc_genes))
rep_genes3 <- keep_genes(gold_genes$Figure3E_representative)
f3e <- make_expr_long(expr_df, prog_clin, rep_genes3, prog_clin$sample_id, "Figure3E_representative_monotonic_protein_trends", group_cols = c("progression3", "group", "major_type", "histology", "n_status"))
p3e <- ggplot(f3e, aes(progression3, z_abundance, fill = progression3)) +
  geom_boxplot(outlier.size = 0.25, linewidth = 0.18) +
  geom_jitter(width = 0.08, size = 0.25, alpha = 0.22) +
  facet_wrap(~factor(genesymbol, levels = rep_genes3), scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(Healthy = gold_pal[["Healthy"]], NSCLC_N0 = gold_pal[["N0"]], NSCLC_Nplus = gold_pal[["Nplus"]]), guide = "none") +
  labs(title = "Representative monotonic protein trends", x = NULL, y = "z abundance") +
  theme_gold(base_size = 5.5) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
	stage_levels3 <- c("Healthy", "NSCLC_N0", "NSCLC_Nplus")
	stage_labels3 <- c(Healthy = "Healthy", NSCLC_N0 = "NSCLC N0", NSCLC_Nplus = "NSCLC N+")
	fig3_expr <- as.matrix(expr_df[, setdiff(names(expr_df), "genesymbol"), drop = FALSE])
	rownames(fig3_expr) <- expr_df$genesymbol
	mode(fig3_expr) <- "numeric"
	fig3_z <- row_z(fig3_expr[, prog_clin$sample_id, drop = FALSE])
	fig3_signature_sets <- list(
	  "Gradual up" = c("SLC13A2", "ACE2", "ANXA5", "ATP1A1", "ATP1A3", "SLC4A4", "AQP1", "STOM"),
	  "Gradual down" = c("S100A8", "TIMP1", "ANXA1"),
	  "Acute phase" = c("LBP", "ORM2", "ORM1", "HP", "SERPINA3", "SAA1", "ITIH4", "A2M", "APCS", "SERPINA1", "SAA2", "CRP"),
	  "Complement/coag." = c("FGB", "SERPINC1", "C4B", "FGG", "C4A", "CFH", "CFB", "C9", "FGA", "C3", "PLG", "CFI", "C1QB", "C1QC", "C2", "C5", "C6", "C7", "C8A", "C8B", "C8G")
	)
	fig3_signature_labels <- c(
	  "Gradual up" = "R-HSA-425407 Tubular transport\nR-HSA-1474244 Adhesion/matrix",
	  "Gradual down" = "R-HSA-1474244 Adhesion/matrix\nR-HSA-6798695 Myeloid/neutrophil",
	  "Acute phase" = "GO:0006953 Acute phase",
	  "Complement/coag." = "R-HSA-166658 Complement/coag."
	)
	fig3_signature_types <- c(
	  "Gradual up" = "Monotonic signal",
	  "Gradual down" = "Monotonic signal",
	  "Acute phase" = "Prior context module",
	  "Complement/coag." = "Prior context module"
	)
	grouped_prog_scores <- bind_rows(lapply(names(fig3_signature_sets), function(sig) {
	  genes <- intersect(fig3_signature_sets[[sig]], rownames(fig3_z))
	  if (length(genes) == 0) stop("No Figure3F proteins found for ", sig)
	  data.frame(
	    sample_id = prog_clin$sample_id,
	    panel_group = sig,
	    panel_type = unname(fig3_signature_types[sig]),
	    signature_label = unname(fig3_signature_labels[sig]),
	    display_group = unname(fig3_signature_labels[sig]),
	    n_proteins = length(genes),
	    proteins = paste(genes, collapse = ";"),
	    module_score = colMeans(fig3_z[genes, prog_clin$sample_id, drop = FALSE], na.rm = TRUE),
	    check.names = FALSE
	  )
	})) %>%
	  left_join(prog_clin %>% dplyr::select(sample_id, progression3), by = "sample_id") %>%
	  mutate(
	    progression_stage = unname(stage_labels3[progression3]),
	    progression_stage = factor(progression_stage, levels = unname(stage_labels3[stage_levels3]))
	  )
	module_group_levels3 <- unname(fig3_signature_labels[names(fig3_signature_sets)])
	module_stage_effects <- grouped_prog_scores %>%
	  group_by(panel_group, panel_type, signature_label, display_group, progression_stage, n_proteins, proteins) %>%
	  summarise(median_score = median(module_score, na.rm = TRUE), .groups = "drop") %>%
	  mutate(
	    stage_index = as.numeric(progression_stage),
	    display_group = factor(display_group, levels = rev(module_group_levels3))
	  )
	f3f <- panel_data("Figure3F_three_stage_protein_group_scores", module_stage_effects)
	p3f <- ggplot(f3f, aes(x = progression_stage, y = display_group, fill = median_score)) +
	  geom_tile(color = "white", linewidth = 0.3, width = 0.72, height = 0.72) +
	  geom_text(aes(label = fmt_num(median_score, 2)), size = 2.0) +
	  scale_fill_gradient2(low = gold_pal[["Healthy"]], mid = "white", high = gold_pal[["LungCancer"]], midpoint = 0) +
	  labs(title = "Three-stage protein-group scores", subtitle = "Protein groups introduced in panel D", x = NULL, y = NULL, fill = "Median\nscore") +
	  theme_gold(base_size = 5.8) +
	  theme(axis.text.x = element_text(angle = 25, hjust = 1), plot.subtitle = element_text(size = 5.0))
	f3g <- panel_data("Figure3G_sample_level_protein_group_scores", grouped_prog_scores %>% mutate(display_group = factor(display_group, levels = module_group_levels3)))
	p3g <- ggplot(f3g, aes(progression_stage, module_score, fill = progression_stage)) +
	  geom_boxplot(outlier.size = 0.25, linewidth = 0.18) +
	  geom_jitter(width = 0.08, size = 0.25, alpha = 0.22) +
	  stat_summary(aes(group = 1), fun = median, geom = "line", color = "black", linewidth = 0.25) +
	  stat_summary(fun = median, geom = "point", color = "black", size = 0.65) +
	  facet_wrap(~display_group, scales = "free_y", ncol = 2) +
	  scale_fill_manual(values = c("Healthy" = gold_pal[["Healthy"]], "NSCLC N0" = gold_pal[["N0"]], "NSCLC N+" = gold_pal[["Nplus"]]), guide = "none") +
	  labs(title = "Sample-level protein-group scores", x = NULL, y = "Mean row z-score group score") +
	  theme_gold(base_size = 5.5) +
	  theme(axis.text.x = element_text(angle = 25, hjust = 1))

	six_module_scores <- prog_scores %>%
	  left_join(module_metadata, by = "module") %>%
	  mutate(
	    progression_stage = unname(stage_labels3[progression3]),
	    progression_stage = factor(progression_stage, levels = unname(stage_labels3[stage_levels3])),
	    module_label = paste0(pathway_id, "\n", module_short)
	  )
	six_module_levels <- paste0(module_metadata$pathway_id, "\n", module_metadata$module_short)
	six_module_stage <- six_module_scores %>%
	  group_by(module, module_short, pathway_source, pathway_id, pathway_name, module_label, progression_stage, genes_used, genes) %>%
	  summarise(median_score = median(module_score, na.rm = TRUE), .groups = "drop") %>%
	  group_by(module) %>%
	  mutate(median_z = as.numeric(scale(median_score))) %>%
	  ungroup() %>%
	  mutate(module_label = factor(module_label, levels = rev(six_module_levels)))
	supp3f_full <- panel_data("Supplementary_Figure_3F_full_six_module_stage_scores", six_module_stage)
	ps3f_full <- ggplot(supp3f_full, aes(x = progression_stage, y = module_label, fill = median_z)) +
	  geom_tile(color = "white", linewidth = 0.3, width = 0.74, height = 0.72) +
	  geom_text(aes(label = fmt_num(median_z, 2)), size = 2.0) +
	  scale_fill_gradient2(low = gold_pal[["Healthy"]], mid = "white", high = gold_pal[["LungCancer"]], midpoint = 0) +
	  labs(title = "Full six curated pathway modules across NSCLC progression", subtitle = "Labels show module-wise scaled stage medians; raw median scores are saved in plotted data", x = NULL, y = NULL, fill = "Median z") +
	  theme_gold(base_size = 6.0) +
	  theme(axis.text.x = element_text(angle = 25, hjust = 1), plot.subtitle = element_text(size = 5.0))
	save_pdf(ps3f_full, file.path(dir_map$figures, "Supplementary_Figure_3F_full_six_module_scores.pdf"), width = 6.4, height = 4.0)

fig3 <- (p3a | p3b | p3c) / p3d / (p3e | p3f) / p3g +
  plot_layout(heights = c(0.9, 1.65, 1.25, 1.4)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 7))
save_gold_pdf(fig3, "Figure3.pdf")

# Figure 4 --------------------------------------------------------------------
luad_de <- read_result("differential_analysis", "LUAD_Nplus_vs_N0_limma_all_proteins.csv")
lusc_de <- read_result("differential_analysis", "LUSC_Nplus_vs_N0_limma_all_proteins.csv")
effect_compare <- read_result("differential_analysis", "subtype_vs_nsclc_effect_size_comparison.csv")
module_effects <- read_result("differential_analysis", "pathway_module_adjusted_effects_by_subtype.csv")
luad_clin <- clinical %>% filter(histology == "LUAD", n_status %in% c("N0", "Nplus"))
lusc_clin <- clinical %>% filter(histology == "LUSC", n_status %in% c("N0", "Nplus"))
f4a <- panel_data("Figure4A_LUAD_nodal_set_ring", ring_data(list(
  LUAD = luad_clin %>% count(category = n_status, name = "n"),
  Status = luad_clin %>% mutate(category = "LUAD") %>% count(category, name = "n")
)))
p4a <- plot_ring(f4a, "LUAD nodal set", center_label = paste0("n=", nrow(luad_clin)))
f4b <- panel_data("Figure4B_LUAD_Nplus_vs_N0_volcano", luad_de)
p4b <- volcano_plot(f4b, "LUAD Nplus vs N0", label_genes = gold_genes$Figure4B_labels)
f4c <- panel_data("Figure4C_LUAD_vs_NSCLC_nodal_effects", effect_compare)
p4c <- ggplot(f4c, aes(nsclc_logFC, luad_logFC, color = selected_luad_nsclc)) +
  geom_hline(yintercept = 0, color = "grey78", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey78", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey65", linewidth = 0.22) +
  geom_point(size = 0.55, alpha = 0.62) +
  scale_color_manual(values = c(`TRUE` = gold_pal[["LungCancer"]], `FALSE` = "grey78"), guide = "none") +
  labs(title = "LUAD vs NSCLC nodal effects", x = "NSCLC logFC", y = "LUAD logFC") +
  theme_gold()
luad_top <- keep_genes(gold_genes$Figure4D_LUAD_heatmap)
f4d <- make_expr_long(expr_df, clinical, luad_top, luad_clin$sample_id, "Figure4D_top_LUAD_nodal_proteins", group_cols = c("histology", "n_status"))
p4d <- patchwork::wrap_elements(full = heatmap_plot(f4d, "Top LUAD nodal proteins", sample_order = unique(f4d$sample_id[order(f4d$n_status)]), feature_order = luad_top))
rep_luad <- keep_genes(gold_genes$Figure4E_LUAD_representative)
f4e <- make_expr_long(expr_df, clinical, rep_luad, luad_clin$sample_id, "Figure4E_LUAD_representative_proteins", group_cols = c("histology", "n_status"))
p4e <- ggplot(f4e, aes(n_status, z_abundance, fill = n_status)) +
  geom_violin(color = "grey35", linewidth = 0.16, alpha = 0.75) +
  geom_boxplot(width = 0.18, outlier.size = 0.2, linewidth = 0.16) +
  facet_wrap(~factor(genesymbol, levels = rep_luad), scales = "free_y", ncol = 3) +
  scale_fill_manual(values = gold_pal, guide = "none") +
  labs(title = "LUAD representative proteins", x = NULL, y = "z abundance") +
  theme_gold(base_size = 5.6)
f4f <- panel_data("Figure4F_LUSC_nodal_set_ring", ring_data(list(
  LUSC = lusc_clin %>% count(category = n_status, name = "n"),
  Status = lusc_clin %>% mutate(category = "LUSC") %>% count(category, name = "n")
)))
p4f <- plot_ring(f4f, "LUSC nodal set", center_label = paste0("n=", nrow(lusc_clin)))
f4g <- panel_data("Figure4G_LUSC_Nplus_vs_N0_volcano", lusc_de)
p4g <- volcano_plot(f4g, "LUSC Nplus vs N0", label_genes = gold_genes$Figure4G_labels)
f4h <- panel_data("Figure4H_LUSC_vs_NSCLC_nodal_effects", effect_compare)
p4h <- ggplot(f4h, aes(nsclc_logFC, lusc_logFC, color = selected_lusc_nsclc)) +
  geom_hline(yintercept = 0, color = "grey78", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey78", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey65", linewidth = 0.22) +
  geom_point(size = 0.55, alpha = 0.62) +
  scale_color_manual(values = c(`TRUE` = gold_pal[["LungCancer"]], `FALSE` = "grey78"), guide = "none") +
  labs(title = "LUSC vs NSCLC nodal effects", x = "NSCLC logFC", y = "LUSC logFC") +
  theme_gold()
lusc_top <- keep_genes(gold_genes$Figure4I_LUSC_heatmap)
f4i <- make_expr_long(expr_df, clinical, lusc_top, lusc_clin$sample_id, "Figure4I_top_LUSC_nodal_proteins", group_cols = c("histology", "n_status"))
p4i <- patchwork::wrap_elements(full = heatmap_plot(f4i, "Top LUSC nodal proteins", sample_order = unique(f4i$sample_id[order(f4i$n_status)]), feature_order = lusc_top))
rep_lusc <- keep_genes(gold_genes$Figure4J_LUSC_representative)
f4j <- make_expr_long(expr_df, clinical, rep_lusc, lusc_clin$sample_id, "Figure4J_LUSC_representative_proteins", group_cols = c("histology", "n_status"))
p4j <- ggplot(f4j, aes(n_status, z_abundance, fill = n_status)) +
  geom_violin(color = "grey35", linewidth = 0.16, alpha = 0.75) +
  geom_boxplot(width = 0.18, outlier.size = 0.2, linewidth = 0.16) +
  facet_wrap(~factor(genesymbol, levels = rep_lusc), scales = "free_y", ncol = 3) +
  scale_fill_manual(values = gold_pal, guide = "none") +
  labs(title = "LUSC representative proteins", x = NULL, y = "z abundance") +
  theme_gold(base_size = 5.6)
f4k <- panel_data("Figure4K_pathway_module_effects_by_subtype", module_effects)
p4k <- ggplot(f4k, aes(adjusted_beta, reorder(module, adjusted_beta), color = subtype, size = -log10(p_value))) +
  geom_vline(xintercept = 0, color = "grey78", linewidth = 0.22) +
  geom_point(alpha = 0.86) +
  scale_color_manual(values = c(NSCLC = gold_pal[["NSCLC"]], LUAD = gold_pal[["LUAD"]], LUSC = gold_pal[["LUSC"]])) +
  labs(title = "Pathway module effects by subtype", x = "Adjusted beta for Nplus vs N0", y = NULL, size = "-log10(P)", color = NULL) +
  theme_gold(base_size = 5.8) +
  theme(legend.position = "bottom")

fig4 <- (p4a | p4b | p4c) / (p4d | p4e) / (p4f | p4g | p4h) / (p4i | p4j) / p4k +
  plot_layout(heights = c(0.9, 1.15, 0.9, 1.1, 0.75)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 7))
save_gold_pdf(fig4, "Figure4.pdf")

# Figure 5 --------------------------------------------------------------------
mb_split <- read_result("modelB", "modelB_sample_split.csv")
mb_perf <- read_result("modelB", "modelB_final_performance.csv")
mb_roc <- read_result("modelB", "modelB_final_roc_coordinates.csv")
mb_conf <- read_result("modelB", "modelB_final_confusion_matrices.csv")
mb_cand <- read_result("modelB", "modelB_candidate_model_performance.csv")
mb_repeat <- read_result("modelB", "modelB_repeated_split_metrics.csv")
mb_disc <- read_result("modelB", "modelB_discovery_oob_predictions.csv") %>% transmute(sample_id, dataset = "Discovery OOB", truth, score = rf_top18_oob_prob, predicted)
mb_val <- read_result("modelB", "modelB_validation_predictions.csv") %>% transmute(sample_id, dataset = "Validation holdout", truth, score = rf_top18_prob, predicted)
mb_score <- bind_rows(mb_disc, mb_val)
mb_heat <- read_result("modelB", "modelB_top18_heatmap_data.csv")
mb_feat <- read_result("modelB", "modelB_feature_table_top18.csv")
modelb_ring <- function(df, view_label) {
  ring_data(list(
    Nodal = df %>% count(category = truth, name = "n"),
    Histology = df %>% count(category = histology, name = "n"),
    Sex = df %>% count(category = sex, name = "n")
  )) %>% mutate(view = view_label)
}
f5a <- panel_data("Figure5A_modelB_discovery_cohort_ring", modelb_ring(mb_split %>% filter(model_set == "Discovery"), "Discovery"))
p5a <- plot_ring(f5a, "Discovery cohort", center_label = paste0("n=", sum(mb_split$model_set == "Discovery")))
f5b <- panel_data("Figure5B_discovery_auc_comparison", mb_cand %>% filter(dataset == "Discovery"))
p5b <- forest_auc(f5b, "Discovery AUC comparison", dataset_filter = "Discovery")
f5c <- panel_data("Figure5C_validation_cohort_ring", modelb_ring(mb_split %>% filter(model_set == "Validation"), "Validation"))
p5c <- plot_ring(f5c, "Validation cohort", center_label = paste0("n=", sum(mb_split$model_set == "Validation")))
f5d <- panel_data("Figure5D_repeated_split_AUC_distribution", mb_repeat)
p5d <- ggplot(f5d, aes(auc)) +
  geom_histogram(bins = 22, fill = gold_pal[["LungCancer"]], color = "white", linewidth = 0.18, alpha = 0.82) +
  geom_vline(xintercept = unique(mb_perf$auc[mb_perf$dataset == "Validation holdout"])[1], linetype = 2, linewidth = 0.32) +
  labs(title = "Repeated split robustness", x = "AUC", y = "Frequency") +
  theme_gold()
f5e <- panel_data("Figure5E_validation_risk_score_distribution", mb_score %>% filter(dataset == "Validation holdout"))
p5e <- ggplot(f5e, aes(truth, score, fill = truth)) +
  geom_hline(yintercept = unique(mb_perf$threshold)[1], linetype = 2, linewidth = 0.22) +
  geom_violin(color = "grey35", linewidth = 0.18, alpha = 0.75) +
  geom_boxplot(width = 0.18, outlier.size = 0.25, linewidth = 0.18) +
  scale_fill_manual(values = gold_pal, guide = "none") +
  labs(title = "Validation risk-score distribution", x = NULL, y = "Nplus probability") +
  theme_gold()
f5f <- panel_data("Figure5F_validation_ROC", mb_roc %>% filter(dataset == "Validation holdout"))
p5f <- roc_plot(f5f, mb_perf %>% filter(dataset == "Validation holdout"), "Validation ROC")
f5g <- panel_data("Figure5G_validation_confusion", mb_conf %>% filter(dataset == "Validation holdout"))
p5g <- confusion_plot(f5g, "Validation confusion")
modelb_gold_features <- keep_genes(gold_genes$Figure5H_ModelB_features)
f5h <- make_expr_long(
  expr_df,
  mb_split,
  modelb_gold_features,
  mb_split$sample_id,
  "Figure5H_discovery_validation_contrast_heatmaps",
  group_cols = c("model_set", "truth")
)
p5h <- (
  patchwork::wrap_elements(full = heatmap_plot(f5h %>% filter(model_set == "Discovery"), "Discovery cohort Model B protein features", sample_order = unique((f5h %>% filter(model_set == "Discovery"))$sample_id[order((f5h %>% filter(model_set == "Discovery"))$truth)]), feature_order = modelb_gold_features)) |
    patchwork::wrap_elements(full = heatmap_plot(f5h %>% filter(model_set == "Validation"), "Validation cohort Model B protein features", sample_order = unique((f5h %>% filter(model_set == "Validation"))$sample_id[order((f5h %>% filter(model_set == "Validation"))$truth)]), feature_order = modelb_gold_features))
)
top9_mb <- keep_genes(gold_genes$Figure5J_representative)
f5j <- make_expr_long(expr_df, clinical, top9_mb, mb_split$sample_id, "Figure5J_representative_modelB_proteins", group_cols = c("n_status", "histology"))
p5j <- ggplot(f5j, aes(n_status, z_abundance, fill = n_status)) +
  geom_violin(color = "grey35", linewidth = 0.15, alpha = 0.75) +
  geom_boxplot(width = 0.18, outlier.size = 0.18, linewidth = 0.15) +
  facet_wrap(~factor(genesymbol, levels = top9_mb), scales = "free_y", ncol = 6) +
  scale_fill_manual(values = gold_pal, guide = "none") +
  labs(title = "Representative Model B proteins", x = NULL, y = "z abundance") +
  theme_gold(base_size = 5.4)
f5k <- panel_data("Figure5K_annotated_modelB_feature_biology", data.frame(
  biology_annotation = c("Vesicle/lysosome", "Acute phase", "Lipid/vitamin", "Iron/heme", "Complement", "Other", "Mitochondrial", "Chromatin"),
  protein_features = c(4, 4, 3, 2, 2, 1, 1, 1),
  genes = c("GLB1, SLC5A8, RTN4, RAB7A", "ORM1, LBP, ORM2, HP", "CETP, APOC2, RBP5", "FTH1, CP", "C4B_2, FGG", "SMIM22", "ATP5L", "HIST1H4K"),
  check.names = FALSE
))
p5k <- ggplot(f5k, aes(protein_features, reorder(biology_annotation, protein_features), fill = biology_annotation)) +
  geom_col(color = "grey35", linewidth = 0.15, width = 0.62, fill = "#6E58A8") +
  geom_text(aes(label = protein_features), hjust = -0.25, size = 2.1) +
  geom_text(aes(x = protein_features + 0.65, label = genes), hjust = 0, size = 2.0) +
  scale_x_continuous(limits = c(0, 11), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Annotated Model B protein-feature biology", x = "Protein features", y = NULL) +
  theme_gold(base_size = 5.5) +
  theme(legend.position = "none")
fig5 <- (p5a | p5b | p5c) / (p5d | p5e | p5f | p5g) / p5h / p5j / p5k +
  plot_layout(heights = c(0.82, 0.78, 1.2, 0.9, 0.55)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 7))
save_gold_pdf(fig5, "Figure5.pdf")

# Figure 6 --------------------------------------------------------------------
mc_split <- read_result("modelC", "modelC_sample_split.csv")
mc_perf <- read_result("modelC", "modelC_final_performance.csv")
mc_roc <- read_result("modelC", "modelC_final_roc_coordinates.csv")
mc_conf <- read_result("modelC", "modelC_final_confusion_matrices.csv")
mc_oof <- read_result("modelC", "modelC_discovery_5fold_oof_auc.csv")
mc_scores <- read_result("modelC", "modelC_final_prediction_scores.csv")
mc_heat <- read_result("modelC", "modelC_top9_heatmap_data.csv")
mc_feat <- read_result("modelC", "modelC_feature_table.csv")
modelc_ring <- function(df, view_label) {
  ring_data(list(
    Nodal = df %>% count(category = truth, name = "n"),
    Sex = df %>% count(category = sex, name = "n"),
    Age = df %>% mutate(category = ifelse(age <= 60, "Age <=60", "Age >60")) %>% count(category, name = "n")
  )) %>% mutate(view = view_label)
}
f6a <- panel_data("Figure6A_SCLC_discovery_cohort_ring", modelc_ring(mc_split %>% filter(model_set == "Discovery"), "Discovery"))
p6a <- plot_ring(f6a, "SCLC discovery cohort", center_label = paste0("n=", sum(mc_split$model_set == "Discovery")))
f6b <- panel_data("Figure6B_discovery_and_model_screen", mc_oof)
p6b <- forest_auc(f6b, "Discovery and model screen", dataset_filter = "Discovery 5-fold out-of-fold")
f6c <- panel_data("Figure6C_SCLC_validation_cohort_ring", modelc_ring(mc_split %>% filter(model_set == "Validation"), "Validation"))
p6c <- plot_ring(f6c, "SCLC validation cohort", center_label = paste0("n=", sum(mc_split$model_set == "Validation")))
modelc_gold_features <- keep_genes(gold_genes$Figure6DE_ModelC_features)
modelc_display_heat <- make_expr_long(
  expr_df,
  mc_split,
  modelc_gold_features,
  mc_split$sample_id,
  "Figure6DE_gold_standard_feature_heatmap_all_sets",
  group_cols = c("model_set", "truth")
)
f6d <- panel_data("Figure6D_discovery_cohort_feature_heatmap", modelc_display_heat %>% filter(model_set == "Discovery"))
p6d <- patchwork::wrap_elements(full = heatmap_plot(f6d, "Discovery cohort feature proteins", sample_order = unique(f6d$sample_id[order(f6d$truth)]), feature_order = modelc_gold_features))
f6e <- panel_data("Figure6E_validation_cohort_feature_heatmap", modelc_display_heat %>% filter(model_set == "Validation"))
p6e <- patchwork::wrap_elements(full = heatmap_plot(f6e, "Validation cohort feature proteins", sample_order = unique(f6e$sample_id[order(f6e$truth)]), feature_order = modelc_gold_features))
f6f <- panel_data("Figure6F_validation_ROC", mc_roc %>% filter(dataset == "Validation holdout"))
p6f <- roc_plot(f6f, mc_perf %>% filter(dataset == "Validation holdout"), "Validation ROC")
f6g <- panel_data("Figure6G_validation_confusion", mc_conf %>% filter(dataset == "Validation holdout"))
p6g <- confusion_plot(f6g, "Validation confusion")
modelc_group_map <- data.frame(
  genesymbol = c("HIST1H4A", "HIST1H4B", "HIST1H4C", "VDAC2", "VDAC3", "TMEM63A", "RAB7A", "NDRG1", "RPSA"),
  feature_group = c(rep("Histone / chromatin", 3), rep("Mitochondrial channel / autophagy", 2), rep("Membrane trafficking", 2), rep("Stress / ribosome / translation", 2)),
  source_proteins = c(rep("HIST1H4A; HIST1H4B; HIST1H4C", 3), rep("VDAC2; VDAC3", 2), rep("TMEM63A; RAB7A", 2), rep("NDRG1; RPSA", 2)),
  check.names = FALSE
)
f6h <- panel_data("Figure6H_feature_group_scores_and_source_proteins", modelc_display_heat %>%
  left_join(modelc_group_map, by = "genesymbol") %>%
  filter(!is.na(feature_group)) %>%
  group_by(sample_id, model_set, truth, feature_group, source_proteins) %>%
  summarise(feature_group_score = mean(z_abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(feature_group = factor(feature_group, levels = c("Histone / chromatin", "Mitochondrial channel / autophagy", "Membrane trafficking", "Stress / ribosome / translation"))))
p6h <- ggplot(f6h, aes(truth, feature_group_score, fill = truth)) +
  geom_boxplot(outlier.size = 0.25, linewidth = 0.18) +
  geom_jitter(width = 0.08, size = 0.28, alpha = 0.32) +
  facet_wrap(~feature_group, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = gold_pal, guide = "none") +
  labs(title = "Feature-group scores and source proteins", x = NULL, y = "Feature-group score") +
  theme_gold(base_size = 5.5)
fig6 <- ((p6a | p6b) + plot_layout(widths = c(0.78, 1.22))) /
  ((p6c | p6d) + plot_layout(widths = c(0.72, 1.28))) /
  ((p6e | (p6f | p6g)) + plot_layout(widths = c(1.05, 0.95))) /
  p6h +
  plot_layout(heights = c(0.82, 1.05, 1.0, 0.72)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 7))
save_gold_pdf(fig6, "Figure6.pdf")

# Supplementary Figure 1 ------------------------------------------------------
s1a <- panel_data("Supplementary_Figure_1A_sample_level_protein_detection", qc_sample)
ps1a <- ggplot(s1a, aes(major_type, detected_fraction, fill = major_type)) +
  geom_violin(color = "grey35", linewidth = 0.18, alpha = 0.78) +
  geom_boxplot(width = 0.16, outlier.size = 0.22, linewidth = 0.16) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = gold_pal, guide = "none") +
  labs(title = "Sample-level protein detection", x = NULL, y = "Fraction of detected proteins") +
  theme_gold()
s1b <- panel_data("Supplementary_Figure_1B_sample_level_QC_scatter", qc_sample)
ps1b <- ggplot(s1b, aes(median_abundance_detected, detected_fraction, color = major_type)) +
  geom_point(size = 0.8, alpha = 0.70) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = gold_pal) +
  labs(title = "Sample-level QC scatter", x = "Median protein abundance", y = "Protein detection rate") +
  theme_gold() +
  theme(legend.position = "bottom")
s1c <- panel_data("Supplementary_Figure_1C_total_protein_abundance", qc_sample)
ps1c <- ggplot(s1c, aes(major_type, total_abundance, fill = major_type)) +
  geom_violin(color = "grey35", linewidth = 0.18, alpha = 0.78) +
  geom_boxplot(width = 0.16, outlier.size = 0.22, linewidth = 0.16) +
  scale_fill_manual(values = gold_pal, guide = "none") +
  labs(title = "Total protein abundance", x = NULL, y = "Total abundance") +
  theme_gold()
supp1 <- (ps1a | ps1b | ps1c) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 7))
save_gold_pdf(supp1, "Supplementary_Figure_1.pdf")

# Supplementary Figure 2 ------------------------------------------------------
lc_de <- read_result("differential_analysis", "LC_vs_Healthy_limma_all_proteins.csv")
lc_modules <- read_result("differential_analysis", "LC_vs_Healthy_curated_module_scores.csv")
module_effect_lc <- read_result("differential_analysis", "LC_vs_Healthy_curated_module_stats.csv") %>%
  transmute(
    module,
    adjusted_beta = adjusted_beta_LungCancer,
    p_value = adjusted_p,
    FDR = adjusted_FDR,
    covariates,
    median_Healthy,
    median_LungCancer,
    delta_median_LungCancer_minus_Healthy,
    adjusted_formula
  )
s2a <- panel_data("Supplementary_Figure_2A_lung_cancer_vs_healthy_volcano", lc_de)
ps2a <- volcano_plot(s2a, "Lung cancer versus healthy", label_genes = gold_genes$Supplementary2A_labels)
s2b_genes <- gold_genes$Supplementary2B_shifts
s2b <- panel_data("Supplementary_Figure_2B_largest_LC_associated_protein_shifts", lc_de %>% filter(gene %in% s2b_genes) %>% mutate(gene = factor(gene, levels = rev(s2b_genes))))
ps2b <- ggplot(s2b, aes(gene, logFC, fill = logFC > 0)) +
  geom_col(color = "grey35", linewidth = 0.14) +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE` = gold_pal[["LungCancer"]], `FALSE` = gold_pal[["Healthy"]]), guide = "none") +
  labs(title = "Largest LC-associated protein shifts", x = NULL, y = "log2 fold change") +
  theme_gold(base_size = 5.8)
lc_top <- keep_genes(gold_genes$Supplementary2C_heatmap)
s2c <- make_expr_long(expr_df, clinical, lc_top, clinical$sample_id[clinical$group %in% c("Healthy", "LungCancer")], "Supplementary_Figure_2C_top_LC_associated_proteins", group_cols = c("group", "major_type", "histology", "n_status"))
ps2c <- patchwork::wrap_elements(full = heatmap_plot(s2c, "Top LC-associated proteins", sample_order = unique(s2c$sample_id[order(s2c$group, s2c$major_type)]), feature_order = lc_top))
s2d <- panel_data("Supplementary_Figure_2D_curated_biological_module_scores", lc_modules)
ps2d <- module_boxplot(s2d, "group", "Curated biological module scores", fill_var = "group", ncol = 3)
s2e <- panel_data("Supplementary_Figure_2E_adjusted_biological_module_effects", module_effect_lc)
ps2e <- ggplot(s2e, aes(adjusted_beta, reorder(module, adjusted_beta), size = -log10(p_value), color = adjusted_beta)) +
  geom_vline(xintercept = 0, color = "grey78", linewidth = 0.2) +
  geom_point(alpha = 0.85) +
  scale_color_gradient2(low = gold_pal[["Healthy"]], mid = "grey85", high = gold_pal[["LungCancer"]], midpoint = 0) +
  labs(title = "Adjusted biological module effects", x = "Adjusted beta", y = NULL, size = "-log10(P)", color = "Beta") +
  theme_gold(base_size = 5.8)
supp2 <- (ps2a | ps2b) / ps2c / ps2d / ps2e +
  plot_layout(heights = c(0.95, 1.2, 1.0, 0.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 7))
save_gold_pdf(supp2, "Supplementary_Figure_2.pdf")

# Supplementary Figure 3 ------------------------------------------------------
s3a <- panel_data("Supplementary_Figure_3A_external_validation_cohorts", bind_rows(
  clinical %>% count(cohort = "Corrected original cohort", category = group, name = "n"),
  zhang_clin %>% count(cohort = "Zhang validation cohort", category = group, name = "n")
) %>%
  group_by(cohort) %>%
  mutate(prop = n / sum(n), label = as.character(n)) %>%
  ungroup() %>%
  mutate(cohort = factor(cohort, levels = c("Zhang validation cohort", "Corrected original cohort"))))
ps3a <- ggplot(s3a, aes(prop, cohort, fill = category)) +
  geom_col(color = "grey35", linewidth = 0.15, width = 0.52) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2.2, color = "black") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = gold_pal, guide = "none") +
  labs(title = "External validation cohorts", x = "Sample proportion", y = NULL) +
  theme_gold()
zhang_expr <- read_result("processed_data", "zhang_expression_ModelA_markers_log2.csv")
marker_consistency <- read_result("modelA", "modelA_marker_direction_consistency_original_vs_Zhang.csv")
s3b <- panel_data("Supplementary_Figure_3B_five_marker_direction_consistency", marker_consistency)
ps3b <- ggplot(s3b, aes(standardized_delta_LC_minus_Healthy, reorder(genesymbol, standardized_delta_LC_minus_Healthy), color = dataset)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.2) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c(`Corrected original cohort` = "grey35", `Zhang validation cohort` = gold_pal[["LungCancer"]])) +
  labs(title = "Five-marker direction consistency", x = "Standardized delta", y = NULL, color = NULL) +
  theme_gold()
supp3 <- (ps3a | ps3b) + plot_layout(widths = c(0.8, 1.4)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 7))
save_gold_pdf(supp3, "Supplementary_Figure_3.pdf")

# Supplementary Figure 4 ------------------------------------------------------
s4a <- panel_data("Supplementary_Figure_4A_NSCLC_nodal_analysis_cohort", bind_rows(
  clinical %>% filter(major_type == "NSCLC", n_status %in% c("N0", "Nplus")) %>% count(row = "NSCLC nodal cohort", category = n_status, name = "n"),
  clinical %>% filter(histology == "LUAD", n_status %in% c("N0", "Nplus")) %>% count(row = "LUAD subset", category = n_status, name = "n"),
  clinical %>% filter(histology == "LUSC", n_status %in% c("N0", "Nplus")) %>% count(row = "LUSC subset", category = n_status, name = "n")
))
ps4a <- ggplot(s4a, aes(n, row, fill = category)) +
  geom_col(color = "grey35", linewidth = 0.15) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 2.2) +
  scale_fill_manual(values = gold_pal) +
  labs(title = "NSCLC nodal analysis cohort", x = "Samples", y = NULL, fill = NULL) +
  theme_gold()
mb_grid <- read_result("modelB", "modelB_oob_threshold_grid.csv")
s4b <- panel_data("Supplementary_Figure_4B_modelB_threshold_sensitivity", mb_grid)
ps4b <- ggplot(s4b, aes(threshold)) +
  geom_line(aes(y = accuracy, color = "Accuracy"), linewidth = 0.45) +
  geom_line(aes(y = balanced_accuracy, color = "Balanced accuracy"), linewidth = 0.45) +
  geom_line(aes(y = sensitivity, color = "Sensitivity"), linewidth = 0.45) +
  geom_line(aes(y = specificity, color = "Specificity"), linewidth = 0.45) +
  geom_vline(xintercept = unique(mb_perf$threshold)[1], linetype = 2, linewidth = 0.25) +
  scale_color_manual(values = c(Accuracy = "grey35", `Balanced accuracy` = "#D1A447", Sensitivity = gold_pal[["LungCancer"]], Specificity = gold_pal[["Healthy"]])) +
  labs(title = "Model B threshold sensitivity", x = "Probability threshold", y = "Metric", color = NULL) +
  theme_gold()
calib <- mb_val %>%
  mutate(risk_quartile = dplyr::ntile(score, 4)) %>%
  group_by(risk_quartile) %>%
  summarise(
    mean_predicted_risk = mean(score),
    observed_Nplus_rate = mean(truth == "Nplus"),
    n = n(),
    .groups = "drop"
  )
s4c <- panel_data("Supplementary_Figure_4C_validation_cohort_calibration", calib)
ps4c <- ggplot(s4c, aes(mean_predicted_risk, observed_Nplus_rate, size = n)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey75", linewidth = 0.22) +
  geom_line(color = gold_pal[["LungCancer"]], linewidth = 0.35) +
  geom_point(color = gold_pal[["LungCancer"]], alpha = 0.82) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(title = "Validation cohort calibration", x = "Mean predicted risk", y = "Observed Nplus rate", size = "n") +
  theme_gold()
supp4 <- (ps4a / (ps4b | ps4c)) + plot_layout(heights = c(0.8, 1.5)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 7))
save_gold_pdf(supp4, "Supplementary_Figure_4.pdf")

# Supplementary Figure 5 ------------------------------------------------------
sclc_de <- read_result("modelC", "modelC_sclc_discovery_limma_all_proteins.csv")
sclc_discovery <- mc_split %>% filter(model_set == "Discovery")
s5a <- panel_data("Supplementary_Figure_5A_SCLC_discovery_cohort_ring", ring_data(list(
  Nodal = sclc_discovery %>% count(category = truth, name = "n"),
  Sex = sclc_discovery %>% count(category = sex, name = "n"),
  Age = sclc_discovery %>% mutate(category = ifelse(age <= 60, "Age <=60", "Age >60")) %>% count(category, name = "n"),
  Smoking = sclc_discovery %>% count(category = smoking_status, name = "n")
)))
ps5a <- plot_ring(s5a, "SCLC discovery cohort", center_label = paste0("n=", nrow(sclc_discovery)))
s5b <- panel_data("Supplementary_Figure_5B_SCLC_discovery_volcano", sclc_de)
ps5b <- volcano_plot(s5b, "SCLC discovery Nplus vs N0")
s5c <- panel_data("Supplementary_Figure_5C_SCLC_discovery_top_pathways", read_result("modelC", "modelC_sclc_top_pathways.csv"))
ps5c <- if (nrow(s5c) > 0) {
  ggplot(s5c, aes(score, reorder(label, score), fill = source)) +
    geom_col(color = "grey35", linewidth = 0.15) +
    scale_fill_manual(values = c(`GO BP ORA` = "#67B7C4", `Reactome ORA` = "#D1A447", `Reactome GSEA` = "#C75A56"), na.value = "grey70") +
    labs(title = "SCLC discovery enriched programs", x = "-log10(adjusted P)", y = NULL, fill = NULL) +
    theme_gold(base_size = 5.6) +
    theme(legend.position = "bottom")
} else {
  empty_plot("No computed SCLC pathway enrichment")
}
supp5 <- (ps5a | ps5b | ps5c) + plot_layout(widths = c(0.8, 1.2, 1.1)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 7))
save_gold_pdf(supp5, "Supplementary_Figure_5.pdf")

standard_figure_files <- c(paste0("Figure", 1:6, ".pdf"), paste0("Supplementary_Figure_", 1:5, ".pdf"))
extra_figure_files <- "Supplementary_Figure_3F_full_six_module_scores.pdf"
figure_index <- bind_rows(
  data.frame(
    figure = standard_figure_files,
    script = "code/04_generate_figures.R",
    reference_structure = standard_figure_files,
    output = file.path("result_file/figures", standard_figure_files),
    check.names = FALSE
  ),
  data.frame(
    figure = extra_figure_files,
    script = "code/04_generate_figures.R",
    reference_structure = NA_character_,
    output = file.path("result_file/figures", extra_figure_files),
    check.names = FALSE
  )
)
write_csv_safe(figure_index, file.path(dir_map$figures, "figure_output_index.csv"))

cat("Gold-standard-style figures completed.\n")
print(figure_index)
