source(file.path(Sys.getenv("UR_FILE0624_PROJECT_DIR", unset = Sys.getenv("UR_PROJECT_DIR", unset = Sys.getenv("UR_FILE0623_PROJECT_DIR", unset = getwd()))), "code", "00_common_functions.R"))

ensure_output_dirs()

dat <- read_project_data()
clinical <- dat$clinical
expr <- dat$expr
expr_df <- dat$expr_df

write_csv_safe(clinical, file.path(dir_map$processed_data, "clinical_processed.csv"))
write_csv_safe(expr_df, file.path(dir_map$processed_data, "expression_processed.csv"))

sample_metrics <- data.frame(
  sample_id = colnames(expr),
  detected_protein_n = colSums(expr > 0, na.rm = TRUE),
  detected_fraction = colMeans(expr > 0, na.rm = TRUE),
  total_abundance = colSums(expr, na.rm = TRUE),
  median_abundance_detected = apply(expr, 2, function(x) median(x[x > 0], na.rm = TRUE)),
  check.names = FALSE
) %>%
  left_join(clinical, by = "sample_id")

protein_metrics <- data.frame(
  genesymbol = rownames(expr),
  detection_fraction = rowMeans(expr > 0, na.rm = TRUE),
  detected_sample_n = rowSums(expr > 0, na.rm = TRUE),
  mean_abundance = rowMeans(expr, na.rm = TRUE),
  median_abundance = apply(expr, 1, median, na.rm = TRUE),
  variance = apply(expr, 1, var, na.rm = TRUE),
  check.names = FALSE
)

top_var <- protein_metrics %>%
  arrange(desc(variance)) %>%
  slice_head(n = min(1000, nrow(.))) %>%
  pull(genesymbol)
pca_input <- t(expr[top_var, , drop = FALSE])
pca_input <- scale(pca_input)
pca_input[is.na(pca_input)] <- 0
pca <- prcomp(pca_input, center = FALSE, scale. = FALSE)
pca_scores <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC1_var = 100 * (pca$sdev[1]^2 / sum(pca$sdev^2)),
  PC2_var = 100 * (pca$sdev[2]^2 / sum(pca$sdev^2)),
  check.names = FALSE
) %>% left_join(clinical, by = "sample_id")

discovery_curve <- function(samples, label, repeats = 100) {
  set.seed(20250600 + length(samples))
  out <- vector("list", repeats)
  detected <- expr[, samples, drop = FALSE] > 0
  for (i in seq_len(repeats)) {
    ord <- sample(seq_along(samples))
    cumulative <- vapply(seq_along(ord), function(k) sum(rowSums(detected[, ord[seq_len(k)], drop = FALSE]) > 0), numeric(1))
    out[[i]] <- data.frame(repeat_id = i, sample_n = seq_along(ord), cumulative_proteins = cumulative, group = label)
  }
  bind_rows(out)
}

healthy_samples <- clinical$sample_id[clinical$group == "Healthy"]
lc_samples <- clinical$sample_id[clinical$group == "LungCancer"]
curve_raw <- bind_rows(
  discovery_curve(healthy_samples, "Healthy"),
  discovery_curve(lc_samples, "Lung cancer")
)
curve_summary <- curve_raw %>%
  group_by(group, sample_n) %>%
  summarise(
    mean_cumulative = mean(cumulative_proteins),
    q10 = quantile(cumulative_proteins, 0.10),
    q90 = quantile(cumulative_proteins, 0.90),
    .groups = "drop"
  )

write_csv_safe(sample_metrics, file.path(dir_map$processed_data, "QC_sample_metrics.csv"))
write_csv_safe(protein_metrics, file.path(dir_map$processed_data, "QC_protein_metrics.csv"))
write_csv_safe(pca_scores, file.path(dir_map$processed_data, "QC_pca_scores.csv"))
write_csv_safe(curve_raw, file.path(dir_map$processed_data, "protein_discovery_curve_repeats.csv"))
write_csv_safe(curve_summary, file.path(dir_map$processed_data, "protein_discovery_curve_summary.csv"))

zhang <- parse_zhang_data(markers = c("RAB1B", "ATP1A3", "NAPSA", "APOA2", "ANXA11"))
zhang_expr_df <- data.frame(genesymbol = rownames(zhang$expr), zhang$expr, check.names = FALSE)
write_csv_safe(zhang$clinical, file.path(dir_map$processed_data, "zhang_clinical_LC_CTL.csv"))
write_csv_safe(zhang_expr_df, file.path(dir_map$processed_data, "zhang_expression_ModelA_markers_log2.csv"))

summary_table <- data.frame(
  item = c(
    "total_samples", "healthy", "lung_cancer", "NSCLC", "SCLC",
    "NSCLC_N0", "NSCLC_Nplus", "LUAD_N0", "LUAD_Nplus",
    "LUSC_N0", "LUSC_Nplus", "SCLC_N0", "SCLC_Nplus",
    "proteins", "protein_detection_fraction_median",
    "PCA_top_variable_proteins", "PCA_PC1_percent", "PCA_PC2_percent",
    "Zhang_LC_CTL_samples"
  ),
  value = c(
    nrow(clinical),
    sum(clinical$group == "Healthy"),
    sum(clinical$group == "LungCancer"),
    sum(clinical$major_type == "NSCLC"),
    sum(clinical$major_type == "SCLC"),
    sum(clinical$major_type == "NSCLC" & clinical$n_status == "N0"),
    sum(clinical$major_type == "NSCLC" & clinical$n_status == "Nplus"),
    sum(clinical$histology == "LUAD" & clinical$n_status == "N0"),
    sum(clinical$histology == "LUAD" & clinical$n_status == "Nplus"),
    sum(clinical$histology == "LUSC" & clinical$n_status == "N0"),
    sum(clinical$histology == "LUSC" & clinical$n_status == "Nplus"),
    sum(clinical$major_type == "SCLC" & clinical$n_status == "N0"),
    sum(clinical$major_type == "SCLC" & clinical$n_status == "Nplus"),
    nrow(expr),
    median(protein_metrics$detection_fraction),
    length(top_var),
    unique(round(pca_scores$PC1_var, 3))[1],
    unique(round(pca_scores$PC2_var, 3))[1],
    nrow(zhang$clinical)
  ),
  check.names = FALSE
)
write_csv_safe(summary_table, file.path(dir_map$tables, "cohort_and_qc_summary.csv"))

cat("Preprocessing completed.\\n")

