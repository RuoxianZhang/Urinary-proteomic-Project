source(file.path(Sys.getenv("UR_FILE0624_PROJECT_DIR", unset = Sys.getenv("UR_PROJECT_DIR", unset = Sys.getenv("UR_FILE0623_PROJECT_DIR", unset = getwd()))), "code", "00_common_functions.R"))

ensure_output_dirs()

read_result <- function(...) read.csv(file.path(result_dir, ...), check.names = FALSE)

md_table <- function(df, digits = 3) {
  if (is.null(df) || nrow(df) == 0) return("_No rows._\n")
  df <- as.data.frame(df, check.names = FALSE)
  df[] <- lapply(df, function(x) {
    if (is.numeric(x)) fmt_num(x, digits) else as.character(x)
  })
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

rel <- function(path) sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", project_dir), "/?"), "", normalizePath(path, mustWork = FALSE))

audit <- read_result("metadata", "input_data_audit.csv")
diff_summary <- read_result("tables", "differential_analysis_summary.csv")
ma_perf <- read_result("modelA", "modelA_final_performance.csv")
mb_perf <- read_result("modelB", "modelB_final_performance.csv")
mc_perf <- read_result("modelC", "modelC_final_performance.csv")
ma_cand <- read_result("modelA", "modelA_candidate_model_performance.csv")
mb_cand <- read_result("modelB", "modelB_candidate_model_performance.csv")
mc_cand <- read_result("modelC", "modelC_discovery_5fold_oof_auc.csv")
mc_feat <- read_result("modelC", "modelC_feature_table.csv")

expected_checks <- data.frame(
  item = c(
    "total n", "Healthy n", "Lung cancer n", "NSCLC n", "SCLC n",
    "NSCLC N0/Nplus", "LUAD N0/Nplus", "LUSC N0/Nplus", "SCLC N0/Nplus",
    "proteins n"
  ),
  expected = c("280", "115", "165", "125", "40", "76/49", "57/34", "19/15", "16/24", "3058"),
  observed = c(
    audit$value[audit$metric == "clinical_rows_total"],
    audit$value[audit$metric == "healthy_n"],
    audit$value[audit$metric == "lung_cancer_n"],
    audit$value[audit$metric == "nsclc_n"],
    audit$value[audit$metric == "sclc_n"],
    paste0(audit$value[audit$metric == "nsclc_n0"], "/", audit$value[audit$metric == "nsclc_nplus"]),
    paste0(audit$value[audit$metric == "luad_n0"], "/", audit$value[audit$metric == "luad_nplus"]),
    paste0(audit$value[audit$metric == "lusc_n0"], "/", audit$value[audit$metric == "lusc_nplus"]),
    paste0(audit$value[audit$metric == "sclc_n0"], "/", audit$value[audit$metric == "sclc_nplus"]),
    audit$value[audit$metric == "proteins_n"]
  ),
  check.names = FALSE
) %>%
  mutate(status = ifelse(expected == observed, "match", "mismatch"))
write_csv_safe(expected_checks, file.path(dir_map$metadata, "required_input_count_checks.csv"))

panel_files <- list.files(dir_map$plotted, pattern = "\\.csv$", full.names = TRUE)
panel_index <- data.frame(
  panel = sub("\\.csv$", "", basename(panel_files)),
  plotted_data = rel(panel_files),
  figure = sub("^(Figure[0-9]+|Supplementary_Figure_[0-9]+).*", "\\1", sub("\\.csv$", "", basename(panel_files))),
  script = "code/04_generate_figures.R",
  check.names = FALSE
) %>%
  mutate(
    primary_result_family = case_when(
      grepl("^Figure2|modelA|ModelA|Zhang", panel) ~ "result_file/modelA and processed Zhang tables",
      grepl("^Figure5|modelB|ModelB", panel) ~ "result_file/modelB",
      grepl("^Figure6|modelC|ModelC|SCLC", panel) ~ "result_file/modelC and SCLC differential outputs",
      grepl("^Figure3|^Figure4|Supplementary_Figure_2|volcano|module|NSCLC|LUAD|LUSC|LC_", panel) ~ "result_file/differential_analysis and result_file/enrichment",
      grepl("QC|PCA|protein|cohort|detection|discovery", panel, ignore.case = TRUE) ~ "result_file/processed_data",
      TRUE ~ "result_file computed tables"
    )
  ) %>%
  arrange(figure, panel)
write_csv_safe(panel_index, file.path(dir_map$tables, "figure_panel_traceability.csv"))

standard_figure_files <- c(paste0("Figure", 1:6, ".pdf"), paste0("Supplementary_Figure_", 1:5, ".pdf"))
extra_figure_files <- "Supplementary_Figure_3F_full_six_module_scores.pdf"
figure_index <- data.frame(
  figure = c(standard_figure_files, extra_figure_files),
  script = "code/04_generate_figures.R",
  output = file.path("result_file/figures", c(standard_figure_files, extra_figure_files)),
  check.names = FALSE
)

observed_nsclc <- diff_summary %>% filter(comparison == "NSCLC_Nplus_vs_N0")
observed_luad <- diff_summary %>% filter(comparison == "LUAD_Nplus_vs_N0")
observed_sclc <- diff_summary %>% filter(comparison == "SCLC_Nplus_vs_N0_all_exploratory")
zhang_auc <- ma_perf$auc[ma_perf$dataset == "Zhang validation"][1]
modelc_min_fdr <- min(mc_feat$FDR, na.rm = TRUE)
mb_val <- mb_perf %>% filter(dataset == "Validation holdout") %>% dplyr::slice(1)
mc_val <- mc_perf %>% filter(dataset == "Validation holdout") %>% dplyr::slice(1)
mc_equal_oof <- mc_cand %>% filter(method == "Equal-direction score") %>% dplyr::slice(1)
mc_candidate_targets <- c(
  `Equal-direction score` = 0.838,
  `Logistic regression` = 0.625,
  `LASSO logistic` = 0.616,
  `Elastic net` = 0.657,
  `Random forest` = 0.773,
  `SVM radial` = 0.759,
  LDA = 0.713
)
mc_candidate_observed <- setNames(mc_cand$auc, mc_cand$method)
mc_candidate_delta <- max(abs(mc_candidate_observed[names(mc_candidate_targets)] - mc_candidate_targets), na.rm = TRUE)

consistency <- data.frame(
  result = c(
    "Input cohort counts",
    "Model A Zhang validation AUC",
    "Model B validation AUC/confusion",
    "Model B OOB threshold",
    "Model C equal-direction OOF and validation AUC",
    "Model C full candidate-model screen",
    "NSCLC Nplus vs N0 nominal/FDR counts",
    "LUAD Nplus vs N0 nominal/FDR counts",
    "SCLC all-sample exploratory nominal/FDR counts",
    "Model C feature FDR"
  ),
  expected_or_reference = c(
    "Counts specified in task",
    "Prepared manuscript/figure target noted as AUC 0.757",
    "Figure5 target AUC 0.859, Acc 0.708, Sens 0.750, Spec 0.688",
    "Figure5 target threshold 0.463",
    "Figure6 target OOF AUC 0.838 and validation AUC 0.833",
    "Figure6 target AUCs: Equal 0.838, Logistic 0.625, LASSO 0.616, Elastic 0.657, RF 0.773, SVM 0.759, LDA 0.713",
    "Prepared manuscript/figure target noted as 352 up / 15 down / 26 FDR",
    "Prepared manuscript/figure target noted as 345 up / 10 down / 4 FDR",
    "Prepared manuscript/figure target noted as approximately 250 up / 15 down",
    "Exploratory; no clinical-grade significance claim"
  ),
  observed_from_raw_reanalysis = c(
    ifelse(all(expected_checks$status == "match"), "All required counts match", "At least one required count mismatches"),
    paste0("AUC=", fmt_num(zhang_auc, 3)),
    paste0("AUC=", fmt_num(mb_val$auc, 3), "; Acc=", fmt_num(mb_val$accuracy, 3), "; Sens=", fmt_num(mb_val$sensitivity, 3), "; Spec=", fmt_num(mb_val$specificity, 3)),
    paste0("OOB Youden threshold=", fmt_num(mb_val$threshold, 3), " using current randomForest seed/version"),
    paste0("OOF AUC=", fmt_num(mc_equal_oof$auc, 3), "; validation AUC=", fmt_num(mc_val$auc, 3), "; Acc=", fmt_num(mc_val$accuracy, 3), "; Sens=", fmt_num(mc_val$sensitivity, 3), "; Spec=", fmt_num(mc_val$specificity, 3)),
    paste0("Max absolute candidate AUC delta vs Figure6 text=", fmt_num(mc_candidate_delta, 3)),
    paste0(observed_nsclc$nominal_positive_up, " up / ", observed_nsclc$nominal_negative_up, " down / ", observed_nsclc$fdr_significant, " FDR"),
    paste0(observed_luad$nominal_positive_up, " up / ", observed_luad$nominal_negative_up, " down / ", observed_luad$fdr_significant, " FDR"),
    paste0(observed_sclc$nominal_positive_up, " up / ", observed_sclc$nominal_negative_up, " down / ", observed_sclc$fdr_significant, " FDR"),
    paste0("Minimum selected-feature FDR=", fmt_num(modelc_min_fdr, 3))
  ),
  status = c(
    ifelse(all(expected_checks$status == "match"), "match", "mismatch"),
    ifelse(abs(zhang_auc - 0.757) < 0.005, "match", "needs confirmation"),
    ifelse(abs(mb_val$auc - 0.859375) < 0.005 && abs(mb_val$accuracy - 0.7083333) < 0.005 && abs(mb_val$sensitivity - 0.75) < 0.005 && abs(mb_val$specificity - 0.6875) < 0.005, "match", "needs confirmation"),
    ifelse(abs(mb_val$threshold - 0.463) < 0.001, "match", "needs confirmation"),
    ifelse(abs(mc_equal_oof$auc - 0.838) < 0.005 && abs(mc_val$auc - 0.8333333) < 0.005, "match", "needs confirmation"),
    ifelse(mc_candidate_delta < 0.02, "match", "needs confirmation"),
    ifelse(observed_nsclc$nominal_positive_up == 352 && observed_nsclc$nominal_negative_up == 15 && observed_nsclc$fdr_significant == 26, "match", "needs confirmation"),
    ifelse(observed_luad$nominal_positive_up == 345 && observed_luad$nominal_negative_up == 10 && observed_luad$fdr_significant == 4, "match", "needs confirmation"),
    "needs confirmation",
    ifelse(modelc_min_fdr < 0.05, "selected features include FDR-significant proteins", "exploratory only")
  ),
  check.names = FALSE
)
write_csv_safe(consistency, file.path(dir_map$tables, "manuscript_consistency_checks.csv"))

model_file_check <- data.frame(
  model = c("Model A", "Model B", "Model C"),
  sample_split = c("result_file/modelA/modelA_sample_split.csv", "result_file/modelB/modelB_sample_split.csv", "result_file/modelC/modelC_sample_split.csv"),
  feature_table = c("result_file/modelA/modelA_feature_table.csv", "result_file/modelB/modelB_feature_table.csv", "result_file/modelC/modelC_feature_table.csv"),
  parameters = c("result_file/modelA/modelA_parameters.json", "result_file/modelB/modelB_parameters.json", "result_file/modelC/modelC_parameters.json"),
  prediction_scores = c("result_file/modelA/modelA_final_prediction_scores.csv", "result_file/modelB/modelB_validation_predictions.csv and modelB_discovery_oob_predictions.csv", "result_file/modelC/modelC_final_prediction_scores.csv"),
  roc_coordinates = c("result_file/modelA/modelA_final_roc_coordinates.csv", "result_file/modelB/modelB_final_roc_coordinates.csv", "result_file/modelC/modelC_final_roc_coordinates.csv"),
  auc_ci = c("result_file/modelA/modelA_final_performance.csv", "result_file/modelB/modelB_final_performance.csv", "result_file/modelC/modelC_final_performance.csv"),
  threshold = c("Discovery Youden", "Discovery OOB Youden", "Discovery 5-fold OOF Youden"),
  confusion_matrix = c("result_file/modelA/modelA_final_confusion_matrices.csv", "result_file/modelB/modelB_final_confusion_matrices.csv", "result_file/modelC/modelC_final_confusion_matrices.csv"),
  candidate_model_comparison = c("result_file/modelA/modelA_candidate_model_performance.csv", "result_file/modelB/modelB_candidate_model_performance.csv", "result_file/modelC/modelC_discovery_5fold_oof_auc.csv"),
  check.names = FALSE
)
write_csv_safe(model_file_check, file.path(dir_map$tables, "model_required_outputs_index.csv"))

report <- c(
  "# Reproducibility report",
  "",
  paste0("Project directory: `", project_dir, "`"),
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope and protection",
  "",
  "- Existing `figure/`, `rawdata/`, and `zhangdata/` files were treated as read-only.",
  "- Figures and tables were regenerated from `rawdata/clinical.csv`, `rawdata/expression.csv`, `zhangdata/mmc7.xlsx`, and `zhangdata/mmc8.xlsx`.",
  "- The prepared manuscript/figures were not used as numeric data sources, and no plotted coordinates were hand-entered.",
  "",
  "## Input audit",
  "",
  md_table(expected_checks),
  "",
  "Full audit file: `result_file/metadata/input_data_audit.csv`.",
  "",
  "## Figure outputs",
  "",
  md_table(figure_index),
  "",
	  "Panel-level plotted data are indexed in `result_file/tables/figure_panel_traceability.csv`. Each PDF panel is drawn from the corresponding CSV in `result_file/figures/plotted_data/`.",
	  "The extra supplemental PDF `Supplementary_Figure_3F_full_six_module_scores.pdf` contains the complete six-module version requested for Figure3F context.",
	  "",
  "## Model outputs",
  "",
  md_table(model_file_check),
  "",
  "### Model A final performance",
  "",
  md_table(ma_perf %>% dplyr::select(model, dataset, threshold_source, auc, ci_low, ci_high, threshold, accuracy, sensitivity, specificity, balanced_accuracy)),
  "",
  "Model A final definition: discovery-derived signed/equal-direction five-protein score. Candidate models remain a screen only.",
  "",
  "### Model B final performance",
  "",
  md_table(mb_perf %>% dplyr::select(model, dataset, threshold_source, auc, ci_low, ci_high, threshold, accuracy, sensitivity, specificity, balanced_accuracy)),
  "",
  "Model B leakage control: feature selection, correlation pruning, model fitting, and threshold selection are performed in discovery/OOB data only; validation is an independent evaluation.",
  "",
  "### Model C final performance",
  "",
  md_table(mc_perf %>% dplyr::select(model, dataset, threshold_source, auc, ci_low, ci_high, threshold, accuracy, sensitivity, specificity, balanced_accuracy)),
  "",
  "Model C is exploratory only. The selected nine proteins had minimum selected-feature FDR shown in the consistency table below.",
  "",
  "## Candidate model comparisons",
  "",
  "Model A candidate screen:",
  "",
  md_table(ma_cand %>% dplyr::select(method, dataset, auc, ci_low, ci_high, accuracy, sensitivity, specificity, balanced_accuracy)),
  "",
  "Model B fixed-feature candidate comparison:",
  "",
  md_table(mb_cand %>% dplyr::select(method, dataset, auc, ci_low, ci_high, accuracy, sensitivity, specificity, balanced_accuracy)),
  "",
  "Model C discovery 5-fold OOF candidate comparison:",
  "",
  md_table(mc_cand %>% dplyr::select(method, dataset, auc, ci_low, ci_high)),
  "",
  "## Differential analysis summary",
  "",
  md_table(diff_summary),
  "",
  "## Database and software provenance",
  "",
  "- Session info: `result_file/metadata/session_info.txt`.",
  "- Database sources: `result_file/metadata/database_sources.csv`.",
  "- Analysis parameters: `result_file/metadata/analysis_parameters.json`.",
  "- GO BP enrichment was computed with `clusterProfiler::enrichGO` using `org.Hs.eg.db` SYMBOL-to-ENTREZ mapping.",
  "- Reactome ORA/GSEA were computed with `ReactomePA`; analysis universe was all 3058 proteins in `rawdata/expression.csv`.",
  "",
  "## Consistency and items needing confirmation",
  "",
  md_table(consistency),
  "",
  "The rows marked `needs confirmation` are not manually corrected in this package. They are reported because the raw-data reanalysis does not reproduce the noted manuscript/reference values exactly.",
  "",
  "## No unresolved hard stops",
  "",
  "All requested PDFs and plotted-data tables were generated from the allowed input data. The remaining issues are numerical discrepancies for manuscript reconciliation, not missing raw inputs for the generated figures."
)

writeLines(report, file.path(result_dir, "reproducibility_report.md"))
cat("Reproducibility report completed.\n")
cat(file.path(result_dir, "reproducibility_report.md"), "\n")
