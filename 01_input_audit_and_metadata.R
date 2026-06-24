source(file.path(Sys.getenv("UR_FILE0624_PROJECT_DIR", unset = Sys.getenv("UR_PROJECT_DIR", unset = Sys.getenv("UR_FILE0623_PROJECT_DIR", unset = getwd()))), "code", "00_common_functions.R"))

ensure_output_dirs()

dat <- read_project_data()
clinical <- dat$clinical
expr <- dat$expr

audit_rows <- list(
  data.frame(metric = "clinical_rows_total", value = nrow(clinical), expected = 280, pass = nrow(clinical) == 280),
  data.frame(metric = "clinical_columns", value = ncol(clinical), expected = NA, pass = NA),
  data.frame(metric = "expression_rows_proteins", value = nrow(expr), expected = 3058, pass = nrow(expr) == 3058),
  data.frame(metric = "expression_sample_columns", value = ncol(expr), expected = 280, pass = ncol(expr) == 280),
  data.frame(metric = "healthy_n", value = sum(clinical$group == "Healthy"), expected = 115, pass = sum(clinical$group == "Healthy") == 115),
  data.frame(metric = "lung_cancer_n", value = sum(clinical$group == "LungCancer"), expected = 165, pass = sum(clinical$group == "LungCancer") == 165),
  data.frame(metric = "nsclc_n", value = sum(clinical$major_type == "NSCLC"), expected = 125, pass = sum(clinical$major_type == "NSCLC") == 125),
  data.frame(metric = "sclc_n", value = sum(clinical$major_type == "SCLC"), expected = 40, pass = sum(clinical$major_type == "SCLC") == 40),
  data.frame(metric = "nsclc_n0", value = sum(clinical$major_type == "NSCLC" & clinical$n_status == "N0"), expected = 76, pass = sum(clinical$major_type == "NSCLC" & clinical$n_status == "N0") == 76),
  data.frame(metric = "nsclc_nplus", value = sum(clinical$major_type == "NSCLC" & clinical$n_status == "Nplus"), expected = 49, pass = sum(clinical$major_type == "NSCLC" & clinical$n_status == "Nplus") == 49),
  data.frame(metric = "luad_n0", value = sum(clinical$histology == "LUAD" & clinical$n_status == "N0"), expected = 57, pass = sum(clinical$histology == "LUAD" & clinical$n_status == "N0") == 57),
  data.frame(metric = "luad_nplus", value = sum(clinical$histology == "LUAD" & clinical$n_status == "Nplus"), expected = 34, pass = sum(clinical$histology == "LUAD" & clinical$n_status == "Nplus") == 34),
  data.frame(metric = "lusc_n0", value = sum(clinical$histology == "LUSC" & clinical$n_status == "N0"), expected = 19, pass = sum(clinical$histology == "LUSC" & clinical$n_status == "N0") == 19),
  data.frame(metric = "lusc_nplus", value = sum(clinical$histology == "LUSC" & clinical$n_status == "Nplus"), expected = 15, pass = sum(clinical$histology == "LUSC" & clinical$n_status == "Nplus") == 15),
  data.frame(metric = "sclc_n0", value = sum(clinical$major_type == "SCLC" & clinical$n_status == "N0"), expected = 16, pass = sum(clinical$major_type == "SCLC" & clinical$n_status == "N0") == 16),
  data.frame(metric = "sclc_nplus", value = sum(clinical$major_type == "SCLC" & clinical$n_status == "Nplus"), expected = 24, pass = sum(clinical$major_type == "SCLC" & clinical$n_status == "Nplus") == 24),
  data.frame(metric = "proteins_n", value = nrow(expr), expected = 3058, pass = nrow(expr) == 3058),
  data.frame(metric = "expression_missing_values", value = sum(is.na(expr)), expected = NA, pass = NA),
  data.frame(metric = "expression_zero_values", value = sum(expr == 0, na.rm = TRUE), expected = NA, pass = NA),
  data.frame(metric = "expression_sample_ids_match_clinical", value = setequal(colnames(expr), clinical$sample_id), expected = TRUE, pass = setequal(colnames(expr), clinical$sample_id))
)

audit <- bind_rows(audit_rows)
write_csv_safe(audit, file.path(dir_map$metadata, "input_data_audit.csv"))

database_sources <- data.frame(
	  resource = c(
	    "GO biological process",
	    "Reactome pathways",
	    "Human gene ID mapping",
	    "Curated Figure3F protein-group signatures",
	    "Curated six pathway-module score definitions",
	    "clusterProfiler",
	    "ReactomePA",
	    "org.Hs.eg.db",
    "limma",
    "pROC",
    "randomForest",
    "glmnet",
    "e1071",
    "MASS"
  ),
	  source_or_package = c(
	    "clusterProfiler::enrichGO ont='BP'",
	    "ReactomePA::enrichPathway and ReactomePA::gsePathway",
	    "clusterProfiler::bitr with org.Hs.eg.db SYMBOL to ENTREZID",
	    "Fixed display protein groups reproduced from Figure3D/3F target structure; scores recomputed from rawdata/expression.csv",
	    "code/00_common_functions.R module_sets with GO/Reactome pathway identifiers",
	    "R package",
	    "R package",
	    "Bioconductor annotation package",
    "R package",
    "R package",
    "R package",
    "R package",
    "R package",
    "R package"
  ),
	  version = c(
	    as.character(utils::packageVersion("clusterProfiler")),
	    as.character(utils::packageVersion("ReactomePA")),
	    as.character(utils::packageVersion("org.Hs.eg.db")),
	    "file0624 publication package code version 2026-06-24",
	    "file0624 publication package code version 2026-06-24",
	    as.character(utils::packageVersion("clusterProfiler")),
	    as.character(utils::packageVersion("ReactomePA")),
	    as.character(utils::packageVersion("org.Hs.eg.db")),
    as.character(utils::packageVersion("limma")),
    as.character(utils::packageVersion("pROC")),
    as.character(utils::packageVersion("randomForest")),
    as.character(utils::packageVersion("glmnet")),
    as.character(utils::packageVersion("e1071")),
    as.character(utils::packageVersion("MASS"))
  ),
	  notes = c(
	    "GO BP over-representation; universe is all 3058 quantified urine proteins mapped to ENTREZID.",
	    "Reactome over-representation and GSEA; universe is all 3058 quantified urine proteins mapped to ENTREZID.",
	    "SYMBOL-to-ENTREZID mapping is recorded in each enrichment output log by mapping rate.",
	    "Main Figure3F uses four protein groups from panel D; sample score is the mean row z-score of listed proteins, and plotted value is the stage median.",
	    "Supplemental six-module heatmap uses the six curated GO/Reactome modules; table saves raw stage median_score and module-wise median_z.",
	    "Used for GO ORA and ID mapping.",
	    "Used for Reactome ORA/GSEA.",
	    "Annotation database version captured from installed package.",
    "Empirical Bayes moderated linear models.",
    "ROC curves, AUC, DeLong CIs, Youden thresholds.",
    "Random forest feature ranking and predictive models.",
    "LASSO and elastic-net logistic candidate models.",
    "Radial SVM candidate models.",
    "Linear discriminant analysis candidate models."
  ),
  check.names = FALSE
)
write_csv_safe(database_sources, file.path(dir_map$metadata, "database_sources.csv"))

analysis_parameters <- list(
  project_dir = project_dir,
  random_seed = list(
    rng_kind = "Mersenne-Twister",
    normal_kind = "Inversion",
    sample_kind = "Rejection",
    preprocessing_discovery_curve_seed = "20250600 + number of samples in the current group",
    candidate_model_rng = "fit_candidate_methods() calls set.seed(seed) before glmnet/randomForest/SVM/LDA fitting",
    modelA_script_seed = 20250613,
    modelA_split = "deterministic: samples arranged by group and sample_id; first 92 Healthy and first 132 LungCancer are Discovery",
    modelA_candidate_model_screen_seed = 20250613,
    modelA_single_marker_rf_cv_seed = "20250611 + length(feature_set)",
    modelB_split = 20250603,
    modelB_screening_rf_seed_source = "inherits the RNG state immediately after stratified_train_indices(seed = 20250603)",
    modelB_candidate_model_screen_seed = 20250604,
    modelB_final_rf = 15,
    modelB_repeated_split_base = 20250603,
    modelB_repeated_split_seed = "20250603 + repeat_id for repeat_id 1..100",
    modelC_split = 20250613,
    modelC_oof = 13,
    modelC_candidate_validation_seed = 20250608
  ),
	  enrichment = list(
	    universe = "all 3058 quantified urine proteins in rawdata/expression.csv, mapped SYMBOL to ENTREZID",
	    pvalue_cutoff_for_selected_proteins = 0.05,
	    absolute_log2fc_cutoff_for_selected_proteins = 0.3,
	    p_adjust_method = "Benjamini-Hochberg",
	    min_gene_set_size = 10,
	    max_gene_set_size = 500
	  ),
	  curated_module_score_definitions = list(
	    score_background = "Rows are z-scored within the samples used for each score set; Figure3 progression scores use Healthy + NSCLC N0 + NSCLC Nplus samples.",
	    main_Figure3F = list(
	      plotted_value = "stage median of sample-level mean row z-score protein-group score",
	      protein_groups = list(
	        gradual_up = c("SLC13A2", "ACE2", "ANXA5", "ATP1A1", "ATP1A3", "SLC4A4", "AQP1", "STOM"),
	        gradual_down = c("S100A8", "TIMP1", "ANXA1"),
	        acute_phase = c("LBP", "ORM2", "ORM1", "HP", "SERPINA3", "SAA1", "ITIH4", "A2M", "APCS", "SERPINA1", "SAA2", "CRP"),
	        complement_coagulation = c("FGB", "SERPINC1", "C4B", "FGG", "C4A", "CFH", "CFB", "C9", "FGA", "C3", "PLG", "CFI", "C1QB", "C1QC", "C2", "C5", "C6", "C7", "C8A", "C8B", "C8G")
	      )
	    ),
	    supplemental_six_modules = list(
	      plotted_value = "module-wise scaled stage medians (median_z); raw median_score is saved in plotted data",
	      module_sets = module_sets
	    )
	  ),
	  model_split_rules = list(
	    modelA = "Healthy/LungCancer stratified discovery n=224 (Healthy=92, LungCancer=132); internal validation n=56; Zhang LC/CTL external validation n=66.",
	    modelB = "NSCLC N0/Nplus stratified discovery n=101 (N0=60, Nplus=41); validation n=24 (N0=16, Nplus=8). Feature selection, pruning, and threshold are discovery-only.",
    modelC = "SCLC stratified discovery n=30 (N0=12, Nplus=18); validation n=10 (N0=4, Nplus=6). Exploratory only."
  ),
  model_random_components = list(
    modelA = list(
      final_score = "No fitted random model; direction weights and Youden threshold are computed from discovery data.",
      candidate_models = c("random forest ntree=500", "SVM radial", "logistic regression", "LASSO logistic", "elastic net", "LDA", "equal-direction score"),
      single_marker_cv = "5-fold CV, 20 repeats, randomForest ntree=300"
    ),
    modelB = list(
      screening_rf = "randomForest ntree=1000 on discovery samples",
      final_rf = "randomForest ntree=1000 on fixed top18 features, set.seed(15)",
      repeated_split = "100 repeats; randomForest ntree=500 for screen and final repeat model"
    ),
    modelC = list(
      oof_candidate_comparison = "5-fold stratified OOF; randomForest ntree=500; glmnet/SVM/LDA run after set.seed(13)",
      candidate_validation = "fit_candidate_methods(..., seed=20250608)",
      final_score = "Equal-direction score; threshold from discovery OOF Youden"
    )
  )
)
write_json_safe(analysis_parameters, file.path(dir_map$metadata, "analysis_parameters.json"))

sink(file.path(dir_map$metadata, "session_info.txt"))
print(sessionInfo())
sink()

failed <- audit %>% filter(!is.na(pass) & pass == FALSE)
if (nrow(failed) > 0) {
  print(failed)
  stop("Input audit failed; stop per user instruction.")
}

cat("Input audit and metadata completed.\\n")
