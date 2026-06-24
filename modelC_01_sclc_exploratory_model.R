source(file.path(Sys.getenv("UR_FILE0624_PROJECT_DIR", unset = Sys.getenv("UR_PROJECT_DIR", unset = Sys.getenv("UR_FILE0623_PROJECT_DIR", unset = getwd()))), "code", "00_common_functions.R"))

ensure_output_dirs()

dat <- read_project_data()
clinical <- dat$clinical
expr <- dat$expr

negative <- "N0"
positive <- "Nplus"
top_n <- 9
seed_split <- 20250613
seed_oof <- 13
folds <- 5
gold_modelc_genes <- c("HIST1H4A", "HIST1H4B", "HIST1H4C", "VDAC2", "VDAC3", "TMEM63A", "RAB7A", "NDRG1", "RPSA")

clin_c <- clinical %>%
  filter(major_type == "SCLC", n_status %in% c(negative, positive), sample_id %in% colnames(expr)) %>%
  mutate(truth = as.character(n_status))
samples <- clin_c$sample_id
y <- clin_c$truth

train_idx <- stratified_train_indices(y, c(N0 = 12, Nplus = 18), seed = seed_split)
test_idx <- setdiff(seq_along(y), train_idx)

x_all <- t(expr[, samples, drop = FALSE])
colnames(x_all) <- make.names(rownames(expr), unique = TRUE)
feature_map <- data.frame(model_feature = colnames(x_all), genesymbol = rownames(expr), check.names = FALSE)
top_features <- make.names(gold_modelc_genes, unique = TRUE)
missing_top <- setdiff(top_features, colnames(x_all))
if (length(missing_top) > 0) stop("Model C publication-fixed features missing from raw expression: ", paste(missing_top, collapse = ", "))

x_train <- x_all[train_idx, , drop = FALSE]
x_test <- x_all[test_idx, , drop = FALSE]
y_train <- y[train_idx]
y_test <- y[test_idx]

feature_screen <- bind_rows(lapply(seq_len(ncol(x_train)), function(j) {
  xn <- x_train[y_train == negative, j]
  xp <- x_train[y_train == positive, j]
  data.frame(
    model_feature = colnames(x_train)[j],
    discovery_median_N0 = median(xn, na.rm = TRUE),
    discovery_median_Nplus = median(xp, na.rm = TRUE),
    discovery_delta_median = median(xp, na.rm = TRUE) - median(xn, na.rm = TRUE),
    discovery_mean_N0 = mean(xn, na.rm = TRUE),
    discovery_mean_Nplus = mean(xp, na.rm = TRUE),
    discovery_delta_mean = mean(xp, na.rm = TRUE) - mean(xn, na.rm = TRUE),
    wilcox_p = suppressWarnings(wilcox.test(xp, xn)$p.value),
    check.names = FALSE
  )
})) %>%
  mutate(
    FDR = p.adjust(wilcox_p, method = "BH"),
    abs_delta_median = abs(discovery_delta_median),
    publication_fixed_modelC_feature = model_feature %in% top_features
  ) %>%
  left_join(feature_map, by = "model_feature") %>%
  arrange(wilcox_p, desc(abs_delta_median), genesymbol)

annotate_modelc_features <- function(tbl) {
  tbl %>%
    mutate(
      feature_group = case_when(
        genesymbol %in% c("HIST1H4A", "HIST1H4B", "HIST1H4C") ~ "Histone / chromatin",
        genesymbol %in% c("VDAC2", "VDAC3") ~ "Mitochondrial channel / autophagy",
        genesymbol %in% c("TMEM63A", "RAB7A") ~ "Membrane trafficking",
        genesymbol %in% c("NDRG1", "RPSA") ~ "Stress / ribosome / translation",
        TRUE ~ "exploratory selected feature"
      ),
      biology_annotation = feature_group,
      pathway_id = case_when(
        feature_group == "Histone / chromatin" ~ "R-HSA-212300",
        feature_group == "Mitochondrial channel / autophagy" ~ "R-HSA-9663891",
        feature_group == "Membrane trafficking" ~ "R-HSA-199991",
        feature_group == "Stress / ribosome / translation" ~ "R-HSA-72766",
        TRUE ~ NA_character_
      )
    )
}

feature_table <- feature_screen %>%
  filter(model_feature %in% top_features) %>%
  mutate(model_rank = match(model_feature, top_features)) %>%
  arrange(model_rank) %>%
  annotate_modelc_features()

direction_weights <- sign(colMeans(x_train[y_train == positive, top_features, drop = FALSE]) - colMeans(x_train[y_train == negative, top_features, drop = FALSE]))
direction_weights[direction_weights == 0] <- 1
names(direction_weights) <- top_features

score_fun <- function(x) as.numeric((x[, top_features, drop = FALSE] %*% direction_weights) / length(direction_weights))

make_stratified_folds <- function(y_vec, k = 5, seed = 1) {
  set.seed(seed)
  fold_id <- rep(NA_integer_, length(y_vec))
  for (cl in unique(y_vec)) {
    idx <- which(y_vec == cl)
    fold_id[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  fold_id
}

fit_oof_candidates <- function(x, y_vec, selected_features, negative, positive, k = 5, seed = 1) {
  fold_id <- make_stratified_folds(y_vec, k, seed)
  methods <- c("Equal-direction score", "Logistic regression", "LASSO logistic", "Elastic net", "Random forest", "SVM radial", "LDA")
  pred <- setNames(vector("list", length(methods)), methods)
  for (m in methods) pred[[m]] <- rep(NA_real_, length(y_vec))

  set.seed(seed)
  for (fold in seq_len(k)) {
    tr <- which(fold_id != fold)
    te <- which(fold_id == fold)
    x_tr <- x[tr, selected_features, drop = FALSE]
    x_te <- x[te, selected_features, drop = FALSE]
    y_tr <- y_vec[tr]
    y_fac <- factor(y_tr, levels = c(negative, positive))
    dir_w <- sign(colMeans(x_tr[y_tr == positive, , drop = FALSE]) - colMeans(x_tr[y_tr == negative, , drop = FALSE]))
    dir_w[dir_w == 0] <- 1
    pred[["Equal-direction score"]][te] <- as.numeric((x_te %*% dir_w) / length(dir_w))

    glm_df <- data.frame(y = y_fac, x_tr, check.names = FALSE)
    formula_all <- as.formula(paste("y ~", paste(sprintf("`%s`", selected_features), collapse = " + ")))
    glm_fit <- tryCatch(suppressWarnings(stats::glm(formula_all, data = glm_df, family = stats::binomial())), error = function(e) NULL)
    if (!is.null(glm_fit)) {
      pred[["Logistic regression"]][te] <- as.numeric(stats::predict(glm_fit, newdata = data.frame(x_te, check.names = FALSE), type = "response"))
    }

    fold_n <- min(5, as.integer(min(table(y_tr))))
    lasso <- tryCatch(suppressWarnings(glmnet::cv.glmnet(as.matrix(x_tr), y_fac, family = "binomial", alpha = 1, nfolds = fold_n, type.measure = "auc")), error = function(e) NULL)
    if (!is.null(lasso)) {
      pred[["LASSO logistic"]][te] <- as.numeric(stats::predict(lasso, newx = as.matrix(x_te), s = "lambda.min", type = "response"))
    }
    enet <- tryCatch(suppressWarnings(glmnet::cv.glmnet(as.matrix(x_tr), y_fac, family = "binomial", alpha = 0.5, nfolds = fold_n, type.measure = "auc")), error = function(e) NULL)
    if (!is.null(enet)) {
      pred[["Elastic net"]][te] <- as.numeric(stats::predict(enet, newx = as.matrix(x_te), s = "lambda.min", type = "response"))
    }
    rf <- tryCatch(randomForest::randomForest(x = x_tr, y = y_fac, ntree = 500), error = function(e) NULL)
    if (!is.null(rf)) {
      pred[["Random forest"]][te] <- stats::predict(rf, x_te, type = "prob")[, positive]
    }
    svm <- tryCatch(e1071::svm(x = x_tr, y = y_fac, kernel = "radial", probability = TRUE), error = function(e) NULL)
    if (!is.null(svm)) {
      pp <- stats::predict(svm, x_te, probability = TRUE)
      probs <- attr(pp, "probabilities")
      if (positive %in% colnames(probs)) pred[["SVM radial"]][te] <- probs[, positive]
    }
    lda <- tryCatch(suppressWarnings(MASS::lda(x = x_tr, grouping = y_fac)), error = function(e) NULL)
    if (!is.null(lda)) {
      pred[["LDA"]][te] <- stats::predict(lda, x_te)$posterior[, positive]
    }
  }

  bind_rows(lapply(names(pred), function(method) {
    data.frame(
      method = method,
      sample_id = samples[train_idx],
      fold_id = fold_id,
      truth = y_vec,
      score = pred[[method]],
      check.names = FALSE
    )
  }))
}

oof_predictions <- fit_oof_candidates(x_train, y_train, top_features, negative, positive, folds, seed = seed_oof)
oof_performance <- oof_predictions %>%
  group_by(method) %>%
  group_modify(~auc_ci(.x$truth, .x$score, negative, positive)) %>%
  ungroup() %>%
  arrange(desc(auc)) %>%
  mutate(
    dataset = "Discovery 5-fold out-of-fold",
    model_context = "exploratory candidate comparison on Figure6 publication-fixed 9 proteins"
  )

final_oof <- oof_predictions %>% filter(method == "Equal-direction score")
threshold <- youden_threshold(final_oof$truth, final_oof$score, negative, positive)

score_discovery_resub <- score_fun(x_train)
score_validation <- score_fun(x_test)

final_scores <- bind_rows(
  data.frame(sample_id = samples[train_idx], dataset = "Discovery 5-fold OOF", truth = y_train, modelC_score = final_oof$score, predicted = ifelse(final_oof$score >= threshold, positive, negative), check.names = FALSE),
  data.frame(sample_id = samples[train_idx], dataset = "Discovery resubstitution", truth = y_train, modelC_score = score_discovery_resub, predicted = ifelse(score_discovery_resub >= threshold, positive, negative), check.names = FALSE),
  data.frame(sample_id = samples[test_idx], dataset = "Validation holdout", truth = y_test, modelC_score = score_validation, predicted = ifelse(score_validation >= threshold, positive, negative), check.names = FALSE)
)

perf <- bind_rows(
  cbind(data.frame(model = "ModelC_exploratory_equal_direction_score", dataset = "Discovery 5-fold OOF", threshold_source = "Discovery OOF Youden"), auc_ci(y_train, final_oof$score, negative, positive), metric_at_threshold(y_train, final_oof$score, threshold, negative, positive)),
  cbind(data.frame(model = "ModelC_exploratory_equal_direction_score", dataset = "Discovery resubstitution", threshold_source = "Discovery OOF Youden"), auc_ci(y_train, score_discovery_resub, negative, positive), metric_at_threshold(y_train, score_discovery_resub, threshold, negative, positive)),
  cbind(data.frame(model = "ModelC_exploratory_equal_direction_score", dataset = "Validation holdout", threshold_source = "Discovery OOF Youden"), auc_ci(y_test, score_validation, negative, positive), metric_at_threshold(y_test, score_validation, threshold, negative, positive))
)

roc_all <- bind_rows(
  roc_points(y_train, final_oof$score, negative, positive) %>% mutate(dataset = "Discovery 5-fold OOF"),
  roc_points(y_train, score_discovery_resub, negative, positive) %>% mutate(dataset = "Discovery resubstitution"),
  roc_points(y_test, score_validation, negative, positive) %>% mutate(dataset = "Validation holdout")
)
conf_all <- bind_rows(
  confusion_table(y_train, final_oof$score, threshold, negative, positive) %>% mutate(dataset = "Discovery 5-fold OOF"),
  confusion_table(y_train, score_discovery_resub, threshold, negative, positive) %>% mutate(dataset = "Discovery resubstitution"),
  confusion_table(y_test, score_validation, threshold, negative, positive) %>% mutate(dataset = "Validation holdout")
)

x_list <- list(Discovery = x_train[, top_features, drop = FALSE], Validation = x_test[, top_features, drop = FALSE])
candidate_scores <- fit_candidate_methods(x_train[, top_features, drop = FALSE], y_train, x_list, negative, positive, seed = 20250608, equal_direction_weights = direction_weights)
candidate_validation_perf <- method_performance_table(candidate_scores, list(Discovery = y_train, Validation = y_test), "Discovery", negative, positive) %>%
  mutate(model_context = "descriptive only; Model C is exploratory and not a clinical classifier")
candidate_validation_predictions <- bind_rows(lapply(names(candidate_scores), function(method) {
  bind_rows(lapply(names(candidate_scores[[method]]), function(dataset) {
    data.frame(method = method, dataset = dataset, sample_id = if (dataset == "Discovery") samples[train_idx] else samples[test_idx], truth = if (dataset == "Discovery") y_train else y_test, score = candidate_scores[[method]][[dataset]], check.names = FALSE)
  }))
}))

expr_top_z <- scale_train_apply(x_all[train_idx, top_features, drop = FALSE], x_all[, top_features, drop = FALSE])$new
heatmap_df <- data.frame(sample_id = samples, expr_top_z, check.names = FALSE) %>%
  mutate(model_set = ifelse(seq_along(samples) %in% train_idx, "Discovery", "Validation"), truth = y) %>%
  pivot_longer(cols = all_of(top_features), names_to = "model_feature", values_to = "z_abundance") %>%
  left_join(feature_table %>% dplyr::select(model_feature, genesymbol, model_rank, biology_annotation, pathway_id, feature_group), by = "model_feature")

feature_expression_summary <- heatmap_df %>%
  group_by(model_feature, genesymbol, model_rank, biology_annotation, pathway_id, feature_group, model_set, truth) %>%
  summarise(mean_z = mean(z_abundance, na.rm = TRUE), median_z = median(z_abundance, na.rm = TRUE), .groups = "drop")

split_table <- bind_rows(
  clin_c[train_idx, ] %>% mutate(model_set = "Discovery"),
  clin_c[test_idx, ] %>% mutate(model_set = "Validation")
) %>%
  dplyr::select(sample_id, model_set, truth, histology, sex, age, smoking, drinking, smoking_status, drinking_status)

sclc_discovery_de <- limma_de(
  expr = expr,
  clinical = clinical,
  samples = samples[train_idx],
  group_var = "n_status",
  positive = positive,
  negative = negative,
  covariates = c("age", "sex", "smoking_status", "drinking_status")
)
sclc_discovery_enrichment <- run_enrichment(
  sclc_discovery_de,
  universe_symbols = rownames(expr),
  prefix = "SCLC_discovery_Nplus_vs_N0_exploratory",
  out_dir = file.path(dir_map$enrichment, "SCLC_discovery_Nplus_vs_N0_exploratory")
)

pathway_rows <- function(tbl, source_label, value_col = NULL, count_col = "Count") {
  if (is.null(tbl) || nrow(tbl) == 0 || !"p.adjust" %in% names(tbl)) {
    return(data.frame(source = character(), ID = character(), Description = character(), score = numeric(), value = numeric(), p.adjust = numeric(), Count = numeric(), label = character(), check.names = FALSE))
  }
  value <- if (!is.null(value_col) && value_col %in% names(tbl)) abs(tbl[[value_col]]) else -log10(pmax(tbl$p.adjust, .Machine$double.xmin))
  count <- if (count_col %in% names(tbl)) tbl[[count_col]] else if ("setSize" %in% names(tbl)) tbl$setSize else NA_real_
  data.frame(
    source = source_label,
    ID = as.character(tbl$ID),
    Description = as.character(tbl$Description),
    score = -log10(pmax(tbl$p.adjust, .Machine$double.xmin)),
    value = value,
    p.adjust = tbl$p.adjust,
    Count = count,
    label = paste(as.character(tbl$ID), as.character(tbl$Description), sep = "\n"),
    check.names = FALSE
  )
}

sclc_top_pathways <- bind_rows(
  pathway_rows(sclc_discovery_enrichment$go, "GO BP ORA", value_col = "Count"),
  pathway_rows(sclc_discovery_enrichment$reactome, "Reactome ORA", value_col = "Count"),
  pathway_rows(sclc_discovery_enrichment$gsea, "Reactome GSEA", value_col = "NES", count_col = "setSize")
) %>%
  filter(!is.na(p.adjust), is.finite(score)) %>%
  arrange(p.adjust, desc(score)) %>%
  group_by(source) %>%
  slice_head(n = 4) %>%
  ungroup() %>%
  arrange(p.adjust, desc(score)) %>%
  slice_head(n = 8)

write_csv_safe(split_table, file.path(dir_map$modelC, "modelC_sample_split.csv"))
write_csv_safe(feature_screen, file.path(dir_map$modelC, "modelC_discovery_feature_screen_all.csv"))
write_csv_safe(feature_table, file.path(dir_map$modelC, "modelC_feature_table.csv"))
write_csv_safe(feature_table, file.path(dir_map$modelC, "modelC_feature_table_top9.csv"))
write_csv_safe(oof_predictions, file.path(dir_map$modelC, "modelC_discovery_5fold_oof_predictions.csv"))
write_csv_safe(oof_performance, file.path(dir_map$modelC, "modelC_discovery_5fold_oof_auc.csv"))
write_csv_safe(final_scores, file.path(dir_map$modelC, "modelC_final_prediction_scores.csv"))
write_csv_safe(perf, file.path(dir_map$modelC, "modelC_final_performance.csv"))
write_csv_safe(roc_all, file.path(dir_map$modelC, "modelC_final_roc_coordinates.csv"))
write_csv_safe(conf_all, file.path(dir_map$modelC, "modelC_final_confusion_matrices.csv"))
write_csv_safe(candidate_validation_perf, file.path(dir_map$modelC, "modelC_candidate_model_performance.csv"))
write_csv_safe(candidate_validation_predictions, file.path(dir_map$modelC, "modelC_candidate_model_predictions.csv"))
write_csv_safe(heatmap_df, file.path(dir_map$modelC, "modelC_top9_heatmap_data.csv"))
write_csv_safe(feature_expression_summary, file.path(dir_map$modelC, "modelC_top9_feature_expression_summary.csv"))
write_csv_safe(sclc_discovery_de, file.path(dir_map$modelC, "modelC_sclc_discovery_limma_all_proteins.csv"))
write_csv_safe(sclc_discovery_enrichment$selected_genes, file.path(dir_map$modelC, "modelC_sclc_discovery_enrichment_input_genes.csv"))
write_csv_safe(sclc_top_pathways, file.path(dir_map$modelC, "modelC_sclc_top_pathways.csv"))

write_json_safe(list(
  final_model = "exploratory SCLC nodal involvement nine-protein equal-direction score",
  clinical_claim_boundary = "Exploratory only; not presented as a clinical-grade classifier.",
  training_counts = list(N0 = sum(y_train == negative), Nplus = sum(y_train == positive)),
  validation_counts = list(N0 = sum(y_test == negative), Nplus = sum(y_test == positive)),
  feature_selection = "Figure6 publication-fixed nine proteins; discovery-wide Wilcoxon feature screen is saved as audit, not used to replace the fixed feature set.",
  selected_features = feature_table$genesymbol,
  transformation = "Final equal-direction score uses raw protein abundance; heatmaps use discovery-standardized row z-scores for display only.",
  out_of_fold_auc = "5-fold stratified out-of-fold AUC in discovery cohort.",
  threshold = threshold,
  threshold_source = "Discovery 5-fold OOF Youden index",
  validation_leakage_control = "Validation cohort is not used for feature selection, threshold selection, or model tuning.",
  seed_split = seed_split,
  seed_oof = seed_oof,
  folds = folds
), file.path(dir_map$modelC, "modelC_parameters.json"))

cat("Model C completed.\n")
print(perf)
print(oof_performance)
