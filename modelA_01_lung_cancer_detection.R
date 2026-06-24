source(file.path(Sys.getenv("UR_FILE0624_PROJECT_DIR", unset = Sys.getenv("UR_PROJECT_DIR", unset = Sys.getenv("UR_FILE0623_PROJECT_DIR", unset = getwd()))), "code", "00_common_functions.R"))

ensure_output_dirs()
set.seed(20250613)

dat <- read_project_data()
clinical <- dat$clinical
expr <- dat$expr

markers <- c("ANXA11", "APOA2", "NAPSA", "ATP1A3", "RAB1B")
missing_markers <- setdiff(markers, rownames(expr))
if (length(missing_markers) > 0) stop("Model A markers missing from raw expression: ", paste(missing_markers, collapse = ", "))

negative <- "Healthy"
positive <- "LungCancer"

clin_a <- clinical %>%
  filter(group %in% c(negative, positive), sample_id %in% colnames(expr)) %>%
  arrange(group, sample_id) %>%
  group_by(group) %>%
  mutate(
    class = as.character(group),
    split_rank = row_number(),
    model_set = case_when(
      group == negative & split_rank <= 92 ~ "Discovery",
      group == positive & split_rank <= 132 ~ "Discovery",
      TRUE ~ "Internal validation"
    )
  ) %>%
  ungroup()

samples <- clin_a$sample_id
y <- clin_a$class
discovery_idx <- which(clin_a$model_set == "Discovery")
validation_idx <- which(clin_a$model_set == "Internal validation")

split_table <- clin_a %>%
  dplyr::select(sample_id, model_set, class, group, major_type, histology, n_status, sex, age, smoking, drinking, split_rank)

x_raw <- t(expr[markers, samples, drop = FALSE])
colnames(x_raw) <- markers
x_rank <- rank_normalize_within(x_raw)
rownames(x_rank) <- samples

w <- direction_weights(
  x_rank[discovery_idx, markers, drop = FALSE],
  y[discovery_idx],
  markers,
  positive,
  negative
)
score <- linear_score(x_rank, w)

score_discovery <- score[discovery_idx]
score_internal <- score[validation_idx]
threshold <- youden_threshold(y[discovery_idx], score_discovery, negative, positive)

zhang <- parse_zhang_data(markers = markers)
zhang_clin <- zhang$clinical %>% filter(group %in% c(negative, positive))
zhang_raw <- t(zhang$expr[markers, zhang_clin$sample_id, drop = FALSE])
colnames(zhang_raw) <- markers
zhang_rank <- rank_normalize_within(zhang_raw)
rownames(zhang_rank) <- zhang_clin$sample_id
score_zhang <- linear_score(zhang_rank, w)

scores <- bind_rows(
  data.frame(sample_id = samples[discovery_idx], dataset = "Discovery", truth = y[discovery_idx], modelA_score = score_discovery, check.names = FALSE),
  data.frame(sample_id = samples[validation_idx], dataset = "Internal validation", truth = y[validation_idx], modelA_score = score_internal, check.names = FALSE),
  data.frame(sample_id = zhang_clin$sample_id, dataset = "Zhang validation", truth = zhang_clin$group, modelA_score = score_zhang, check.names = FALSE)
) %>%
  mutate(predicted = ifelse(modelA_score >= threshold, positive, negative))

perf <- bind_rows(
  cbind(data.frame(model = "ModelA_final_equal_direction_score", dataset = "Discovery", threshold_source = "Discovery Youden"), auc_ci(y[discovery_idx], score_discovery, negative, positive), metric_at_threshold(y[discovery_idx], score_discovery, threshold, negative, positive)),
  cbind(data.frame(model = "ModelA_final_equal_direction_score", dataset = "Internal validation", threshold_source = "Discovery Youden"), auc_ci(y[validation_idx], score_internal, negative, positive), metric_at_threshold(y[validation_idx], score_internal, threshold, negative, positive)),
  cbind(data.frame(model = "ModelA_final_equal_direction_score", dataset = "Corrected original cohort", threshold_source = "Discovery Youden"), auc_ci(y, score, negative, positive), metric_at_threshold(y, score, threshold, negative, positive)),
  cbind(data.frame(model = "ModelA_final_equal_direction_score", dataset = "Zhang validation", threshold_source = "Discovery Youden"), auc_ci(zhang_clin$group, score_zhang, negative, positive), metric_at_threshold(zhang_clin$group, score_zhang, threshold, negative, positive))
)

roc_all <- bind_rows(
  roc_points(y[discovery_idx], score_discovery, negative, positive) %>% mutate(dataset = "Discovery"),
  roc_points(y[validation_idx], score_internal, negative, positive) %>% mutate(dataset = "Internal validation"),
  roc_points(zhang_clin$group, score_zhang, negative, positive) %>% mutate(dataset = "Zhang validation")
)

conf_all <- bind_rows(
  confusion_table(y[discovery_idx], score_discovery, threshold, negative, positive) %>% mutate(dataset = "Discovery"),
  confusion_table(y[validation_idx], score_internal, threshold, negative, positive) %>% mutate(dataset = "Internal validation"),
  confusion_table(zhang_clin$group, score_zhang, threshold, negative, positive) %>% mutate(dataset = "Zhang validation")
)

marker_direction <- data.frame(
  genesymbol = markers,
  direction_weight = as.numeric(w[markers]),
  direction_label = ifelse(w[markers] > 0, "Higher in lung cancer", "Lower in lung cancer"),
  discovery_mean_healthy = colMeans(x_rank[discovery_idx, markers, drop = FALSE][y[discovery_idx] == negative, , drop = FALSE]),
  discovery_mean_lung_cancer = colMeans(x_rank[discovery_idx, markers, drop = FALSE][y[discovery_idx] == positive, , drop = FALSE]),
  check.names = FALSE
)

marker_delta_table <- function(mat, clin, dataset_label) {
  clin_use <- clin %>%
    filter(group %in% c(negative, positive), sample_id %in% colnames(mat)) %>%
    mutate(group = factor(as.character(group), levels = c(negative, positive)))
  z <- row_z(mat[markers, clin_use$sample_id, drop = FALSE])
  long <- data.frame(gene = rownames(z), z, check.names = FALSE) %>%
    tidyr::pivot_longer(cols = all_of(clin_use$sample_id), names_to = "sample_id", values_to = "z_abundance") %>%
    left_join(clin_use %>% dplyr::select(sample_id, group), by = "sample_id")
  long %>%
    group_by(gene) %>%
    summarise(
      dataset = dataset_label,
      mean_Healthy_z = mean(z_abundance[group == negative], na.rm = TRUE),
      mean_LungCancer_z = mean(z_abundance[group == positive], na.rm = TRUE),
      standardized_delta_LC_minus_Healthy = mean_LungCancer_z - mean_Healthy_z,
      wilcox_p = suppressWarnings(wilcox.test(z_abundance ~ group)$p.value),
      .groups = "drop"
    ) %>%
    mutate(
      FDR = p.adjust(wilcox_p, method = "BH"),
      genesymbol = gene
    ) %>%
    dplyr::select(dataset, gene, genesymbol, mean_Healthy_z, mean_LungCancer_z, standardized_delta_LC_minus_Healthy, wilcox_p, FDR)
}

marker_consistency <- bind_rows(
  marker_delta_table(expr, clinical, "Corrected original cohort"),
  marker_delta_table(zhang$expr, zhang_clin, "Zhang validation cohort")
) %>%
  left_join(marker_direction %>% transmute(gene = genesymbol, direction_weight, discovery_direction_label = direction_label), by = "gene") %>%
  mutate(
    observed_direction = ifelse(standardized_delta_LC_minus_Healthy >= 0, "Higher in lung cancer", "Lower in lung cancer"),
    concordant_with_discovery_direction = sign(standardized_delta_LC_minus_Healthy) == sign(direction_weight)
  ) %>%
  arrange(match(gene, markers), dataset)

truth_list <- list(
  Discovery = y[discovery_idx],
  `Internal validation` = y[validation_idx],
  `Zhang validation` = zhang_clin$group
)
x_list <- list(
  Discovery = x_rank[discovery_idx, markers, drop = FALSE],
  `Internal validation` = x_rank[validation_idx, markers, drop = FALSE],
  `Zhang validation` = zhang_rank[, markers, drop = FALSE]
)
candidate_scores <- fit_candidate_methods(
  x_rank[discovery_idx, markers, drop = FALSE],
  y[discovery_idx],
  x_list,
  negative,
  positive,
  seed = 20250613,
  equal_direction_weights = w
)
candidate_perf <- method_performance_table(candidate_scores, truth_list, "Discovery", negative, positive) %>%
  mutate(model_context = "candidate model screen; final Model A remains discovery-derived signed/equal-direction score")

candidate_pred <- bind_rows(lapply(names(candidate_scores), function(method) {
  bind_rows(lapply(names(candidate_scores[[method]]), function(dataset) {
    data.frame(
      method = method,
      dataset = dataset,
      sample_id = if (dataset == "Discovery") samples[discovery_idx] else if (dataset == "Internal validation") samples[validation_idx] else zhang_clin$sample_id,
      truth = truth_list[[dataset]],
      score = candidate_scores[[method]][[dataset]],
      check.names = FALSE
    )
  }))
}))

single_marker_cv <- function(feature_set, name, folds = 5, repeats = 20) {
  set.seed(20250611 + length(feature_set))
  out <- vector("list", repeats)
  y_train <- y[discovery_idx]
  x <- x_rank[discovery_idx, feature_set, drop = FALSE]
  for (r in seq_len(repeats)) {
    fold_id <- rep(NA_integer_, length(y_train))
    for (cl in unique(y_train)) {
      idx <- which(y_train == cl)
      fold_id[idx] <- sample(rep(seq_len(folds), length.out = length(idx)))
    }
    pred <- rep(NA_real_, length(y_train))
    for (k in seq_len(folds)) {
      tr <- which(fold_id != k)
      te <- which(fold_id == k)
      rf <- randomForest::randomForest(
        x = x[tr, , drop = FALSE],
        y = factor(y_train[tr], levels = c(negative, positive)),
        ntree = 300
      )
      pred[te] <- predict(rf, x[te, , drop = FALSE], type = "prob")[, positive]
    }
    out[[r]] <- cbind(data.frame(feature_set = name, repeat_id = r), auc_ci(y_train, pred, negative, positive))
  }
  bind_rows(out)
}

feature_sets <- c(list(`Five-protein panel` = markers), setNames(as.list(markers), markers))
feature_set_cv <- bind_rows(lapply(names(feature_sets), function(nm) single_marker_cv(feature_sets[[nm]], nm))) %>%
  group_by(feature_set) %>%
  summarise(
    mean_cv_auc = mean(auc),
    low_2.5 = quantile(auc, 0.025),
    high_97.5 = quantile(auc, 0.975),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_cv_auc))

heatmap_df <- bind_rows(
  data.frame(sample_id = samples[discovery_idx], dataset = "Discovery", x_rank[discovery_idx, markers, drop = FALSE], check.names = FALSE),
  data.frame(sample_id = samples[validation_idx], dataset = "Internal validation", x_rank[validation_idx, markers, drop = FALSE], check.names = FALSE),
  data.frame(sample_id = zhang_clin$sample_id, dataset = "Zhang validation", zhang_rank[, markers, drop = FALSE], check.names = FALSE)
) %>%
  pivot_longer(cols = all_of(markers), names_to = "genesymbol", values_to = "rank_normalized_abundance") %>%
  left_join(marker_direction %>% dplyr::select(genesymbol, direction_weight), by = "genesymbol") %>%
  mutate(direction_adjusted_z = rank_normalized_abundance * direction_weight)

write_csv_safe(split_table, file.path(dir_map$modelA, "modelA_sample_split.csv"))
write_csv_safe(marker_direction, file.path(dir_map$modelA, "modelA_feature_table.csv"))
write_csv_safe(marker_direction, file.path(dir_map$modelA, "modelA_feature_table_direction_weights.csv"))
write_csv_safe(scores, file.path(dir_map$modelA, "modelA_final_prediction_scores.csv"))
write_csv_safe(perf, file.path(dir_map$modelA, "modelA_final_performance.csv"))
write_csv_safe(roc_all, file.path(dir_map$modelA, "modelA_final_roc_coordinates.csv"))
write_csv_safe(conf_all, file.path(dir_map$modelA, "modelA_final_confusion_matrices.csv"))
write_csv_safe(candidate_perf, file.path(dir_map$modelA, "modelA_candidate_model_performance.csv"))
write_csv_safe(candidate_pred, file.path(dir_map$modelA, "modelA_candidate_model_predictions.csv"))
write_csv_safe(feature_set_cv, file.path(dir_map$modelA, "modelA_rf_feature_set_cross_validation.csv"))
write_csv_safe(heatmap_df, file.path(dir_map$modelA, "modelA_direction_adjusted_heatmap_data.csv"))
write_csv_safe(marker_consistency, file.path(dir_map$modelA, "modelA_marker_direction_consistency_original_vs_Zhang.csv"))

write_json_safe(list(
  final_model = "discovery-derived signed/equal-direction five-protein score",
  features = markers,
  transformation = "rank-normalization within the corrected original LC/healthy cohort; Zhang markers rank-normalized within the independent LC/CTL cohort",
  score = "linear score with L2-normalized direction weights",
  direction_weights = as.list(as.numeric(w[markers])),
  direction_weight_names = markers,
  threshold = threshold,
  threshold_source = "Discovery set Youden index",
  split_rule = "Clinical samples arranged by group and sample_id; first 92 Healthy and first 132 LungCancer assigned to discovery, remaining 23 Healthy and 33 LungCancer to internal validation.",
  candidate_model_screen = "random forest, radial SVM, logistic regression, LASSO logistic regression, elastic net, LDA, and equal-direction score; not used to redefine final Model A"
), file.path(dir_map$modelA, "modelA_parameters.json"))

cat("Model A completed.\n")
print(perf)
