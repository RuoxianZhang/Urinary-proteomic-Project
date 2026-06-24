source(file.path(Sys.getenv("UR_FILE0624_PROJECT_DIR", unset = Sys.getenv("UR_PROJECT_DIR", unset = Sys.getenv("UR_FILE0623_PROJECT_DIR", unset = getwd()))), "code", "00_common_functions.R"))

ensure_output_dirs()

dat <- read_project_data()
clinical <- dat$clinical
expr <- dat$expr
universe <- rownames(expr)

run_and_save_de <- function(name, samples, group_var, positive, negative, covariates, enrich = TRUE) {
  res <- limma_de(expr, clinical, samples, group_var, positive, negative, covariates)
  out_path <- file.path(dir_map$differential, paste0(name, "_limma_all_proteins.csv"))
  write_csv_safe(res, out_path)
  nominal <- res %>% filter(direction_nominal != "not_selected")
  write_csv_safe(nominal, file.path(dir_map$differential, paste0(name, "_nominal_candidates.csv")))
  covariate_df <- clinical[match(samples, clinical$sample_id), , drop = FALSE]
  used_covariates <- clean_design_covariates(covariate_df, covariates)
  summary <- data.frame(
    comparison = name,
    samples = length(samples),
    negative = negative,
    positive = positive,
    proteins_tested = nrow(res),
    nominal_positive_up = sum(res$direction_nominal == "positive_up"),
    nominal_negative_up = sum(res$direction_nominal == "negative_up"),
    fdr_significant = sum(res$direction_fdr != "not_fdr_significant", na.rm = TRUE),
    covariates = paste(used_covariates, collapse = ";"),
    check.names = FALSE
  )
  write_csv_safe(summary, file.path(dir_map$differential, paste0(name, "_summary.csv")))
  if (enrich) {
    run_enrichment(
      res_table = res,
      universe_symbols = universe,
      prefix = name,
      out_dir = file.path(dir_map$enrichment, name),
      p_cutoff = 0.05,
      lfc_cutoff = 0.3,
      min_gs = 10,
      max_gs = 500
    )
  }
  list(res = res, nominal = nominal, summary = summary)
}

lc_samples <- clinical$sample_id[clinical$group %in% c("Healthy", "LungCancer")]
lc <- clinical %>%
  mutate(lc_status = ifelse(group == "LungCancer", "LungCancer", "Healthy"))
lc_clinical_backup <- clinical
clinical$lc_status <- lc$lc_status
lc_res <- run_and_save_de(
  "LC_vs_Healthy",
  lc_samples,
  group_var = "lc_status",
  positive = "LungCancer",
  negative = "Healthy",
  covariates = c("age", "sex", "smoking_status", "drinking_status")
)

progression_samples <- clinical$sample_id[clinical$group == "Healthy" | (clinical$major_type == "NSCLC" & clinical$n_status %in% c("N0", "Nplus"))]
clinical$progression3 <- case_when(
  clinical$group == "Healthy" ~ "Healthy",
  clinical$major_type == "NSCLC" & clinical$n_status == "N0" ~ "NSCLC_N0",
  clinical$major_type == "NSCLC" & clinical$n_status == "Nplus" ~ "NSCLC_Nplus",
  TRUE ~ NA_character_
)

hn0_samples <- clinical$sample_id[clinical$progression3 %in% c("Healthy", "NSCLC_N0")]
hn0 <- run_and_save_de(
  "NSCLC_N0_vs_Healthy",
  hn0_samples,
  group_var = "progression3",
  positive = "NSCLC_N0",
  negative = "Healthy",
  covariates = c("age", "sex")
)

nsclc_samples <- clinical$sample_id[clinical$major_type == "NSCLC" & clinical$n_status %in% c("N0", "Nplus")]
nsclc <- run_and_save_de(
  "NSCLC_Nplus_vs_N0",
  nsclc_samples,
  group_var = "n_status",
  positive = "Nplus",
  negative = "N0",
  covariates = c("histology", "age", "sex", "smoking_status", "drinking_status")
)

luad_samples <- clinical$sample_id[clinical$histology == "LUAD" & clinical$n_status %in% c("N0", "Nplus")]
luad <- run_and_save_de(
  "LUAD_Nplus_vs_N0",
  luad_samples,
  group_var = "n_status",
  positive = "Nplus",
  negative = "N0",
  covariates = c("age", "sex", "smoking_status", "drinking_status")
)

lusc_samples <- clinical$sample_id[clinical$histology == "LUSC" & clinical$n_status %in% c("N0", "Nplus")]
lusc <- run_and_save_de(
  "LUSC_Nplus_vs_N0",
  lusc_samples,
  group_var = "n_status",
  positive = "Nplus",
  negative = "N0",
  covariates = c("age", "sex", "smoking_status", "drinking_status")
)

sclc_samples <- clinical$sample_id[clinical$major_type == "SCLC" & clinical$n_status %in% c("N0", "Nplus")]
sclc_all <- run_and_save_de(
  "SCLC_Nplus_vs_N0_all_exploratory",
  sclc_samples,
  group_var = "n_status",
  positive = "Nplus",
  negative = "N0",
  covariates = c("age", "sex", "smoking_status", "drinking_status")
)

write_csv_safe(bind_rows(lc_res$summary, hn0$summary, nsclc$summary, luad$summary, lusc$summary, sclc_all$summary), file.path(dir_map$tables, "differential_analysis_summary.csv"))

score_progression <- module_scores(expr, progression_samples) %>%
  left_join(clinical %>% dplyr::select(sample_id, progression3, group, major_type, histology, n_status, sex, age, smoking, drinking, smoking_status, drinking_status), by = "sample_id")
score_nsclc <- module_scores(expr, nsclc_samples) %>%
  left_join(clinical %>% dplyr::select(sample_id, n_status, histology, sex, age, smoking, drinking, smoking_status, drinking_status), by = "sample_id")
score_luad <- module_scores(expr, luad_samples) %>%
  left_join(clinical %>% dplyr::select(sample_id, n_status, histology, sex, age, smoking, drinking, smoking_status, drinking_status), by = "sample_id")
score_lusc <- module_scores(expr, lusc_samples) %>%
  left_join(clinical %>% dplyr::select(sample_id, n_status, histology, sex, age, smoking, drinking, smoking_status, drinking_status), by = "sample_id")
score_lc <- module_scores(expr, lc_samples) %>%
  left_join(clinical %>% dplyr::select(sample_id, group, major_type, histology, n_status, sex, age, smoking, drinking, smoking_status, drinking_status), by = "sample_id")

write_csv_safe(score_progression, file.path(dir_map$differential, "progression_curated_module_scores.csv"))
write_csv_safe(score_nsclc, file.path(dir_map$differential, "NSCLC_curated_module_scores.csv"))
write_csv_safe(score_luad, file.path(dir_map$differential, "LUAD_curated_module_scores.csv"))
write_csv_safe(score_lusc, file.path(dir_map$differential, "LUSC_curated_module_scores.csv"))
write_csv_safe(score_lc, file.path(dir_map$differential, "LC_vs_Healthy_curated_module_scores.csv"))

stats_nsclc <- module_stats(score_nsclc, clinical, "n_status", "Nplus", "N0")
stats_luad <- module_stats(score_luad, clinical, "n_status", "Nplus", "N0")
stats_lusc <- module_stats(score_lusc, clinical, "n_status", "Nplus", "N0")
write_csv_safe(stats_nsclc, file.path(dir_map$differential, "NSCLC_curated_module_statistics.csv"))
write_csv_safe(stats_luad, file.path(dir_map$differential, "LUAD_curated_module_statistics.csv"))
write_csv_safe(stats_lusc, file.path(dir_map$differential, "LUSC_curated_module_statistics.csv"))

fit_module_effect <- function(df, term, covariates) {
  covs <- clean_design_covariates(df, covariates)
  formula_text <- paste("module_score ~", term, if (length(covs)) paste("+", paste(covs, collapse = "+")) else "")
  fit <- lm(as.formula(formula_text), data = df)
  sm <- summary(fit)$coefficients
  coef_name <- grep(paste0("^", term), rownames(sm), value = TRUE)[1]
  data.frame(
    adjusted_beta = sm[coef_name, "Estimate"],
    adjusted_p = sm[coef_name, "Pr(>|t|)"],
    adjusted_formula = formula_text,
    covariates = paste(covs, collapse = ";"),
    check.names = FALSE
  )
}

adjust_lc_module <- function(scores, covariates) {
  bind_rows(lapply(unique(scores$module), function(m) {
    df <- scores %>%
      filter(module == m, group %in% c("Healthy", "LungCancer")) %>%
      mutate(group_factor = factor(as.character(group), levels = c("Healthy", "LungCancer")))
    effect <- fit_module_effect(df, "group_factor", covariates)
    data.frame(
      module = m,
      adjusted_beta_LungCancer = effect$adjusted_beta,
      adjusted_p = effect$adjusted_p,
      adjusted_formula = effect$adjusted_formula,
      covariates = effect$covariates,
      median_Healthy = median(df$module_score[df$group_factor == "Healthy"], na.rm = TRUE),
      median_LungCancer = median(df$module_score[df$group_factor == "LungCancer"], na.rm = TRUE),
      delta_median_LungCancer_minus_Healthy = median(df$module_score[df$group_factor == "LungCancer"], na.rm = TRUE) - median(df$module_score[df$group_factor == "Healthy"], na.rm = TRUE),
      wilcox_p = suppressWarnings(wilcox.test(module_score ~ group_factor, data = df)$p.value),
      check.names = FALSE
    )
  })) %>%
    mutate(
      adjusted_FDR = p.adjust(adjusted_p, method = "BH"),
      wilcox_FDR = p.adjust(wilcox_p, method = "BH")
    ) %>%
    arrange(adjusted_FDR)
}

adjust_module <- function(scores, label, covariates) {
  bind_rows(lapply(unique(scores$module), function(m) {
    df <- scores %>%
      filter(module == m, n_status %in% c("N0", "Nplus")) %>%
      mutate(n_status = factor(as.character(n_status), levels = c("N0", "Nplus")))
    effect <- fit_module_effect(df, "n_status", covariates)
    data.frame(
      subtype = label,
      subset = label,
      module = m,
      adjusted_beta = effect$adjusted_beta,
      adjusted_beta_Nplus = effect$adjusted_beta,
      p_value = effect$adjusted_p,
      adjusted_p = effect$adjusted_p,
      adjusted_formula = effect$adjusted_formula,
      covariates = effect$covariates,
      check.names = FALSE
    )
  })) %>%
    group_by(subtype) %>%
    mutate(FDR = p.adjust(p_value, method = "BH"), adjusted_FDR = FDR) %>%
    ungroup()
}

module_effect_lc <- adjust_lc_module(score_lc, c("age", "sex", "smoking_status", "drinking_status"))
write_csv_safe(module_effect_lc, file.path(dir_map$differential, "LC_vs_Healthy_curated_module_stats.csv"))

module_effects <- bind_rows(
  adjust_module(score_progression %>% filter(major_type == "NSCLC"), "NSCLC", c("histology", "age", "sex", "smoking_status", "drinking_status")),
  adjust_module(score_progression %>% filter(histology == "LUAD"), "LUAD", c("age", "sex", "smoking_status", "drinking_status")),
  adjust_module(score_progression %>% filter(histology == "LUSC"), "LUSC", c("age", "sex", "smoking_status", "drinking_status"))
)
write_csv_safe(module_effects, file.path(dir_map$differential, "pathway_module_adjusted_effects_by_subtype.csv"))

effect_compare <- nsclc$res %>%
  transmute(gene, nsclc_logFC = logFC, nsclc_P.Value = P.Value, nsclc_adj.P.Val = adj.P.Val) %>%
  inner_join(luad$res %>% transmute(gene, luad_logFC = logFC, luad_P.Value = P.Value, luad_adj.P.Val = adj.P.Val), by = "gene") %>%
  inner_join(lusc$res %>% transmute(gene, lusc_logFC = logFC, lusc_P.Value = P.Value, lusc_adj.P.Val = adj.P.Val), by = "gene") %>%
  mutate(
    nsclc_selected = nsclc_P.Value < 0.05 & abs(nsclc_logFC) > 0.3,
    luad_selected = luad_P.Value < 0.05 & abs(luad_logFC) > 0.3,
    lusc_selected = lusc_P.Value < 0.05 & abs(lusc_logFC) > 0.3,
    selected_luad_nsclc = nsclc_selected | luad_selected,
    selected_lusc_nsclc = nsclc_selected | lusc_selected,
    luad_same_direction = sign(nsclc_logFC) == sign(luad_logFC),
    lusc_same_direction = sign(nsclc_logFC) == sign(lusc_logFC),
    luad_evidence = case_when(
      nsclc_selected & luad_selected ~ "both_nominal",
      nsclc_selected ~ "NSCLC_nominal",
      luad_selected ~ "LUAD_nominal",
      TRUE ~ "not_selected"
    ),
    lusc_evidence = case_when(
      nsclc_selected & lusc_selected ~ "both_nominal",
      nsclc_selected ~ "NSCLC_nominal",
      lusc_selected ~ "LUSC_nominal",
      TRUE ~ "not_selected"
    )
  )
write_csv_safe(effect_compare, file.path(dir_map$differential, "subtype_vs_nsclc_effect_size_comparison.csv"))

cat("Differential, module, and enrichment analyses completed.\\n")
print(read.csv(file.path(dir_map$tables, "differential_analysis_summary.csv"), check.names = FALSE))
