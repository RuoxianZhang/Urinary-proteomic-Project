project_dir <- normalizePath(Sys.getenv("UR_FILE0624_PROJECT_DIR", unset = Sys.getenv("UR_PROJECT_DIR", unset = Sys.getenv("UR_FILE0623_PROJECT_DIR", unset = getwd()))), mustWork = TRUE)
result_dir <- file.path(project_dir, "result_file")
fig_dir <- file.path(result_dir, "figures")
table_dir <- file.path(result_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

read_result <- function(path) read.csv(file.path(result_dir, path), check.names = FALSE)
fmt <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, sprintf(paste0("%.", digits, "f"), as.numeric(x)))
}

standard_perf <- function(tbl, model_label) {
  data.frame(
    model = model_label,
    dataset = tbl$dataset,
    AUC = fmt(tbl$auc),
    CI95 = paste0(fmt(tbl$ci_low), "-", fmt(tbl$ci_high)),
    threshold = fmt(tbl$threshold),
    TP = tbl$tp,
    FP = tbl$fp,
    FN = tbl$fn,
    TN = tbl$tn,
    accuracy = fmt(tbl$accuracy),
    sensitivity = fmt(tbl$sensitivity),
    specificity = fmt(tbl$specificity),
    balanced_accuracy = fmt(tbl$balanced_accuracy),
    stringsAsFactors = FALSE
  )
}

modelA_perf <- read_result("modelA/modelA_final_performance.csv")
modelB_perf <- read_result("modelB/modelB_final_performance.csv")
modelC_perf <- read_result("modelC/modelC_final_performance.csv")

final_perf <- rbind(
  standard_perf(modelA_perf, "Model A final five-protein score"),
  standard_perf(modelB_perf, "Model B final RF top18"),
  standard_perf(modelC_perf, "Model C exploratory nine-protein score")
)
write.csv(final_perf, file.path(table_dir, "modelABC_final_performance_readable.csv"), row.names = FALSE)

candidate_auc <- rbind(
  transform(read_result("modelA/modelA_candidate_model_performance.csv")[, c("method", "dataset", "auc", "ci_low", "ci_high", "threshold")], model = "Model A"),
  transform(read_result("modelB/modelB_candidate_model_performance.csv")[, c("method", "dataset", "auc", "ci_low", "ci_high", "threshold")], model = "Model B"),
  transform(read_result("modelC/modelC_discovery_5fold_oof_auc.csv")[, c("method", "dataset", "auc", "ci_low", "ci_high")], threshold = NA_real_, model = "Model C")
)
candidate_auc <- candidate_auc[, c("model", "method", "dataset", "auc", "ci_low", "ci_high", "threshold")]
candidate_auc$AUC <- fmt(candidate_auc$auc)
candidate_auc$CI95 <- paste0(fmt(candidate_auc$ci_low), "-", fmt(candidate_auc$ci_high))
candidate_auc$threshold_display <- fmt(candidate_auc$threshold)
write.csv(candidate_auc, file.path(table_dir, "modelABC_candidate_auc_readable.csv"), row.names = FALSE)

score_summary_one <- function(path, model, score_col, set_col) {
  x <- read_result(path)
  split_vals <- split(x[[score_col]], x[[set_col]])
  do.call(rbind, lapply(names(split_vals), function(dataset) {
    v <- split_vals[[dataset]]
    data.frame(
      model = model,
      dataset = dataset,
      n = length(v),
      min = fmt(min(v, na.rm = TRUE)),
      q1 = fmt(as.numeric(quantile(v, 0.25, na.rm = TRUE))),
      median = fmt(median(v, na.rm = TRUE)),
      mean = fmt(mean(v, na.rm = TRUE)),
      q3 = fmt(as.numeric(quantile(v, 0.75, na.rm = TRUE))),
      max = fmt(max(v, na.rm = TRUE)),
      score_column = score_col,
      source_file = path,
      stringsAsFactors = FALSE
    )
  }))
}

score_summary <- rbind(
  score_summary_one("modelA/modelA_final_prediction_scores.csv", "Model A", "modelA_score", "dataset"),
  score_summary_one("modelB/modelB_discovery_oob_predictions.csv", "Model B discovery", "rf_top18_oob_prob", "model_set"),
  score_summary_one("modelB/modelB_validation_predictions.csv", "Model B validation", "rf_top18_prob", "model_set"),
  score_summary_one("modelC/modelC_final_prediction_scores.csv", "Model C", "modelC_score", "dataset")
)
write.csv(score_summary, file.path(table_dir, "modelABC_score_distribution_summary.csv"), row.names = FALSE)

feature_lines <- c(
  "Model A features: ANXA11(+), APOA2(-), NAPSA(-), ATP1A3(+), RAB1B(+)",
  paste("Model B top18 features:", paste(read_result("modelB/modelB_feature_table.csv")$genesymbol, collapse = ", ")),
  paste("Model C top9 features:", paste(read_result("modelC/modelC_feature_table.csv")$genesymbol, collapse = ", "))
)
writeLines(feature_lines, file.path(table_dir, "modelABC_feature_sets_readable.txt"))

wrap_text <- function(x, width = 118) {
  unlist(strwrap(x, width = width), use.names = FALSE)
}

draw_page <- function(lines, title = NULL, cex = 0.64) {
  grid::grid.newpage()
  if (!is.null(title)) {
    grid::grid.text(title, x = grid::unit(0.04, "npc"), y = grid::unit(0.965, "npc"),
                    just = c("left", "top"), gp = grid::gpar(fontface = "bold", fontsize = 13))
    y0 <- 0.925
  } else {
    y0 <- 0.965
  }
  line_h <- 0.024
  max_lines <- floor((y0 - 0.035) / line_h)
  for (i in seq_len(min(length(lines), max_lines))) {
    grid::grid.text(lines[i], x = grid::unit(0.04, "npc"),
                    y = grid::unit(y0 - (i - 1) * line_h, "npc"),
                    just = c("left", "top"),
                    gp = grid::gpar(fontfamily = "mono", fontsize = 10.5 * cex))
  }
}

paginate_lines <- function(lines, n = 34) {
  split(lines, ceiling(seq_along(lines) / n))
}

final_perf_lines <- function(df) {
  out <- c()
  for (i in seq_len(nrow(df))) {
    out <- c(
      out,
      paste0(df$model[i], " | ", df$dataset[i]),
      paste0("  AUC=", df$AUC[i], " (95% CI ", df$CI95[i], "); threshold=", df$threshold[i]),
      paste0("  Confusion: TP=", df$TP[i], ", FP=", df$FP[i], ", FN=", df$FN[i], ", TN=", df$TN[i]),
      paste0("  Accuracy=", df$accuracy[i], "; Sensitivity=", df$sensitivity[i],
             "; Specificity=", df$specificity[i], "; Balanced accuracy=", df$balanced_accuracy[i]),
      ""
    )
  }
  out
}

candidate_lines <- function(df) {
  out <- c()
  for (model_name in unique(df$model)) {
    out <- c(out, paste0(model_name, ":"), "")
    sub <- df[df$model == model_name, , drop = FALSE]
    for (i in seq_len(nrow(sub))) {
      th <- ifelse(is.na(sub$threshold_display[i]), "NA", sub$threshold_display[i])
      out <- c(out, paste0(
        "  ", sub$dataset[i], " | ", sub$method[i],
        " | AUC=", sub$AUC[i], " (95% CI ", sub$CI95[i], ")",
        " | threshold=", th
      ))
    }
    out <- c(out, "")
  }
  out
}

score_summary_lines <- function(df) {
  out <- c()
  for (i in seq_len(nrow(df))) {
    out <- c(out, paste0(
      df$model[i], " | ", df$dataset[i],
      " | n=", df$n[i],
      " | min=", df$min[i],
      " | Q1=", df$q1[i],
      " | median=", df$median[i],
      " | mean=", df$mean[i],
      " | Q3=", df$q3[i],
      " | max=", df$max[i]
    ))
  }
  out
}

write_table_pages <- function(pdf_path) {
  grDevices::pdf(pdf_path, width = 8.27, height = 11.69, onefile = TRUE, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)

  intro <- c(
    "Model A/B/C readable score report",
    "",
    "All numeric values are read from file0624/result_file model tables and plotted_data outputs.",
    "This PDF is intended for human reading; CSV copies are saved under result_file/tables/.",
    "",
    feature_lines,
    "",
    "Key source files:",
    "  result_file/modelA/modelA_final_performance.csv",
    "  result_file/modelB/modelB_final_performance.csv",
    "  result_file/modelC/modelC_final_performance.csv",
    "  result_file/modelA/modelA_final_prediction_scores.csv",
    "  result_file/modelB/modelB_discovery_oob_predictions.csv",
    "  result_file/modelB/modelB_validation_predictions.csv",
    "  result_file/modelC/modelC_final_prediction_scores.csv"
  )
  draw_page(wrap_text(intro, 100), title = "Readable Model Scores", cex = 0.78)

  final_lines <- final_perf_lines(final_perf)
  chunks <- paginate_lines(final_lines, 32)
  for (chunk in chunks) draw_page(chunk, title = "Final Model A/B/C Performance", cex = 0.76)

  cand_view <- candidate_auc[, c("model", "method", "dataset", "AUC", "CI95", "threshold_display")]
  cand_lines <- candidate_lines(cand_view)
  chunks <- paginate_lines(wrap_text(cand_lines, 112), 36)
  for (chunk in chunks) draw_page(chunk, title = "Candidate Model Screen", cex = 0.70)

  score_lines <- score_summary_lines(score_summary)
  chunks <- paginate_lines(wrap_text(score_lines, 112), 36)
  for (chunk in chunks) draw_page(chunk, title = "Score Distribution Summary", cex = 0.70)
}

pdf_path <- file.path(fig_dir, "ModelABC_scores_readable_A4.pdf")
write_table_pages(pdf_path)

cat("Readable Model A/B/C score report written to: ", pdf_path, "\n", sep = "")
cat("Readable tables written to result_file/tables/modelABC_* files.\n")
