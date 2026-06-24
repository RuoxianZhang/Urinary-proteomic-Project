options(stringsAsFactors = FALSE, warn = 1)

set_reproducible_rng <- function() {
  suppressWarnings({
    tryCatch(
      RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection"),
      error = function(e) RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion")
    )
  })
}

set_reproducible_rng()

required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork", "scales",
  "limma", "pROC", "randomForest", "glmnet", "e1071", "MASS",
  "clusterProfiler", "ReactomePA", "org.Hs.eg.db", "openxlsx",
  "jsonlite", "stringr"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(limma)
  library(pROC)
  library(randomForest)
  library(glmnet)
  library(e1071)
  library(MASS)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(openxlsx)
  library(jsonlite)
  library(stringr)
})

project_dir <- normalizePath(Sys.getenv("UR_FILE0624_PROJECT_DIR", unset = Sys.getenv("UR_PROJECT_DIR", unset = Sys.getenv("UR_FILE0623_PROJECT_DIR", unset = getwd()))), mustWork = TRUE)
result_dir <- file.path(project_dir, "result_file")

dir_map <- list(
  processed_data = file.path(result_dir, "processed_data"),
  modelA = file.path(result_dir, "modelA"),
  modelB = file.path(result_dir, "modelB"),
  modelC = file.path(result_dir, "modelC"),
  differential = file.path(result_dir, "differential_analysis"),
  enrichment = file.path(result_dir, "enrichment"),
  tables = file.path(result_dir, "tables"),
  figures = file.path(result_dir, "figures"),
  plotted = file.path(result_dir, "figures", "plotted_data"),
  metadata = file.path(result_dir, "metadata"),
  logs = file.path(result_dir, "logs")
)

ensure_output_dirs <- function() {
  invisible(lapply(dir_map, dir.create, recursive = TRUE, showWarnings = FALSE))
}

write_csv_safe <- function(x, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = row.names)
  invisible(path)
}

write_json_safe <- function(x, path, pretty = TRUE, auto_unbox = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(x, path, pretty = pretty, auto_unbox = auto_unbox, null = "null")
  invisible(path)
}

clean_history <- function(x, yes_label) {
  raw <- as.character(x)
  raw_trim <- trimws(raw)
  missing_like <- is.na(raw_trim) | raw_trim == "" | toupper(raw_trim) %in% c("NA", "N/A", "UNKNOWN", "NOT RECORDED", "NONE")
  yes_like <- raw_trim %in% yes_label | tolower(raw_trim) %in% c("yes", "y", "true", "1", tolower(yes_label))
  out <- ifelse(!missing_like & yes_like, yes_label, "No_or_not_recorded")
  factor(out, levels = c("No_or_not_recorded", yes_label))
}

read_project_data <- function() {
  clinical <- read.csv(file.path(project_dir, "rawdata", "clinical.csv"), check.names = FALSE) %>%
    mutate(
      sample_id = as.character(sample_id),
      group = factor(group, levels = c("Healthy", "LungCancer")),
      major_type = factor(major_type, levels = c("Healthy", "NSCLC", "SCLC")),
      n_status = factor(n_status, levels = c("N0", "Nplus")),
      histology = factor(histology),
      sex = factor(sex),
      smoking = factor(smoking),
      drinking = factor(drinking),
      smoking_status = clean_history(smoking, "Smoker"),
      drinking_status = clean_history(drinking, "Drinker"),
      age_group = ifelse(age <= 60, "Age <=60", "Age >60")
    )
  expr_df <- read.csv(file.path(project_dir, "rawdata", "expression.csv"), check.names = FALSE)
  stopifnot("genesymbol" %in% names(expr_df))
  expr <- as.matrix(expr_df[, setdiff(names(expr_df), "genesymbol"), drop = FALSE])
  rownames(expr) <- make.unique(as.character(expr_df$genesymbol))
  mode(expr) <- "numeric"
  missing_samples <- setdiff(clinical$sample_id, colnames(expr))
  if (length(missing_samples) > 0) stop("Clinical samples missing from expression matrix: ", paste(missing_samples, collapse = ", "))
  expr <- expr[, clinical$sample_id, drop = FALSE]
  list(clinical = clinical, expr = expr, expr_df = data.frame(genesymbol = rownames(expr), expr, check.names = FALSE))
}

is_usable_covariate <- function(x) {
  if (is.numeric(x)) return(length(unique(stats::na.omit(x))) > 1)
  length(unique(stats::na.omit(as.character(x)))) > 1
}

clean_design_covariates <- function(df, covariates) {
  covariates[vapply(covariates, function(v) {
    v %in% names(df) && is_usable_covariate(df[[v]])
  }, logical(1))]
}

stratified_train_indices <- function(y, counts, seed) {
  set.seed(seed)
  idx <- integer()
  for (cl in names(counts)) {
    candidates <- which(y == cl)
    if (length(candidates) < counts[[cl]]) stop("Not enough samples for class ", cl)
    idx <- c(idx, sample(candidates, counts[[cl]]))
  }
  sort(idx)
}

stratified_train_legacy <- function(y, counts) {
  idx <- seq_along(y)
  train <- integer()
  for (cl in names(counts)) {
    candidates <- idx[y == cl]
    if (length(candidates) < counts[[cl]]) stop("Not enough samples for class ", cl)
    train <- c(train, sample(candidates, counts[[cl]]))
  }
  sort(train)
}

normalize_weights <- function(w) {
  w[is.na(w) | !is.finite(w)] <- 0
  denom <- sqrt(sum(w^2))
  if (denom == 0) denom <- 1
  w / denom
}

linear_score <- function(x, w) {
  as.numeric(as.matrix(x[, names(w), drop = FALSE]) %*% normalize_weights(w))
}

direction_weights <- function(x, y, markers, positive, negative) {
  w <- sapply(markers, function(g) {
    delta <- mean(x[y == positive, g], na.rm = TRUE) - mean(x[y == negative, g], na.rm = TRUE)
    ifelse(delta >= 0, 1, -1)
  })
  names(w) <- markers
  w
}

scale_train_apply <- function(x_train, x_new = NULL) {
  center <- colMeans(x_train, na.rm = TRUE)
  scale <- apply(x_train, 2, stats::sd, na.rm = TRUE)
  scale[is.na(scale) | scale == 0] <- 1
  train_scaled <- sweep(sweep(x_train, 2, center, "-"), 2, scale, "/")
  out <- list(train = train_scaled, center = center, scale = scale)
  if (!is.null(x_new)) {
    out$new <- sweep(sweep(x_new, 2, center, "-"), 2, scale, "/")
  }
  out
}

rank_normalize_train_apply <- function(x_train, x_new = NULL) {
  normalize_train_col <- function(v) {
    r <- rank(v, ties.method = "average", na.last = "keep")
    n <- sum(!is.na(v))
    qnorm((r - 0.5) / n)
  }
  train_norm <- apply(x_train, 2, normalize_train_col)
  if (is.vector(train_norm)) train_norm <- matrix(train_norm, ncol = ncol(x_train), dimnames = dimnames(x_train))
  out <- list(train = train_norm)
  if (!is.null(x_new)) {
    new_norm <- matrix(NA_real_, nrow = nrow(x_new), ncol = ncol(x_new), dimnames = dimnames(x_new))
    for (j in seq_len(ncol(x_train))) {
      v <- x_train[, j]
      n <- sum(!is.na(v))
      for (i in seq_len(nrow(x_new))) {
        p <- (sum(v < x_new[i, j], na.rm = TRUE) + 0.5 * sum(v == x_new[i, j], na.rm = TRUE) + 0.5) / (n + 1)
        p <- min(max(p, 1 / (n + 1)), n / (n + 1))
        new_norm[i, j] <- qnorm(p)
      }
    }
    out$new <- new_norm
  }
  out
}

rank_normalize_within <- function(x) {
  out <- apply(x, 2, function(v) {
    r <- rank(v, ties.method = "average", na.last = "keep")
    n <- sum(!is.na(v))
    qnorm((r - 0.5) / n)
  })
  if (is.vector(out)) out <- matrix(out, ncol = ncol(x), dimnames = dimnames(x))
  out
}

row_z <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

cap_value <- function(x, limit = 2.5) pmax(pmin(x, limit), -limit)

auc_ci <- function(truth, score, negative, positive) {
  truth <- factor(truth, levels = c(negative, positive))
  roc_obj <- pROC::roc(truth, score, levels = c(negative, positive), direction = "<", quiet = TRUE)
  ci <- suppressWarnings(pROC::ci.auc(roc_obj))
  data.frame(
    auc = as.numeric(pROC::auc(roc_obj)),
    ci_low = as.numeric(ci[1]),
    ci_high = as.numeric(ci[3]),
    check.names = FALSE
  )
}

roc_points <- function(truth, score, negative, positive) {
  truth <- factor(truth, levels = c(negative, positive))
  roc_obj <- pROC::roc(truth, score, levels = c(negative, positive), direction = "<", quiet = TRUE)
  data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    threshold = roc_obj$thresholds,
    auc = as.numeric(pROC::auc(roc_obj)),
    check.names = FALSE
  ) %>% arrange(fpr, tpr)
}

youden_threshold <- function(truth, score, negative, positive) {
  roc_obj <- pROC::roc(factor(truth, levels = c(negative, positive)), score, levels = c(negative, positive), direction = "<", quiet = TRUE)
  as.numeric(pROC::coords(roc_obj, x = "best", best.method = "youden", ret = "threshold", transpose = FALSE)[1, 1])
}

metric_at_threshold <- function(truth, score, threshold, negative, positive) {
  pred <- ifelse(score >= threshold, positive, negative)
  tp <- sum(pred == positive & truth == positive, na.rm = TRUE)
  fp <- sum(pred == positive & truth == negative, na.rm = TRUE)
  fn <- sum(pred == negative & truth == positive, na.rm = TRUE)
  tn <- sum(pred == negative & truth == negative, na.rm = TRUE)
  sensitivity <- ifelse(tp + fn > 0, tp / (tp + fn), NA_real_)
  specificity <- ifelse(tn + fp > 0, tn / (tn + fp), NA_real_)
  data.frame(
    threshold = threshold,
    tp = tp, fp = fp, fn = fn, tn = tn,
    accuracy = (tp + tn) / (tp + fp + fn + tn),
    sensitivity = sensitivity,
    specificity = specificity,
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    check.names = FALSE
  )
}

confusion_table <- function(truth, score, threshold, negative, positive) {
  pred <- ifelse(score >= threshold, positive, negative)
  expand.grid(
    actual = c(positive, negative),
    predicted = c(positive, negative),
    stringsAsFactors = FALSE
  ) %>%
    rowwise() %>%
    mutate(n = sum(truth == actual & pred == predicted, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      outcome = case_when(
        actual == positive & predicted == positive ~ "TP",
        actual == positive & predicted == negative ~ "FN",
        actual == negative & predicted == positive ~ "FP",
        TRUE ~ "TN"
      )
    )
}

threshold_grid_metrics <- function(truth, score, negative, positive, n = 201) {
  grid <- seq(min(score, na.rm = TRUE), max(score, na.rm = TRUE), length.out = n)
  bind_rows(lapply(grid, function(th) {
    cbind(data.frame(threshold = th), metric_at_threshold(truth, score, th, negative, positive)[, -1, drop = FALSE])
  })) %>% arrange(desc(balanced_accuracy), desc(accuracy))
}

limma_de <- function(expr, clinical, samples, group_var, positive, negative, covariates = character()) {
  df <- clinical %>% filter(sample_id %in% samples)
  df <- df[match(samples, df$sample_id), , drop = FALSE]
  rownames(df) <- df$sample_id
  df[[group_var]] <- factor(as.character(df[[group_var]]), levels = c(negative, positive))
  df <- df[!is.na(df[[group_var]]), , drop = FALSE]
  covariates <- clean_design_covariates(df, covariates)
  formula_text <- paste("~", paste(c(group_var, covariates), collapse = " + "))
  design_formula <- as.formula(formula_text)
  model_frame <- model.frame(design_formula, data = df, na.action = stats::na.omit)
  design <- model.matrix(design_formula, data = model_frame)
  used_samples <- rownames(model_frame)
  fit <- limma::lmFit(expr[, used_samples, drop = FALSE], design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  coef_name <- grep(paste0("^", group_var), colnames(design), value = TRUE)[1]
  tab <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  tab$gene <- rownames(tab)
  group_values <- as.character(df[used_samples, group_var])
  negative_samples <- used_samples[group_values == negative]
  positive_samples <- used_samples[group_values == positive]
  tab %>%
    dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val, B) %>%
    mutate(
      mean_negative = rowMeans(expr[gene, negative_samples, drop = FALSE], na.rm = TRUE),
      mean_positive = rowMeans(expr[gene, positive_samples, drop = FALSE], na.rm = TRUE),
      direction_nominal = case_when(
        P.Value < 0.05 & logFC > 0.3 ~ "positive_up",
        P.Value < 0.05 & logFC < -0.3 ~ "negative_up",
        TRUE ~ "not_selected"
      ),
      direction_fdr = case_when(
        adj.P.Val < 0.05 & logFC > 0.3 ~ "positive_up",
        adj.P.Val < 0.05 & logFC < -0.3 ~ "negative_up",
        TRUE ~ "not_fdr_significant"
      )
    )
}

module_sets <- list(
  "Acute-phase response" = c("SAA1", "SAA2", "SAA4", "CRP", "SERPINA1", "SERPINA3", "ORM1", "ORM2", "HP", "LBP", "ITIH4", "A2M", "APCS"),
  "Complement / coagulation" = c("C1QA", "C1QB", "C1QC", "C2", "C3", "C4A", "C4B", "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9", "CFB", "CFH", "CFI", "FGA", "FGB", "FGG", "PLG", "SERPINC1"),
  "Myeloid / neutrophil activation" = c("LCN2", "LRG1", "MPO", "ELANE", "CTSC", "CTSZ", "LYZ", "OLFM4", "S100A8", "S100A9", "CD14", "CAMP"),
  "Iron / heme handling" = c("FTH1", "FTL", "TF", "HP", "HBA1", "HBA2", "HBB", "HBD", "CP", "HMOX1"),
  "Adhesion / matrix remodeling" = c("FN1", "VCAN", "VIM", "LGALS3", "ANXA1", "ANXA5", "MMP9", "TIMP1", "ICAM1", "ITGA6", "ITGB1", "SPP1"),
  "Tubular / epithelial transport" = c("ACE2", "SLC13A2", "SLC4A4", "ATP1A1", "ATP1A3", "CUBN", "LRP2", "AQP1", "SLC2A1", "VAMP7", "STOM")
)

module_metadata <- data.frame(
  module = names(module_sets),
  module_short = c("Acute phase", "Complement/coag.", "Myeloid/neutrophil", "Iron/heme", "Adhesion/matrix", "Tubular transport"),
  pathway_source = c("GO biological process ORA", rep("Reactome GSEA", 5)),
  pathway_id = c("GO:0006953", "R-HSA-166658", "R-HSA-6798695", "R-HSA-2168880", "R-HSA-1474244", "R-HSA-425407"),
  pathway_name = c(
    "acute-phase response",
    "Complement cascade",
    "Neutrophil degranulation",
    "Scavenging of heme from plasma",
    "Extracellular matrix organization",
    "SLC-mediated transmembrane transport"
  ),
  check.names = FALSE
)

module_scores <- function(expr, samples, sets = module_sets) {
  z <- row_z(expr[, samples, drop = FALSE])
  bind_rows(lapply(names(sets), function(module) {
    genes <- intersect(sets[[module]], rownames(z))
    if (length(genes) == 0) return(data.frame())
    data.frame(
      sample_id = samples,
      module = module,
      module_score = colMeans(z[genes, , drop = FALSE], na.rm = TRUE),
      genes_used = length(genes),
      genes = paste(genes, collapse = ";"),
      check.names = FALSE
    )
  }))
}

module_stats <- function(score_df, clinical, group_var = "n_status", positive = "Nplus", negative = "N0") {
  df <- score_df
  if (!group_var %in% names(df)) {
    df <- df %>% left_join(clinical, by = "sample_id")
  }
  df %>%
    filter(.data[[group_var]] %in% c(negative, positive)) %>%
    group_by(module, genes_used, genes) %>%
    summarise(
      median_negative = median(module_score[.data[[group_var]] == negative], na.rm = TRUE),
      median_positive = median(module_score[.data[[group_var]] == positive], na.rm = TRUE),
      delta_median = median_positive - median_negative,
      wilcox_p = suppressWarnings(wilcox.test(module_score ~ .data[[group_var]])$p.value),
      .groups = "drop"
    ) %>%
    mutate(FDR = p.adjust(wilcox_p, method = "BH")) %>%
    arrange(FDR)
}

run_enrichment <- function(res_table, universe_symbols, prefix, out_dir, p_cutoff = 0.05, lfc_cutoff = 0.3, min_gs = 10, max_gs = 500) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  selected <- res_table %>%
    filter(P.Value < p_cutoff, abs(logFC) > lfc_cutoff) %>%
    pull(gene) %>%
    unique()
  selected_genes <- data.frame(gene = selected, check.names = FALSE)
  universe_map <- suppressWarnings(clusterProfiler::bitr(universe_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  selected_map <- suppressWarnings(clusterProfiler::bitr(selected, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  go_all <- data.frame()
  reactome_all <- data.frame()
  if (nrow(selected_map) > 0) {
    ego <- tryCatch(
      clusterProfiler::enrichGO(
        gene = unique(selected_map$ENTREZID),
        universe = unique(universe_map$ENTREZID),
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = min_gs,
        maxGSSize = max_gs,
        readable = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(ego)) go_all <- as.data.frame(ego)
    er <- tryCatch(
      ReactomePA::enrichPathway(
        gene = unique(selected_map$ENTREZID),
        universe = unique(universe_map$ENTREZID),
        organism = "human",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        qvalueCutoff = 1,
        minGSSize = min_gs,
        maxGSSize = max_gs,
        readable = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(er)) reactome_all <- as.data.frame(er)
  }
  ranked <- res_table$logFC
  names(ranked) <- res_table$gene
  ranked <- ranked[!is.na(ranked)]
  ranked <- sort(ranked, decreasing = TRUE)
  ranked_map <- suppressWarnings(clusterProfiler::bitr(names(ranked), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  ranked_df <- data.frame(SYMBOL = names(ranked), stat = ranked, check.names = FALSE) %>%
    inner_join(ranked_map, by = "SYMBOL") %>%
    group_by(ENTREZID) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()
  gsea <- data.frame()
  if (nrow(ranked_df) > 0) {
    gene_list <- ranked_df$stat
    names(gene_list) <- ranked_df$ENTREZID
    gene_list <- sort(gene_list, decreasing = TRUE)
    gp <- tryCatch(
      ReactomePA::gsePathway(
        geneList = gene_list,
        organism = "human",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        minGSSize = min_gs,
        maxGSSize = max_gs,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(gp)) gsea <- as.data.frame(gp)
  }
  write_csv_safe(go_all, file.path(out_dir, paste0(prefix, "_GO_BP_ORA.csv")))
  write_csv_safe(reactome_all, file.path(out_dir, paste0(prefix, "_Reactome_ORA.csv")))
  write_csv_safe(gsea, file.path(out_dir, paste0(prefix, "_Reactome_GSEA.csv")))
  write_csv_safe(selected_genes, file.path(out_dir, paste0(prefix, "_selected_nominal_genes.csv")))
  membership <- bind_rows(
    if (nrow(go_all) > 0) go_all %>% transmute(source = "GO_BP_ORA", ID, Description, p.adjust, geneID) else data.frame(),
    if (nrow(reactome_all) > 0) reactome_all %>% transmute(source = "Reactome_ORA", ID, Description, p.adjust, geneID) else data.frame()
  )
  if (nrow(membership) > 0) {
    long <- membership %>%
      tidyr::separate_rows(geneID, sep = "/") %>%
      dplyr::rename(gene = geneID)
  } else {
    long <- data.frame()
  }
  write_csv_safe(long, file.path(out_dir, paste0(prefix, "_pathway_protein_membership_long.csv")))
  write_csv_safe(membership, file.path(out_dir, paste0(prefix, "_pathway_protein_membership_wide.csv")))
  list(go = go_all, reactome = reactome_all, gsea = gsea, membership = membership, membership_long = long, selected_genes = selected_genes)
}

parse_zhang_data <- function(markers = NULL) {
  meta_raw <- openxlsx::read.xlsx(file.path(project_dir, "zhangdata", "mmc7.xlsx"), sheet = "Sheet1", colNames = FALSE)
  header <- as.character(unlist(meta_raw[4, ]))
  header <- make.names(header, unique = TRUE)
  meta <- meta_raw[-(1:4), , drop = FALSE]
  names(meta) <- header
  meta <- meta %>%
    transmute(
      sample_id = as.character(Sample.ID),
      label = as.character(Label),
      dataset = as.character(Type.of.data.set),
      sex = as.character(Gender),
      age = suppressWarnings(as.numeric(Age)),
      group = ifelse(label == "LC", "LungCancer", ifelse(label == "CTL", "Healthy", "Other"))
    )
  all_raw <- openxlsx::read.xlsx(file.path(project_dir, "zhangdata", "mmc8.xlsx"), sheet = "iFOT of all proteins", colNames = FALSE)
  matrix_labels <- as.character(unlist(all_raw[1, -(1:2)]))
  mat <- all_raw[-1, , drop = FALSE]
  gene <- make.unique(as.character(mat[[1]]))
  x <- as.matrix(mat[, -(1:2), drop = FALSE])
  mode(x) <- "numeric"
  rownames(x) <- gene
  sample_ids <- character(length(matrix_labels))
  label_seen <- setNames(integer(length(unique(matrix_labels))), unique(matrix_labels))
  for (i in seq_along(matrix_labels)) {
    lab <- matrix_labels[i]
    label_seen[lab] <- label_seen[lab] + 1L
    candidates <- meta$sample_id[meta$label == lab]
    if (length(candidates) < label_seen[lab]) {
      stop("Zhang matrix has more columns for label ", lab, " than metadata rows.")
    }
    sample_ids[i] <- candidates[label_seen[lab]]
  }
  colnames(x) <- sample_ids
  if (!is.null(markers)) {
    missing <- setdiff(markers, rownames(x))
    if (length(missing) > 0) stop("Zhang data missing markers: ", paste(missing, collapse = ", "))
    x <- x[markers, , drop = FALSE]
  }
  keep <- meta$sample_id %in% colnames(x) & meta$group %in% c("Healthy", "LungCancer")
  meta <- meta[keep, , drop = FALSE]
  x <- x[, meta$sample_id, drop = FALSE]
  expr_log <- log2(x + 1)
  list(clinical = meta, expr = expr_log)
}

fit_candidate_methods <- function(x_train, y_train, x_list, negative, positive, seed = 1, equal_direction_weights = NULL) {
  set.seed(seed)
  y_factor <- factor(y_train, levels = c(negative, positive))
  out <- list()
  if (is.null(equal_direction_weights)) {
    equal_direction_weights <- sign(colMeans(x_train[y_train == positive, , drop = FALSE]) - colMeans(x_train[y_train == negative, , drop = FALSE]))
    equal_direction_weights[equal_direction_weights == 0] <- 1
  }
  score_equal <- function(x) as.numeric((x %*% equal_direction_weights) / length(equal_direction_weights))
  out[["Equal-direction score"]] <- lapply(x_list, score_equal)
  fit_glm <- tryCatch(glm(y_factor ~ ., data = data.frame(y_factor = y_factor, x_train, check.names = FALSE), family = binomial()), error = function(e) NULL)
  if (!is.null(fit_glm)) {
    out[["Logistic regression"]] <- lapply(x_list, function(x) as.numeric(predict(fit_glm, newdata = data.frame(x, check.names = FALSE), type = "response")))
  }
  y_bin <- ifelse(y_train == positive, 1, 0)
  fit_lasso <- tryCatch(glmnet::cv.glmnet(x_train, y_bin, family = "binomial", alpha = 1, nfolds = min(5, table(y_train))), error = function(e) NULL)
  if (!is.null(fit_lasso)) {
    out[["LASSO logistic"]] <- lapply(x_list, function(x) as.numeric(predict(fit_lasso, newx = x, s = "lambda.min", type = "response")))
  }
  fit_enet <- tryCatch(glmnet::cv.glmnet(x_train, y_bin, family = "binomial", alpha = 0.5, nfolds = min(5, table(y_train))), error = function(e) NULL)
  if (!is.null(fit_enet)) {
    out[["Elastic net"]] <- lapply(x_list, function(x) as.numeric(predict(fit_enet, newx = x, s = "lambda.min", type = "response")))
  }
  fit_rf <- tryCatch(randomForest::randomForest(x = x_train, y = y_factor, ntree = 500, importance = TRUE), error = function(e) NULL)
  if (!is.null(fit_rf)) {
    out[["Random forest"]] <- lapply(x_list, function(x) as.numeric(predict(fit_rf, x, type = "prob")[, positive]))
  }
  fit_svm <- tryCatch(e1071::svm(x = x_train, y = y_factor, kernel = "radial", probability = TRUE, scale = FALSE), error = function(e) NULL)
  if (!is.null(fit_svm)) {
    out[["SVM radial"]] <- lapply(x_list, function(x) {
      pred <- predict(fit_svm, x, probability = TRUE)
      probs <- attr(pred, "probabilities")
      as.numeric(probs[, positive])
    })
  }
  fit_lda <- tryCatch(MASS::lda(x = x_train, grouping = y_factor), error = function(e) NULL)
  if (!is.null(fit_lda)) {
    out[["LDA"]] <- lapply(x_list, function(x) as.numeric(predict(fit_lda, x)$posterior[, positive]))
  }
  out
}

method_performance_table <- function(score_list, truth_list, train_dataset_name, negative, positive) {
  bind_rows(lapply(names(score_list), function(method) {
    bind_rows(lapply(names(score_list[[method]]), function(dataset) {
      score <- score_list[[method]][[dataset]]
      truth <- truth_list[[dataset]]
      train_score <- score_list[[method]][[train_dataset_name]]
      train_truth <- truth_list[[train_dataset_name]]
      th <- youden_threshold(train_truth, train_score, negative, positive)
      cbind(
        data.frame(method = method, dataset = dataset, threshold_source = paste0(train_dataset_name, " Youden"), check.names = FALSE),
        auc_ci(truth, score, negative, positive),
        metric_at_threshold(truth, score, th, negative, positive)
      )
    }))
  }))
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "f", digits = digits))
}

theme_ur <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.22),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      plot.subtitle = element_text(color = "grey25", size = base_size - 0.2),
      strip.background = element_rect(fill = "grey95", color = "grey50", linewidth = 0.25),
      strip.text = element_text(face = "bold", size = base_size - 0.3),
      legend.title = element_text(face = "bold"),
      legend.key.size = unit(0.30, "cm"),
      panel.border = element_rect(fill = NA, color = "grey40", linewidth = 0.20)
    )
}

pal <- c(
  Healthy = "#5BA7B4",
  LungCancer = "#C65A4A",
  NSCLC = "#3F8F78",
  SCLC = "#596C9D",
  LUAD = "#3F8F78",
  LUSC = "#9D7A55",
  N0 = "#8BA1B2",
  Nplus = "#B74F50",
  Discovery = "#68727D",
  Validation = "#1F2933",
  Zhang = "#C65A4A"
)

save_pdf <- function(plot, filename, width = 8.27, height = 11.69) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename, plot, width = width, height = height, units = "in", device = "pdf", bg = "white", useDingbats = FALSE)
  invisible(filename)
}
