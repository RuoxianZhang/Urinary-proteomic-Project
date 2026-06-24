project_dir <- normalizePath(Sys.getenv("UR_FILE0624_PROJECT_DIR", unset = Sys.getenv("UR_PROJECT_DIR", unset = Sys.getenv("UR_FILE0623_PROJECT_DIR", unset = getwd()))), mustWork = TRUE)
report_dir <- file.path(project_dir, "reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

all_paths <- list.files(project_dir, recursive = TRUE, all.files = TRUE, full.names = TRUE, no.. = TRUE)
info <- file.info(all_paths)
file_paths <- all_paths[!is.na(info$isdir) & !info$isdir]
rel_paths <- sub(paste0("^", gsub("([\\^\\$\\.\\|\\(\\)\\[\\]\\*\\+\\?\\{\\}\\\\])", "\\\\\\1", project_dir), "/?"), "", normalizePath(file_paths, mustWork = TRUE))

exclude <- grepl("(^|/)\\.DS_Store$", rel_paths) |
  grepl("(^|/)__pycache__/", rel_paths) |
  rel_paths %in% c(
    "reports/file0624_publication_manifest.csv",
    "reports/file0624_package_summary.txt",
    "result_file/logs/11_make_publication_manifest.log"
  )

file_paths <- file_paths[!exclude]
rel_paths <- rel_paths[!exclude]
info <- file.info(file_paths)

classify_file <- function(path) {
  if (grepl("^rawdata/", path)) return("primary_raw_data")
  if (grepl("^zhangdata/", path)) return("external_validation_data")
  if (grepl("^code/", path)) return("analysis_code")
  if (grepl("^result_file/modelA/", path)) return("modelA_outputs")
  if (grepl("^result_file/modelB/", path)) return("modelB_outputs")
  if (grepl("^result_file/modelC/", path)) return("modelC_outputs")
  if (grepl("^result_file/figures/", path)) return("code_generated_figures_and_plotted_data")
  if (grepl("^result_file/differential_analysis/", path)) return("differential_analysis_outputs")
  if (grepl("^result_file/enrichment/", path)) return("enrichment_outputs")
  if (grepl("^result_file/processed_data/", path)) return("processed_data")
  if (grepl("^result_file/metadata/", path)) return("metadata")
  if (grepl("^result_file/tables/", path)) return("tables")
  if (grepl("^result_file/logs/", path)) return("pipeline_logs")
  if (grepl("^reports/", path)) return("audit_reports")
  if (grepl("^(Figure|Supplementary_Figure|ModelABC).*\\.pdf$", path)) return("publication_ready_root_pdfs")
  if (path %in% c("README.md", "PUBLICATION_PACKAGE_INDEX.md")) return("package_documentation")
  "other"
}

manifest <- data.frame(
  path = rel_paths,
  category = vapply(rel_paths, classify_file, character(1)),
  size_bytes = as.numeric(info$size),
  modified_time = format(as.POSIXct(info$mtime), "%Y-%m-%d %H:%M:%S %z"),
  md5 = unname(tools::md5sum(file_paths)),
  stringsAsFactors = FALSE
)
manifest <- manifest[order(manifest$category, manifest$path), ]

manifest_path <- file.path(report_dir, "file0624_publication_manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)

summary_df <- aggregate(size_bytes ~ category, manifest, function(x) c(files = length(x), bytes = sum(x)))
summary_lines <- c(
  "file0624 publication reproducibility package summary",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("Project directory: ", project_dir),
  paste0("Total tracked files: ", nrow(manifest)),
  paste0("Total tracked size bytes: ", sum(manifest$size_bytes)),
  "",
  "Category summary:"
)
for (i in seq_len(nrow(summary_df))) {
  vals <- summary_df$size_bytes[i, ]
  summary_lines <- c(
    summary_lines,
    paste0("- ", summary_df$category[i], ": ", vals[["files"]], " files; ", vals[["bytes"]], " bytes")
  )
}
writeLines(summary_lines, file.path(report_dir, "file0624_package_summary.txt"))

cat("Publication manifest written: ", manifest_path, "\n", sep = "")
