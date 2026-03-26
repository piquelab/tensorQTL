#!/usr/bin/env Rscript
# =============================================================================
# 04_tQTL_fdr_analysis.R
# Post-processing of tensorQTL cis (permutation) results
#
# Input:  tensorQTL cis permutation output — one file per condition:
#           <input_dir>/<CONDITION>_cis.cis_qtl.txt.gz
#
# What this script does:
#   1. Loads per-condition cis permutation results
#   2. Selects the best available p-value (pval_beta > pval_perm > pval_nominal)
#   3. Applies cross-gene FDR correction via p.adjust()
#   4. Calls eGenes at a user-defined FDR threshold
#   5. Writes:
#        - QQ plot per condition
#        - Bar chart of eGene counts across conditions
#        - Per-condition eGene lists
#        - Plain-text summary report
#
# Run:
#   Rscript 04_tQTL_fdr_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

# =============================================================================
# CONFIGURATION — edit these variables before running
# =============================================================================

input_dir  <- "/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/eQTL_mapping/tensor/tensorqtl_output_cis_SV15_100kb"
output_dir <- file.path(input_dir, "fdr_analysis")

# FDR threshold for eGene calling
fdr_threshold <- 0.10

# =============================================================================
# SETUP
# =============================================================================

dir.create(file.path(output_dir, "figures"),     showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "eGene_lists"), showWarnings = FALSE, recursive = TRUE)

today <- format(Sys.Date(), "%m%d%Y")

# Auto-detect conditions from files matching <CONDITION>_cis.cis_qtl.txt.gz
qtl_files  <- list.files(input_dir, pattern = "_cis\\.cis_qtl\\.txt\\.gz$", full.names = FALSE)
conditions <- sub("_cis\\.cis_qtl\\.txt\\.gz$", "", qtl_files)

if (length(conditions) == 0) {
  stop("No files matching *_cis.cis_qtl.txt.gz found in: ", input_dir)
}

cat("=== tensorQTL FDR Analysis ===\n")
cat("Input dir    :", input_dir, "\n")
cat("Conditions   :", paste(conditions, collapse = ", "), "\n")
cat("FDR threshold:", fdr_threshold, "\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Extract pval_beta; stop with a clear message if the column is absent
select_pval <- function(df, condition) {
  if (!"pval_beta" %in% colnames(df)) {
    stop("pval_beta column not found for condition: ", condition,
         "\n  Available columns: ", paste(colnames(df), collapse = ", "))
  }
  list(pvals = df$pval_beta, method = "pval_beta")
}

# QQ plot for a single condition
make_qq_plot <- function(pvals, condition, method_label, fdr_threshold, out_path) {
  pvals <- sort(pvals[!is.na(pvals) & pvals > 0])
  n     <- length(pvals)
  if (n == 0) {
    message("  No valid p-values for QQ plot: ", condition)
    return(invisible(NULL))
  }

  df <- data.frame(
    expected = -log10((seq_len(n)) / n),
    observed = -log10(pvals),
    ci_lo    = -log10(qbeta(0.975, seq_len(n), n - seq_len(n) + 1)),
    ci_hi    = -log10(qbeta(0.025, seq_len(n), n - seq_len(n) + 1))
  )

  p <- ggplot(df, aes(x = expected, y = observed)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = "grey80", alpha = 0.5) +
    geom_point(size = 0.8, alpha = 0.6, colour = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
    labs(
      title    = condition,
      subtitle = paste0("Base p-value: ", method_label,
                        "  |  FDR threshold: ", fdr_threshold),
      x = expression(Expected ~ -log[10](p)),
      y = expression(Observed ~ -log[10](p))
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = 9)
    )

  ggsave(out_path, p, width = 7, height = 7, dpi = 300)
  invisible(p)
}

# =============================================================================
# MAIN LOOP — load, correct, call eGenes
# =============================================================================

cat("--- Loading results and applying FDR correction ---\n")

results      <- list()  # full corrected data frames
egene_lists  <- list()  # significant phenotype IDs per condition
summary_rows <- list()  # one-row summary per condition

for (cond in conditions) {
  cat("\nCondition:", cond, "\n")

  qtl_file <- file.path(input_dir, paste0(cond, "_cis.cis_qtl.txt.gz"))
  if (!file.exists(qtl_file)) {
    message("  WARNING: File not found — skipping: ", qtl_file)
    next
  }

  res <- fread(qtl_file)
  cat("  Genes loaded:", nrow(res), "\n")

  # Extract pval_beta (required)
  pv            <- select_pval(res, cond)
  res$pval_base <- pv$pvals
  method_label  <- pv$method
  cat("  P-value column:", method_label, "\n")

  # Cross-gene FDR correction
  res$pval_fdr <- p.adjust(res$pval_base, method = "fdr")

  # Call eGenes
  sig      <- !is.na(res$pval_fdr) & res$pval_fdr < fdr_threshold
  egenes   <- unique(res$phenotype_id[sig])
  n_egenes <- length(egenes)
  cat("  eGenes (FDR <", fdr_threshold, "):", n_egenes, "\n")

  # Store
  results[[cond]]     <- res
  egene_lists[[cond]] <- egenes

  summary_rows[[cond]] <- data.frame(
    condition     = cond,
    n_genes       = nrow(res),
    pval_method   = method_label,
    n_egenes      = n_egenes,
    fdr_threshold = fdr_threshold,
    stringsAsFactors = FALSE
  )

  # QQ plot
  qq_path <- file.path(output_dir, "figures",
                       paste0(today, "_QQ_", cond, ".png"))
  make_qq_plot(res$pval_base, cond, method_label, fdr_threshold, qq_path)
  cat("  QQ plot saved:", basename(qq_path), "\n")

  # eGene list
  egene_path <- file.path(output_dir, "eGene_lists",
                          paste0(today, "_eGenes_", cond, ".txt"))
  writeLines(egenes, egene_path)
  cat("  eGene list saved:", basename(egene_path), "\n")
}

# =============================================================================
# SUMMARY TABLE
# =============================================================================

summary_df   <- do.call(rbind, summary_rows)
summary_path <- file.path(output_dir, paste0(today, "_eGene_summary.txt"))
fwrite(summary_df, summary_path, sep = "\t", quote = FALSE)

cat("\n--- Summary table ---\n")
print(summary_df, row.names = FALSE)

# =============================================================================
# BAR CHART — eGene counts per condition
# =============================================================================

p_bar <- ggplot(summary_df, aes(x = condition, y = n_egenes)) +
  geom_col(fill = "steelblue", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = n_egenes), vjust = -0.4, size = 4) +
  labs(
    title = paste0("eGenes per condition  (FDR < ", fdr_threshold, ")"),
    x     = NULL,
    y     = "Number of eGenes"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

bar_path <- file.path(output_dir, "figures",
                      paste0(today, "_eGene_counts.png"))
ggsave(bar_path, p_bar,
       width  = max(5, length(conditions) * 1.5),
       height = 5, dpi = 300)
cat("\nBar chart saved:", basename(bar_path), "\n")

# =============================================================================
# PLAIN-TEXT SUMMARY REPORT
# =============================================================================

unique_egenes <- unique(unlist(egene_lists))
writeLines(unique_egenes, file.path(input_dir, "unique_significant_egenes.txt"))
cat("Unique eGene list saved: unique_significant_egenes.txt\n")

report <- c(
  "tensorQTL FDR Analysis Report",
  strrep("=", 50),
  paste("Date             :", today),
  paste("Input directory  :", input_dir),
  paste("FDR threshold    :", fdr_threshold),
  paste("Conditions       :", length(summary_rows), "processed,",
        length(conditions) - length(summary_rows), "missing"),
  "",
  "--- eGene counts ---",
  apply(summary_df, 1, function(r)
    sprintf("  %-25s  %s eGenes  [%s]",
            r["condition"], r["n_egenes"], r["pval_method"])
  ),
  "",
  paste("Unique eGenes (any condition):", length(unique_egenes)),
  "",
  "--- Output files ---",
  paste("  Summary table :", basename(summary_path)),
  paste("  Bar chart     :", basename(bar_path)),
  paste("  QQ plots      : figures/"),
  paste("  eGene lists   : eGene_lists/")
)

report_path <- file.path(output_dir, paste0(today, "_SUMMARY.txt"))
writeLines(report, report_path)
cat("\nSummary report saved:", basename(report_path), "\n\n")
cat(paste(report, collapse = "\n"), "\n")

cat("\n=== Done ===\n")
