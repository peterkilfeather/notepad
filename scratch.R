# Sex mismatch QC for bulk RNA-seq (R)
# - Input: a raw count matrix (genes x samples) and sample metadata with reported sex
# - Output: inferred sex calls, mismatch flags, and QC plots

suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mclust)
  library(pheatmap)
})

# -----------------------------
# 0) USER INPUTS
# -----------------------------
# counts: matrix/data.frame with genes as rownames and samples as colnames
# meta:   data.frame with one row per sample, and a sample id column that matches colnames(counts)
#
# Example:
# counts <- readRDS("counts.rds")  # genes x samples
# meta   <- read.csv("metadata.csv")

sample_id_col <- "sample_id"          # <-- change if needed
reported_sex_col <- "reported_sex"     # <-- change if needed (values like M/F/Male/Female)

# -----------------------------
# 1) Harmonize reported sex
# -----------------------------
harmonize_sex <- function(x) {
  x <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x %in% c("f", "female", "woman", "w") ~ "F",
    x %in% c("m", "male", "man") ~ "M",
    x %in% c("", "na", "n/a", "unknown", "unk", "u", "other") ~ "Unknown",
    is.na(x) ~ "Unknown",
    TRUE ~ "Unknown"
  )
}

# -----------------------------
# 2) Compute log2(CPM+1)
# -----------------------------
compute_logcpm <- function(counts_mat) {
  stopifnot(!is.null(rownames(counts_mat)))
  y <- DGEList(counts = counts_mat)
  # minimal filtering to stabilize CPM; keep markers regardless later
  keep <- filterByExpr(y)
  y2 <- y[keep, , keep.lib.sizes = FALSE]
  y2 <- calcNormFactors(y2)
  logcpm <- cpm(y2, log = TRUE, prior.count = 1)

  # Re-add any genes that were filtered out but might be needed (e.g., XIST / Y markers)
  # We'll handle missing genes downstream gracefully; no need to force re-add here.
  logcpm
}

# -----------------------------
# 3) Marker panel
# -----------------------------
markers_y <- c("DDX3Y","KDM5D","UTY","RPS4Y1","ZFY","EIF1AY","USP9Y")
marker_xist <- "XIST"

# If your counts are Ensembl IDs, map them to symbols first (e.g., via biomaRt/org.Hs.eg.db).
# For now, this assumes gene symbols as rownames(counts).

# -----------------------------
# 4) Score samples
# -----------------------------
score_sex_markers <- function(logcpm, markers_y, marker_xist) {
  genes <- rownames(logcpm)

  y_present <- intersect(markers_y, genes)
  x_present <- intersect(marker_xist, genes)

  # Y-score as mean across available Y markers
  y_score <- if (length(y_present) > 0) colMeans(logcpm[y_present, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(logcpm))
  names(y_score) <- colnames(logcpm)

  # XIST score (single gene)
  xist_score <- if (length(x_present) == 1) as.numeric(logcpm[x_present, ]) else rep(NA_real_, ncol(logcpm))
  names(xist_score) <- colnames(logcpm)

  # Count Y markers above a small threshold
  y_thresh <- 1  # on log2(CPM+1) scale; adjust if needed
  y_n_above <- if (length(y_present) > 0) colSums(logcpm[y_present, , drop = FALSE] > y_thresh, na.rm = TRUE) else rep(NA_integer_, ncol(logcpm))
  names(y_n_above) <- colnames(logcpm)

  tibble(
    sample = colnames(logcpm),
    y_score = y_score,
    xist_score = xist_score,
    y_n_above = y_n_above,
    y_markers_used = length(y_present),
    has_xist = length(x_present) == 1
  )
}

# -----------------------------
# 5) Infer sex via 2D clustering (Mclust or kmeans fallback)
# -----------------------------
infer_sex <- function(scores_df) {
  df <- scores_df %>%
    mutate(
      ok = is.finite(y_score) & is.finite(xist_score)
    )

  df_ok <- df %>% filter(ok)

  if (nrow(df_ok) < 10) {
    warning("Too few samples with finite scores to infer sex robustly.")
    return(df %>% mutate(inferred_sex = NA_character_, inferred_conf = NA_real_))
  }

  X <- as.matrix(df_ok[, c("xist_score","y_score")])

  # Prefer model-based clustering; if it fails, use kmeans
  fit <- tryCatch(Mclust(X, G = 2, verbose = FALSE), error = function(e) NULL)

  if (!is.null(fit)) {
    cl <- fit$classification
    # confidence: max posterior
    conf <- apply(fit$z, 1, max)

    df_ok <- df_ok %>%
      mutate(cluster = cl, inferred_conf = conf)
  } else {
    set.seed(1)
    km <- kmeans(X, centers = 2)
    df_ok <- df_ok %>%
      mutate(cluster = km$cluster, inferred_conf = NA_real_)
  }

  # Label clusters: cluster with higher median y_score is "M"; higher xist_score is "F"
  cluster_summary <- df_ok %>%
    group_by(cluster) %>%
    summarise(
      med_y = median(y_score, na.rm = TRUE),
      med_xist = median(xist_score, na.rm = TRUE),
      .groups = "drop"
    )

  male_cluster <- cluster_summary %>% arrange(desc(med_y)) %>% slice(1) %>% pull(cluster)
  df_ok <- df_ok %>%
    mutate(inferred_sex = if_else(cluster == male_cluster, "M", "F"))

  # Merge back
  df %>%
    left_join(df_ok %>% select(sample, inferred_sex, inferred_conf), by = "sample")
}

# -----------------------------
# 6) Flag mismatches + ambiguity
# -----------------------------
flag_mismatches <- function(scores_df, meta_df, sample_id_col, reported_sex_col) {
  meta2 <- meta_df %>%
    mutate(
      sample = .data[[sample_id_col]],
      reported_sex = harmonize_sex(.data[[reported_sex_col]])
    ) %>%
    select(sample, reported_sex, everything())

  out <- scores_df %>%
    left_join(meta2, by = "sample") %>%
    mutate(
      inferred_sex = as.character(inferred_sex),
      concordance = case_when(
        reported_sex %in% c("M","F") & inferred_sex %in% c("M","F") & reported_sex == inferred_sex ~ "Concordant",
        reported_sex %in% c("M","F") & inferred_sex %in% c("M","F") & reported_sex != inferred_sex ~ "Discordant",
        TRUE ~ "Unknown/Ambiguous"
      ),
      # optional heuristic for ambiguity: low Y markers and low XIST, or both high
      # tune these cutoffs per dataset
      ambiguity_flag = case_when(
        !is.finite(xist_score) | !is.finite(y_score) ~ TRUE,
        (xist_score < 1 & y_score < 1) ~ TRUE,
        (xist_score > 5 & y_score > 3) ~ TRUE,
        TRUE ~ FALSE
      )
    )

  out
}

# -----------------------------
# 7) Plots
# -----------------------------
plot_scatter <- function(df) {
  ggplot(df, aes(x = xist_score, y = y_score)) +
    geom_point(aes(color = reported_sex, shape = inferred_sex), alpha = 0.85, size = 2.2) +
    geom_text(
      data = df %>% filter(concordance == "Discordant"),
      aes(label = sample),
      vjust = -0.6, size = 2.6, check_overlap = TRUE
    ) +
    labs(
      title = "Sex marker QC: XIST vs Y-marker score",
      x = "XIST (log2(CPM+1))",
      y = "Mean Y-marker score (log2(CPM+1))",
      color = "Reported sex",
      shape = "Inferred sex"
    ) +
    theme_bw()
}

plot_marker_heatmap <- function(logcpm, df, markers_y, marker_xist, top_n = 200) {
  genes <- rownames(logcpm)
  panel <- unique(c(marker_xist, markers_y))
  panel_present <- intersect(panel, genes)

  if (length(panel_present) < 2) {
    warning("Too few marker genes found for heatmap.")
    return(invisible(NULL))
  }

  # Order samples: inferred sex then discordance
  samp_order <- df %>%
    arrange(inferred_sex, desc(concordance == "Discordant"), reported_sex) %>%
    pull(sample)

  mat <- logcpm[panel_present, samp_order, drop = FALSE]

  ann <- df %>%
    select(sample, reported_sex, inferred_sex, concordance) %>%
    distinct() %>%
    tibble::column_to_rownames("sample") %>%
    .[samp_order, , drop = FALSE]

  pheatmap(
    mat,
    annotation_col = ann,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_colnames = FALSE,
    fontsize_row = 10,
    main = "Sex marker expression (log2(CPM+1))"
  )
}

# -----------------------------
# 8) Run end-to-end
# -----------------------------
run_sex_qc <- function(counts, meta, sample_id_col, reported_sex_col,
                       markers_y = markers_y, marker_xist = marker_xist) {

  # Align metadata to counts
  stopifnot(all(colnames(counts) %in% meta[[sample_id_col]]))
  meta <- meta %>% filter(.data[[sample_id_col]] %in% colnames(counts))

  logcpm <- compute_logcpm(counts)
  scores <- score_sex_markers(logcpm, markers_y, marker_xist)
  scores2 <- infer_sex(scores)
  results <- flag_mismatches(scores2, meta, sample_id_col, reported_sex_col)

  list(logcpm = logcpm, scores = results)
}

# -----------------------------
# 9) Example usage
# -----------------------------
# res <- run_sex_qc(counts, meta, sample_id_col = "sample_id", reported_sex_col = "sex")
# head(res$scores %>% arrange(desc(concordance == "Discordant")))

# p1 <- plot_scatter(res$scores)
# print(p1)

# plot_marker_heatmap(res$logcpm, res$scores, markers_y, marker_xist)

# Write a discordant table
# res$scores %>%
#   filter(concordance == "Discordant") %>%
#   arrange(desc(inferred_conf), desc(y_n_above)) %>%
#   write.csv("sex_mismatch_candidates.csv", row.names = FALSE)
