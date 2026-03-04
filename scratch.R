# ============================================================
# Bulk RNA-seq count-matrix QC (TCGA-style) — RStudio script
# ============================================================
# What this does:
# - Loads raw counts (genes x samples) + sample metadata
# - Aligns metadata to counts
# - Basic library QC metrics
# - Gene filtering (edgeR::filterByExpr)
# - logCPM transform (TMM)
# - Sample-sample correlation + low median-correlation outliers
# - WGCNA-style sample connectivity (adjacency from correlation)
# - PCA distance outliers (top variable genes)
# - OPTIONAL: within-patient concordance summary (longitudinal)
# - Key plot: Sample connectivity (Z.k) vs PCA distance scatterplot
#
# Outputs (written to outdir):
#   sample_qc_metrics.csv
#   candidate_outliers.csv
#   qc_thresholds.txt
#   patient_within_correlation_summary.csv (if patient_id present)
#   plots/*.png
#   logCPM_filtered.rds
# ============================================================

# ----------------------------
# 0) User inputs (EDIT THESE)
# ----------------------------
counts_path <- "counts.rds"      # .rds (preferred) OR .csv/.tsv with genes in rows, samples in columns, gene IDs in first column
meta_path   <- "metadata.csv"    # .csv/.tsv containing at least sample_id

outdir <- "qc_out"

# Metadata column names (edit if your file differs)
sample_col   <- "sample_id"
patient_col  <- "patient_id"     # optional
timepoint_col <- "timepoint"     # optional
treatment_col <- "treatment"     # optional

# Parameters (reasonable defaults for ~1500 samples)
top_var_genes <- 3000
n_pcs <- 5
corr_bottom_quantile <- 0.01     # bottom 1% by median correlation
z_connectivity_cut <- -2.0       # Z.k cutoff
pca_sd_cut <- 3.0                # mean + 3*sd
remove_rule_min_flags <- 2       # candidate outlier if flagged by >=2 methods
prior_count <- 1.0               # logCPM prior.count
adj_power <- 2.0                 # adjacency = |cor|^power (2 is common/simple)

# ----------------------------
# 1) Libraries
# ----------------------------
suppressPackageStartupMessages({
  library(edgeR)
  library(WGCNA)         # optional; loaded for familiarity (not strictly required)
  library(matrixStats)   # rowVars
  library(tidyverse)
})

# ----------------------------
# 2) Create output folders
# ----------------------------
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
plots_dir <- file.path(outdir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 3) Helper functions (I/O)
# ----------------------------
read_counts <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    x <- readRDS(path)
    if (is.data.frame(x)) x <- as.matrix(x)
    return(x)
  }
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    x <- read.csv(path, row.names = 1, check.names = FALSE)
    return(as.matrix(x))
  }
  if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) {
    x <- read.delim(path, row.names = 1, check.names = FALSE)
    return(as.matrix(x))
  }
  stop("Unsupported counts format. Use .rds, .csv, .tsv, or .txt")
}

read_meta <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) return(read.csv(path, stringsAsFactors = FALSE))
  if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) return(read.delim(path, stringsAsFactors = FALSE))
  stop("Unsupported metadata format. Use .csv, .tsv, or .txt")
}

save_png <- function(filename, width=1200, height=900, res=150) {
  png(filename, width=width, height=height, res=res)
}

# ----------------------------
# 4) Load data
# ----------------------------
counts <- read_counts(counts_path)
meta   <- read_meta(meta_path)

stopifnot(!is.null(colnames(counts)))
stopifnot(!is.null(rownames(counts)))
stopifnot(sample_col %in% colnames(meta))

meta <- meta %>% mutate(.sample = .data[[sample_col]])

# Drop metadata rows not present in counts
meta_in_counts <- meta$.sample %in% colnames(counts)
if (!all(meta_in_counts)) {
  missing <- meta$.sample[!meta_in_counts]
  warning(length(missing), " metadata samples not found in counts and will be dropped.\n",
          "Example: ", paste(head(missing, 5), collapse=", "))
  meta <- meta %>% filter(.sample %in% colnames(counts))
}

# Reorder counts to match metadata order
counts <- counts[, meta$.sample, drop=FALSE]

# Drop any all-zero samples (guard)
allzero <- colSums(counts) == 0
if (any(allzero)) {
  warning("Dropping all-zero samples: ", paste(colnames(counts)[allzero], collapse=", "))
  counts <- counts[, !allzero, drop=FALSE]
  meta <- meta[!allzero, , drop=FALSE]
}

# ----------------------------
# 5) Basic library QC metrics
# ----------------------------
lib_size       <- colSums(counts)
genes_detected <- colSums(counts > 0)
zero_fraction  <- colMeans(counts == 0)

qc_basic <- tibble(
  sample = colnames(counts),
  library_size = as.numeric(lib_size),
  genes_detected = as.integer(genes_detected),
  zero_fraction = as.numeric(zero_fraction)
)

# Attach optional metadata columns if present
if (patient_col %in% colnames(meta))    qc_basic$patient_id <- meta[[patient_col]]
if (timepoint_col %in% colnames(meta))  qc_basic$timepoint  <- meta[[timepoint_col]]
if (treatment_col %in% colnames(meta))  qc_basic$treatment  <- meta[[treatment_col]]

# ----------------------------
# 6) Gene filtering + logCPM transform
# ----------------------------
dge <- DGEList(counts = counts)
keep_genes <- filterByExpr(dge)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

logCPM <- cpm(dge, log = TRUE, prior.count = prior_count)
saveRDS(logCPM, file.path(outdir, "logCPM_filtered.rds"))

# ----------------------------
# 7) Correlation QC (median correlation)
# ----------------------------
cor_matrix <- cor(logCPM, use = "pairwise.complete.obs", method = "pearson")
median_corr <- apply(cor_matrix, 2, median, na.rm = TRUE)

corr_threshold <- as.numeric(quantile(median_corr, probs = corr_bottom_quantile, na.rm = TRUE))
corr_outliers <- names(median_corr)[median_corr < corr_threshold]

qc_corr <- tibble(
  sample = names(median_corr),
  median_corr = as.numeric(median_corr),
  corr_outlier = sample %in% corr_outliers
)

# ----------------------------
# 8) Sample connectivity (WGCNA-style)
# ----------------------------
adj <- abs(cor_matrix)^adj_power
diag(adj) <- 0
connectivity <- rowSums(adj, na.rm = TRUE)

z_connectivity <- as.numeric(scale(connectivity))
names(z_connectivity) <- names(connectivity)

connectivity_outliers <- names(connectivity)[z_connectivity < z_connectivity_cut]

qc_conn <- tibble(
  sample = names(connectivity),
  connectivity = as.numeric(connectivity),
  z_connectivity = as.numeric(z_connectivity),
  connectivity_outlier = sample %in% connectivity_outliers
)

# ----------------------------
# 9) PCA + PCA distance
# ----------------------------
top_var_genes <- min(top_var_genes, nrow(logCPM))
vars <- rowVars(logCPM, na.rm = TRUE)
var_idx <- order(vars, decreasing = TRUE)[seq_len(top_var_genes)]

pca <- prcomp(t(logCPM[var_idx, , drop=FALSE]), center = TRUE, scale. = TRUE)

n_pcs <- min(n_pcs, ncol(pca$x))
pc_scores <- pca$x[, seq_len(n_pcs), drop=FALSE]
centroid <- colMeans(pc_scores)

pca_distance <- sqrt(rowSums((pc_scores - matrix(centroid, nrow=nrow(pc_scores), ncol=n_pcs, byrow=TRUE))^2))

pca_threshold <- mean(pca_distance, na.rm = TRUE) + pca_sd_cut * sd(pca_distance, na.rm = TRUE)
pca_outliers <- rownames(pc_scores)[pca_distance > pca_threshold]

qc_pca <- tibble(
  sample = rownames(pc_scores),
  pca_distance = as.numeric(pca_distance),
  pca_outlier = sample %in% pca_outliers
)

# ----------------------------
# 10) Combine QC metrics + define candidate outliers
# ----------------------------
qc <- qc_basic %>%
  left_join(qc_corr, by="sample") %>%
  left_join(qc_conn, by="sample") %>%
  left_join(qc_pca,  by="sample") %>%
  mutate(
    flag_corr = as.integer(corr_outlier),
    flag_connectivity = as.integer(connectivity_outlier),
    flag_pca = as.integer(pca_outlier),
    total_flags = flag_corr + flag_connectivity + flag_pca,
    candidate_outlier = total_flags >= remove_rule_min_flags
  )

candidate_outliers <- qc %>% filter(candidate_outlier)

# ----------------------------
# 11) Optional: within-patient concordance summary
# ----------------------------
if (patient_col %in% colnames(meta)) {
  tmp <- tibble(sample = meta[[sample_col]], patient_id = meta[[patient_col]]) %>%
    filter(!is.na(patient_id))

  within_patient_mean <- rep(NA_real_, nrow(tmp))
  names(within_patient_mean) <- tmp$sample

  for (pid in unique(tmp$patient_id)) {
    s <- tmp$sample[tmp$patient_id == pid]
    if (length(s) <= 1) next
    sub <- cor_matrix[s, s, drop=FALSE]
    diag(sub) <- NA
    within_patient_mean[s] <- rowMeans(sub, na.rm = TRUE)
  }

  qc <- qc %>%
    mutate(within_patient_mean_corr = within_patient_mean[sample])

  patient_summary <- qc %>%
    filter(!is.na(patient_id)) %>%
    group_by(patient_id) %>%
    summarise(
      n_samples = n(),
      mean_within_patient_corr = mean(within_patient_mean_corr, na.rm = TRUE),
      min_within_patient_corr  = min(within_patient_mean_corr, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv(patient_summary, file.path(outdir, "patient_within_correlation_summary.csv"), row.names = FALSE)
}

# ----------------------------
# 12) Save tables + thresholds
# ----------------------------
write.csv(qc, file.path(outdir, "sample_qc_metrics.csv"), row.names = FALSE)
write.csv(candidate_outliers, file.path(outdir, "candidate_outliers.csv"), row.names = FALSE)

thr_path <- file.path(outdir, "qc_thresholds.txt")
writeLines(c(
  paste0("Date: ", Sys.time()),
  paste0("Counts file: ", counts_path),
  paste0("Meta file:   ", meta_path),
  "",
  "Thresholds:",
  paste0("Median-correlation outliers: < quantile(", corr_bottom_quantile, ") = ", signif(corr_threshold, 4)),
  paste0("Connectivity outliers: Z.k < ", z_connectivity_cut),
  paste0("PCA-distance outliers: > mean + ", pca_sd_cut, "*sd = ", signif(pca_threshold, 4)),
  paste0("Candidate outlier rule: total_flags >= ", remove_rule_min_flags),
  "",
  "Notes:",
  paste0("Adjacency power for connectivity: ", adj_power),
  paste0("Top variable genes for PCA: ", top_var_genes),
  paste0("PCs used for PCA-distance: ", n_pcs),
  paste0("Filtered genes: ", nrow(counts), " -> ", nrow(logCPM)),
  paste0("Samples: ", ncol(logCPM))
), thr_path)

# ----------------------------
# 13) Plots (PNG)
# ----------------------------
# 13.1 Library size histogram (log10)
save_png(file.path(plots_dir, "library_size_log10_hist.png"))
hist(log10(qc$library_size), main="Library size (log10)", xlab="log10(total counts)", breaks=50)
dev.off()

# 13.2 Genes detected
save_png(file.path(plots_dir, "genes_detected_hist.png"))
hist(qc$genes_detected, main="Genes detected per sample", xlab="# genes with >0 counts", breaks=50)
dev.off()

# 13.3 Zero fraction
save_png(file.path(plots_dir, "zero_fraction_hist.png"))
hist(qc$zero_fraction, main="Zero-count fraction per sample", xlab="Fraction of genes with 0 counts", breaks=50)
dev.off()

# 13.4 Median correlation distribution
save_png(file.path(plots_dir, "median_correlation_hist.png"))
hist(qc$median_corr, main="Median sample correlation", xlab="Median correlation to all samples", breaks=50)
abline(v=corr_threshold, lwd=2)
dev.off()

# 13.5 PCA: PC1 vs PC2 (candidate outliers highlighted)
pca_df <- as_tibble(pca$x[, 1:3, drop=FALSE], rownames="sample") %>%
  left_join(qc %>% select(sample, candidate_outlier, total_flags, any_of(c("treatment","timepoint","patient_id"))),
            by="sample")

save_png(file.path(plots_dir, "pca_pc1_pc2.png"))
plot(pca_df$PC1, pca_df$PC2, xlab="PC1", ylab="PC2", main="PCA (top variable genes)")
points(pca_df$PC1[pca_df$candidate_outlier], pca_df$PC2[pca_df$candidate_outlier], pch=19)
dev.off()

# 13.6 KEY PLOT: Sample connectivity (Z.k) vs PCA distance
save_png(file.path(plots_dir, "connectivity_zk_vs_pca_distance.png"))
plot(qc$z_connectivity, qc$pca_distance,
     xlab="Sample connectivity (Z.k)", ylab="PCA distance (top PCs)",
     main="Sample connectivity vs PCA distance")
abline(v=z_connectivity_cut, lwd=2)
abline(h=pca_threshold, lwd=2)
points(qc$z_connectivity[qc$candidate_outlier], qc$pca_distance[qc$candidate_outlier], pch=19)
dev.off()

# 13.7 Optional: within-patient mean correlation distribution
if ("within_patient_mean_corr" %in% colnames(qc)) {
  save_png(file.path(plots_dir, "within_patient_mean_correlation_hist.png"))
  hist(qc$within_patient_mean_corr,
       main="Within-patient mean correlation (per sample)",
       xlab="Mean correlation to other samples from same patient", breaks=50)
  dev.off()
}

# ----------------------------
# 14) Quick console summary + return objects in session
# ----------------------------
message("QC complete. Outputs in: ", normalizePath(outdir))
message("Filtered genes: ", nrow(counts), " -> ", nrow(logCPM))
message("Samples: ", ncol(logCPM))
message("Outliers flagged:")
message("  Median-corr (bottom quantile): ", length(corr_outliers))
message("  Connectivity (Z.k < cut):      ", length(connectivity_outliers))
message("  PCA distance (> threshold):     ", length(pca_outliers))
message("  Candidate outliers (>= flags):  ", nrow(candidate_outliers))

# Useful objects left in the environment:
#   qc, candidate_outliers, logCPM, dge, cor_matrix, pca



library(ggplot2)
library(dplyr)

# Extract PCA scores
pca_df <- as.data.frame(pca$x[,1:2])
pca_df$sample <- rownames(pca_df)

# Compute centroid
centroid <- colMeans(pca_df[,c("PC1","PC2")])

# Distance from centroid
pca_df$dist <- sqrt(
  (pca_df$PC1 - centroid["PC1"])^2 +
  (pca_df$PC2 - centroid["PC2"])^2
)

# Standard deviation of distances
sd_dist <- sd(pca_df$dist)

# Function to create circle coordinates
circle_df <- function(radius, center, n=200){
  tibble(
    x = center["PC1"] + radius * cos(seq(0,2*pi,length.out=n)),
    y = center["PC2"] + radius * sin(seq(0,2*pi,length.out=n)),
    sd = radius/sd_dist
  )
}

# Circles for 1,2,3 SD
circles <- bind_rows(
  circle_df(sd_dist, centroid),
  circle_df(2*sd_dist, centroid),
  circle_df(3*sd_dist, centroid)
)

# Plot
ggplot(pca_df, aes(PC1, PC2)) +

  geom_point(alpha=0.7, size=2) +

  geom_path(
    data=circles,
    aes(x,y,group=sd,linetype=factor(sd)),
    colour="red",
    linewidth=1
  ) +

  geom_point(
    aes(x=centroid["PC1"], y=centroid["PC2"]),
    colour="blue",
    size=4
  ) +

  labs(
    title="PCA Sample Distribution",
    subtitle="Contours represent 1, 2 and 3 SD from centroid",
    linetype="SD from centroid"
  ) +

  theme_bw()

------------------------------------------------------------------------------------


# ============================================================
# PCA outlier calling from prcomp() scores (10 PCs)
# - Implements "more than 2*SD away from the centroid"
# - Implements a more scientifically appropriate method:
#     robust Mahalanobis distance in PC space (MCD) with chi-square cutoff
# - Visualises both on:
#     1) PC1 vs PC2 scatter (centroid highlighted)
#     2) Distance-to-centroid distribution (Euclidean; client + sensible variant)
#     3) Mahalanobis distance diagnostic (QQ-style) + threshold
# ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)

# Robust covariance for Mahalanobis (recommended)
# install.packages("rrcov") if needed
library(rrcov)

# -----------------------------
# Inputs assumed to exist:
#   pca  <- prcomp(t(logCPM[var_idx,]), center=TRUE, scale.=TRUE)
# Optional:
#   meta data frame with a sample_id column to join (edit names below)
# -----------------------------
n_pcs <- 10
scores <- as.data.frame(pca$x[, 1:n_pcs, drop = FALSE])
scores$sample <- rownames(scores)

# If you have metadata you want to merge (optional)
# meta must have sample_id that matches rownames(pca$x)
# scores <- scores %>% left_join(meta, by = c("sample" = "sample_id"))

# -----------------------------
# 1) Euclidean distance to centroid in 10D PC space
# -----------------------------
pc_cols <- paste0("PC", 1:n_pcs)
centroid <- colMeans(scores[, pc_cols, drop = FALSE])

scores <- scores %>%
  mutate(
    dist_euclid = sqrt(rowSums((as.matrix(across(all_of(pc_cols))) -
                                 matrix(centroid, nrow = n(), ncol = n_pcs, byrow = TRUE))^2))
  )

sd_dist  <- sd(scores$dist_euclid)
mean_dist <- mean(scores$dist_euclid)

# Client interpretation A (literal): distance > 2*SD (from 0)
# This is uncommon because distances are not mean-zero.
thr_client_literal <- 2 * sd_dist

# Client interpretation B (more common "2 SD rule"): distance > mean + 2*SD
thr_client_meanplus <- mean_dist + 2 * sd_dist

scores <- scores %>%
  mutate(
    out_client_literal  = dist_euclid > thr_client_literal,
    out_client_meanplus = dist_euclid > thr_client_meanplus
  )

# -----------------------------
# 2) Scientifically preferred: robust Mahalanobis distance in 10D PC space
#    - Use MCD robust covariance to reduce influence of outliers
#    - Under multivariate normality, MD^2 ~ Chi-square(df = n_pcs)
# -----------------------------
X <- as.matrix(scores[, pc_cols, drop = FALSE])

mcd <- CovMcd(X)  # robust location + covariance
md2 <- mahalanobis(X, center = mcd@center, cov = mcd@cov)

# Choose a cutoff. Common QC choices:
#  - 0.975: "more inclusive" (flags more)
#  - 0.99 : "stricter"
#  - 0.999: "very strict"
alpha <- 0.99
thr_md2 <- qchisq(alpha, df = n_pcs)

scores <- scores %>%
  mutate(
    md2 = md2,
    out_robust_md = md2 > thr_md2
  )

# -----------------------------
# 3) Visualisation A: PC1 vs PC2 scatter (centroid shown)
#    Note: outlier calls are in 10D; points may not look extreme in PC1/PC2.
# -----------------------------
centroid_12 <- tibble(PC1 = centroid["PC1"], PC2 = centroid["PC2"])

p_pc12_client <- ggplot(scores, aes(PC1, PC2)) +
  geom_point(aes(color = out_client_meanplus), alpha = 0.7, size = 2) +
  geom_point(data = centroid_12, aes(PC1, PC2), inherit.aes = FALSE, shape = 4, size = 5, stroke = 1.2) +
  labs(
    title = "PCA (PC1 vs PC2) with client-style outliers",
    subtitle = sprintf("Outliers: Euclidean distance in PC1–PC%d > mean + 2*SD", n_pcs),
    color = "Outlier"
  ) +
  theme_bw()

p_pc12_robust <- ggplot(scores, aes(PC1, PC2)) +
  geom_point(aes(color = out_robust_md), alpha = 0.7, size = 2) +
  geom_point(data = centroid_12, aes(PC1, PC2), inherit.aes = FALSE, shape = 4, size = 5, stroke = 1.2) +
  labs(
    title = "PCA (PC1 vs PC2) with robust multivariate outliers",
    subtitle = sprintf("Outliers: robust Mahalanobis distance in PC1–PC%d (MCD), cutoff Chi-square(DF=%d) at %.3f", n_pcs, n_pcs, alpha),
    color = "Outlier"
  ) +
  theme_bw()

p_pc12_client
p_pc12_robust

# -----------------------------
# 4) Visualisation B: Distance-to-centroid distribution with thresholds
# -----------------------------
dist_long <- scores %>%
  select(sample, dist_euclid) %>%
  mutate(dummy = 1)

p_dist <- ggplot(scores, aes(x = dist_euclid)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = thr_client_literal, linewidth = 1, linetype = 2) +
  geom_vline(xintercept = thr_client_meanplus, linewidth = 1, linetype = 1) +
  labs(
    title = "Euclidean distance to PCA centroid (PC1–PC10)",
    subtitle = paste(
      sprintf("Client literal: dist > 2*SD = %.3f (dashed)", thr_client_literal),
      sprintf("2SD rule: dist > mean + 2*SD = %.3f (solid)", thr_client_meanplus),
      sep = " | "
    ),
    x = sprintf("Euclidean distance in PC1–PC%d space", n_pcs),
    y = "Count"
  ) +
  theme_bw()

p_dist

# -----------------------------
# 5) Visualisation C: Robust Mahalanobis diagnostic plot (QQ-style)
#    If points deviate strongly above the line, there are heavy-tail/outliers.
# -----------------------------
qq_df <- tibble(
  md2_sorted = sort(scores$md2),
  chi_sorted = qchisq(ppoints(nrow(scores)), df = n_pcs)
)

p_md_qq <- ggplot(qq_df, aes(x = chi_sorted, y = md2_sorted)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_abline(intercept = 0, slope = 1, linewidth = 1) +
  geom_hline(yintercept = thr_md2, linetype = 1, linewidth = 1) +
  labs(
    title = "Robust Mahalanobis distance diagnostic (PC1–PC10)",
    subtitle = sprintf("Horizontal line = Chi-square cutoff at alpha=%.3f (%.3f)", alpha, thr_md2),
    x = sprintf("Theoretical Chi-square quantiles (df=%d)", n_pcs),
    y = "Observed robust MD^2 (sorted)"
  ) +
  theme_bw()

p_md_qq






# distance in PC1–PC10
scores$dist_pc10 <- sqrt(
  rowSums((scores[,1:10] - 
           matrix(centroid, nrow=nrow(scores), ncol=10, byrow=TRUE))^2)
)

# flag the Mahalanobis outliers you already calculated
scores$outlier <- scores$md2 > qchisq(0.999, df=10)

ggplot(scores, aes(dist_pc12, dist_pc10)) +

  geom_point(aes(color=outlier), size=2, alpha=0.8) +

  labs(
    x="Distance from centroid (PC1–PC2)",
    y="Distance from centroid (PC1–PC10)",
    title="Comparison of PCA distances",
    subtitle="Outliers are defined in the full 10-PC space"
  ) +

  theme_bw()























df <- qc %>%
  select(sample, z_connectivity, median_corr) %>%
  left_join(scores %>% select(sample, md2), by = "sample")

# thresholds (edit if you used different values)
n_pcs <- 10
alpha <- 0.999
thr_md2 <- qchisq(alpha, df = n_pcs)
thr_zk  <- -2

# optional median-corr threshold (e.g. bottom 1%)
thr_medcorr <- quantile(df$median_corr, 0.01, na.rm = TRUE)

df <- df %>%
  mutate(
    flag_md2 = md2 > thr_md2,
    flag_zk  = z_connectivity < thr_zk,
    flag_cor = median_corr < thr_medcorr,
    n_flags  = flag_md2 + flag_zk + flag_cor,
    candidate_outlier = n_flags >= 2
  )

p_combo <- ggplot(df, aes(x = z_connectivity, y = md2)) +
  geom_point(aes(color = candidate_outlier, size = 1 - median_corr), alpha = 0.75) +
  geom_vline(xintercept = thr_zk, linewidth = 1) +
  geom_hline(yintercept = thr_md2, linewidth = 1) +
  labs(
    title = "Combined sample QC: connectivity vs robust PCA outlier score",
    subtitle = paste0(
      "Horizontal line: robust MD² cutoff (Chi-square df=", n_pcs, ", alpha=", alpha, "). ",
      "Vertical line: connectivity cutoff (Z.k=", thr_zk, "). ",
      "Point size increases as median correlation decreases."
    ),
    x = "Sample connectivity (Z.k)  (lower = more dissimilar)",
    y = "Robust Mahalanobis distance² (PC1–PC10)  (higher = more outlying)",
    color = "Candidate outlier",
    size = "1 - median corr"
  ) +
  theme_bw()

