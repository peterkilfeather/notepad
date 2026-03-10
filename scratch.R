library(tidyverse)
library(lme4)
library(performance)
library(rlang)

#------------------------------------------------------------------------------
# Helper: safe factor conversion
#------------------------------------------------------------------------------

.convert_metadata_types <- function(df, vars) {
  for (v in vars) {
    if (v %in% names(df)) {
      if (is.character(df[[v]]) || is.logical(df[[v]])) {
        df[[v]] <- factor(df[[v]])
      }
    }
  }
  droplevels(df)
}

#------------------------------------------------------------------------------
# Helper: single model comparison
#------------------------------------------------------------------------------

.fit_pc_lrt_single <- function(df, pc, var, subject_col, adjust_vars = NULL) {
  needed <- unique(c(pc, var, adjust_vars, subject_col))
  
  dsub <- df %>%
    select(all_of(needed)) %>%
    filter(if_all(everything(), ~ !is.na(.))) %>%
    droplevels()
  
  n_samples  <- nrow(dsub)
  n_subjects <- dplyr::n_distinct(dsub[[subject_col]])
  
  if (n_samples < 10 || n_subjects < 2) {
    return(tibble(
      pc = pc,
      variable = var,
      adjust_vars = paste(adjust_vars, collapse = " + "),
      variable_type = NA_character_,
      n_samples = n_samples,
      n_subjects = n_subjects,
      n_levels = NA_integer_,
      logLik_null = NA_real_,
      logLik_full = NA_real_,
      LRT_stat = NA_real_,
      df_diff = NA_real_,
      p_value = NA_real_,
      r2_marginal = NA_real_,
      r2_conditional = NA_real_,
      note = "Too few samples/subjects"
    ))
  }
  
  if (is.character(dsub[[var]]) || is.logical(dsub[[var]])) {
    dsub[[var]] <- factor(dsub[[var]])
  }
  
  variable_type <- if (is.factor(dsub[[var]])) "categorical" else "continuous"
  n_levels <- if (is.factor(dsub[[var]])) nlevels(dsub[[var]]) else NA_integer_
  
  if (is.factor(dsub[[var]]) && nlevels(dsub[[var]]) < 2) {
    return(tibble(
      pc = pc,
      variable = var,
      adjust_vars = paste(adjust_vars, collapse = " + "),
      variable_type = variable_type,
      n_samples = n_samples,
      n_subjects = n_subjects,
      n_levels = n_levels,
      logLik_null = NA_real_,
      logLik_full = NA_real_,
      LRT_stat = NA_real_,
      df_diff = NA_real_,
      p_value = NA_real_,
      r2_marginal = NA_real_,
      r2_conditional = NA_real_,
      note = "Variable has <2 levels after filtering"
    ))
  }
  
  rhs_null <- c(adjust_vars, sprintf("(1|%s)", subject_col))
  rhs_full <- c(adjust_vars, var, sprintf("(1|%s)", subject_col))
  
  if (length(adjust_vars) == 0) {
    null_formula <- as.formula(sprintf("%s ~ 1 + (1|%s)", pc, subject_col))
  } else {
    null_formula <- as.formula(sprintf("%s ~ %s", pc, paste(rhs_null, collapse = " + ")))
  }
  
  full_formula <- as.formula(sprintf("%s ~ %s", pc, paste(rhs_full, collapse = " + ")))
  
  null_fit <- tryCatch(
    lmer(null_formula, data = dsub, REML = FALSE),
    error = function(e) e
  )
  
  if (inherits(null_fit, "error")) {
    return(tibble(
      pc = pc,
      variable = var,
      adjust_vars = paste(adjust_vars, collapse = " + "),
      variable_type = variable_type,
      n_samples = n_samples,
      n_subjects = n_subjects,
      n_levels = n_levels,
      logLik_null = NA_real_,
      logLik_full = NA_real_,
      LRT_stat = NA_real_,
      df_diff = NA_real_,
      p_value = NA_real_,
      r2_marginal = NA_real_,
      r2_conditional = NA_real_,
      note = paste("Null model failed:", null_fit$message)
    ))
  }
  
  full_fit <- tryCatch(
    lmer(full_formula, data = dsub, REML = FALSE),
    error = function(e) e
  )
  
  if (inherits(full_fit, "error")) {
    return(tibble(
      pc = pc,
      variable = var,
      adjust_vars = paste(adjust_vars, collapse = " + "),
      variable_type = variable_type,
      n_samples = n_samples,
      n_subjects = n_subjects,
      n_levels = n_levels,
      logLik_null = as.numeric(logLik(null_fit)),
      logLik_full = NA_real_,
      LRT_stat = NA_real_,
      df_diff = NA_real_,
      p_value = NA_real_,
      r2_marginal = NA_real_,
      r2_conditional = NA_real_,
      note = paste("Full model failed:", full_fit$message)
    ))
  }
  
  an <- tryCatch(
    anova(null_fit, full_fit),
    error = function(e) e
  )
  
  if (inherits(an, "error")) {
    return(tibble(
      pc = pc,
      variable = var,
      adjust_vars = paste(adjust_vars, collapse = " + "),
      variable_type = variable_type,
      n_samples = n_samples,
      n_subjects = n_subjects,
      n_levels = n_levels,
      logLik_null = as.numeric(logLik(null_fit)),
      logLik_full = as.numeric(logLik(full_fit)),
      LRT_stat = NA_real_,
      df_diff = NA_real_,
      p_value = NA_real_,
      r2_marginal = NA_real_,
      r2_conditional = NA_real_,
      note = paste("anova() failed:", an$message)
    ))
  }
  
  r2vals <- tryCatch(
    performance::r2_nakagawa(full_fit),
    error = function(e) NULL
  )
  
  tibble(
    pc = pc,
    variable = var,
    adjust_vars = paste(adjust_vars, collapse = " + "),
    variable_type = variable_type,
    n_samples = n_samples,
    n_subjects = n_subjects,
    n_levels = n_levels,
    logLik_null = as.numeric(logLik(null_fit)),
    logLik_full = as.numeric(logLik(full_fit)),
    LRT_stat = an$Chisq[2],
    df_diff = an$`Chi Df`[2],
    p_value = an$`Pr(>Chisq)`[2],
    r2_marginal = if (!is.null(r2vals)) unname(r2vals$R2_marginal) else NA_real_,
    r2_conditional = if (!is.null(r2vals)) unname(r2vals$R2_conditional) else NA_real_,
    note = NA_character_
  )
}

#------------------------------------------------------------------------------
# Main wrapper
#------------------------------------------------------------------------------

run_pc_metadata_lrt <- function(
  pca,
  metadata,
  subject_col,
  vars_to_test,
  pcs = 1:10,
  sample_id_col = NULL,
  adjust_vars = NULL,
  fdr_scope = c("global", "within_pc", "within_variable")
) {
  fdr_scope <- match.arg(fdr_scope)
  
  if (!inherits(pca, "prcomp")) {
    stop("pca must be a prcomp object.")
  }
  
  if (is.null(rownames(pca$x))) {
    stop("pca$x must have rownames corresponding to sample IDs.")
  }
  
  pc_names <- paste0("PC", pcs)
  missing_pcs <- setdiff(pc_names, colnames(pca$x))
  if (length(missing_pcs) > 0) {
    stop("These PCs are not present in pca$x: ", paste(missing_pcs, collapse = ", "))
  }
  
  # Build PC score table
  pc_df <- as.data.frame(pca$x[, pc_names, drop = FALSE]) %>%
    rownames_to_column("sample_id")
  
  # Build metadata table with sample_id
  meta_df <- metadata
  
  if (is.null(sample_id_col)) {
    if (is.null(rownames(meta_df))) {
      stop("Provide sample_id_col or set rownames(metadata) to sample IDs.")
    }
    meta_df <- meta_df %>%
      rownames_to_column("sample_id")
  } else {
    if (!sample_id_col %in% names(meta_df)) {
      stop("sample_id_col not found in metadata.")
    }
    meta_df <- meta_df %>%
      rename(sample_id = all_of(sample_id_col))
  }
  
  # Checks
  required_cols <- unique(c("sample_id", subject_col, vars_to_test, adjust_vars))
  missing_cols <- setdiff(required_cols, names(meta_df))
  if (length(missing_cols) > 0) {
    stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
  }
  
  dat <- pc_df %>%
    inner_join(meta_df, by = "sample_id")
  
  if (nrow(dat) == 0) {
    stop("No overlapping sample IDs between pca$x rownames and metadata.")
  }
  
  dat[[subject_col]] <- factor(dat[[subject_col]])
  dat <- .convert_metadata_types(dat, unique(c(vars_to_test, adjust_vars)))
  
  # Avoid testing a variable while also adjusting for it
  vars_run <- setdiff(vars_to_test, adjust_vars)
  
  results <- purrr::map_dfr(
    pc_names,
    function(pc) {
      purrr::map_dfr(
        vars_run,
        ~ .fit_pc_lrt_single(
          df = dat,
          pc = pc,
          var = .x,
          subject_col = subject_col,
          adjust_vars = adjust_vars
        )
      )
    }
  )
  
  # FDR
  if (fdr_scope == "global") {
    results <- results %>%
      mutate(fdr = p.adjust(p_value, method = "fdr"))
  } else if (fdr_scope == "within_pc") {
    results <- results %>%
      group_by(pc) %>%
      mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
      ungroup()
  } else if (fdr_scope == "within_variable") {
    results <- results %>%
      group_by(variable) %>%
      mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
      ungroup()
  }
  
  results <- results %>%
    arrange(fdr, p_value)
  
  # Variance explained by each PC from prcomp
  pc_var_df <- tibble(
    pc = paste0("PC", seq_along(pca$sdev)),
    pc_variance_fraction = (pca$sdev^2) / sum(pca$sdev^2)
  ) %>%
    filter(pc %in% pc_names)
  
  results <- results %>%
    left_join(pc_var_df, by = "pc")
  
  # Heatmap-ready
  heatmap_df <- results %>%
    mutate(
      pc = factor(pc, levels = pc_names),
      variable = factor(variable, levels = rev(unique(vars_run))),
      minus_log10_fdr = ifelse(is.na(fdr) | fdr <= 0, NA_real_, -log10(fdr))
    )
  
  # Plots
  p_fdr <- ggplot(heatmap_df, aes(x = pc, y = variable, fill = minus_log10_fdr)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(na.value = "grey90") +
    labs(
      x = "Principal component",
      y = "Metadata variable",
      fill = "-log10(FDR)",
      title = "Metadata associations with principal components"
    ) +
    theme_bw(base_size = 12)
  
  p_r2 <- ggplot(heatmap_df, aes(x = pc, y = variable, fill = r2_marginal)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(na.value = "grey90") +
    labs(
      x = "Principal component",
      y = "Metadata variable",
      fill = "Marginal R²",
      title = "Effect-size summary for metadata-PC associations"
    ) +
    theme_bw(base_size = 12)
  
  list(
    input_data = dat,
    results = results,
    top_hits = results %>% filter(!is.na(fdr)) %>% arrange(fdr, p_value),
    heatmap_data = heatmap_df,
    plot_fdr = p_fdr,
    plot_r2 = p_r2
  )
}











vars_to_test <- c(
  "treatment",
  "timepoint",
  "plate",
  "seq_pool",
  "RIN",
  "library_size",
  "percent_mapped",
  "sex",
  "age"
)

res_unadjusted <- run_pc_metadata_lrt(
  pca = pca_obj,
  metadata = meta,
  subject_col = "subject_id",
  vars_to_test = vars_to_test,
  pcs = 1:10,
  sample_id_col = NULL,
  adjust_vars = NULL,
  fdr_scope = "global"
)

res_unadjusted$results
res_unadjusted$plot_fdr
res_unadjusted$plot_r2










plot_pc_vs_variable <- function(result_obj, pc, variable) {
  df <- result_obj$input_data
  
  if (is.character(df[[variable]]) || is.logical(df[[variable]])) {
    df[[variable]] <- factor(df[[variable]])
  }
  
  if (is.factor(df[[variable]])) {
    ggplot(df, aes(x = .data[[variable]], y = .data[[pc]])) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, height = 0, alpha = 0.5) +
      labs(
        x = variable,
        y = pc,
        title = paste(pc, "vs", variable)
      ) +
      theme_bw(base_size = 12)
  } else {
    ggplot(df, aes(x = .data[[variable]], y = .data[[pc]])) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(
        x = variable,
        y = pc,
        title = paste(pc, "vs", variable)
      ) +
      theme_bw(base_size = 12)
  }
}











res_adjusted <- run_pc_metadata_lrt(
  pca = pca_obj,
  metadata = meta,
  subject_col = "subject_id",
  vars_to_test = c(
    "plate",
    "seq_pool",
    "RIN",
    "library_size",
    "percent_mapped",
    "sex",
    "age",
    "treatment",
    "timepoint"
  ),
  pcs = 1:10,
  adjust_vars = c("treatment", "timepoint"),
  fdr_scope = "global"
)

res_adjusted$results
res_adjusted$plot_fdr






    
