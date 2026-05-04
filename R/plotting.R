# R/plotting.R
# Dependencies: dplyr, ggplot2, tidyr, purrr

# Consistent color palette for all methods
ics_method_palette <- c(
  "ICS_bins=1"                     = "#E41A1C",
  "ICS_bins=2"                     = "#FF7F00",
  "ICS_bins=4"                     = "#4DAF4A",
  "ICS_bins=6"                     = "#377EB8",
  "ICS_bins=10"                    = "#984EA3",
  "ICS (norm; null bins=1)"        = "#E41A1C",
  "ICS (norm; null bins=2)"        = "#FF7F00",
  "ICS (norm; null bins=4)"        = "#4DAF4A",
  "ICS (norm; null bins=6)"        = "#377EB8",
  "ICS (norm; null bins=10)"       = "#984EA3",
  "GlobalMoran"                    = "#FF7F00",
  "Global Moran's I (uniform null)" = "#FF7F00",
  "Joincount"                      = "#A65628",
  "Join count (uniform null)"      = "#A65628"
)

plot_pvals_clean <- function(res_long, pathway_title = NULL) {

  df <- res_long %>%
    mutate(
      stat = case_when(
        grepl("^ICS", method)         ~ "ICS",
        grepl("^GlobalMoran", method) ~ "GlobalMoran",
        grepl("[Jj]oin", method)      ~ "Joincount",
        TRUE ~ "Other"
      ),
      # Color key: ICS gets bin-specific labels; classical methods get their method name
      color_key  = case_when(
        stat == "ICS" ~ paste0("ICS bins=", ics_bins),
        TRUE ~ method
      ),
      line_group = ifelse(stat == "ICS", paste0("ICS_", ics_bins), stat)
    )

  # Build a palette that covers all color_key values present
  palette <- ics_method_palette
  # Add any missing keys from the data
  missing_keys <- setdiff(unique(df$color_key), names(palette))
  if (length(missing_keys) > 0) {
    extra_cols <- setNames(rep("grey60", length(missing_keys)), missing_keys)
    palette <- c(palette, extra_cols)
  }

  ggplot(df, aes(x = pct_random, y = p_value, color = color_key,
                 linetype = stat, shape = stat, group = line_group)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = palette) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = sort(unique(df$pct_random))) +
    labs(
      title    = pathway_title %||% unique(df$pathway)[1],
      x        = "% Random genes (out of 100)",
      y        = "One-sided p-value",
      color    = "Method",
      linetype = "Statistic type",
      shape    = "Statistic type"
    ) +
    theme_minimal()
}


plot_effect_by_bins <- function(res_long) {
  df <- res_long %>%
    mutate(
      stat = case_when(
        grepl("^ICS", method)         ~ "ICS",
        grepl("^GlobalMoran", method) ~ "GlobalMoran",
        grepl("[Jj]oin", method)      ~ "Joincount",
        TRUE ~ "Other"
      ),
      bins_label = case_when(
        stat == "ICS" ~ paste0("bins=", ics_bins),
        TRUE ~ "uniform"
      )
    )

  ggplot(df, aes(x = pct_random, y = effect_directional, color = bins_label, linetype = stat)) +
    geom_hline(yintercept = 0, linewidth = 0.6) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = sort(unique(df$pct_random))) +
    labs(
      title    = paste0(unique(df$pathway)[1], " (Directional Effect)"),
      x        = "% Random genes (out of 100)",
      y        = "Directional effect size",
      color    = "Null (bins)",
      linetype = "Statistic"
    ) +
    theme_minimal()
}


plot_roc_full <- function(auc_tbl, title = "ROC curves (full)") {
  roc_df <- auc_tbl %>%
    select(method, auc, roc) %>%
    unnest(roc) %>%
    mutate(method = paste0(method, " (AUC=", round(auc, 3), ")"))

  ggplot(roc_df, aes(x = FPR, y = TPR, color = method)) +
    geom_line(linewidth = 1.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_equal() +
    theme_minimal(base_size = 13) +
    labs(title = title, x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", color = "Method")
}


plot_roc_zoom <- function(auc_tbl, fpr_max = 0.1, alpha = 0.05,
                          title = "ROC curves (zoom: FPR 0-0.1)") {
  roc_df <- auc_tbl %>%
    select(method, auc, roc) %>%
    unnest(roc) %>%
    mutate(method = paste0(method, " (AUC=", round(auc, 3), ")"))

  ggplot(roc_df, aes(x = FPR, y = TPR, color = method)) +
    geom_line(linewidth = 1.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_vline(xintercept = alpha, linetype = "dotted", linewidth = 1) +
    coord_cartesian(xlim = c(0, fpr_max), ylim = c(0, 1)) +
    theme_minimal(base_size = 13) +
    labs(
      title    = title,
      subtitle = paste0("Dashed = random classifier; dotted vertical = FPR=", alpha),
      x        = "False Positive Rate (FPR)",
      y        = "True Positive Rate (TPR)",
      color    = "Method"
    )
}


compute_tpr_at_fpr <- function(auc_tbl, fpr_target = 0.05) {
  auc_tbl %>%
    transmute(
      method,
      auc,
      tpr_at = purrr::map_dbl(roc, function(df) {
        # take best TPR achievable with FPR <= target
        max(df$TPR[df$FPR <= fpr_target], na.rm = TRUE)
      })
    ) %>%
    arrange(desc(tpr_at))
}


plot_auc_vs_signal <- function(auc_by_signal, title = "AUC vs signal strength (mixture simulation)") {
  auc_by_signal %>%
    select(signal_frac, method, auc) %>%
    ggplot(aes(x = signal_frac, y = auc, group = method, color = method)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = sort(unique(auc_by_signal$signal_frac))) +
    labs(
      title = title,
      x     = "Signal fraction",
      y     = "AUC",
      color = "Method"
    ) +
    theme_minimal(base_size = 13)
}


# --- LOO diagnostic plots ---

#' Plot gene-level LOO results: delta vs subgraph degree
plot_loo_importance_vs_contribution <- function(loo_df, title = "LOO: Importance vs Contribution") {
  cor_val <- cor(loo_df$subgraph_degree, loo_df$delta,
                 use = "complete.obs", method = "spearman")

  ggplot(loo_df, aes(x = subgraph_degree, y = delta)) +
    geom_point(aes(color = is_articulation), alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
    labs(
      title    = title,
      subtitle = paste0("Spearman r = ", round(cor_val, 3)),
      x        = "Subgraph degree (importance)",
      y        = "LOO delta (contribution to ICS signal)",
      color    = "Articulation point"
    ) +
    theme_minimal(base_size = 13)
}

#' Plot component-level LOO results
plot_loo_components <- function(comp_loo_df, title = "Component-level LOO") {
  comp_loo_df <- comp_loo_df %>%
    arrange(desc(abs(delta))) %>%
    mutate(comp_label = paste0("C", comp_id, " (n=", comp_size, ")"))

  ggplot(comp_loo_df, aes(x = reorder(comp_label, -delta), y = delta)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    labs(
      title = title,
      x     = "Component",
      y     = "LOO delta (ICS_full - ICS_without)"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Plot degree distribution: GWAS genes vs network background
plot_degree_diagnostic <- function(graph, gwas_nodes, title = "Degree distribution: GWAS vs background") {
  deg <- degree(graph)
  node_names <- V(graph)$name

  df <- data.frame(
    gene   = node_names,
    degree = deg,
    group  = ifelse(node_names %in% gwas_nodes, "GWAS genes", "Background")
  )

  ggplot(df, aes(x = degree, fill = group)) +
    geom_density(alpha = 0.4) +
    scale_x_log10() +
    labs(
      title = title,
      x     = "Degree (log scale)",
      y     = "Density",
      fill  = ""
    ) +
    theme_minimal(base_size = 13)
}
