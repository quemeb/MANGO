# R/simulation.R
# Dependencies: igraph, dplyr
# Requires: prep_graph_fast(), ics_norm(), global_moran_I(),
#           fast_joincount(), sample_degree_matched() (from stats.R, sampling.R)

run_pathway_noise_sweep <- function(pathway_name,
                                    n_total = 100,
                                    noise_steps = seq(0, 1, by = 0.2),
                                    n_iterations = 1000,
                                    num_bins_ics = 4,
                                    graph,
                                    pathway_data,
                                    bin_cache,
                                    seed = 42,
                                    replace_null = TRUE,
                                    run_uniform_null = FALSE) {

  prep <- prep_graph_fast(graph)
  g <- prep$g; A <- prep$A; nodes <- prep$nodes

  row <- pathway_data %>% filter(PATHWAY_NAMES == pathway_name)
  if (nrow(row) == 0) stop("Pathway name not found: ", pathway_name)

  pathway_genes    <- unique(unlist(row$Genes))
  background_genes <- setdiff(nodes, pathway_genes)

  out <- vector("list", length(noise_steps))

  for (k in seq_along(noise_steps)) {
    frac_rand <- noise_steps[k]
    n_rand    <- as.integer(round(n_total * frac_rand))
    n_path    <- n_total - n_rand

    set.seed(as.integer(seed) + k)

    if (n_path > length(pathway_genes))    stop("Not enough pathway genes for n_path=", n_path)
    if (n_rand > length(background_genes)) stop("Not enough background genes for n_rand=", n_rand)

    target_nodes <- c(
      if (n_path > 0) sample(pathway_genes,    n_path) else character(0),
      if (n_rand > 0) sample(background_genes, n_rand) else character(0)
    )

    # --- observed stats ---
    idx_target <- match(target_nodes, nodes)
    x_obs <- integer(length(nodes)); x_obs[idx_target] <- 1L

    obs_ics    <- ics_norm(induced_subgraph(g, vids = target_nodes))
    obs_gmoran <- global_moran_I(A, x_obs)
    obs_jc     <- fast_joincount(A, x_obs)

    # --- degree-matched null (primary) — shared sample for ALL methods ---
    rand_ics    <- numeric(n_iterations)
    rand_gmoran <- numeric(n_iterations)
    rand_jc     <- numeric(n_iterations)

    for (i in seq_len(n_iterations)) {
      s_dm <- sample_degree_matched(target_nodes, bin_cache,
                                    num_bins = num_bins_ics, replace = replace_null)
      idx_s <- match(s_dm, nodes)
      x_s   <- integer(length(nodes)); x_s[idx_s] <- 1L

      rand_ics[i]    <- ics_norm(induced_subgraph(g, vids = s_dm))
      rand_gmoran[i] <- global_moran_I(A, x_s)
      rand_jc[i]     <- fast_joincount(A, x_s)
    }

    p_ics    <- mean(rand_ics    <= obs_ics)
    p_gmoran <- mean(rand_gmoran >= obs_gmoran)
    p_jc     <- mean(rand_jc     >= obs_jc)

    row_out <- data.frame(
      pathway    = pathway_name,
      pct_random = 100 * frac_rand,
      n_total    = n_total,
      n_pathway  = n_path,
      n_random   = n_rand,
      null_type  = "degree_matched",
      ics_bins   = num_bins_ics,

      obs_ics_norm       = obs_ics,
      null_mean_ics_norm = mean(rand_ics),
      null_sd_ics_norm   = sd(rand_ics),
      p_ics_norm         = p_ics,

      obs_global_moran          = obs_gmoran,
      null_mean_global_moran    = mean(rand_gmoran),
      null_sd_global_moran      = sd(rand_gmoran),
      p_global_moran            = p_gmoran,

      obs_joincount        = obs_jc,
      null_mean_joincount  = mean(rand_jc),
      null_sd_joincount    = sd(rand_jc),
      p_joincount          = p_jc
    )

    # --- optional uniform null (sensitivity analysis) ---
    if (run_uniform_null) {
      rand_gmoran_u <- numeric(n_iterations)
      rand_jc_u     <- numeric(n_iterations)
      rand_ics_u    <- numeric(n_iterations)

      for (i in seq_len(n_iterations)) {
        u     <- sample(nodes, n_total, replace = replace_null)
        idx_u <- match(u, nodes)
        x_u   <- integer(length(nodes)); x_u[idx_u] <- 1L

        rand_ics_u[i]    <- ics_norm(induced_subgraph(g, vids = u))
        rand_gmoran_u[i] <- global_moran_I(A, x_u)
        rand_jc_u[i]     <- fast_joincount(A, x_u)
      }

      row_uniform <- data.frame(
        pathway    = pathway_name,
        pct_random = 100 * frac_rand,
        n_total    = n_total,
        n_pathway  = n_path,
        n_random   = n_rand,
        null_type  = "uniform",
        ics_bins   = NA_integer_,

        obs_ics_norm       = obs_ics,
        null_mean_ics_norm = mean(rand_ics_u),
        null_sd_ics_norm   = sd(rand_ics_u),
        p_ics_norm         = mean(rand_ics_u <= obs_ics),

        obs_global_moran          = obs_gmoran,
        null_mean_global_moran    = mean(rand_gmoran_u),
        null_sd_global_moran      = sd(rand_gmoran_u),
        p_global_moran            = mean(rand_gmoran_u >= obs_gmoran),

        obs_joincount        = obs_jc,
        null_mean_joincount  = mean(rand_jc_u),
        null_sd_joincount    = sd(rand_jc_u),
        p_joincount          = mean(rand_jc_u >= obs_jc)
      )
      row_out <- bind_rows(row_out, row_uniform)
    }

    out[[k]] <- row_out
  }

  bind_rows(out)
}


run_pathway_noise_sweep_compare_ics_bins <- function(pathway_name,
                                                     n_total = 100,
                                                     noise_steps = seq(0, 1, by = 0.2),
                                                     n_iterations = 1000,
                                                     ics_bins_set = c(1, 2, 4, 6, 10),
                                                     graph,
                                                     pathway_data,
                                                     bin_cache,
                                                     seed = 42,
                                                     replace_null = TRUE,
                                                     run_uniform_null = FALSE) {

  prep <- prep_graph_fast(graph)
  g <- prep$g; A <- prep$A; nodes <- prep$nodes

  row <- pathway_data %>% filter(PATHWAY_NAMES == pathway_name)
  if (nrow(row) == 0) stop("Pathway name not found: ", pathway_name)

  pathway_genes    <- unique(unlist(row$Genes))
  background_genes <- setdiff(nodes, pathway_genes)

  rows <- list()

  for (k in seq_along(noise_steps)) {
    frac_rand <- noise_steps[k]
    n_rand    <- as.integer(round(n_total * frac_rand))
    n_path    <- n_total - n_rand

    set.seed(as.integer(seed) + k)

    if (n_path > length(pathway_genes))    stop("Not enough pathway genes for n_path=", n_path)
    if (n_rand > length(background_genes)) stop("Not enough background genes for n_rand=", n_rand)

    target_nodes <- c(
      if (n_path > 0) sample(pathway_genes,    n_path) else character(0),
      if (n_rand > 0) sample(background_genes, n_rand) else character(0)
    )

    # --- observed ---
    idx_target <- match(target_nodes, nodes)
    x_obs <- integer(length(nodes)); x_obs[idx_target] <- 1L

    obs_ics    <- ics_norm(induced_subgraph(g, vids = target_nodes))
    obs_gmoran <- global_moran_I(A, x_obs)
    obs_jc     <- fast_joincount(A, x_obs)

    # ---- Degree-matched null per bin scheme ----
    for (b in ics_bins_set) {
      rand_ics    <- numeric(n_iterations)
      rand_gmoran <- numeric(n_iterations)
      rand_jc     <- numeric(n_iterations)

      for (i in seq_len(n_iterations)) {
        s_dm <- if (b > 1) {
          sample_degree_matched(target_nodes, bin_cache, num_bins = b, replace = replace_null)
        } else {
          sample(nodes, n_total, replace = replace_null)
        }
        idx_s <- match(s_dm, nodes)
        x_s   <- integer(length(nodes)); x_s[idx_s] <- 1L

        rand_ics[i]    <- ics_norm(induced_subgraph(g, vids = s_dm))
        rand_gmoran[i] <- global_moran_I(A, x_s)
        rand_jc[i]     <- fast_joincount(A, x_s)
      }

      null_label <- if (b > 1) paste0("degree-matched (bins=", b, ")") else "uniform"

      rows[[length(rows) + 1]] <- data.frame(
        pathway = pathway_name, pct_random = 100 * frac_rand,
        method = paste0("ICS (norm; null bins=", b, ")"), ics_bins = b,
        null_type = null_label,
        observed = obs_ics, null_mean = mean(rand_ics), null_sd = sd(rand_ics),
        p_value = mean(rand_ics <= obs_ics),
        effect_directional = mean(rand_ics) - obs_ics
      )

      rows[[length(rows) + 1]] <- data.frame(
        pathway = pathway_name, pct_random = 100 * frac_rand,
        method = paste0("Global Moran's I (", null_label, ")"), ics_bins = b,
        null_type = null_label,
        observed = obs_gmoran, null_mean = mean(rand_gmoran), null_sd = sd(rand_gmoran),
        p_value = mean(rand_gmoran >= obs_gmoran),
        effect_directional = obs_gmoran - mean(rand_gmoran)
      )

      rows[[length(rows) + 1]] <- data.frame(
        pathway = pathway_name, pct_random = 100 * frac_rand,
        method = paste0("Join count (", null_label, ")"), ics_bins = b,
        null_type = null_label,
        observed = obs_jc, null_mean = mean(rand_jc), null_sd = sd(rand_jc),
        p_value = mean(rand_jc >= obs_jc),
        effect_directional = obs_jc - mean(rand_jc)
      )
    }

    # ---- Optional uniform null sensitivity ----
    if (run_uniform_null) {
      rand_gmoran_u <- numeric(n_iterations)
      rand_jc_u     <- numeric(n_iterations)
      rand_ics_u    <- numeric(n_iterations)

      for (i in seq_len(n_iterations)) {
        u     <- sample(nodes, n_total, replace = replace_null)
        idx_u <- match(u, nodes)
        x_u   <- integer(length(nodes)); x_u[idx_u] <- 1L

        rand_ics_u[i]    <- ics_norm(induced_subgraph(g, vids = u))
        rand_gmoran_u[i] <- global_moran_I(A, x_u)
        rand_jc_u[i]     <- fast_joincount(A, x_u)
      }

      rows[[length(rows) + 1]] <- data.frame(
        pathway = pathway_name, pct_random = 100 * frac_rand,
        method = "ICS (norm; uniform null)", ics_bins = NA_integer_,
        null_type = "uniform_sensitivity",
        observed = obs_ics, null_mean = mean(rand_ics_u), null_sd = sd(rand_ics_u),
        p_value = mean(rand_ics_u <= obs_ics),
        effect_directional = mean(rand_ics_u) - obs_ics
      )
      rows[[length(rows) + 1]] <- data.frame(
        pathway = pathway_name, pct_random = 100 * frac_rand,
        method = "Global Moran's I (uniform sensitivity)", ics_bins = NA_integer_,
        null_type = "uniform_sensitivity",
        observed = obs_gmoran, null_mean = mean(rand_gmoran_u), null_sd = sd(rand_gmoran_u),
        p_value = mean(rand_gmoran_u >= obs_gmoran),
        effect_directional = obs_gmoran - mean(rand_gmoran_u)
      )
      rows[[length(rows) + 1]] <- data.frame(
        pathway = pathway_name, pct_random = 100 * frac_rand,
        method = "Join count (uniform sensitivity)", ics_bins = NA_integer_,
        null_type = "uniform_sensitivity",
        observed = obs_jc, null_mean = mean(rand_jc_u), null_sd = sd(rand_jc_u),
        p_value = mean(rand_jc_u >= obs_jc),
        effect_directional = obs_jc - mean(rand_jc_u)
      )
    }
  }

  bind_rows(rows)
}


false_positive_compare_ics_bins <- function(
    graph,
    bin_cache,
    n_total      = 100,
    ics_bins_set = c(1, 2, 4, 6, 10),
    n_outer      = 500,
    n_iterations = 500,
    alpha        = 0.05,
    seed         = 123,
    replace_null = TRUE,
    n_cores      = parallel::detectCores() - 1,
    return_pvals = FALSE,
    run_uniform_null = FALSE
) {
  suppressPackageStartupMessages({
    library(foreach)
    library(doParallel)
  })

  prep    <- prep_graph_fast(graph)
  g       <- prep$g
  A       <- prep$A
  nodes   <- prep$nodes
  Bset    <- sort(unique(ics_bins_set))
  deg_all <- degree(g)

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)

  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(igraph)
      library(Matrix)
      library(dplyr)
    })
    NULL
  })

  parallel::clusterExport(
    cl,
    varlist = c(
      "g", "A", "nodes", "deg_all",
      "n_total", "n_iterations", "Bset", "seed", "replace_null", "bin_cache",
      "run_uniform_null",
      "ics_norm", "fast_joincount",
      "global_moran_I", "sample_degree_matched"
    ),
    envir = environment()
  )

  res_list <- foreach(
    k         = seq_len(n_outer),
    .combine  = dplyr::bind_rows,
    .packages = c("igraph", "Matrix", "dplyr")
  ) %dopar% {

    set.seed(seed + k)

    target_nodes <- sample(nodes, n_total)
    idx_target   <- match(target_nodes, nodes)
    x_obs        <- integer(length(nodes)); x_obs[idx_target] <- 1L

    mean_deg   <- mean(deg_all[target_nodes])
    median_deg <- median(deg_all[target_nodes])

    obs_gmoran    <- global_moran_I(A, x_obs)
    obs_joincount <- fast_joincount(A, x_obs)
    obs_ics       <- ics_norm(induced_subgraph(g, vids = target_nodes))

    # --- degree-matched null per bin scheme (primary) ---
    p_ics_row    <- setNames(numeric(length(Bset)), paste0("bins_", Bset))
    p_gmoran_row <- setNames(numeric(length(Bset)), paste0("gmoran_bins_", Bset))
    p_jc_row     <- setNames(numeric(length(Bset)), paste0("jc_bins_", Bset))

    for (b in Bset) {
      rand_ics    <- numeric(n_iterations)
      rand_gmoran <- numeric(n_iterations)
      rand_jc     <- numeric(n_iterations)

      for (i in seq_len(n_iterations)) {
        s_dm <- if (b > 1) {
          sample_degree_matched(target_nodes, bin_cache,
                                num_bins = b, replace = replace_null)
        } else {
          sample(nodes, n_total, replace = replace_null)
        }
        idx_s <- match(s_dm, nodes)
        x_s   <- integer(length(nodes)); x_s[idx_s] <- 1L

        rand_ics[i]    <- ics_norm(induced_subgraph(g, vids = s_dm))
        rand_gmoran[i] <- global_moran_I(A, x_s)
        rand_jc[i]     <- fast_joincount(A, x_s)
      }

      p_ics_row[paste0("bins_", b)]       <- mean(rand_ics    <= obs_ics)
      p_gmoran_row[paste0("gmoran_bins_", b)] <- mean(rand_gmoran >= obs_gmoran)
      p_jc_row[paste0("jc_bins_", b)]     <- mean(rand_jc     >= obs_joincount)
    }

    # --- optional uniform null ---
    p_gmoran_u <- NA_real_
    p_jc_u     <- NA_real_
    p_ics_u    <- NA_real_

    if (run_uniform_null) {
      rand_gm_u <- numeric(n_iterations)
      rand_j_u  <- numeric(n_iterations)
      rand_i_u  <- numeric(n_iterations)

      for (i in seq_len(n_iterations)) {
        u     <- sample(nodes, n_total, replace = replace_null)
        idx_u <- match(u, nodes)
        x_u   <- integer(length(nodes)); x_u[idx_u] <- 1L

        rand_gm_u[i] <- global_moran_I(A, x_u)
        rand_j_u[i]  <- fast_joincount(A, x_u)
        rand_i_u[i]  <- ics_norm(induced_subgraph(g, vids = u))
      }

      p_gmoran_u <- mean(rand_gm_u >= obs_gmoran)
      p_jc_u     <- mean(rand_j_u  >= obs_joincount)
      p_ics_u    <- mean(rand_i_u  <= obs_ics)
    }

    data.frame(
      k             = k,
      mean_degree   = mean_deg,
      median_degree = median_deg,
      p_gmoran_uniform = p_gmoran_u,
      p_jc_uniform     = p_jc_u,
      p_ics_uniform    = p_ics_u,
      as.data.frame(t(p_ics_row)),
      as.data.frame(t(p_gmoran_row)),
      as.data.frame(t(p_jc_row)),
      stringsAsFactors = FALSE
    )
  }

  # --- summarize false positive rates ---
  ics_cols    <- paste0("bins_", Bset)
  gmoran_cols <- paste0("gmoran_bins_", Bset)
  jc_cols     <- paste0("jc_bins_", Bset)

  fpr_vals <- c(
    sapply(ics_cols,    function(nm) setNames(mean(res_list[[nm]] < alpha), paste0("ICS_", nm))),
    sapply(gmoran_cols, function(nm) setNames(mean(res_list[[nm]] < alpha), paste0("GMoran_", nm))),
    sapply(jc_cols,     function(nm) setNames(mean(res_list[[nm]] < alpha), paste0("JC_", nm)))
  )

  if (run_uniform_null) {
    fpr_vals <- c(
      fpr_vals,
      GMoran_uniform = mean(res_list$p_gmoran_uniform < alpha, na.rm = TRUE),
      JC_uniform     = mean(res_list$p_jc_uniform     < alpha, na.rm = TRUE),
      ICS_uniform    = mean(res_list$p_ics_uniform    < alpha, na.rm = TRUE)
    )
  }

  out_summary <- data.frame(
    method              = names(fpr_vals),
    alpha               = alpha,
    n_outer             = n_outer,
    n_iterations        = n_iterations,
    false_positive_rate = as.numeric(fpr_vals),
    row.names           = NULL
  )

  if (return_pvals) {
    return(list(
      summary = out_summary,
      pvals   = res_list,
      degree  = res_list[, c("k", "mean_degree", "median_degree")]
    ))
  } else {
    return(out_summary)
  }
}
