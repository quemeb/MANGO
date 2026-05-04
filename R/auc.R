# R/auc.R
# Dependencies: igraph, Matrix, dplyr, foreach, doParallel
# Requires: prep_graph_fast(), ics_norm(), global_moran_I(),
#           fast_joincount() (from stats.R)
#           sample_degree_matched() (from sampling.R)

roc_auc_from_p <- function(df, score_col = "score", label_col = "label") {
  df <- df %>% filter(is.finite(.data[[score_col]]), !is.na(.data[[label_col]]))

  # sort by decreasing score (strongest first)
  df <- df %>% arrange(desc(.data[[score_col]]))

  n_pos <- sum(df[[label_col]] == 1)
  n_neg <- sum(df[[label_col]] == 0)
  stopifnot(n_pos > 0, n_neg > 0)

  TP <- cumsum(df[[label_col]] == 1)
  FP <- cumsum(df[[label_col]] == 0)

  roc <- data.frame(
    FPR = FP / n_neg,
    TPR = TP / n_pos
  )
  roc <- rbind(data.frame(FPR = 0, TPR = 0), roc, data.frame(FPR = 1, TPR = 1))

  auc <- sum(diff(roc$FPR) * (head(roc$TPR, -1) + tail(roc$TPR, -1)) / 2)
  list(roc = roc, auc = auc)
}


run_single_pathway_auc_parallel <- function(
  pathway_name,
  graph,
  pathway_data,
  bin_cache,
  n_total = 100,
  reps_per_level = 25,
  signal_noise_levels = c(0, 20, 40),
  null_noise_level = 100,
  n_perm = 200,
  ics_bins_set = c(1, 4, 10),
  seed = 1,
  replace_null = TRUE,
  n_cores = 7,
  run_uniform_null = FALSE
) {
  prep  <- prep_graph_fast(graph)
  g     <- prep$g
  A     <- prep$A
  nodes <- prep$nodes

  # pathway genes
  row <- pathway_data %>% filter(PATHWAY_NAMES == pathway_name)
  if (nrow(row) == 0) stop("Pathway not found: ", pathway_name)
  pathway_genes <- intersect(unique(unlist(row$Genes)), nodes)
  if (length(pathway_genes) == 0) stop("No pathway genes in graph: ", pathway_name)
  background_genes <- setdiff(nodes, pathway_genes)

  # tasks
  task_df <- bind_rows(
    tibble(pct_noise = signal_noise_levels, label = 1L),
    tibble(pct_noise = null_noise_level,    label = 0L)
  ) %>%
    slice(rep(1:n(), each = reps_per_level)) %>%
    group_by(pct_noise, label) %>%
    mutate(rep = row_number()) %>%
    ungroup() %>%
    mutate(task_id = row_number())

  # cluster
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)

  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(igraph); library(Matrix)
    })
    NULL
  })

  parallel::clusterExport(
    cl,
    varlist = c(
      "g", "A", "nodes", "pathway_genes", "background_genes",
      "n_total", "n_perm", "ics_bins_set", "seed", "replace_null",
      "run_uniform_null",
      "fast_joincount", "global_moran_I", "ics_norm",
      "sample_degree_matched", "bin_cache"
    ),
    envir = environment()
  )

  res_long <- foreach(
    ii = seq_len(nrow(task_df)),
    .combine = dplyr::bind_rows,
    .packages = c("igraph", "Matrix", "dplyr")
  ) %dopar% {

    pct_noise <- task_df$pct_noise[ii]
    label     <- task_df$label[ii]
    rep_id    <- task_df$rep[ii]

    set.seed(seed + ii)

    frac   <- pct_noise / 100
    n_rand <- as.integer(round(n_total * frac))
    n_path <- n_total - n_rand

    target_nodes <- c(
      if (n_path > 0) sample(pathway_genes,    n_path) else character(0),
      if (n_rand > 0) sample(background_genes, n_rand) else character(0)
    )

    idx_target <- match(target_nodes, nodes)

    # observed:
    x_obs   <- integer(length(nodes)); x_obs[idx_target] <- 1L
    obs_gm  <- global_moran_I(A, x_obs)
    obs_jc  <- fast_joincount(A, x_obs)
    obs_ics <- ics_norm(induced_subgraph(g, vids = target_nodes))

    # --- degree-matched null per bin scheme (primary) ---
    all_out <- list()

    for (b in ics_bins_set) {
      rand_ics <- numeric(n_perm)
      rand_gm  <- numeric(n_perm)
      rand_jc  <- numeric(n_perm)

      for (p in seq_len(n_perm)) {
        sset <- if (b > 1) {
          sample_degree_matched(target_nodes, bin_cache, num_bins = b, replace = replace_null)
        } else {
          sample(nodes, n_total, replace = replace_null)
        }
        idx_s <- match(sset, nodes)
        x_s   <- integer(length(nodes)); x_s[idx_s] <- 1L

        rand_ics[p] <- ics_norm(induced_subgraph(g, vids = sset))
        rand_gm[p]  <- global_moran_I(A, x_s)
        rand_jc[p]  <- fast_joincount(A, x_s)
      }

      null_label <- if (b > 1) paste0("dm_bins=", b) else "uniform"

      all_out[[length(all_out) + 1]] <- c(
        setNames(mean(rand_ics <= obs_ics),  paste0("ICS_bins=", b)),
        setNames(mean(rand_gm >= obs_gm),    paste0("GlobalMoran_", null_label)),
        setNames(mean(rand_jc >= obs_jc),    paste0("Joincount_", null_label))
      )
    }

    out <- unlist(all_out)

    # --- optional uniform null sensitivity ---
    if (run_uniform_null) {
      rand_gm_u <- numeric(n_perm)
      rand_jc_u <- numeric(n_perm)
      rand_i_u  <- numeric(n_perm)

      for (p in seq_len(n_perm)) {
        u     <- sample(nodes, n_total, replace = replace_null)
        idx_u <- match(u, nodes)
        x_u   <- integer(length(nodes)); x_u[idx_u] <- 1L

        rand_gm_u[p] <- global_moran_I(A, x_u)
        rand_jc_u[p] <- fast_joincount(A, x_u)
        rand_i_u[p]  <- ics_norm(induced_subgraph(g, vids = u))
      }

      out <- c(out,
        GlobalMoran_uniform_sens = mean(rand_gm_u >= obs_gm),
        Joincount_uniform_sens   = mean(rand_jc_u >= obs_jc),
        ICS_uniform_sens         = mean(rand_i_u  <= obs_ics)
      )
    }

    data.frame(
      pathway   = pathway_name,
      pct_noise = pct_noise,
      label     = label,
      rep       = rep_id,
      method    = names(out),
      p_value   = as.numeric(out),
      stringsAsFactors = FALSE
    )
  }

  res_long
}


.get_pathway_gene_list <- function(pathway_data, pathway_name, nodes_in_graph) {
  row <- pathway_data %>% filter(PATHWAY_NAMES == pathway_name)
  if (nrow(row) == 0) return(character(0))
  intersect(unique(unlist(row$Genes)), nodes_in_graph)
}


run_fixed_mixture_auc_parallel <- function(
  graph,
  pathway_data,
  bin_cache,
  fixed_pathways,
  K = 200,
  signal_frac_vec = c(0.05, 0.1, 0.2, 0.4),
  reps = 40,
  n_perm = 200,
  ics_bins_set = c(1, 2, 4, 6, 10),
  seed = 1,
  replace_null = TRUE,
  n_cores = 7,
  run_uniform_null = FALSE
) {
  prep      <- prep_graph_fast(graph)
  g         <- prep$g
  A         <- prep$A
  nodes     <- prep$nodes
  all_genes <- nodes

  fixed_pathways <- intersect(fixed_pathways, unique(pathway_data$PATHWAY_NAMES))
  if (length(fixed_pathways) == 0) stop("fixed_pathways not found in pathway_data.")

  pw_gene_map <- lapply(fixed_pathways, function(pw) {
    row <- pathway_data %>% filter(PATHWAY_NAMES == pw)
    intersect(unique(unlist(row$Genes)), nodes)
  })
  names(pw_gene_map) <- fixed_pathways

  union_pool <- unique(unlist(pw_gene_map))
  union_pool <- intersect(union_pool, all_genes)

  if (length(union_pool) < 20) stop("Union pool too small; pick different pathways.")
  max_s <- max(1L, as.integer(round(K * max(signal_frac_vec))))
  if (length(union_pool) < max_s) {
    stop("Union pool (", length(union_pool),
         ") smaller than max signal needed (", max_s,
         "). Reduce K/signal_frac or choose pathways with more genes.")
  }

  fixed_pw_str <- paste(fixed_pathways, collapse = " | ")

  task_df <- expand.grid(
    rep         = seq_len(reps),
    signal_frac = signal_frac_vec,
    label       = c(1L, 0L),
    stringsAsFactors = FALSE
  ) %>%
    mutate(task_id = row_number())

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)

  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({library(igraph); library(Matrix); library(dplyr)})
    NULL
  })

  parallel::clusterExport(
    cl,
    varlist = c(
      "g", "A", "nodes", "all_genes", "union_pool", "fixed_pw_str",
      "K", "n_perm", "ics_bins_set", "seed", "replace_null", "bin_cache",
      "run_uniform_null",
      "fast_joincount", "global_moran_I", "ics_norm",
      "sample_degree_matched"
    ),
    envir = environment()
  )

  res_long <- foreach(
    ii = seq_len(nrow(task_df)),
    .combine = dplyr::bind_rows,
    .packages = c("igraph", "Matrix", "dplyr")
  ) %dopar% {

    rep_id <- task_df$rep[ii]
    frac   <- task_df$signal_frac[ii]
    label  <- task_df$label[ii]

    set.seed(seed + ii)

    s <- as.integer(round(K * frac))
    if (s < 1) s <- 1
    if (s > K) s <- K

    if (label == 1L) {
      sig_core <- sample(union_pool, s, replace = FALSE)
      fillers  <- if (K - s > 0) sample(setdiff(all_genes, sig_core), K - s, replace = FALSE) else character(0)
      target_nodes <- c(sig_core, fillers)
    } else {
      target_nodes <- sample(all_genes, K, replace = FALSE)
    }

    idx_target <- match(target_nodes, nodes)

    # observed
    x_obs   <- integer(length(nodes)); x_obs[idx_target] <- 1L
    obs_gm  <- global_moran_I(A, x_obs)
    obs_jc  <- fast_joincount(A, x_obs)
    obs_ics <- ics_norm(induced_subgraph(g, vids = target_nodes))

    # --- degree-matched null per bin scheme (primary) ---
    all_out <- list()

    for (b in ics_bins_set) {
      rand_ics <- numeric(n_perm)
      rand_gm  <- numeric(n_perm)
      rand_jc  <- numeric(n_perm)

      for (p in seq_len(n_perm)) {
        sset <- if (b > 1) {
          sample_degree_matched(target_nodes, bin_cache, num_bins = b, replace = replace_null)
        } else {
          sample(all_genes, K, replace = replace_null)
        }
        idx_s <- match(sset, nodes)
        x_s   <- integer(length(nodes)); x_s[idx_s] <- 1L

        rand_ics[p] <- ics_norm(induced_subgraph(g, vids = sset))
        rand_gm[p]  <- global_moran_I(A, x_s)
        rand_jc[p]  <- fast_joincount(A, x_s)
      }

      null_label <- if (b > 1) paste0("dm_bins=", b) else "uniform"

      all_out[[length(all_out) + 1]] <- c(
        setNames(mean(rand_ics <= obs_ics),  paste0("ICS_bins=", b)),
        setNames(mean(rand_gm >= obs_gm),    paste0("GlobalMoran_", null_label)),
        setNames(mean(rand_jc >= obs_jc),    paste0("Joincount_", null_label))
      )
    }

    out <- unlist(all_out)

    # --- optional uniform null sensitivity ---
    if (run_uniform_null) {
      rand_gm_u <- numeric(n_perm)
      rand_jc_u <- numeric(n_perm)
      rand_i_u  <- numeric(n_perm)

      for (p in seq_len(n_perm)) {
        u     <- sample(all_genes, K, replace = replace_null)
        idx_u <- match(u, nodes)
        x_u   <- integer(length(nodes)); x_u[idx_u] <- 1L

        rand_gm_u[p] <- global_moran_I(A, x_u)
        rand_jc_u[p] <- fast_joincount(A, x_u)
        rand_i_u[p]  <- ics_norm(induced_subgraph(g, vids = u))
      }

      out <- c(out,
        GlobalMoran_uniform_sens = mean(rand_gm_u >= obs_gm),
        Joincount_uniform_sens   = mean(rand_jc_u >= obs_jc),
        ICS_uniform_sens         = mean(rand_i_u  <= obs_ics)
      )
    }

    data.frame(
      rep             = rep_id,
      signal_frac     = frac,
      label           = label,
      chosen_pathways = fixed_pw_str,
      method          = names(out),
      p_value         = as.numeric(out),
      stringsAsFactors = FALSE
    )
  }

  res_long %>%
    mutate(score = -log10(pmax(p_value, 1e-300)))
}
