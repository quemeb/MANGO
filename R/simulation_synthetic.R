# R/simulation_synthetic.R
# Dependencies: igraph, dplyr
# Requires: functions from stats.R, sampling.R

#' Plant a single BFS cluster on the network
#'
#' @param graph Undirected igraph
#' @param seed_node Starting node (character). If NULL, picks a random node.
#' @param cluster_size Number of genes in the cluster
#' @return Character vector of genes in the planted cluster
plant_bfs_cluster <- function(graph, seed_node = NULL, cluster_size = 50) {
  nodes <- V(graph)$name
  if (is.null(seed_node)) seed_node <- sample(nodes, 1)

  # BFS from seed
  bfs_result <- bfs(graph, root = seed_node)
  bfs_order <- V(graph)$name[bfs_result$order]
  bfs_order <- bfs_order[!is.na(bfs_order)]

  if (length(bfs_order) < cluster_size) {
    warning("BFS neighborhood smaller than cluster_size. Using all ",
            length(bfs_order), " reachable nodes.")
    return(bfs_order)
  }

  bfs_order[1:cluster_size]
}

#' Plant multiple dispersed BFS clusters
#'
#' Seeds are chosen to be far apart (by shortest path distance).
#'
#' @param graph Undirected igraph
#' @param n_clusters Number of clusters to plant
#' @param cluster_size Size of each cluster
#' @param min_separation Minimum shortest-path distance between seeds
#' @return List with: signal_genes (character vector), cluster_assignments (named vector)
plant_dispersed_clusters <- function(graph, n_clusters = 3, cluster_size = 20,
                                     min_separation = 5) {
  nodes <- V(graph)$name
  seeds <- character(0)

  # Greedy seed selection: pick seeds that are far from all existing seeds
  attempts <- 0
  while (length(seeds) < n_clusters && attempts < 1000) {
    candidate <- sample(nodes, 1)
    if (length(seeds) == 0) {
      seeds <- candidate
    } else {
      dists <- distances(graph, v = candidate, to = seeds)
      if (all(dists >= min_separation)) {
        seeds <- c(seeds, candidate)
      }
    }
    attempts <- attempts + 1
  }

  if (length(seeds) < n_clusters) {
    warning("Could only place ", length(seeds), " of ", n_clusters,
            " clusters with min_separation=", min_separation)
  }

  all_signal <- character(0)
  cluster_map <- list()

  for (i in seq_along(seeds)) {
    cluster_genes <- plant_bfs_cluster(graph, seed_node = seeds[i],
                                       cluster_size = cluster_size)
    # Remove overlap with previously planted clusters
    cluster_genes <- setdiff(cluster_genes, all_signal)
    all_signal <- c(all_signal, cluster_genes)
    cluster_map[[paste0("cluster_", i)]] <- cluster_genes
  }

  list(signal_genes = all_signal, clusters = cluster_map, seeds = seeds)
}

#' Select hub genes that are NOT topologically clustered
#'
#' For Scenario 3 (hub-driven apparent clustering validation).
#' Picks high-degree genes that are dispersed across the network.
#'
#' @param graph Undirected igraph
#' @param n_hubs Number of hub genes to select
#' @param degree_percentile Only consider genes above this degree percentile (default 0.9)
#' @return Character vector of hub gene names
select_dispersed_hubs <- function(graph, n_hubs = 50, degree_percentile = 0.9) {
  deg <- degree(graph)
  threshold <- quantile(deg, probs = degree_percentile)
  hub_candidates <- V(graph)$name[deg >= threshold]

  if (length(hub_candidates) <= n_hubs) return(hub_candidates)

  # Random sample of hubs (they are spread across the network because we're
  # sampling from the top-degree nodes without any neighborhood constraint)
  sample(hub_candidates, n_hubs)
}

#' Select signal genes with degree distribution matched to background
#'
#' For Scenario 7 (degree-stratified confound validation).
#' Plants a BFS cluster but constrains seed selection to start from a node
#' whose degree is near the network median, and checks that the resulting
#' cluster's degree distribution is not significantly different from the
#' full network's degree distribution.
#'
#' @param graph Undirected igraph
#' @param cluster_size Target cluster size
#' @param max_attempts Number of random seeds to try
#' @return List with signal_genes and degree comparison stats
plant_degree_matched_cluster <- function(graph, cluster_size = 50, max_attempts = 100) {
  deg <- degree(graph)
  median_deg <- median(deg)
  nodes <- V(graph)$name

  best_pval <- 0
  best_cluster <- NULL

  for (a in seq_len(max_attempts)) {
    # Pick a seed near median degree
    candidates <- nodes[abs(deg - median_deg) <= 2]
    if (length(candidates) == 0) candidates <- nodes  # fallback
    seed <- sample(candidates, 1)

    cluster <- plant_bfs_cluster(graph, seed_node = seed, cluster_size = cluster_size)
    cluster_degs <- deg[match(cluster, nodes)]

    # KS test: is this cluster's degree distribution different from network?
    ks <- ks.test(cluster_degs, deg)

    if (ks$p.value > best_pval) {
      best_pval <- ks$p.value
      best_cluster <- cluster
    }

    # If p > 0.3, good enough — cluster degree looks like background
    if (ks$p.value > 0.3) break
  }

  list(
    signal_genes = best_cluster,
    ks_pvalue    = best_pval,
    cluster_degs = deg[match(best_cluster, nodes)],
    network_degs = deg
  )
}

#' Run a single synthetic simulation scenario and compute p-values for all methods
#'
#' @param graph Undirected igraph
#' @param signal_genes Character vector of planted signal genes
#' @param n_total Total selected set size (signal + noise)
#' @param n_perm Number of permutations
#' @param bin_cache Degree bin cache
#' @param ics_bins_set Vector of bin counts to test for ICS
#' @param replace_null Sampling with replacement
#' @return data.frame with method, p_value, observed, null_mean, null_sd
run_synthetic_scenario <- function(graph, signal_genes, n_total, n_perm = 500,
                                   bin_cache, ics_bins_set = c(1, 4, 10),
                                   replace_null = TRUE, seed = 42,
                                   run_uniform_null = FALSE) {
  prep <- prep_graph_fast(graph)
  g <- prep$g; A <- prep$A; nodes <- prep$nodes

  signal_genes <- intersect(signal_genes, nodes)
  background <- setdiff(nodes, signal_genes)

  n_signal <- min(length(signal_genes), n_total)
  n_noise <- n_total - n_signal

  set.seed(seed)
  target_nodes <- c(
    sample(signal_genes, n_signal),
    if (n_noise > 0) sample(background, n_noise) else character(0)
  )

  idx_target <- match(target_nodes, nodes)
  x_obs <- integer(length(nodes)); x_obs[idx_target] <- 1L

  obs_ics <- ics_norm(induced_subgraph(g, vids = target_nodes))
  obs_gm  <- global_moran_I(A, x_obs)
  obs_jc  <- fast_joincount(A, x_obs)

  rows <- list()

  # --- degree-matched null per bin scheme (primary) ---
  # All methods computed on the SAME sample per iteration
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

    rows[[length(rows) + 1]] <- data.frame(
      method = paste0("ICS_bins=", b), null_type = null_label,
      observed = obs_ics, null_mean = mean(rand_ics), null_sd = sd(rand_ics),
      p_value = mean(rand_ics <= obs_ics))
    rows[[length(rows) + 1]] <- data.frame(
      method = paste0("GlobalMoran_", null_label), null_type = null_label,
      observed = obs_gm, null_mean = mean(rand_gm), null_sd = sd(rand_gm),
      p_value = mean(rand_gm >= obs_gm))
    rows[[length(rows) + 1]] <- data.frame(
      method = paste0("Joincount_", null_label), null_type = null_label,
      observed = obs_jc, null_mean = mean(rand_jc), null_sd = sd(rand_jc),
      p_value = mean(rand_jc >= obs_jc))
  }

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

    rows[[length(rows) + 1]] <- data.frame(
      method = "GlobalMoran_uniform_sens", null_type = "uniform_sensitivity",
      observed = obs_gm, null_mean = mean(rand_gm_u), null_sd = sd(rand_gm_u),
      p_value = mean(rand_gm_u >= obs_gm))
    rows[[length(rows) + 1]] <- data.frame(
      method = "Joincount_uniform_sens", null_type = "uniform_sensitivity",
      observed = obs_jc, null_mean = mean(rand_jc_u), null_sd = sd(rand_jc_u),
      p_value = mean(rand_jc_u >= obs_jc))
    rows[[length(rows) + 1]] <- data.frame(
      method = "ICS_uniform_sens", null_type = "uniform_sensitivity",
      observed = obs_ics, null_mean = mean(rand_i_u), null_sd = sd(rand_i_u),
      p_value = mean(rand_i_u <= obs_ics))
  }

  do.call(rbind, rows)
}


#' Plant degree-typical dispersed signal genes (GWAS-like pattern)
#'
#' Creates a signal set that mimics real GWAS hits: a few small BFS clusters
#' plus isolated singletons, all with degree distribution indistinguishable
#' from the network background (KS p > ks_alpha).
#'
#' @param graph Undirected igraph
#' @param n_signal Total number of signal genes
#' @param n_small_clusters Number of small BFS clusters to plant
#' @param cluster_size Integer vector of length 2: min and max cluster size (sampled uniformly)
#' @param ks_alpha Reject candidates whose KS p-value is below this threshold
#' @param max_attempts Maximum rejection-sampling attempts
#' @param seed Optional RNG seed
#' @return List with signal_genes, clusters, singletons, ks_pvalue, n_attempts
plant_degree_typical_dispersed <- function(graph, n_signal = 50,
                                           n_small_clusters = 3,
                                           cluster_size = c(5, 8),
                                           ks_alpha = 0.05,
                                           max_attempts = 500,
                                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  nodes <- V(graph)$name
  net_degs <- degree(graph)

  best <- NULL
  best_ks <- 0


  for (a in seq_len(max_attempts)) {
    all_signal <- character(0)
    cluster_list <- list()

    # Plant small BFS clusters with well-separated seeds
    seeds_used <- character(0)
    for (i in seq_len(n_small_clusters)) {
      cs <- sample(cluster_size[1]:cluster_size[2], 1)
      # Pick a seed far from existing seeds
      candidate_ok <- FALSE
      for (try in 1:50) {
        cand <- sample(setdiff(nodes, all_signal), 1)
        if (length(seeds_used) == 0 ||
            all(distances(graph, v = cand, to = seeds_used) >= 3)) {
          candidate_ok <- TRUE
          break
        }
      }
      if (!candidate_ok) cand <- sample(setdiff(nodes, all_signal), 1)

      cl_genes <- plant_bfs_cluster(graph, seed_node = cand, cluster_size = cs)
      cl_genes <- setdiff(cl_genes, all_signal)
      all_signal <- c(all_signal, cl_genes)
      cluster_list[[paste0("cluster_", i)]] <- cl_genes
      seeds_used <- c(seeds_used, cand)
    }

    # Fill remaining slots with isolated singletons
    n_remaining <- n_signal - length(all_signal)
    if (n_remaining > 0) {
      pool <- setdiff(nodes, all_signal)
      singletons <- sample(pool, min(n_remaining, length(pool)))
      all_signal <- c(all_signal, singletons)
    } else {
      singletons <- character(0)
    }

    # KS test: degree-typicality check
    signal_degs <- net_degs[match(all_signal, nodes)]
    ks_p <- ks.test(signal_degs, net_degs)$p.value

    if (ks_p > best_ks) {
      best_ks <- ks_p
      best <- list(
        signal_genes = all_signal,
        clusters     = cluster_list,
        singletons   = singletons,
        ks_pvalue    = ks_p,
        n_attempts   = a
      )
    }

    if (ks_p > ks_alpha) break
  }

  if (best_ks <= ks_alpha) {
    warning("Could not achieve KS p > ", ks_alpha, " after ", max_attempts,
            " attempts. Best KS p = ", round(best_ks, 4))
  }

  best
}


#' Select degree-typical dispersed genes (no planted clustering)
#'
#' Null control for S3b: random genes whose degree distribution matches the
#' network background (KS p > ks_alpha). No BFS clusters planted.
#'
#' @param graph Undirected igraph
#' @param n_genes Number of genes to select
#' @param ks_alpha Reject candidates whose KS p-value is below this threshold
#' @param max_attempts Maximum rejection-sampling attempts
#' @param seed Optional RNG seed
#' @return List with genes, ks_pvalue, n_attempts
select_degree_typical_dispersed <- function(graph, n_genes = 50,
                                            ks_alpha = 0.05,
                                            max_attempts = 500,
                                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  nodes <- V(graph)$name
  net_degs <- degree(graph)

  best <- NULL
  best_ks <- 0

  for (a in seq_len(max_attempts)) {
    selected <- sample(nodes, n_genes)
    sel_degs <- net_degs[match(selected, nodes)]
    ks_p <- ks.test(sel_degs, net_degs)$p.value

    if (ks_p > best_ks) {
      best_ks <- ks_p
      best <- list(
        genes      = selected,
        ks_pvalue  = ks_p,
        n_attempts = a
      )
    }

    if (ks_p > ks_alpha) break
  }

  if (best_ks <= ks_alpha) {
    warning("Could not achieve KS p > ", ks_alpha, " after ", max_attempts,
            " attempts. Best KS p = ", round(best_ks, 4))
  }

  best
}
