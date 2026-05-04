# R/sampling.R
# Dependencies: base R only

build_bin_cache <- function(deg, node_names, max_bins = 10) {

  bin_cache <- list()

  for (b in 1:max_bins) {

    if (b == 1) {
      breaks <- c(min(deg), max(deg))
    } else {
      probs  <- seq(0, 1, length.out = b + 1)
      breaks <- as.numeric(quantile(deg, probs = probs, na.rm = TRUE))
      breaks <- unique(breaks)
    }

    bin_id <- cut(deg, breaks = breaks, include.lowest = TRUE)

    bin_df <- data.frame(
      node   = node_names,
      degree = deg,
      bin    = bin_id
    )

    bin_cache[[paste0("bins_", b)]] <- list(
      breaks      = breaks,
      node_table  = bin_df,
      bin_members = split(node_names, bin_id)
    )
  }

  bin_cache
}

sample_degree_matched <- function(target_nodes, bin_cache, num_bins, replace = TRUE) {
  key <- paste0("bins_", num_bins)
  bins_obj <- bin_cache[[key]]
  if (is.null(bins_obj)) stop("bin_cache missing ", key)

  node_table  <- bins_obj$node_table
  bin_members <- bins_obj$bin_members

  tgt_bins <- node_table$bin[match(target_nodes, node_table$node)]
  if (any(is.na(tgt_bins))) stop("Some target_nodes not found in bin_cache node_table.")

  target_counts <- table(tgt_bins)
  bin_levels    <- names(target_counts)

  unlist(lapply(bin_levels, function(b) {
    pool   <- bin_members[[b]]
    n_need <- as.integer(target_counts[[b]])
    if (!replace && length(pool) < n_need) {
      stop("Bin ", b, " too small for replace=FALSE; use replace=TRUE or fewer bins.")
    }
    sample(pool, size = n_need, replace = replace)
  }), use.names = FALSE)
}
