# R/loo.R
# Dependencies: igraph
# Requires: ics_norm() from stats.R

#' Component-level leave-one-out for ICS
#'
#' For each connected component of the GWAS subgraph, remove it and recompute ICS.
#' Returns a data.frame with component ID, size, and delta.
#'
#' @param graph The full undirected network (igraph object)
#' @param target_nodes Character vector of GWAS gene names
#' @return data.frame with columns: comp_id, comp_size, comp_genes (list-column),
#'         ics_full, ics_without, delta (= ics_full - ics_without)
loo_component <- function(graph, target_nodes) {
  target_nodes <- intersect(target_nodes, V(graph)$name)
  subg <- induced_subgraph(graph, vids = target_nodes)

  ics_full <- ics_norm(subg)

  comps <- components(subg, mode = "weak")
  comp_membership <- comps$membership
  comp_ids <- sort(unique(comp_membership))

  results <- lapply(comp_ids, function(cid) {
    comp_genes <- V(subg)$name[comp_membership == cid]
    remaining <- setdiff(target_nodes, comp_genes)

    if (length(remaining) == 0) {
      ics_without <- NA_real_
    } else {
      subg_without <- induced_subgraph(graph, vids = remaining)
      ics_without <- ics_norm(subg_without)
    }

    data.frame(
      comp_id     = cid,
      comp_size   = length(comp_genes),
      ics_full    = ics_full,
      ics_without = ics_without,
      delta       = ics_full - ics_without,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, results)
  # Attach gene lists as a list-column
  out$comp_genes <- lapply(comp_ids, function(cid) {
    V(subg)$name[comp_membership == cid]
  })

  out
}


#' Gene-level leave-one-out for ICS with articulation point optimization
#'
#' For each gene in the target set, compute ICS after removing that gene.
#' Uses Tarjan's articulation point detection to avoid redundant recomputation:
#' - Non-articulation-point removal: component count unchanged, only size changes by 1.
#'   ICS can be updated analytically.
#' - Articulation point removal: must recompute (component may split).
#'
#' @param graph The full undirected network (igraph object)
#' @param target_nodes Character vector of GWAS gene names
#' @return data.frame with columns: gene, subgraph_degree, is_articulation,
#'         ics_full, ics_without, delta, comp_id
loo_gene <- function(graph, target_nodes) {
  target_nodes <- intersect(target_nodes, V(graph)$name)
  subg <- induced_subgraph(graph, vids = target_nodes)
  n <- vcount(subg)

  ics_full <- ics_norm(subg)

  # Component membership
  comps <- components(subg, mode = "weak")
  comp_membership <- comps$membership
  names(comp_membership) <- V(subg)$name
  comp_sizes <- comps$csize

  # Subgraph degree (this is the "importance" dimension)
  sub_deg <- degree(subg)
  names(sub_deg) <- V(subg)$name

  # Articulation points
  art_points <- articulation_points(subg)
  art_names <- V(subg)$name[art_points]
  is_art <- V(subg)$name %in% art_names
  names(is_art) <- V(subg)$name

  # For non-articulation points, ICS update is analytical:
  # Removing a non-AP gene from component c of size s:
  #   - Component count stays the same
  #   - That component's size goes from s to s-1
  #   - ICS_raw changes from sum(1/c_j) to sum(1/c_j) - 1/s + 1/(s-1)
  #   - ICS_norm divides by (n-1) instead of n
  #
  # For articulation points, we must recompute via induced_subgraph.

  ics_raw_full <- sum(1 / comp_sizes)

  results <- lapply(V(subg)$name, function(gene) {
    cid <- comp_membership[gene]
    s <- comp_sizes[cid]

    if (!is_art[gene]) {
      # Analytical update
      if (s == 1) {
        # Isolated node — removing it removes a component
        new_raw <- ics_raw_full - 1  # remove the 1/1 term
      } else {
        new_raw <- ics_raw_full - (1/s) + (1/(s - 1))
      }
      ics_without <- new_raw / (n - 1)
    } else {
      # Must recompute
      remaining <- setdiff(target_nodes, gene)
      if (length(remaining) == 0) {
        ics_without <- NA_real_
      } else {
        subg_without <- induced_subgraph(graph, vids = remaining)
        ics_without <- ics_norm(subg_without)
      }
    }

    data.frame(
      gene            = gene,
      subgraph_degree = sub_deg[gene],
      is_articulation = is_art[gene],
      comp_id         = cid,
      comp_size       = s,
      ics_full        = ics_full,
      ics_without     = ics_without,
      delta           = ics_full - ics_without,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}


#' Brute-force gene-level LOO (validation reference)
#'
#' For each gene, removes it, induces the subgraph from scratch, and
#' recomputes ics_norm. No analytical shortcuts. Use to validate loo_gene().
#'
#' @param graph The full undirected network (igraph object)
#' @param target_nodes Character vector of gene names
#' @return data.frame with columns: gene, ics_full, ics_without, delta
loo_gene_naive <- function(graph, target_nodes) {
  target_nodes <- intersect(target_nodes, V(graph)$name)
  subg <- induced_subgraph(graph, vids = target_nodes)
  ics_full <- ics_norm(subg)

  results <- lapply(target_nodes, function(gene) {
    remaining <- setdiff(target_nodes, gene)
    if (length(remaining) == 0) {
      ics_without <- NA_real_
    } else {
      subg_without <- induced_subgraph(graph, vids = remaining)
      ics_without <- ics_norm(subg_without)
    }
    data.frame(
      gene        = gene,
      ics_full    = ics_full,
      ics_without = ics_without,
      delta       = ics_full - ics_without,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}
