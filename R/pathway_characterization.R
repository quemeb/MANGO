# R/pathway_characterization.R
# Dependencies: igraph, dplyr

#' Characterize pathway gene connectivity on the integrated network
#'
#' @param graph The full undirected network
#' @param pathway_data The pathway_stats data.frame (with Genes list-column)
#' @param pathway_name Name of the pathway
#' @return A list with connectivity summary
characterize_pathway <- function(graph, pathway_data, pathway_name) {
  row <- pathway_data %>% filter(PATHWAY_NAMES == pathway_name)
  if (nrow(row) == 0) stop("Pathway not found: ", pathway_name)

  all_nodes <- V(graph)$name
  pw_genes <- intersect(unique(unlist(row$Genes)), all_nodes)

  if (length(pw_genes) == 0) return(NULL)

  # Induced subgraph of pathway genes
  subg <- induced_subgraph(graph, vids = pw_genes)
  comps <- components(subg, mode = "weak")

  # Subgraph degree
  sub_deg <- degree(subg)

  # Full network degree of pathway genes vs all genes
  full_deg_pw  <- degree(graph, v = pw_genes)
  full_deg_all <- degree(graph)

  list(
    pathway_name        = pathway_name,
    n_genes_on_network  = length(pw_genes),
    n_components        = comps$no,
    largest_component   = max(comps$csize),
    component_sizes     = sort(comps$csize, decreasing = TRUE),
    n_isolated          = sum(sub_deg == 0),
    frac_isolated       = mean(sub_deg == 0),
    n_edges_subgraph    = ecount(subg),
    mean_subgraph_deg   = mean(sub_deg),
    mean_full_deg_pw    = mean(full_deg_pw),
    mean_full_deg_all   = mean(full_deg_all),
    subgraph_degrees    = sub_deg,
    full_degrees_pw     = full_deg_pw
  )
}

#' Characterize all pathways in batch
#'
#' @param graph The full undirected network
#' @param pathway_data The pathway_stats data.frame
#' @param min_genes Minimum genes on network to include (default 10)
#' @return data.frame summary
characterize_all_pathways <- function(graph, pathway_data, min_genes = 10) {
  pw_names <- unique(pathway_data$PATHWAY_NAMES)

  results <- lapply(pw_names, function(pw) {
    tryCatch({
      ch <- characterize_pathway(graph, pathway_data, pw)
      if (is.null(ch) || ch$n_genes_on_network < min_genes) return(NULL)
      data.frame(
        pathway           = ch$pathway_name,
        n_genes           = ch$n_genes_on_network,
        n_components      = ch$n_components,
        largest_component = ch$largest_component,
        n_isolated        = ch$n_isolated,
        frac_isolated     = ch$frac_isolated,
        n_edges           = ch$n_edges_subgraph,
        mean_subgraph_deg = round(ch$mean_subgraph_deg, 2),
        mean_full_deg     = round(ch$mean_full_deg_pw, 2),
        stringsAsFactors  = FALSE
      )
    }, error = function(e) NULL)
  })

  do.call(rbind, Filter(Negate(is.null), results))
}
