# R/gwas_mapping.R
# Dependencies: dplyr, stringr, tidyr, igraph

#' Build combined UniProt-to-PARTICIPANT reference from BioPAX node tables
#'
#' Extracts UniProt accessions from the UNIFICATION_XREF field and pairs them
#' with the PARTICIPANT identifier used as network node names.
#'
#' @param node_tables A list of node data.frames from readSifnx()$nodes
#' @return A data.frame with columns: PARTICIPANT, Uniprot
build_uniprot_ref <- function(node_tables) {
  ref_list <- lapply(node_tables, function(nodes) {
    ref <- nodes[, c("PARTICIPANT", "UNIFICATION_XREF")]
    ref <- ref %>% tidyr::separate_rows(UNIFICATION_XREF, sep = ";")
    ref$Uniprot <- stringr::str_extract(ref$UNIFICATION_XREF, "(?<=knowledgebase:)[A-Z0-9]+")
    ref[!is.na(ref$Uniprot), c("PARTICIPANT", "Uniprot")]
  })

  combined <- do.call(rbind, ref_list)
  combined[!duplicated(combined), ]
}

#' Map UniProt IDs to network node names
#'
#' Goes directly from UniProt → PARTICIPANT → network vertex, without any
#' intermediate HGNC translation step.
#'
#' @param uniprot_ids Character vector of UniProt accessions (your GWAS gene list)
#' @param uniprot_ref Data.frame from build_uniprot_ref()
#' @param graph igraph object (the biological network)
#' @return Character vector of PARTICIPANT names present in the network
map_uniprot_to_network <- function(uniprot_ids, uniprot_ref, graph) {
  unique_ids <- unique(uniprot_ids)

  # Match UniProt IDs to PARTICIPANT names via the reference table
  matched_idx <- match(unique_ids, uniprot_ref$Uniprot)
  participant_names <- uniprot_ref$PARTICIPANT[na.omit(matched_idx)]

  # Filter to those actually present in the network
  on_network <- participant_names[participant_names %in% V(graph)$name]

  # Report mapping stats
  cat("Input UniProt IDs:", length(unique_ids), "\n")
  cat("Matched to PARTICIPANT:", length(participant_names), "\n")
  cat("Present on network:", length(on_network), "\n")
  cat("Lost at UniProt match:", length(unique_ids) - length(na.omit(matched_idx)), "\n")
  cat("Lost at network match:", length(participant_names) - length(on_network), "\n")

  unique(on_network)
}
