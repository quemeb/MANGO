# R/stats.R
# Dependencies: igraph, Matrix

prep_graph_fast <- function(graph) {
  if (is.null(V(graph)$name)) stop("graph vertices need $name")
  g <- graph
  if (is_directed(g)) g <- as.undirected(g, mode = "collapse")
  A <- as_adj(g, sparse = TRUE)   # dgCMatrix
  diag(A) <- 0                    # neighbor-only
  list(g = g, A = A, nodes = V(g)$name)
}

ics_norm <- function(subg) {
  n <- vcount(subg)
  if (n == 0) return(NA_real_)
  comps <- components(subg, mode = "weak")$csize
  (sum(1 / comps)) / n
}

# Using Moran definition: Ii = z_i * (A z)_i / m2; averaged over target nodes.
fast_moran_mean <- function(A, x, idx_target) {
  n <- length(x)
  xb <- mean(x)
  z <- x - xb
  m2 <- sum(z^2) / n
  if (m2 == 0) return(0)
  lag_z <- as.numeric(A %*% z)
  Ii <- (z * lag_z) / m2
  mean(Ii[idx_target])
}

# Gi* with Astar = A + I (self included), all weights 1.
fast_gistar_mean <- function(A, x, idx_target) {
  n <- length(x)
  xb <- mean(x)
  s <- sqrt((sum(x^2) / n) - xb^2)
  if (!is.finite(s) || s == 0) return(0)

  # w* x  = (A + I) %*% x  = A %*% x + x
  wx <- as.numeric(A %*% x) + x

  # w1 = (A + I) %*% 1  = degree + 1
  deg <- Matrix::rowSums(A)
  w1 <- deg + 1

  # w2 = sum_j w_ij^2 ; with 0/1 weights still equals w1
  w2 <- w1

  denom <- s * sqrt((n * w2 - (w1^2)) / (n - 1))
  Gi <- (wx - xb * w1) / denom
  Gi[!is.finite(Gi)] <- 0
  mean(Gi[idx_target])
}

# Global Moran's I: I = (n / W) * (z' A z) / (z' z)
# where z = x - mean(x), W = sum(A).
# This is the proper global spatial autocorrelation statistic.
global_moran_I <- function(A, x) {
  n <- length(x)
  z <- x - mean(x)
  ztz <- sum(z^2)
  if (ztz == 0) return(0)
  W <- sum(A)
  if (W == 0) return(0)
  ztAz <- as.numeric(Matrix::crossprod(z, A %*% z))
  (n / W) * (ztAz / ztz)
}

# Join count statistic for binary labels on an undirected graph.
# Counts selected-selected edges (1-1 joins) relative to total possible.
# x'Ax / 2 gives the number of edges where both endpoints are 1
# (divide by 2 because A is symmetric and each edge is counted twice)
fast_joincount <- function(A, x) {
  as.numeric(Matrix::crossprod(x, A %*% x)) / 2
}

# Fast version of compute_stats(): uses matrix-based helpers instead of
# local_moran_I() / local_getis_ord_gistar(). Same statistical output.
compute_stats_fast <- function(graph, target_nodes) {
  prep <- prep_graph_fast(graph)
  g <- prep$g; A <- prep$A; nodes <- prep$nodes

  subg <- induced_subgraph(g, vids = target_nodes)
  obs_ics <- ics_norm(subg)

  idx_target <- match(target_nodes, nodes)
  x <- integer(length(nodes)); x[idx_target] <- 1L

  list(
    ics_norm      = obs_ics,
    moran_mean    = fast_moran_mean(A, x, idx_target),
    global_moran  = global_moran_I(A, x),
    gistar_mean   = fast_gistar_mean(A, x, idx_target),
    joincount_11  = fast_joincount(A, x)
  )
}
