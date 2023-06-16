

#' Deduced Zeros
#'
#' Computing the zeros that can be deduced from a partial observation (row-wise) of a transitively closed graph
#'
#' @param graph Graph where zeros need to be induced.
#' @param observations A matrix with observed rows and unobserved rows (filled with NA's).
#' @param interventions Indices of the observed rows.
#' @param reduced.output If \code{TRUE} only the graph is returned. Default is \code{FALSE}.
#' 
#' @return A list including the resulting \code{$graph} and \code{$enforced.zeros} the matrix indicating by \code{0} all deduced zeros (other entries are \code{1}).
#'


impossible_edges = function(graph,
                            observations,
                            interventions,
                            reduced.output = FALSE){
  G = graph
  # Number of nodes
  nodes = base::nrow(graph)
  enforced.zeros = base::matrix(1, ncol = nodes, nrow = nodes)
  for(k in interventions){
    causal = base::which(observations[k,] == 1)
    non.causal = base::which(observations[k,] == 0)
    enforced.zeros[causal,non.causal] = 0
  }
  # Replacing the deduced zeros in G
  replace = base::which(enforced.zeros == 0)
  g = base::as.vector(G)
  g[replace] = 0
  G = base::matrix(g, ncol = nodes, nrow = nodes)
  if(reduced.output){
    output = G
  } else{
    output = list(graph = G,
                  enforced.zeros = enforced.zeros)
  }
  return(output)
}
