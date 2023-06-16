
#' Simulate Spanning Trees to Ensure One-constraints
#'
#' Simulating spanning trees rooted in the interventional nodes so that one-constraints are satisfied.
#'
#' @param graph Graph where ones need to be ensured.
#' @param enforced.zeros Indices of the vectorized edge matrix indicating which entries are deduced to be zero.
#' @param observations A matrix with observed rows and unobserved rows (filled with NA's).
#' @param interventions Indices of the observed rows.
#' @param reduced.output If \code{TRUE} only the graph is returned. Default is \code{FALSE}.
#'
#' This function depends on simulate_spanning_tree().
#'
#' @return A list including the resulting \code{$graph} and \code{$enforced.ones} the matrix indicating with a \code{1} where ones are simulated (other entries are \code{0})
#'


enforce_ones = function(graph,
                        enforced.zeros,
                        observations,
                        interventions,
                        reduced.output = FALSE){
  G = graph
  # Number of nodes
  nodes = base::nrow(graph)
  enforced.ones = base::matrix(0, ncol = nodes, nrow = nodes)
  for(k in interventions){
    # Compute Decendants Set
    causal = base::setdiff(base::which(observations[k,] == 1), k)
    if(length(causal) > 0){
      causal.submatrix = enforced.zeros[c(k,causal),c(k,causal)]
      # Simulate the spanning tree
      spanning.tree = simulate_spanning_tree(graph = causal.submatrix)
      # Simulate the spanning tree into the enforced ones matrix
      enforced.ones[c(k,causal),c(k,causal)] = base::pmax(enforced.ones[c(k,causal),c(k,causal)],
                                                          spanning.tree)
    }
  }
  # Replace the simulated spanning trees in G
  replace = base::which(enforced.ones == 1)
  g = base::as.vector(G)
  g[replace] = 1
  G = base::matrix(g, ncol = nodes, nrow = nodes)
  if(reduced.output){
    output = G
  } else{
    output = list(graph = G,
                  enforced.ones = enforced.ones)
  }
  return(output)
}

