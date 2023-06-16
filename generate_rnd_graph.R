
#' Drawing a Random Graph with Known Interventions
#'
#' Generates an Erdös-Rènyi graph with given number of nodes
#' and sparsity and draws a random number of known rows. It returns
#' the graph and the index set of known interventions.
#'
#' @param nodes Number of nodes.
#' @param i.nodes Number of known interventions. Default to 0.
#' @param sparsity Probability of an edge between two nodes. Default to 0.5.
#' @param reduced.output If \code{TRUE} only the graph is returned. Default is \code{FALSE}.
#' @return A list with entry \code{$graph} for the random graphs edge matrix (diagonal is set to 1) and with entry \code{$interventions} for the vector of indices for the known interventions.
#'

generate_rnd_graph = function(nodes,
                              i.nodes = 0,
                              sparsity = 0.5,
                              reduced.output = FALSE){

  # Sample a random edge matrix with given sparsity
  rnd = base::sample(c(0,1), size = nodes^2 , replace = TRUE,
                     prob = c(1-sparsity, sparsity))
  R.V = base::matrix(rnd, ncol = nodes, nrow = nodes)
  # Override the diagonal with 1's
  diag(R.V) = base::rep(c(1), times = nodes)
  # Sample at random the known edge matrix rows
  interventions = base::sample(1:nodes, replace = FALSE, size = i.nodes)
  if(reduced.output){
    output = R.V
  } else{
    output = list(graph = R.V,
                  interventions = interventions)
  }
  return(output)
}

