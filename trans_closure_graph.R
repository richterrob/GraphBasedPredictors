
#' Transitive Closure of a Graph
#'
#' Computing the Transitive Closure of a Graph with known Edgematrix-Rows.
#' Note: This code is by far not the most efficient way to compute the transitive closure.
#'
#' @param G edge matrix of the input graph with diagonal equal to 1.
#' @param V.i the indices of known rows. Default is NULL.
#' @param k.closed If NULL the transitive closure is computed, if a number between 1 and number of nodes the corresponding k-reachability graph is computed. Default is NULL.
#' @param reduced.output If \code{TRUE} only the graph is returned. Default is \code{FALSE}.
#' 
#' @return A list with entry \code{$ancestral.graph} for the transitive closure (of k-reachability graph) of \code{G} and \code{$partial.observations}, the same matrix with all non-known rows masked as \code{NA}.
#' @export
#'

trans_closure_graph = function(G,
                               V.i = NULL,
                               k.closed = NULL,
                               reduced.output = FALSE){
  # Number of nodes
  nodes = base::nrow(G)
  G.A = G
  zeros = base::which(G == 0)
  new.zeros = zeros
  
  # Transitive Closure:
  if(base::is.null(k.closed)){
    t= 1
    finished = 0
    while((t <= (nodes-2)) & (finished == 0))  {
      zeros = new.zeros
      G.old = G.A
      G.A = G.A %*% G.A
      t = t + 1
      new.zeros = base::which(G.A == 0)
      if(base::length(zeros) == base::length(new.zeros)){
        finished = 1
      }
    }
    # k-reachability Graph:
  } else{
    # k.closed = 1
    if(k != 1){
      power = base::floor(log2(k.closed))
      rest = k.closed - (2^base::floor(log2(k.closed)))
      for(go in 1:power){
        G.old = G.A
        G.A = G.A %*% G.A
      }
      for(goon in 1:rest){
        G.A = G %*% G.A
      }
    }
  }
  # Getting a binary matrix
  G.A = base::pmin(G.A, 1)
  # Constructing the partial observations matrix
  G.A.po = base::matrix(NA, nrow = nodes, ncol = nodes)
  G.A.po[V.i,] = G.A[V.i,]
  if(reduced.output){
    output = G.A
  } else{
    output = list(ancestral.graph = G.A,
                  partial.observation = G.A.po)
  }
  return(output)
}


