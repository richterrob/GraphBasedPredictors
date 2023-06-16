
#' Observed Indegree Predictor (OIP)
#'
#' Computing the OIP
#'
#' @param observations A binary edge matrix with observed and unobserved (filled with \code{NA}) rows.
#' @param interventions A vector with the indices of the observed rows.
#' @param verbose Outputs the iteration number. Default is \code{FALSE}.
#' @param transitivity.assumed A logical. If \code{TRUE} the T-OIP is computed. Default is \code{FALSE}.
#' @param transitivity.preprocessed A logical. If \code{TRUE} partial observations are transitively
#' closed before computing the OIP, or, T-OIP. Default is \code{FALSE}.
#' @return Prediction matrix given by the OIP, T-OIP, P-OIP or P-T-OIP (entries withing 0 and 1).
#'
#' @export


oip = function(observations,
               interventions,
               verbose = FALSE,
               transitivity.assumed = FALSE,
               transitivity.preprocessed = FALSE){
  G.A.po = observations
  V.i = interventions
  if(transitivity.preprocessed){
    G.A.po.old = G.A.po
    G.A.po[which(!(1:nrow(G.A.po) %in% V.i)),] = 0
    G.A.po = trans_closure_graph(G = G.A.po, V.i = V.i)
    G.A.po = G.A.po$ancestral.graph
    G.A.po[which(!(1:nrow(G.A.po) %in% V.i)),] = NA
  }

  nodes = base::nrow(G.A.po)
  # Number of interventions
  Int = base::length(interventions)
  P = G.A.po
  # Set of non-interventions
  not.interventions = base::setdiff(1:nodes, interventions)

  for(l in 1:nodes){
    # If l itself is in interventions do not count this intervention in the denominator
    one.less = if(l %in% interventions){1} else{0}
    add.one = if(transitivity.assumed){1} else{0}
    # Computing the OIP
    P[not.interventions,l] = (base::sum(G.A.po[base::setdiff(interventions, l),l]) + add.one)/
      (Int - one.less + add.one)
    # Set the diagonal to 1
    P[l,l] = 1
    if(verbose){
      print(paste("Computed incoming edges for node ", l, " of ", nodes, sep = ""))
    }
  }
  if(transitivity.assumed){
    # Compute logically deduced zeros
    G.list.temp = generate_rnd_graph(nodes = nodes,
                                     sparsity = 0.5)
    G.temp = G.list.temp$graph
    G.list.temp = impossible_edges(graph = G.temp,
                                   observations = G.A.po,
                                   interventions = interventions)
    enforced.zeros = G.list.temp$enforced.zeros
    # Replace the logically deduced zeros in P
    replace = base::which(enforced.zeros == 0)
    g = as.vector(P)
    g[replace] = 0
    pred = base::matrix(g, ncol = nodes, nrow = nodes)
  } else{
    pred = P
  }
  return(pred)
}
