
#' Transitivity Assuming Predictor (TAP)
#'
#' Simulating the graph-based predictor, as a baseline for causal structure learning, TAP (TAP-q) by rejection sampling as detailed in Algorithm 1 of
#' "Improved baselines for causal structure learning on interventional data".
#'
#' @param observations A binary edge matrix with observed and unobserved (filled with \code{NA}) rows.
#' @param interventions A vector with the indices of the observed rows.
#' @param nbr.simulations Number of consistent simulations. Default is \code{100}.
#' @param q Sparsity used to generate random graphs. Default is \code{0.5}.
#' @param verbose Outputs the iteration number. Default is \code{FALSE}.
#'
#' This function depends on trans_closure_graph(), generate_rnd_graph(), impossible_edges().
#' 
#' @return Prediction matrix approximating the SRB (entries withing 0 and 1).
#' @export

tap = function(observations,
               interventions,
               nbr.simulations = 100,
               q = 0.5,
               verbose = FALSE,
               more.info = FALSE){
  
  # k-reachibility extention is not workable.
  k = NULL
  
  G.A.po = observations
  V.i = interventions
  transitivity.preprocessed = TRUE
  if(transitivity.preprocessed){
    G.A.po.old = G.A.po
    G.A.po[which(!(1:nrow(G.A.po) %in% V.i)),] = 0
    G.A.po = trans_closure_graph(G = G.A.po, V.i = V.i)
    G.A.po = G.A.po$ancestral.graph
    G.A.po[which(!(1:nrow(G.A.po) %in% V.i)),] = NA
  }
  # Number of nodes
  nodes = base::nrow(G.A.po)
  # Counter how many matrices t.c.'s are consistent with the observations
  P = base::matrix(0,ncol = nodes, nrow = nodes)
  # Counter of all matrices whose t.c. are consistent with the observations
  Denom = base::matrix(0,ncol = nodes, nrow = nodes)
  # Compute logically deduced zeros
  G.list = generate_rnd_graph(nodes = nodes,
                              sparsity = q)
  G = G.list$graph
  if(is.null(k)){
    G.list = impossible_edges(graph = G,
                              observations = G.A.po,
                              interventions = interventions)
  } else{
    G.list = enforce_zeros_k_reachability(graph = G,
                                          observations = G.A.po,
                                          interventions = interventions,
                                          k = k)
  }
  G = G.list$graph
  enforced.zeros = G.list$enforced.zeros

  while(nbr.simulations > 0){
    # Making sure that impossible and possible zeros have distinct values in P
    if(nbr.simulations == 1){
      G.list = generate_rnd_graph(nodes = nodes,
                                  sparsity = 1)
    } else{
      # Sample random graph with given sparsity
      G.list = generate_rnd_graph(nodes = nodes,
                                  sparsity = q)
    }
    # Enforce deduced zeros
    G = G.list$graph
    replace = base::which(enforced.zeros == 0)
    g = base::as.vector(G)
    g[replace] = 0
    G = base::matrix(g, ncol = nodes, nrow = nodes)
    # Compute transitive closure (k-reachability graph)
    G.A.list = trans_closure_graph(G,
                                   interventions,
                                   k.closed = k)
    G.A = G.A.list$ancestral.graph
    G.A.newpo = G.A.list$partial.observation
    # Checking consistency
    of.interest = (G.A.newpo == G.A.po)
    if(base::min(of.interest[V.i,]) == 1){
      Denom = Denom + base::matrix(1,ncol = nodes, nrow = nodes)
      P = P + G.A
      nbr.simulations = nbr.simulations - 1
      if(verbose){
        print(paste("Simulations remaining",nbr.simulations, sep = " "))
      }
    }
  }
  # Compute SRB approximation
  pred = P/Denom

  if(more.info){
    output = list(pred = pred,
                  enforced.zeros = enforced.zeros)
    return(output)
  } else{
    return(pred)
  }

}
