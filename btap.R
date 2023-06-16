
#' Biased Transitivity Assuming Predictor (B-TAP)
#'
#' Simulating the graph-based predictor, as a baseline for causal structure learning, biased TAP (B-TAP-q) by rejection sampling as detailed in Algorithm 2 of
#' "Improved baselines for causal structure learning on interventional data".
#'
#' @param observations A binary edge matrix with observed and unobserved (filled with \code{NA}) rows.
#' @param interventions A vector with the indices of the observed rows.
#' @param nbr.simulations Number of simulations. Default is \code{100}.
#' @param q Sparsity used to generate random graphs. Default is \code{0.5}.
#' @param verbose Outputs the iteration number. Default is \code{FALSE}.
#'
#' This function depends on trans_closure_graph(), generate_rnd_graph(), impossible_edges() and enforce_ones().
#'
#' @return Prediction matrix approximating the B-TAP (entries withing 0 and 1).
#' @export


btap = function(observations,
                interventions,
                nbr.simulations = 100,
                q = 0.5,
                verbose = FALSE){

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
  # Counter for all simulated matrices
  P = base::matrix(0, ncol = nodes, nrow = nodes)
  # Compute logically deduced zeros
  G.list.temp = generate_rnd_graph(nodes = nodes,
                                   sparsity = q)
  G.temp = G.list.temp$graph
  G.list.temp = impossible_edges(graph = G.temp,
                              observations = G.A.po,
                              interventions = interventions)
  G.temp = G.list.temp$graph
  enforced.zeros = G.list.temp$enforced.zeros
  for(k in 1:nbr.simulations){
    # Making sure that impossible and possible zeros have distinct values in P
    if(k == 1){
      G.list.temp = generate_rnd_graph(nodes = nodes,
                                       sparsity = 1)
    } else{
      # Sample random graph with given sparsity
      G.list.temp = generate_rnd_graph(nodes = nodes,
                                       sparsity = q)
    }
    G.temp = G.list.temp$graph
    # Enforce the given zeros
    replace = base::which(enforced.zeros == 0)
    g = base::as.vector(G.temp)
    g[replace] = 0
    G.temp = base::matrix(g, ncol = nodes, nrow = nodes)
    # Sample logically induced spanning trees to ensure the t.c. includes the given ones
    G.temp.list = enforce_ones(graph = G.temp,
                               enforced.zeros = enforced.zeros,
                               observations = G.A.po,
                               interventions = interventions)
    G.temp = G.temp.list$graph
    # Compute the transitive closure
    G.A.list.temp = trans_closure_graph(G.temp, interventions)
    G.A.temp = G.A.list.temp$ancestral.graph
    G.A.newpo = G.A.list.temp$partial.observation
    # Checking consistency
    of.interest = (G.A.newpo == G.A.po)
    if(base::min(of.interest[V.i,]) == 1){
      if(verbose){
        print(paste("Simulations remaining",(nbr.simulations - k), sep = " "))
      }
      P = P + G.A.temp
    }
  }
  if(transitivity.preprocessed){
    P[V.i,] = G.A.po.old[V.i,]
  }
  # Compute B-SRB approximation
  pred = (1/nbr.simulations)*P
  return(pred)
}

