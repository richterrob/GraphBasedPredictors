
#' Simulate a spanning tree
#'
#' Simulates a spanning tree in a given directed graph.
#'
#' @param graph Binary graph where a spanning tree needs to be simulated on.
#' @param start.node Root node of the spanning tree.
#' @return spanning tree.
#'

simulate_spanning_tree = function(graph,
                                  start.node = 1){
  # Number of nodes
  nodes = base::nrow(graph)
  set.of.nodes = 1:nodes
  # List for every node of possible parents
  possible.parent = list()
  for(nod in set.of.nodes){
    possible.parent[[nod]] = base::setdiff(base::which(graph[,nod] == 1), nod)
  }
  # In "remaining" the not-yet visited nodes are denoted
  remaining = base::setdiff(set.of.nodes, start.node)
  # In "visited" the visited nodes are denoted
  visited = start.node
  # In "spanning.tree" the edge will be added to compose the spanning tree
  spanning.tree = base::matrix(0, ncol = nodes, nrow = nodes)
  # Start the random walk at the start.node
  random.walk = start.node
  broken = FALSE
  while(!(base::length(remaining) == 0)){
    # Pick a random sample from the set of nodes
    random.seq = base::sample(set.of.nodes, replace = TRUE, size = 100)
    # Paste the root/start node at the start
    random.seq = c(random.walk[length(random.walk)], random.seq)
    for(kau in 1:100){
      # Check if the node has been visited
      if(!(random.seq[kau + 1] %in% visited)){
        # If yes, check if its predecessor is a possible parent
        if(random.seq[kau] %in% possible.parent[[random.seq[kau + 1]]]){
          # If yes, include the node in the set of visited nodes and mark the edge in the spanning tree
          visited = c(visited, random.seq[kau + 1])
          spanning.tree[random.seq[kau], random.seq[kau + 1]] = 1
        } else{
          # If the predecessor is not a possible parent stop the for-loop.
          if(kau < 100){
          random.seq[(kau + 2):101] = base::rep( start.node, times = (100 - kau))
          }
          broken = TRUE
        }
      }
    }
    # If the chain was not broken start from random.seq[101] the next random sample, otherwise go back to the start.node
    if(!broken){
      random.walk = random.seq[101]
    } else{
      random.walk = start.node
    }
    # Compute the remaining nodes
    remaining = base::setdiff(remaining, visited)
  }
  return(spanning.tree)
}
