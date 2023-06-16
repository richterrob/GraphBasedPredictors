





check_non_trivial = function(ground_truth,
                             observations){
  for(k in 1:nrow(observations)){
    observations[k,k] = 1
  }
  gt.vec = as.vector(ground_truth)
  ob.vec = as.vector(observations)
  interest.vec = gt.vec[is.na(ob.vec)]
  if(min(as.numeric(interest.vec == 0)) == 1){
    return(TRUE)
  } else{
    if(min(as.numeric(interest.vec == 1)) == 1){
      return(TRUE)
    } else{
      return(FALSE)
    }
  }
}
