offsetfun <-
function(Offset=0, targets){

  if(Offset > 0){
    ranges(targets) <- ranges(targets) + Offset
    targets <- reduce(targets)
  }
  return(targets)
}

