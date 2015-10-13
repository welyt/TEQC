readsPerTarget <-
function(reads, targets, Offset=0){

  # add Offset if required
  targets <- offsetfun(Offset=Offset, targets=targets)

  # overlaps between reads and targets (of at least 1 base)
  OL <- findOverlaps(targets, reads)

  # numbers of reads per target
  n.reads <- table(queryHits(OL))
  
  targets$nReads <- 0
  targets$nReads[as.numeric(names(n.reads))] <- n.reads
  return(targets)
}

