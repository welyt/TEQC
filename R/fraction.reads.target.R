fraction.reads.target <-
function(reads, targets, Offset=0, mappingReads=FALSE){

  # in case 'reads' is output of 'reads2pairs' and contains also 'singleReads'
  if(is.list(reads) & ("readpairs" %in% names(reads)))
    reads <- reads$readpairs

  # add Offset if required
  targets <- offsetfun(Offset=Offset, targets=targets)

  # overlaps between reads and targets (of at least 1 base)
  OL <- overlapsAny(reads, targets)

  # fraction of reads / read pairs mapping to targets
  res <- sum(sum(OL)) / nrow(reads)

  # return reads that map to targets
  if(mappingReads){
    mr <- reads[OL,]
    res <- c(onTargetFraction=res, list(mappingReads=mr))
  }
  return(res)
}
