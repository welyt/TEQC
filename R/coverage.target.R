coverage.target <-
function(reads, targets, Offset=0, perTarget=TRUE, perBase=TRUE){

  # add Offset if required
  targets <- offsetfun(Offset=Offset, targets=targets)

  # get largest position per chromosome - either target or reads position
  max.targ <- end(range(targets))
  chr.targ <- names(max.targ)
  max.read <- end(range(reads))
  noread <- setdiff(chr.targ, names(max.read))
  sel <- setdiff(chr.targ, noread)
  tmp <- max.read[sel] < max.targ[sel]
  max.read[sel][tmp] <- max.targ[sel][tmp]
  
  # coverage for each base
  covercounts.all <- coverage(reads, width=max.read)

  # in case on some chromosome(s) there are targets but no reads, add 0 coverages for that chromosome
  if(length(noread) > 0){
    N <- names(covercounts.all)
    for(chr in noread)
      covercounts.all <- c(covercounts.all, RleList(Rle(0, max.targ[[chr]])))
    names(covercounts.all) <- c(N, noread)
    covercounts.all <- covercounts.all[sort(names(covercounts.all))]
  }

  # restrict to target bases
  covercounts.target <- RleList()
  targetcov <- NULL
  for(chr in names(covercounts.all)){
    if(chr %in% chr.targ){
      cov.chr <- covercounts.all[[chr]]
      ir.chr <- ranges(targets)[[chr]]
      tmp <- lapply(ir.chr, function(x) seqselect(cov.chr, x))

      # average coverage per target
      if(perTarget){
        avgcov <- sapply(tmp, mean)
        targetcov <- c(targetcov, avgcov)
      }
      
      # coverage per base
      cov.chr <- do.call(c, tmp)
      covercounts.target <- c(covercounts.target, RleList(cov.chr))
    }
  }
  names(covercounts.target) <- chr.targ
  
  # average coverage for targeted bases
  S <- sum(sum(covercounts.target))
  n <- sum(sapply(covercounts.target, length))
  res <- S / n
  names(res) <- "avgTargetCoverage"

  if(perTarget){
    targets$avgCoverage <- targetcov
    res <- c(res, list(targetCoverages=targets))
  }

  if(perBase)
    res <- c(res, list(coverageAll=covercounts.all, coverageTarget=covercounts.target))

  return(res)
}

