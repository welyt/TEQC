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
  targetcov <- targetSD <- NULL
  for(chr in names(covercounts.all)){
    if(chr %in% chr.targ){
      cov.chr <- covercounts.all[[chr]]
      ir.chr <- ranges(targets)[[chr]]
      tmp <- lapply(ir.chr, function(x) seqselect(cov.chr, x))

      # coverage average and SD per target
      if(perTarget){
        avgcov <- sapply(tmp, mean)
        sdcov <- sapply(tmp, sd)
        targetcov <- c(targetcov, avgcov)
        targetSD <- c(targetSD, sdcov)
      }
      
      # coverage per base
      cov.chr <- do.call(c, tmp)
      covercounts.target <- c(covercounts.target, RleList(cov.chr))
    }
  }
  names(covercounts.target) <- chr.targ
  
  # coverage average, SD and quartiles for all targeted bases
  tmp <- as.integer(unlist(covercounts.target))
  avg <- mean(tmp)
  std <- sd(tmp)
  qu <- quantile(tmp)

  res <- list(avgTargetCoverage=avg, targetCoverageSD=std, targetCoverageQuantiles=qu)

  if(perTarget){
    targets$avgCoverage <- targetcov
    targets$coverageSD <- targetSD
    res <- c(res, list(targetCoverages=targets))
  }

  if(perBase)
    res <- c(res, list(coverageAll=covercounts.all, coverageTarget=covercounts.target))

  return(res)
}

