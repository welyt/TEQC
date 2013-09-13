coverage.target <-
function(reads, targets, Offset=0, perTarget=TRUE, perBase=TRUE){

  # add Offset if required
  targets <- offsetfun(Offset=Offset, targets=targets)

  # get largest position per chromosome - either target or reads position
  max.targ <- end(range(targets))
  chr.targ <- names(max.targ)
  max.read <- end(range(reads))

  w <- sapply(max.read, length) == 0
  reads <- reads[!w]
  max.read <- max.read[!w]

  # if on a chromosome there are targets but no reads
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

#!! BUG fix
#  for(chr in names(covercounts.all)){  -> chr's don't come in order of chr.targ (when reads come from bam file)
#    if(chr %in% chr.targ){             -> avgCoverage, coverageSD and coverageTarget outputs are wrong !!
   for(chr in chr.targ){
      cov.chr <- covercounts.all[[chr]]
      ir.chr <- ranges(targets)[[chr]]
      tmp <- lapply(ir.chr, function(x) cov.chr[x]) # use [ instead of
                                                    # deprecated seqselect

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
#    }
  }

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

  if(perBase){
    names(covercounts.target) <- chr.targ
    res <- c(res, list(coverageAll=covercounts.all, coverageTarget=covercounts.target))
  }
  
  return(res)
}

