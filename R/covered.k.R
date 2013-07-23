covered.k <-
function(coverageTarget, k=c(1,2,3,5,10,20)){

  # coverages
  covercounts <- unlist(coverageTarget, use.names=FALSE)
  tab <- as.matrix(table(covercounts))

  # add evtl. missing numbers of covering reads to 'tab'
  miss <- setdiff(k, rownames(tab))
  if(length(miss) > 0){
    add <- matrix(0, length(miss))
    rownames(add) <- miss
    tab <- rbind(tab, add)
  }
  tab <- tab[order(as.numeric(rownames(tab))),]

  # cumulative coverages
  cs <- cumsum(rev(tab))

  # get coverages for k-values
  cs.k <- cs[as.character(k)]

  # fractions of bases with coverage >k
  cs.k / sum(tab)
}

