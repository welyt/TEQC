coverage.hist <-
function(coverageTarget, col.hist="lightblue", col.line="orange", covthreshold, outline=FALSE, breaks="Sturges", xlab, ylab, main, lwd, ...){

  # graphical parameters for the histogram and line
  par(mar=c(5,4,4,4))
  if(missing(xlab)) xlab <- "Coverage"
  if(missing(ylab)) ylab <- "Fraction of target bases"
  if(missing(main)) main <- "Coverage Distribution"
  if(missing(lwd)) lwd <- 2

  covercounts <- as.numeric(unlist(coverageTarget, use.names=FALSE))
  
  # remove outliers
  if(!outline){
    x.out <- boxplot.stats(covercounts)$out
    m <- mean(covercounts)
    x.out <- min(x.out[x.out > m])

    # !! use this reduced version only for plotting - but all calculations still have to be based on the complete coverage data !!
    covercounts2 <- covercounts[covercounts < x.out]
  }
  else
    covercounts2 <- covercounts
  # !!
  
  # histogram of per-target-base coverages (with relative frequencies)
  H <- hist(covercounts2, breaks=breaks, plot=FALSE)
  
  # !!
  #H$counts <- H$counts / sum(H$counts)
  H$counts <- H$counts / length(covercounts)
  # !!
  
  plot(H, freq=TRUE, xlab=xlab, ylab=ylab, main=main, col=col.hist, ...)

  # cumulative coverages
  tab <- table(covercounts)
  cs <- cumsum(rev(tab))

  # standardize cumulative coverages in order to fit them into the histogram plot
  m <- max(H$counts)
  cs.s <- cs * (m / sum(tab))

  # line plot
  lines(x=names(cs.s), y=cs.s, col=col.line, lwd=lwd)
  axis(side=4, at=seq(0, m, length.out=11), labels=seq(0, 1, by=0.1))
  mtext("Cumulative fraction of target bases", side=4, line=2.2)
  
  # lines indicating which fraction of bases is covered by at least 'covthreshold' reads
  if(!missing(covthreshold)){
    a <- abs(as.numeric(names(cs.s)) - covthreshold)
    b <- which(a == min(a))
    lines(x=c(covthreshold, covthreshold), y=c(0, cs.s[b]), lty=2)
    lines(x=c(covthreshold, par("usr")[2]), y=rep(cs.s[b], 2), lty=2)
    legend("topright", paste(covthreshold, "X coverage", sep=""), lty=2)
  }
}



# internal function - same as coverage.hist, just the cumulative sum of target
#    base fraction with coverage x is returned -> needed for multiTEQCreport.R
.coverage.hist <-
function(coverageTarget, col.hist="lightblue", col.line="orange", covthreshold, outline=FALSE, breaks="Sturges", xlab, ylab, main, lwd, ...){

  # graphical parameters for the histogram and line
  par(mar=c(5,4,4,4))
  if(missing(xlab)) xlab <- "Coverage"
  if(missing(ylab)) ylab <- "Fraction of target bases"
  if(missing(main)) main <- "Coverage Distribution"
  if(missing(lwd)) lwd <- 2

  covercounts <- as.numeric(unlist(coverageTarget, use.names=FALSE))

  if(!outline){
    x.out <- boxplot.stats(covercounts)$out
    m <- mean(covercounts)
    x.out <- min(x.out[x.out > m])
    
    # !! use this reduced version only for plotting - but all calculations still have to be based on the complete coverage data !!
    covercounts2 <- covercounts[covercounts < x.out]
  }
  else
    covercounts2 <- covercounts
  # !!
  
  # histogram of per-target-base coverages (with relative frequencies)
  H <- hist(covercounts2, breaks=breaks, plot=FALSE)
  
  # !!
  #H$counts <- H$counts / sum(H$counts)
  H$counts <- H$counts / length(covercounts)
  # !!

  plot(H, freq=TRUE, xlab=xlab, ylab=ylab, main=main, col=col.hist, ...)

  # cumulative coverages
  tab <- table(covercounts)
  cs <- cumsum(rev(tab))

  # standardize cumulative coverages in order to fit them into the histogram plot
  m <- max(H$counts)
  cs.s <- cs * (m / sum(tab))

  # line plot
  lines(x=names(cs.s), y=cs.s, col=col.line, lwd=lwd)
  axis(side=4, at=seq(0, m, length.out=11), labels=seq(0, 1, by=0.1))
  mtext("Cumulative fraction of target bases", side=4, line=2.2)

  # lines indicating which fraction of bases is covered by at least 'covthreshold' reads
  if(!missing(covthreshold)){
    a <- abs(as.numeric(names(cs.s)) - covthreshold)
    b <- which(a == min(a))
    lines(x=c(covthreshold, covthreshold), y=c(0, cs.s[b]), lty=2)
    lines(x=c(covthreshold, par("usr")[2]), y=rep(cs.s[b], 2), lty=2)
    legend("topright", paste(covthreshold, "X coverage", sep=""), lty=2)
  }

  # cumulative coverage fractions
#!! wrong
  #res <- cumsum(cs / sum(cs))
  res <- cs / sum(tab)
#!!

  res <- data.frame(coverage=rev(names(res)), fractionTargetBases=rev(res))
  return(res)
}
