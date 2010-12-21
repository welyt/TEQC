coverage.hist <-
function(coverageTarget, col.hist="lightblue", col.line="orange", covthreshold, breaks="Sturges", xlab, ylab, main, lwd, ...){

  # graphical parameters for the histogram and line
  par(mar=c(5,4,4,4))
  if(missing(xlab)) xlab <- "Coverage"
  if(missing(ylab)) ylab <- "Fraction of bases"
  if(missing(main)) main <- "Coverage Distribution"
  if(missing(lwd)) lwd <- 2

  # histogram of per-base coverages (with relative frequencies)
  covercounts <- as.numeric(unlist(coverageTarget))
  H <- hist(covercounts, breaks=breaks, plot=FALSE)
  H$counts <- H$counts / sum(H$counts)
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
  mtext("Cumulative fraction of bases", side=4, line=2.2)
  
  # lines indicating which fraction of bases is covered by at least 'covthreshold' reads
  if(!missing(covthreshold)){
    a <- abs(as.numeric(names(cs.s)) - covthreshold)
    b <- which(a == min(a))
    lines(x=c(covthreshold, covthreshold), y=c(0, cs.s[b]), lty=2)
    lines(x=c(covthreshold, par("usr")[2]), y=rep(cs.s[b], 2), lty=2)
    legend("topright", paste(covthreshold, "X coverage", sep=""), lty=2)
  }
}

