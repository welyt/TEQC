coverage.uniformity <-
function(coveragelist, addlines=TRUE, add=FALSE, xlab, ylab, xlim, ylim, col, lwd, ...){

  if(!("coverageTarget" %in% names(coveragelist)))
    stop("input list needs element 'coverageTarget' (run 'coverage.target()' with option 'perBase=TRUE')")

  # normalize per-base coverages by mean coverage
  covercounts <- coveragelist$coverageTarget
  covercounts <- unlist(covercounts)
  covercounts.norm <- covercounts / coveragelist$avgTargetCoverage

  # cumulative fraction of bases per normalized coverage
  frac <- table(covercounts.norm) / length(covercounts.norm)
  cf <- cumsum(rev(frac))

  # graphical parameters
  if(missing(xlab)) xlab <- "Normalized coverage"
  if(missing(ylab)) ylab <- "Fraction of bases"
  if(missing(xlim)) xlim <- 0:1
  if(missing(ylim)) ylim <- 0:1
  if(missing(col)) col <- "darkred"
  if(missing(lwd)) lwd <- 2

  # coverage distribution
  if(!add)
    plot(as.numeric(names(cf)), cf, type="l", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=col, lwd=lwd, ...)
  else
    lines(as.numeric(names(cf)), cf, col=col, lwd=lwd, ...)

  # indicate with lines the fraction of bases achieving at least mean coverage or half the mean coverage
  if(addlines){
    a <- cbind(abs(as.numeric(names(cf)) - 1), abs(as.numeric(names(cf)) - 0.5))
    b <- apply(a, 2, function(x) which(x == min(x)))
    lines(x=c(1, 1), y=c(0, cf[b[1]]), lty=2)
    lines(x=c(0, 1), y=rep(cf[b[1]], 2), lty=2)
    lines(x=c(0.5, 0.5), y=c(0, cf[b[2]]), lty=2)
    lines(x=c(0, 0.5), y=rep(cf[b[2]], 2), lty=2)
  }
}
