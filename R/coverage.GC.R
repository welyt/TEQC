coverage.GC <-
function(coverageAll, baits, returnBaitValues=FALSE, linecol="darkred", lwd, xlab, ylab, pch, col, cex, ...){

  # get GC content
  Seq <- unlist(values(baits))[,1]
  seqlen <- nchar(Seq)
  Cn <- gregexpr("C", Seq)
  Cn <- sapply(Cn, length)
  Gn <- gregexpr("G", Seq)
  Gn <- sapply(Gn, length)
  GC <- (Gn + Cn) / seqlen

  # get per-bait coverage
  covercounts.baits <- RleList()
  baitcov <- NULL
  for(chr in unique(space(baits))){
    cov.chr <- coverageAll[[chr]]
    ir.chr <- ranges(baits)[[chr]]
    tmp <- lapply(ir.chr, function(x) seqselect(cov.chr, x))
    avgcov <- sapply(tmp, mean)
    baitcov <- c(baitcov, avgcov)

    # per-base coverage
    # seqselect has to be done again on 'reduced' ranges, since baits might be overlapping!
    cov.chr <- seqselect(cov.chr, reduce(ir.chr))
    covercounts.baits <- c(covercounts.baits, RleList(cov.chr))
  }

  # average coverage for bait-covered bases
  avgcov <- mean(as.numeric(unlist(covercounts.baits)))

  # normalized per-bait coverage
  baitcov.norm <- baitcov / avgcov

  # set graphical parameters
  if(missing(xlab)) xlab <- "Bait GC content"
  if(missing(ylab)) ylab <- "Normalized coverage"
  if(missing(lwd)) lwd <- 3
  if(missing(col)) col <- "grey"
  if(missing(pch)) pch <- "."
  if(missing(cex)) cex <- 3

  # scatter plot
  plot(GC, baitcov.norm, xlab=xlab, ylab=ylab, col=col, pch=pch, cex=cex, ...)

  # loess curve
  if(length(unique(GC)) > 3)  # otherwise smooth.spline doesn't work
    lines(smooth.spline(x=GC, y=baitcov.norm), col=linecol, lwd=lwd)

  # return GC content and per-bait coverage within the 'baits' table if required
  if(returnBaitValues){
    baits$GCcontent <- GC
    baits$Coverage <- baitcov
    baits$NormalizedCoverage <- baitcov.norm
    return(baits)
  }
}
