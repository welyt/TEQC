coverage.density <-
function(coveragelist, normalized=TRUE, legend, main, xlab, col, lwd, lty, xlim, ylim, ...){

  # in case coverages of just one sample are given
  if(!is.null(names(coveragelist)) & all(names(coveragelist) %in% c("avgTargetCoverage", "targetCoverageSD", "targetCoverageQuantiles", "targetCoverages", "coverageAll", "coverageTarget")))
    coveragelist <- list(coveragelist)

  # check input
  n.samples <- length(coveragelist)
  listnames <- lapply(coveragelist, names)
  if(!all(sapply(listnames, function(x) "coverageTarget" %in% x)))
    stop("elements of input list need element 'coverageTarget' (run 'coverage.target()' with option 'perBase=TRUE')")

  coverages <- lapply(coveragelist, function(x) x$coverageTarget)

  chroms <- names(coverages[[1]])
  L <- sapply(coverages[[1]], length)
  if(!all(sapply(coverages, function(x) identical(names(x), chroms))))
    stop("elements of 'covlist' do not seem to be coverages for the same targets (there are positions on different chromosomes for different samples)")
  if(!all(sapply(coverages, function(x) identical(sapply(x, length), L))))
    stop("elements of 'covlist' do not seem to be coverages for the same targets (there are different numbers of target positions for different samples)")

  coverages <- lapply(coverages, unlist, use.names=FALSE)

  # normalized coverages
  if(normalized){
    for(i in 1:n.samples)
      coverages[[i]] <- coverages[[i]] / coveragelist[[i]]$avgTargetCoverage
    if(missing(xlab)) xlab <- "Normalized coverage"
  }
  else if(missing(xlab)) xlab <- "Coverage"

  # densities
  dens <- lapply(coverages, function(x) density(as.numeric(x)))

  # graphical parameters
  if(missing(legend)){
    if(length(dens) == 1)
      legend <- NULL
    else if(is.null(names(coverages)))
      legend <- paste("sample", 1:n.samples)
    else
      legend <- names(coverages)
  }
  if(missing(main)) main <- ""
  if(missing(lwd)) lwd <- rep(2, n.samples)
  else if(length(lwd) < n.samples) lwd <- rep(lwd, length.out=n.samples)
  if(missing(lty)) lty <- rep(1, n.samples)
  else if(length(lty) < n.samples) lty <- rep(lty, length.out=n.samples)
  if(missing(col)) col <- 1:n.samples
  else if(length(col) < n.samples) col <- rep(col, length.out=n.samples)
  if(missing(xlim)){
    ma <- max(sapply(dens, function(d) max(d$x)))
    mi <- min(sapply(dens, function(d) min(d$x)))
    xlim <- c(mi, ma)
  }
  if(missing(ylim)) ylim <- c(0, max(sapply(dens, function(d) max(d$y))))

  # density plot
  plot(dens[[1]], col=col[1], lwd=lwd[1], lty=lty[1], xlab=xlab, main=main, xlim=xlim, ylim=ylim, ...)
  if(length(dens) > 1){
    for(i in 2:n.samples)
      lines(dens[[i]], col=col[i], lwd=lwd[i], lty=lty[i])
    if(!is.null(legend))
      legend("topright", legend=legend, col=col, lwd=lwd, lty=lty)
  }
}

