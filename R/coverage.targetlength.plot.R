coverage.targetlength.plot <-
function(targets, plotcolumn, linecol=2, xlab, ylab, lwd, pch, cex, ...){

  if(ncol(targets) == 0)
    stop("'targets' does not include values to plot")

  if(is.character(plotcolumn) & !(plotcolumn %in% colnames(targets)))
    stop("selected 'plotcolumn' does not exist")

  targetlen <- width(targets)
  y <- targets[[plotcolumn]]

  # set graphical parameters
  if(is.numeric(plotcolumn))
    plotcolumn <- colnames(targets)[plotcolumn]
  if(missing(xlab)) xlab <- "Target length (bp)"
  if(missing(ylab)) ylab <- plotcolumn
  if(missing(lwd)) lwd <- 3
  if(missing(cex)) cex <- 2
  if(missing(pch)) pch <- "."

  # scatter plot
  plot(targetlen, y, xlab=xlab, ylab=ylab, pch=pch, cex=cex, ...)

  # loess curve
  if(length(unique(targetlen)) > 3)  # otherwise smooth.spline doesn't work
    lines(smooth.spline(x=targetlen, y=y), col=linecol, lwd=lwd)
}

