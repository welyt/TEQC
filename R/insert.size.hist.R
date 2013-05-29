insert.size.hist <-

# !!
#function(readpairs, returnInserts=FALSE, legendpos="topleft", main, xlab, ylab, breaks, col, ...){
function(readpairs, returnInserts=FALSE, legendpos="topleft", outline=FALSE, main, xlab, ylab, breaks, col, ...){
# !!

  # in case 'readpairs' contains also 'singleReads'
  if(is.list(readpairs) & ("readpairs" %in% names(readpairs)))
    readpairs <- readpairs$readpairs

  # insert sizes
  inserts <- width(readpairs)

  # graphical parameters
  if(missing(main)) main <- ""
  if(missing(xlab)) xlab <- "Insert size"
  if(missing(ylab)) ylab <- "Frequency"
  if(missing(breaks)) breaks <- 100
  if(missing(col)) col <- "lightblue"

  # average and median insert size
  m <- mean(inserts)
  me <- median(inserts)
  std <- sd(inserts)

# !!
  # if outline=F, remove "outliers" (according to boxplot.stats) before plotting
  if(!outline){
    x.out <- boxplot.stats(inserts)$out
    x.out <- min(x.out[x.out > m])
    inserts2 <- inserts[inserts < x.out]
  }
  else {
    inserts2 <- inserts
  }

  hist(inserts2, freq=TRUE, xlab=xlab, ylab=ylab, breaks=breaks, col=col, main=main, ...)
# !!

  abline(v=c(m - std, m, m + std, me), lty=2, col=c(1,2,1,3), lwd=2)
  legend(legendpos, c(paste("average (", round(m, 2), ")", sep=""),
                  paste("average +- SD (", round(std, 2), ")", sep=""),
                  paste("median (", round(me, 2), ")", sep="")), lty=2, col=c(2,1,3), lwd=2)

  if(returnInserts){
    names(inserts) <- values(readpairs)[[1L]][,1]
    return(inserts)
  }
}


