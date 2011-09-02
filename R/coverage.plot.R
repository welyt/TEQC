coverage.plot <-
function(coverageAll, targets, chr, Start, End, Offset=0, add=FALSE,
                          col.line=1, col.target="orange", col.offset="yellow", xlab, ylab, ylim, ...){

  covercounts <- coverageAll[[chr]]
  chrom <- substr(chr, 4, nchar(chr))

  # stop if all reads lie "left" of the selected Start position
  L <- length(covercounts)
  if(L < Start)
    stop(paste("no reads falls into the selected region on chromosome", chrom))

  # add 0's when 'End' is "right" of largest read position
  if(L < End)
    covercounts <- c(covercounts, Rle(rep(0, End-L)))

  ir <- IRanges(start=Start, end=End)
  covsel <- seqselect(covercounts, ir)

  # also stop if coverage is 0 for all bases in selected region
  if(all(covsel == 0))
    stop(paste("no reads falls into the selected region on chromosome", chrom))

  # graphical parameters
  if(missing(xlab)) xlab <- paste("Chromosome", chrom)
  if(missing(ylab)) ylab <- "Coverage"

  # plot coverages along the selected chromosomal region
  if(!add){
    if(missing(ylim)){
      ma <- max(covsel)
      mi <- .04 * ma
      ylim <- c(-mi, ma)
    }
    plot(Start:End, covsel, ylim=ylim, type="l", col=col.line, xlab=xlab, ylab=ylab, ...)
    abline(h=0, lty=3)
  }
  else
    lines(Start:End, covsel, col=col.line, ...)

  # add bars showing the location of targets [+ Offset]
  if(!missing(targets)){
    tar <- intersect(ir, ranges(targets)[[chr]])
    if(length(tar) > 0){   # ... only when there are targets in the selected region
      mi <- par("usr")[3] / 2
      if(Offset > 0){
        tar2 <- tar + Offset
        rect(start(tar2), mi*1.2, end(tar2), mi*.6, col=col.offset)
      }
      rect(start(tar), mi*1.2, end(tar), mi*.6, col=col.target)
    }
  }
}


