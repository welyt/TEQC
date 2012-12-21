# get on- and off-target read multiplicities per chromosome
multfun <- function(x){
  if(nrow(x) > 0){
    r <- sort(x$ranges)  # just in case reads are not sorted
    dups <- Rle(duplicated(r))
    m <- runLength(dups)[runValue(dups)] + 1  # -> reads with > 1 copies
    m1 <- nrow(x) - sum(m)                    # -> unique reads
    c("1"=m1, table(m))
  }
}


duplicates.barplot <-
function(reads, targets, returnDups=FALSE, truncateX, col=c("red","lightblue"), xlab, ylab, ylim, ...){

  # in case 'reads' is output of 'reads2pairs' and contains also 'singleReads'
  if(is.list(reads) & ("readpairs" %in% names(reads)))
    reads <- reads$readpairs

  # which reads are on target
  on.target <- overlapsAny(reads, targets)
  reads.on <- reads[on.target,]
  reads.off <- reads[!on.target,]

  params.on <- RDApplyParams(rangedData=reads.on, applyFun=multfun)
  multi.on <- unlist(rdapply(params.on))
  params.off <- RDApplyParams(rangedData=reads.off, applyFun=multfun)
  multi.off <- unlist(rdapply(params.off))

  # summarize over chromosomes
  m.on <- sapply(strsplit(names(multi.on), "\\."), function(x) x[2])
  T1 <- tapply(multi.on, m.on, sum)
  T1 <- T1[order(as.numeric(names(T1)))]
  m.off <- sapply(strsplit(names(multi.off), "\\."), function(x) x[2])
  T2 <- tapply(multi.off, m.off, sum)
  T2 <- T2[order(as.numeric(names(T2)))]

  # prepare data for barplot
  N <- union(names(T1), names(T2))
  l <- length(N)
  barmat0 <- matrix(0, nrow=2, ncol=l)
  colnames(barmat0) <- N
  barmat0[1, names(T1)] <- T1
  barmat0[2, names(T2)] <- T2

  # plot fractions instead of absolute numbers (separate for on- and off-target)
  barmat <- barmat0 / rowSums(barmat0)

  # graphical parameters
  if(missing(xlab)) xlab <- "Read multiplicity"
  if(missing(ylab)) ylab <- "Fraction of reads"
  if(missing(ylim)) ylim <- c(0, max(barmat) + 0.05)

  if(!missing(truncateX))
    l <- truncateX
  B <- barplot(barmat[,1:l], beside=TRUE, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
  legend("topright", c("on target", "off target"), fill=col)

  # add absolute numbers (in millions) on top of the bars
  nr <- round(as.vector(barmat0[,1:l]) / 1e6, 2)
  text(x=as.vector(B), y=as.vector(barmat[,1:l])+ 0.02, labels=nr, cex=.8)

  if(returnDups){
    rownames(barmat) <- rownames(barmat0) <- c("on target", "off target")
    return(list(absolute=barmat0, relative=barmat))
  }
}
