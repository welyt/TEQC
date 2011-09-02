get.reads <-
function(readsfile, filetype=c("bed", "bam"), chrcol=1, startcol=2, endcol=3, idcol, zerobased=TRUE, sep="\t", skip=1, header=FALSE, ...){

  filetype <- match.arg(filetype)
  if(filetype == "bam"){
    # read BAM file
    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE), what=c("qname", "pos", "qwidth", "rname"))
    aln <- scanBam(readsfile, param=param)[[1]]

    # create RangedData object
    rd <- with(aln, RangedData(IRanges(pos, width=qwidth), ID=qname,  space=rname))
  }

  else {
    # read only the required columns
    n <- max(count.fields(readsfile, sep=sep))
    colclasses <- rep("NULL", n)
    colclasses[chrcol] <- "character"
    colclasses[c(startcol, endcol)] <- "integer"
    if(!missing(idcol))
      colclasses[idcol] <- "character"

    dat <- read.delim(readsfile, colClasses=colclasses, sep=sep, skip=skip, header=header, ...)

    # sort reads (better for example when calling 'findOverlaps()' several times)
    o <- order(dat[,chrcol], dat[,startcol])
    if(!identical(o, 1:nrow(dat)))
      dat <- dat[o,]

    # make IRanges object
    ir <- IRanges(start=dat[,startcol], end=dat[,endcol])

    # shift start position forward by 1 to go from 0-based to 1-based system
    if(zerobased)
      start(ir) <- start(ir) + 1

    # make RangedData object
    if(missing(idcol))
      rd <- RangedData(ranges=ir, space=dat[,chrcol])
    else
      rd <- RangedData(ranges=ir, space=dat[,chrcol], ID=dat[,idcol])
  }

  # check for Illumina read pair IDs - #0/1 and #0/2 have to be removed
  if(length(colnames(rd)) > 0){
    illu <- grep("#0/", rd$ID[1])
    if(length(illu) > 0)
      rd$ID <- gsub("#0/[[:digit:]]", "", rd$ID)
  }
  
  print(paste("read", nrow(rd), "sequenced reads"))
  return(rd)
}


