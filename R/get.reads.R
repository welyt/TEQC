get.reads <-
function(readsfile, chrcol=1, startcol=2, endcol=3, idcol, zerobased=TRUE, sep="\t", skip=1, header=FALSE, ...){

  dat <- read.delim(readsfile, sep=sep, as.is=TRUE, skip=skip, header=header, ...)

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

  print(paste("read", nrow(rd), "sequenced reads"))
  return(rd)
}

