reads2pairs <-
function(reads){

  if(ncol(reads) == 0)
    stop("a 'values' column (first column) with read pair IDs is needed in the 'reads' RangedData table")

  # check if there are as many 1st as 2nd reads
  id.all <- sapply(values(reads), nrow)
  id.pair <- sapply(values(reads), function(x) length(unique(x[,1])))
  singlereads <- any(id.all / id.pair != 2)

  # if there are single reads, give them back separately
  if(singlereads){
    ID <- unlist(values(reads))[,1]
    dups <- duplicated(ID)
    nondups <- !(ID %in% ID[dups])
    print(paste("there were", sum(nondups), "single reads found without matching second read"))
    res <- list(singleReads=reads[nondups,,drop=TRUE])
    reads <- reads[!nondups,,drop=TRUE]
  }

  # merge reads of each pair
  mergefun <- function(x){
               # sort by 1) pair ID and 2) position
               o <- order(values(x)[[1L]][,1], start(x$ranges))
               r <- x$ranges[o,]

               # split ranges table into one for 1st and one for 2nd read
               x.split <- split(r, rep(1:2, length.out=nrow(r)))

               # merge the tables such that each line starts with start position of first read
               #   of the pair and ends with end position of second read
               ir <- IRanges(start=start(x.split[[1]]), end=end(x.split[[2]]))

               # sort again according to position
               o2 <- order(ir)
               ir <- ir[o2,]
               sel <- seq(1, length(r)-1, by=2)
               ID <- values(x)[[1L]][o[sel][o2], 1]

               RangedData(ir, space=space(x)[sel], ID=ID)
              }
  res2 <- endoapply(reads, mergefun)

  if(singlereads)
    res <- c(res, readpairs=res2)
  else
    res <- res2
  return(res)
}
