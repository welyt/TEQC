# function to merge reads of each pair
mergefun <-
function(x){
  # sort by 1) pair ID and 2) position
  o <- order(values(x)[[1L]][,1], start(x$ranges))
  r <- x$ranges[o,]

  # split ranges table into one for 1st and one for 2nd read
  x.split <- split(r, rep(1:2, length.out=length(r)))  

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



reads2pairs <-
#function(reads){
function(reads, max.distance){

  if(ncol(reads) == 0)
    stop("a 'values' column (first column) with read pair IDs is needed in the 'reads' RangedData table")
  
  # check if there are as many 1st as 2nd reads on each chromosome
  # -> single reads or reads where pairs matching different chromosomes will be removed
  ID <- values(reads)[,1]
  dups <- sapply(ID, function(x) duplicated(x))
  singlereads <- any(sapply(dups, function(x) sum(x) != length(x)/2))   

  # give back single reads separately
  if(singlereads){
    ID <- unlist(ID)
    ID2 <- ID[unlist(dups)]     

    # check if some IDs are not unique to one pair  
    if(any(duplicated(ID2)))                        
      stop("read pair IDs do not seem to be unique") 
    
    nondups <- Rle(!(ID %in% ID2))     
    print(paste("there were", sum(nondups), "reads found without matching second read, or whose second read matches to a different chromosome"))
    res <- list(singleReads=reads[which(nondups),,drop=TRUE])
    reads <- reads[which(!nondups),,drop=TRUE]
  }
  
  # merge reads of each pair
  res2 <- endoapply(reads, mergefun)
  
  # remove read pairs that are more than max.distance apart (within same chromosome)
  if(!missing(max.distance)){
    toofar <- which(width(res2) > max.distance)
    if(length(toofar) > 0){
      print(paste("there were", length(toofar), "read pairs found with distance between the two reads exceeding max.distance =", max.distance))
      id.toofar <- res2$ID[toofar]
      reads.toofar <- reads[values(reads)[,1] %in% id.toofar,,drop=T]
      res2 <- res2[-toofar,,drop=T]
      
      # add them to 'singleReads' output
      if(singlereads)
        res$singleReads <- rbind(res$singleReads, reads.toofar)
      else
        res <- list(singleReads=reads.toofar)
    }
  }

  if(singlereads | !missing(max.distance))
    res <- c(res, readpairs=res2)
  else
    res <- res2
  return(res)
}
