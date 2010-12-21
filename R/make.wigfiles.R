make.wigfiles <-
function(coverageAll, chroms, trackname="Coverage", filename="Coverage"){

  if(missing(chroms))
    chroms <- names(coverageAll)
  else if(!all(chroms %in% names(coverageAll)))
    stop("chromosome names are wrong, or there are no reads/targets on selected chromosomes")

  for(chr in chroms){
    print(paste("preparing wiggle for chromsome", substr(chr, 4, nchar(chr))))

    # wiggle track definition line
    wig <- NULL
    wig[1] <- paste("track type=wiggle_0 name=\"", chr, trackname, "_WIG\" description=\"", chr, trackname, "_WIG\" visibility=full color=255,0,0 altColor=0,128,0 priority=20 graphType=points yLineMark=0 yLineOnOff=on windowingFunction=mean", sep="")

    # wiggle format ('variableStep' format because of irregular intervals between data points)
    wig[2] <- paste("variableStep chrom=", chr, sep="")
    outfile <- paste(filename, "_", chr, ".wig", sep="")
    write(wig, outfile)

    # coverage values
    covercounts.chr <- coverageAll[[chr]]
    
    # keep only non-0 coverages
    ind <- which(covercounts.chr > 0)
    covercounts.chr <- as.integer(covercounts.chr[ind])
    names(covercounts.chr) <- ind
    write.table(covercounts.chr, file=outfile, sep="\t", quote=FALSE, col.names=FALSE, append=TRUE)
  }
}

