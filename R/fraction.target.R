fraction.target <-
function(targets, Offset=0, genome=c(NA, "hg19", "hg18"), genomesize){

  # add Offset if required
  targets <- offsetfun(Offset=Offset, targets=targets)

  if(missing(genomesize)){
    genome <- match.arg(genome)
    if(is.na(genome))
      stop("either 'genome' or 'genomesize' has to be specified")
    genomesize <- switch(genome,
                        hg18 = 3107677273,
                        hg19 = 3137161264)
  }
  regionsize <- sum(width(targets))
  regionsize / genomesize
}

