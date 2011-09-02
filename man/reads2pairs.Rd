\name{reads2pairs}
\alias{reads2pairs}
\title{Merges reads to read pairs}
\description{Combines the two reads of a read pair (in case of paired-end data) to a new 'range'
    starting at the first reads's start position and ending at the second read's end position.}

\usage{
reads2pairs(reads, max.distance)
}

\arguments{
  \item{reads}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of sequenced reads, i.e. output of \code{\link{get.reads}}. The first 'values' column has
      to contain read pair identifiers, i.e. when \code{reads} was created by \code{\link{get.reads}},
      the option \code{idcol} had to be specified. The input can also contain single reads without
      'read mate' (e.g. when the first read of a pair did not align to the reference genome, however the second one did
      align and was still kept). Those single reads will be returned in a separate table \code{singleReads}.
      When the two reads in a pair align to different chromosomes, they will also be returned in table \code{singleReads}.}
  \item{max.distance}{Integer value defining the maximum allowed distance between two reads of a pair
      (from start position of first read to end position of second read). Reads exceeding this
      distance will be returned in the separate table \code{singleReads}. If \code{max.distance} is
      not specified, reads will be joined to pairs regardless of their distance. Only when the two reads in a pair
      align to different chromosomes, they will be removed in any case and added to table \code{singleReads}.}
}
\details{The function puts together the two reads of each pair and creates new ranges spanning both
    reads and everything in between. Those ranges correspond to the extent of the actual DNA molecules
    for which both ends were sequenced. The output of the function can be used by several other functions,
    whenever calculations should be based on read pairs rather than on single reads, e.g.
    \code{\link{fraction.reads.target}}, \code{\link{readsPerTarget}}, \code{\link{duplicates.barplot}}}

\value{If \code{reads} only contains complete read pairs and for all pairs the respective reads
    align to the same chromosome and their distances do not exceed \code{max.distance} (if specified),
    a \code{\link[IRanges:RangedData-class]{RangedData}}
    object is returned containing positions of the merged reads per pair, ranging from start
    position of the first read to end position of the second read.

    If \code{reads} also contains single reads, or if reads within a pair are further apart than
    \code{max.distance} (if specified) or align to different chromosome, a list is returned with elements
    \item{singleReads}{\code{RangedData} object containing original positions of single reads without
          'read mates' and/or read pairs aligning too far apart or on different chromosomes}
    \item{readpairs}{\code{RangedData} object containing positions of the merged reads per pair, ranging from start
    position of the first read to end position of the second read}}

%\references{}
\author{Manuela Hummel \email{manuela.hummel@crg.es}}

%\note{}

\seealso{\code{\link{get.reads}}, \code{\link{fraction.reads.target}},
          \code{\link{readsPerTarget}}, \code{\link{duplicates.barplot}}, \code{\link{insert.size.hist}}}

\examples{
exptPath <- system.file("extdata", package="TEQC")
readsfile <- file.path(exptPath, "ExampleSet_Reads.bed")
reads <- get.reads(readsfile, idcol=4, skip=0)
readpairs <- reads2pairs(reads)
}

\keyword{ math }
