\name{get.reads}
\alias{get.reads}
\title{Read genomic positions of sequencing data}
\description{Reads a bedfile containing positions of sequenced read aligned to a reference genome
      and creates a \code{RangedData} object.}
      
\usage{get.reads(readsfile, filetype = c("bed", "bam"), chrcol = 1, startcol = 2, endcol = 3, idcol, zerobased = TRUE, sep = "\t", skip = 1, header = FALSE, ...)}

\arguments{
  \item{readsfile}{name of bedfile giving the positions of aligned reads}
  \item{filetype}{Input file type. If \code{"bam"}, the .bam file is read using
  \code{\link[Rsamtools]{scanBam}}, where flag option \code{isUnmappedQuery=FALSE}
  is used. Defaults to \code{"bed"}}
  \item{chrcol}{In which column in the reads bedfile there is the chromosome information
    (chromosome information in the file should be in string format, e.g. "chrX").
    Ignored if \code{filetype = "bam"}.}
  \item{startcol}{In which column there are the starting positions of the reads.
  Ignored if \code{filetype = "bam"}.}
  \item{endcol}{In which column there are the end positions of the reads.
  Ignored if \code{filetype = "bam"}.}
  \item{idcol}{In which column there are read identifiers. For single-end data it is optionally.
    For paired-end data it is required for some functionalities. The two reads of one pair need to have the same ID.
    Ignored if \code{filetype = "bam"} (the ID column is automatically included then).
    If read IDs include "#0/1" and "#0/2"
        in the end (indicating read 1 and read 2 of a pair), those characters will
        be removed from the IDs.}
  \item{zerobased}{if \code{TRUE}, start coordinates in \code{readsfile} are assumed to be 0-based
    and are then converted to 1-based system by adding 1. If \code{FALSE}, coordinates are not
    shifted. In this case they should already be 1-based in \code{readsfile}.
    Ignored if \code{filetype = "bam"}, since \code{\link[Rsamtools]{scanBam}} converts 0-based
    to 1-based coordinates.}
  \item{sep}{Column separator character, defaults to tabs. Ignored if \code{filetype = "bam"}.}
  \item{skip}{Number of lines of the bedfile to skip before beginning to read data; defaults to 1.
    Ignored if \code{filetype = "bam"}.}
  \item{header}{A logical value indicating whether the file contains the names of the variables as its first line;
    defaults to FALSE. Ignored if \code{filetype = "bam"}.}
  \item{\dots}{Further arguments passed to \code{read.delim}. Ignored if \code{filetype = "bam"}.}
}

%\details{}

\value{A \code{\link[IRanges:RangedData-class]{RangedData}} table holding the read positions}

%\references{}
\author{Manuela Hummel \email{manuela.hummel@crg.es}}
%\note{}

\seealso{\code{\link{get.targets}}}

\examples{
exptPath <- system.file("extdata", package="TEQC")
readsfile <- file.path(exptPath, "ExampleSet_Reads.bed")
reads <- get.reads(readsfile, idcol=4, skip=0)
}

\keyword{ file }
