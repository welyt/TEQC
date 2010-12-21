\name{get.reads}
\alias{get.reads}
\title{Read genomic positions of sequencing data}
\description{Reads a bedfile containing positions of sequenced read aligned to a reference genome
      and creates a \code{RangedData} object.}
      
\usage{get.reads(readsfile, chrcol = 1, startcol = 2, endcol = 3, idcol, zerobased = TRUE, sep = "\t", skip = 1, header = FALSE, ...)}

\arguments{
  \item{readsfile}{name of bedfile giving the positions of aligned reads}
  \item{chrcol}{in which column in the reads bedfile there is the chromosome information
    (chromosome information in the file should be in string format, e.g. "chrX")}
  \item{startcol}{in which column there are the starting positions of the reads}
  \item{endcol}{in which column there are the end positions of the reads}
  \item{idcol}{in which column there are read identifiers. For single-end data it is optionally.
    For paired-end data it is required for some functionalities. The two reads of one pair need to have the same ID.}
  \item{zerobased}{if \code{TRUE}, start coordinates in \code{readsfile} are assumed to be 0-based
    and are then converted to 1-based system by adding 1. If \code{FALSE}, coordinates are not
    shifted. In this case they should already be 1-based in \code{readsfile}.}
  \item{sep}{column separator character, defaults to tabs}
  \item{skip}{number of lines of the bedfile to skip before beginning to read data; defaults to 1}
  \item{header}{a logical value indicating whether the file contains the names of the variables as its first line; defaults to FALSE}
  \item{\dots}{further arguments passed to \code{read.delim}}
}
%\details{}

\value{A \code{\link[IRanges:RangedData-class]{RangedData}} table holding the read positions}

%\references{}
\author{Manuela Hummel \email{manuela.hummel@crg.es}}
%\note{}

\seealso{\code{\link{get.targets}}}

\examples{
exptPath <- system.file("extdata", package="TEQC")
readsfile <- paste(exptPath, "ExampleSet_Reads.bed", sep="/")
reads <- get.reads(readsfile, idcol=4, skip=0)
}

\keyword{ file }
