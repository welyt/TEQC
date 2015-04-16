\name{get.baits}
\alias{get.baits}
\title{Read capture hybridization probe positions}
\description{Reads a file containing positions and sequences of the capture hybridization probes
  and creates a \code{RangedData} object.}

\usage{get.baits(baitsfile, chrcol = 1, startcol = 2, endcol = 3, seqcol = 4, zerobased = TRUE, sep = "\t", header = TRUE, ...)}

\arguments{
  \item{baitsfile}{name of file giving the positions and sequences of each hybridization probe ("bait")}
  \item{chrcol}{in which column in \code{baitsfile} there is the chromosome information
    (chromosome information in the file should be in string format, e.g. "chrX")}
  \item{startcol}{in which column there are the starting positions of the baits}
  \item{endcol}{in which column there are the end positions of the baits}
  \item{seqcol}{in which column there are the sequences of the baits}
  \item{zerobased}{if \code{TRUE}, start coordinates in \code{baitsfile} are assumed to be 0-based
    and are then converted to 1-based system by adding 1. If \code{FALSE}, coordinates are not
    shifted. In this case they should already be 1-based in \code{baitsfile}.}
  \item{sep}{column separator character, defaults to tabs}
  \item{header}{a logical value indicating whether the file contains the names of the variables as its first line; defaults to FALSE}
  \item{\dots}{further arguments passed to \code{read.delim}}
}
\details{The \code{baitsfile} containing positions and sequences of hybridization probes
    has to be created beforehand, in many cases manually. (The function was made like this in order to
    keep things as general and platform independent as possible.) E.g. with baits designed by
    Agilent's eArray tool, the \code{baitsfile} can be created by merging the files
    '..._D_BED_...bed' and '..._D_DNAFront_BCBottom_...txt'.}

\value{A \code{\link[IRanges:RangedData-class]{RangedData}} table holding the hybridization probe
  ("bait") positions and sequences. Overlapping or adjacent baits are not collapsed.}

%\references{}
\author{Manuela Hummel \email{m.hummel@dkfz.de}}

%\note{}

\seealso{\code{\link{get.reads}}, \code{\link{get.targets}}}

\examples{
exptPath <- system.file("extdata", package="TEQC")
baitsfile <- file.path(exptPath, "ExampleSet_Baits.txt")
baits <- get.baits(baitsfile, chrcol=3, startcol=4, endcol=5, seqcol=2)
}

\keyword{ file }
