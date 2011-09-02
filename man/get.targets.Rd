\name{get.targets}
\alias{get.targets}
\title{Read capture target positions}
\description{Reads a bedfile containing positions of the capture targets and creates a \code{RangedData} object.}

\usage{get.targets(targetsfile, chrcol = 1, startcol = 2, endcol = 3, zerobased = TRUE, sep = "\t", skip = 1, header = FALSE, ...)}

\arguments{
  \item{targetsfile}{name of bedfile giving the positions of each target region}
  \item{chrcol}{in which column in the targets bedfile there is the chromosome information
    (chromosome information in the file should be in string format, e.g. "chrX")}
  \item{startcol}{in which column there are the starting positions of the targeted regions}
  \item{endcol}{in which column there are the end positions of the targeted regions}
  \item{zerobased}{if \code{TRUE}, start coordinates in \code{targetsfile} are assumed to be 0-based
    and are then converted to 1-based system by adding 1. If \code{FALSE}, coordinates are not
    shifted. In this case they should already be 1-based in \code{targetsfile}.}
  \item{sep}{column separator character, defaults to tabs}
  \item{skip}{number of lines of the bedfile to skip before beginning to read data; defaults to 1}
  \item{header}{a logical value indicating whether the file contains the names of the variables as its first line; defaults to FALSE}
  \item{\dots}{further arguments passed to \code{read.delim}}
}
%\details{}

\value{A \code{\link[IRanges:RangedData-class]{RangedData}} table holding the target region positions. Note that overlapping or
  adjacent regions are collapsed to one region.}

%\references{}
\author{Manuela Hummel \email{manuela.hummel@crg.es}}

\note{Since overlapping regions are collapsed, the input bedfile can also contain positions
  of the (in most cases overlapping) hybridization probes used for the target capture.}

\seealso{\code{\link{get.reads}}}

\examples{
exptPath <- system.file("extdata", package="TEQC")
targetsfile <- file.path(exptPath, "ExampleSet_Targets.bed")
targets <- get.targets(targetsfile, skip=0)
}

\keyword{ file }
