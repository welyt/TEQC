\name{duplicates.barplot}
\alias{duplicates.barplot}
\title{Read duplicates barplot}
\description{Barplot showing fractions of reads / read pairs which are unique and for which there are two, three, ... copies.
    Separate bars are made for on- and off-target reads / read pairs}

\usage{
duplicates.barplot(reads, targets, returnDups=FALSE, truncateX, col=c("red","lightblue"), xlab, ylab, ylim, ...)
}

\arguments{
  \item{reads}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of sequenced reads, i.e. output from \code{\link{get.reads}}. Alternatively, for paired-end data, it can be the output
      of \code{\link{reads2pairs}} when multiplicities of read \emph{pairs} instead
      of fraction of single reads shall be visualized.}
  \item{targets}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of target regions, i.e. output from \code{\link{get.targets}}}
  \item{returnDups}{if \code{TRUE}, on- and off-target read / read pair multiplicities are returned}
  \item{truncateX}{integer; show bars only up to a read / read pair multiplicity of \code{truncateX} (x-axis)}
  \item{col}{vector specifying the two colors of bars and legend for on- and off-target read multiplicities}
  \item{xlab, ylab}{x- and y-axis labels}
  \item{ylim}{y-axis coordinate ranges}
  \item{\dots}{further graphical parameters passed to \code{barplot}}
}
\details{Single-end reads are considered as duplicates if they have same start end end position.
    Paired-end read pairs are considered as duplicates if start and end positions of both reads of
    the pairs are identical. Usually, duplicates are removed before further analyses (e.g. SNP detection),
    because they could represent PCR artefacts. However, in target capture experiments it is likely
    to have also many "real" duplicates (actual different molecules that happen to start at same position)
    due to the "enrichment" of the target regions.
    The separation in the barplot between on- and off-target reads / read pairs gives an impression on
    whether on-target there are more reads with higher multiplicites, which hence might indicate a reasonable amount
    of "real" duplication. A paired-end read pair is considered on-target if at least one of its reads
    overlaps with a target.}

\value{Barplot where the bar heights correspond to fractions of reads / read pairs which are present in the data
    with the respective number of copies (x-axis). Fractions are calculated separately for on- and off-target
    reads / read pairs. A read pair is considered on-target if at least one of its reads
    overlaps with a target. Absolute numbers (in millions) are additionally written on top of the bars.

    If \code{returnDups} equals \code{TRUE}, a list with two elements \code{absolute} and
    \code{relative} is returned. The former is a matrix that contains the absolute numbers of reads / read pairs
    for each multiplicity (columns), for both on- and off-target reads / read pairs (rows).
    The latter gives row-based fractions which correspond to the bar heights.}

%\references{}
\author{Manuela Hummel \email{manuela.hummel@crg.es}}
%\note{}

\seealso{\code{\link{get.reads}}, \code{\link{reads2pairs}}, \code{\link{get.targets}}}

\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- paste(exptPath, "ExampleSet_Reads.bed", sep="/")
reads <- get.reads(readsfile, idcol=4, skip=0)
targetsfile <- paste(exptPath, "ExampleSet_Targets.bed", sep="/")
targets <- get.targets(targetsfile, skip=0)

## duplicates barplot for single reads
duplicates.barplot(reads, targets, returnDups=TRUE)

## duplicates barplot for read pairs
readpairs <- reads2pairs(reads)
duplicates.barplot(readpairs, targets, returnDups=TRUE)
}

\keyword{ hplot }
