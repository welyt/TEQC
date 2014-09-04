\name{insert.size.hist}
\alias{insert.size.hist}
\title{Insert sizes histogram}
\description{Computes read pair insert sizes, i.e. distance from first base of first read
    to last base of second read of a read pair, and plots a histogram for all insert sizes.}

\usage{
insert.size.hist(readpairs, returnInserts = FALSE, legendpos="topleft", outline=FALSE, main, xlab, ylab, breaks, col, ...)
}

\arguments{
  \item{readpairs}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of read pairs, i.e. output of \code{\link{reads2pairs}} (or the element \code{readpairs} from the
      \code{\link{reads2pairs}} output in case single reads without matching pair were found).}
  \item{returnInserts}{if \code{TRUE}, the vector of read pair insert sizes is returned}
  \item{legendpos}{position of the legend, e.g. 'topleft' or 'topright'}
  \item{outline}{if \code{FALSE}, outliers (according to \code{\link{boxplot.stats}}) are removed before plotting.
  If \code{returnInserts=TRUE}, those outliers are still included in the output.}
  \item{main}{plot title}
  \item{xlab, ylab}{x- and y-axis labels}
  \item{breaks}{e.g. integer specifying the number of cells for the histogram, see \code{?hist}}
  \item{col}{histogram color}
  \item{\dots}{further graphical parameters passed to \code{hist}}
}
%\details{}

\value{Histogram of read pair insert sizes. Average, standard deviation and median insert size
    are given in the legend and indicated by lines.

    If \code{returnInserts = TRUE}, a named vector of insert sizes is returned.}

%\references{}
\author{Manuela Hummel \email{m.hummel@dkfz.de}}

%\note{}

\seealso{\code{\link{get.reads}}, \code{\link{reads2pairs}}, \code{\link{duplicates.barplot}}}

\examples{
## get reads
exptPath <- system.file("extdata", package="TEQC")
readsfile <- file.path(exptPath, "ExampleSet_Reads.bed")
reads <- get.reads(readsfile, idcol=4, skip=0)

## merge to read pairs
readpairs <- reads2pairs(reads)

## insert size histogram
insert.size.hist(readpairs, breaks=10)
}

\keyword{ hplot }
