\name{coverage.plot}
\alias{coverage.plot}
\title{Coverage versus base position plot}
\description{Line plot of per-base coverages along a genomic region. Position of target regions can be shown.}
\usage{
coverage.plot(coverageAll, targets, chr, Start, End, Offset = 0, add = FALSE, col.line = 1, col.target = "orange", col.offset = "yellow", xlab, ylab, ylim, ...)
}

\arguments{
  \item{coverageAll}{\code{\link[IRanges:AtomicList-class]{RleList}} containing \code{\link[S4Vectors:Rle-class]{Rle}}
        vectors of per-base coverages for each chromosome, i.e. \code{coverageAll} output from \code{\link{coverage.target}}}
  \item{targets}{optional; \code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of target regions, i.e. output from \code{\link{get.targets}}; if missing no genomic regions are highlighted}
  \item{chr}{on which chromosome the region to plot is located (string, e.g. "chr1")}
  \item{Start}{genomic position where to start the plot}
  \item{End}{genomic position where to end the plot}
  \item{Offset}{integer; highlight \code{Offset} bases on both sides of each targeted region; defaults to 0}
  \item{add}{if \code{TRUE}, the coverage line of a new sample is added to an already existing plot}
  \item{col.line}{color of the coverage line}
  \item{col.target}{color of the bar indicating target regions}
  \item{col.offset}{color for highlighting \code{Offset} on the sides of target regions}
  \item{xlab, ylab}{x- and y-axis labels}
  \item{ylim}{y-axis coordinate ranges}
  \item{\dots}{further graphical parameters, passed to \code{plot}}
}
\details{If coverage of a new sample is added to an existing plot with \code{add = TRUE}, parameters
  \code{chr}, \code{Start}, \code{End} still have to be specified and should be the same as
  in the previous call in order to make sense. Parameters \code{targets} and \code{Offset} can but
  do not have to be given again. They can also differ from the previous ones, if for the additional sample
  a different target was captured.}

\value{Line plot showing per-base read coverages for a specified genomic region.
      When positions of target regions are provided, a bar on the bottom indicates their location
      such that coverage can be related to the captured targets.}

%\references{}
\author{Manuela Hummel \email{m.hummel@dkfz.de}}
%\note{}

\seealso{\code{\link{coverage.target}}, \code{\link{make.wigfiles}}, \code{\link{covered.k}},
    \code{\link{coverage.hist}}, \code{\link{coverage.uniformity}}, \code{\link{coverage.targetlength.plot}}}

\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- file.path(exptPath, "ExampleSet_Reads.bed")
reads <- get.reads(readsfile, idcol=4, skip=0)
targetsfile <- file.path(exptPath, "ExampleSet_Targets.bed")
targets <- get.targets(targetsfile, skip=0)

## calculate per-base coverages
Coverage <- coverage.target(reads, targets, perBase=TRUE)

## coverage plot
coverage.plot(Coverage$coverageAll, targets, Offset=100, chr="chr1", Start=11157524, End=11158764)
}

\keyword{ hplot }
