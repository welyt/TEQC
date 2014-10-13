\name{coverage.hist}
\alias{coverage.hist}
\title{Coverage histogram}
\description{Histogram and cumulative density of target base coverages}

\usage{
coverage.hist(coverageTarget, col.hist = "lightblue", col.line = "orange", covthreshold, outline = FALSE, breaks = "Sturges", xlab, ylab, main, lwd, ...)
}

\arguments{
  \item{coverageTarget}{\code{\link[IRanges:Rle-class]{RleList}} containing \code{\link[IRanges:Rle-class]{Rle}}
        vectors of per-target-base coverages for each chromosome, i.e. \code{coverageTarget} output from \code{\link{coverage.target}}}
  \item{col.hist}{histogram color}
  \item{col.line}{color of the cumulative density line}
  \item{covthreshold}{indicates with dashed vertical and horizontal lines, which fraction of bases
         has a coverage of at least \code{covthreshold}; if missing, no dashed lines are drawn}
  \item{outline}{if \code{FALSE}, outliers (according to \code{\link{boxplot.stats}}) are removed before plotting.}
  \item{breaks}{number of cells for the histogram, or string naming an algorithm to compute
         the number of cells, or function to compute the number of cells,
         or vector giving the breakpoints between histogram cells (see \code{?hist}) but the
         latter option only with equidistant breakpoints}
  \item{xlab, ylab}{x- and y-axis labels}
  \item{main}{plot title}
  \item{lwd}{line width}
  \item{\dots}{further graphical parameters, passed to \code{plot(histogram)}}
}
%\details{}

\value{Histogram of read coverages for bases within the target. Additionally, a line and the right
      axis indicate the cumulative fraction of target bases with coverage of at least x.
      If option \code{covthreshold} is specified, red dashed lines highlight the cumulative fraction
      of target bases with at least the specified coverage.}

%\references{}
\author{Manuela Hummel \email{m.hummel@dkfz.de}}
%\note{}

\seealso{\code{\link{coverage.target}}, \code{\link{coverage.uniformity}}, \code{\link{coverage.density}}, \code{\link{coverage.plot}},
    \code{\link{coverage.targetlength.plot}}}
\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- file.path(exptPath, "ExampleSet_Reads.bed")
reads <- get.reads(readsfile, idcol=4, skip=0)
targetsfile <- file.path(exptPath, "ExampleSet_Targets.bed")
targets <- get.targets(targetsfile, skip=0)

## calculate per-base coverages
Coverage <- coverage.target(reads, targets, perBase=TRUE)

## coverage histogram
coverage.hist(Coverage$coverageTarget, covthreshold=8)
}

\keyword{ hplot }
