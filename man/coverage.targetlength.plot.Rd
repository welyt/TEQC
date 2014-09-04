\name{coverage.targetlength.plot}
\alias{coverage.targetlength.plot}
\title{Number of reads or average coverage versus target length scatter plot}
\description{Plots either numbers of on-target reads or average per-target coverage
    (or potentially other per-target values) against respective target lengths.
    A smoothing spline is added to the scatter plot.}

\usage{
coverage.targetlength.plot(targets, plotcolumn, linecol = 2, xlab, ylab, lwd, pch, cex, ...)
}

\arguments{
  \item{targets}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of target regions and further 'values' columns that should be plotted,
      i.e. output from \code{\link{coverage.target}} or \code{\link{readsPerTarget}}}
  \item{plotcolumn}{name or index of column to plot (of the 'values' DataFrame within \code{targets})}
  \item{linecol}{color of spline curve}
  \item{xlab, ylab}{x- and y-axis labels}
  \item{lwd}{line width of spline curve}
  \item{pch}{plotting character}
  \item{cex}{size of plotting character}
  \item{\dots}{further graphical parameters, passed to \code{plot}}
}
\details{\code{\link{coverage.target}} and \code{\link{readsPerTarget}} can be used to calculate
    average per-target coverages and numbers of reads overlapping each target. The values are
    added to the \code{RangedData} table containing the target positions. Such \code{RangedData}
    table can then be used for plotting the calculated values against the respecitve target lengths.}

\value{A scatterplot with the given per-target values on the y-axis and corresponding target lengths on the x-axis.
    A smoothing spline is added to the plot.}

%\references{}
\author{Manuela Hummel \email{m.hummel@dkfz.de}}
%\note{}

\seealso{\code{\link{coverage.target}}, \code{\link{readsPerTarget}}, \code{\link{covered.k}}, \code{\link{coverage.hist}},
    \code{\link{coverage.uniformity}}, \code{\link{coverage.plot}}}

\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- file.path(exptPath, "ExampleSet_Reads.bed")
reads <- get.reads(readsfile, idcol=4, skip=0)
targetsfile <- file.path(exptPath, "ExampleSet_Targets.bed")
targets <- get.targets(targetsfile, skip=0)

## get average per-target coverage
Coverage <- coverage.target(reads, targets, perTarget=TRUE)
targets2 <- Coverage$targetCoverages

## get numbers of reads per target
targets2 <- readsPerTarget(reads, targets2)

## coverage vs target length
coverage.targetlength.plot(targets2, plotcolumn="avgCoverage", pch="o")

## coverage vs number of reads per target
coverage.targetlength.plot(targets2, plotcolumn="nReads", pch="o")
}

\keyword{ hplot }
