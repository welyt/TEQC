\name{coverage.density}
\alias{coverage.density}
\title{Coverage density plot}
\description{Visualization of target coverage density for one or more samples.}

\usage{coverage.density(coveragelist, normalized = TRUE, legend, main, xlab, col, lwd, lty, xlim, ylim, ...)}

\arguments{
  \item{coveragelist}{Output of function \code{\link{coverage.target}}, where option \code{perBase}
      had to be set to \code{TRUE}, i.e. a list with elements \code{coverageTarget} and \code{avgTargetCoverage}.
      Or, when density of several samples shall be visualized, a list with respective outputs of \code{\link{coverage.target}}.}
  \item{normalized}{if \code{TRUE}, densities of normalized coverages will be shown; original coverages otherwise}
  \item{legend}{legend text. If missing, names of \code{coveragelist} will be taken.
      If \code{NULL}, no legend will be drawn.}
  \item{main}{main title}
  \item{xlab}{x-axis label}
  \item{col}{line color(s)}
  \item{lwd}{line width(s)}
  \item{lty}{line style(s)}
  \item{xlim, ylim}{x- and y-axis coordinate ranges}
  \item{\dots}{further graphical parameters passed to \code{plot}}
}
\details{If \code{normalized = TRUE}, the function calculates normalized coverages: per-base coverages divided by
    average coverage over all targeted bases. Normalized coverages are not dependent
    on the absolute quantity of reads and are hence better comparable between different samples
    or even different experiments.}

\value{Line plot(s) showing densities.}

%\references{}

\author{Manuela Hummel \email{m.hummel@dkfz.de}}
%\note{}

\seealso{\code{\link{coverage.target}}, \code{\link{covered.k}}, \code{\link{coverage.hist}},
    \code{\link{coverage.uniformity}}, \code{\link{coverage.correlation}}, \code{\link{coverage.plot}}}

\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- file.path(exptPath, "ExampleSet_Reads.bed")
reads <- get.reads(readsfile, idcol=4, skip=0)
targetsfile <- file.path(exptPath, "ExampleSet_Targets.bed")
targets <- get.targets(targetsfile, skip=0)

## calculate per-base coverages
Coverage <- coverage.target(reads, targets, perBase=TRUE)

## coverage density
coverage.density(Coverage)
}

\keyword{ hplot }
