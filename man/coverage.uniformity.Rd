\name{coverage.uniformity}
\alias{coverage.uniformity}
\title{Coverage uniformity plot}
\description{Visualization of target coverage uniformity. A line shows the cumulative fraction
    of targeted bases that reach at least a certain normalized coverage.}

\usage{coverage.uniformity(coveragelist, addlines = TRUE, add = FALSE, xlab, ylab, xlim, ylim, col, lwd, ...)}

\arguments{
  \item{coveragelist}{output of function \code{\link{coverage.target}}, where option \code{perBase}
      had to be set to \code{TRUE}, i.e. a list with elements \code{coverageTarget} and \code{avgTargetCoverage}}
  \item{addlines}{if \code{TRUE}, dashed lines are added to the plot that indicate the fractions of bases
    achieving at least the average or at least half the average coverage}
  \item{add}{if \code{TRUE}, the coverage uniformity line of a new sample is added to an already existing plot}
  \item{xlab, ylab}{x- and y-axis labels}
  \item{xlim, ylim}{x- and y-axis coordinate ranges}
  \item{col}{line color}
  \item{lwd}{line width}
  \item{\dots}{further graphical parameters passed to \code{plot}}
}
\details{The function calculates normalized coverages: per-base coverages divided by
    average coverage over all targeted bases. Normalized coverages are not dependent
    on the absolute quantity of reads and are hence better comparable between different samples
    or even different experiments.}

\value{Line plot showing the fraction of targeted bases (y-axis) achieving a normalized
    coverage of at least x. The x-axis by default is truncated at 1, which corresponds to the average
    normalized coverage. The steeper the curve is falling, the less uniform is the coverage.
    If \code{addlines = TRUE}, dashed lines indicate the fractions of bases achieving at
    least the average (=1) or at least half (=0.5) the average coverage.}

\references{Gnirke A, Melnikov A, Maguire J, Rogov P, LeProust EM, Brockman W, Fennell T, Giannoukos G, Fisher S, Russ C, Gabriel S, Jaffe DB, Lander ES, Nusbaum C.
Solution hybrid selection with ultra-long oligonucleotides for massively parallel targeted sequencing. Nat Biotechnol. 2009; 27(2): 182-9.
}

\author{Manuela Hummel \email{manuela.hummel@crg.es}}
%\note{}

\seealso{\code{\link{coverage.target}}, \code{\link{covered.k}}, \code{\link{coverage.hist}}, \code{\link{coverage.density}},
    \code{\link{coverage.plot}}, \code{\link{coverage.targetlength.plot}}}

\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- paste(exptPath, "ExampleSet_Reads.bed", sep="/")
reads <- get.reads(readsfile, idcol=4, skip=0)
targetsfile <- paste(exptPath, "ExampleSet_Targets.bed", sep="/")
targets <- get.targets(targetsfile, skip=0)

## calculate per-base coverages
Coverage <- coverage.target(reads, targets, perBase=TRUE)

## coverage uniformity plot
coverage.uniformity(Coverage)
}

\keyword{ hplot }
