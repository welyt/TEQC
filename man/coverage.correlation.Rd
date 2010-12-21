\name{coverage.correlation}
\alias{coverage.correlation}
\title{Coverage correlation plot}
\description{Visualization of target coverage correlations between pairs of samples.}

\usage{coverage.correlation(coveragelist, normalized = TRUE, plotfrac = 0.001, seed = 123, labels, main, pch = ".", cex.labels,
                  cex.pch = 2, cex.main = 1.2, cex.corr, font.labels = 1, font.main = 2, ...)}

\arguments{
  \item{coveragelist}{List where each element is the output of function \code{\link{coverage.target}}, where option \code{perBase}
      had to be set to \code{TRUE}.}
  \item{normalized}{if \code{TRUE}, correlation of normalized target coverages will be shown; original coverages otherwise}
  \item{plotfrac}{numeric value between 0 and 1. Coverages for a fraction of \code{plotfrac}
      of all target bases are shown.}
  \item{seed}{seed for random selection of \code{plotfrac} bases}
  \item{labels}{sample names that are written in the diagonal panels; if missing, names of \code{coveragelist}
      are taken; if those are \code{NULL}, "sample 'i'" is shown}
  \item{main}{main title}
  \item{pch}{plot symbol for the scatter plots}
  \item{cex.labels, cex.pch, cex.main}{sizes of sample labels, plot symbols, main title}
  \item{cex.corr}{size of the correlation values; if missing, sizes are made proportionally to the values
      of (positive) correlation.}
  \item{font.labels, font.main}{fonts for sample labels and main title}
  \item{\dots}{further graphical parameters, e.g. limits and symbol color for the scatter plots}
}
\details{If \code{normalized = TRUE}, the function calculates normalized coverages: per-base coverages divided by
    average coverage over all targeted bases. Normalized coverages are not dependent
    on the absolute quantity of reads and are hence better comparable between different samples
    or even different experiments.}

\value{'pairs'-style plot where upper panels show scatter plot of (a randomly chosen fraction of)
    coverage values for pairs of samples. The lower panels show the respective Pearson correlation
    coefficients, calculated using all coverage values (even if not all of them are shown in the
    scatter plot).}

%\references{}

\author{Manuela Hummel \email{manuela.hummel@crg.es}}
%\note{}

\seealso{\code{\link{coverage.target}}, \code{\link{covered.k}}, \code{\link{coverage.hist}},
    \code{\link{coverage.density}}, \code{\link{coverage.uniformity}}, \code{\link{coverage.plot}}}

\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- paste(exptPath, "ExampleSet_Reads.bed", sep="/")
reads <- get.reads(readsfile, idcol=4, skip=0)
targetsfile <- paste(exptPath, "ExampleSet_Targets.bed", sep="/")
targets <- get.targets(targetsfile, skip=0)

## calculate per-base coverages
Coverage <- coverage.target(reads, targets, perBase=TRUE)

## simulate another sample
r <- sample(nrow(reads), 0.1 * nrow(reads))
reads2 <- reads[-r,,drop=TRUE]
Coverage2 <- coverage.target(reads2, targets, perBase=TRUE)

## coverage uniformity plot
covlist <- list(Coverage, Coverage2)
coverage.correlation(covlist, plotfrac=0.1)
}

\keyword{ hplot }
