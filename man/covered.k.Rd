\name{covered.k}
\alias{covered.k}
\title{Target capture sensitivity}
\description{Calculates fraction of target bases covered by at least k reads}

\usage{
covered.k(coverageTarget, k = c(1, 2, 3, 5, 10, 20))
}

\arguments{
  \item{coverageTarget}{\code{\link[IRanges:Rle-class]{RleList}} containing \code{\link[IRanges:Rle-class]{Rle}}
        vectors of per-target-base coverages for each chromosome, i.e. \code{coverageTarget} output from \code{\link{coverage.target}}}
  \item{k}{integer vector of \code{k}-values for which to show fraction of target bases with coverage >= \code{k}}
}
%\details{}

\value{Named vector of same length as \code{k} giving the corresponding fractions of target bases
      achieving coverages >= \code{k}}

%\references{}
\author{Manuela Hummel \email{manuela.hummel@crg.es}}
%\note{}

\seealso{\code{\link{coverage.target}}, \code{\link{coverage.hist}}, \code{\link{coverage.uniformity}},
    \code{\link{coverage.plot}}, \code{\link{coverage.targetlength.plot}}}

\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- file.path(exptPath, "ExampleSet_Reads.bed")
reads <- get.reads(readsfile, idcol=4, skip=0)
targetsfile <- file.path(exptPath, "ExampleSet_Targets.bed")
targets <- get.targets(targetsfile, skip=0)

## calculate per-base coverages
Coverage <- coverage.target(reads, targets, perBase=TRUE)

covered.k(Coverage$coverageTarget, k=c(1,10,20))
}

\keyword{ math }
