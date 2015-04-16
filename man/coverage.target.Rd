\name{coverage.target}
\alias{coverage.target}
\title{Calculates read coverage}
\description{Calculates average coverage over all target bases, average coverage for each target separately,
      and per-base coverage for all and for targeted bases}
\usage{
coverage.target(reads, targets, Offset = 0, perTarget = TRUE, perBase = TRUE)
}

\arguments{
  \item{reads}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of sequenced reads, i.e. output from \code{\link{get.reads}}}
  \item{targets}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of target regions, i.e. output from \code{\link{get.targets}}}
  \item{Offset}{integer; add \code{Offset} bases on both sides to targeted regions and potentially
      collapse resulting overlapping target regions}
  \item{perTarget}{if TRUE, coverage average and standard deviation per target are calculated and returned}
  \item{perBase}{if TRUE, the per-base coverages i) only for targeted bases and
     ii) for all sequenced and/or targeted bases, are returned}
}
%\details{}

\value{A list is returned with elements
  \item{avgTargetCoverage}{average coverage over all target bases}
  \item{targetCoverageSD}{standard deviation of coverage of all target bases}
  \item{targetCoverageQuantiles}{0\% (minium), 25\%, 50\% (median), 75\% and 100\% (maximum) quantiles of coverage of all target bases}
  \item{targetCoverages}{Input \code{RangedData} table \code{targets} with two additional 'values'
      columns \code{avgCoverage} and \code{coverageSD}. The former contains the average coverage for each target,
      the latter the respective coverage standard deviation. Only returned if \code{perTarget} equals \code{TRUE}.}
  \item{coverageAll}{\code{\link[IRanges:AtomicList-class]{RleList}} containing a \code{\link[S4Vectors:Rle-class]{Rle}}
        vector for each chromosome with coverages for all bases that are sequenced and/or within a targeted;
        only returned if \code{perBase} equals \code{TRUE}}
  \item{coverageTarget}{\code{\link[IRanges:AtomicList-class]{RleList}} containing a \code{\link[S4Vectors:Rle-class]{Rle}}
        vector for each chromosome with coverages for target bases only;
        only returned if \code{perBase} equals \code{TRUE}}
}

%\references{}

\author{Manuela Hummel \email{m.hummel@dkfz.de}}
%\note{}

\seealso{\code{\link{covered.k}}, \code{\link{coverage.hist}}, \code{\link{coverage.uniformity}}, \code{\link{coverage.plot}},
    \code{\link{coverage.targetlength.plot}}}

\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- file.path(exptPath, "ExampleSet_Reads.bed")
reads <- get.reads(readsfile, idcol=4, skip=0)
targetsfile <- file.path(exptPath, "ExampleSet_Targets.bed")
targets <- get.targets(targetsfile, skip=0)

## total average, per-base and per-target coverages
Coverage <- coverage.target(reads, targets)
}

\keyword{ math }
