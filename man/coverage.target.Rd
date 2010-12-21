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
  \item{perTarget}{if TRUE, average coverages per target are calculated and returned}
  \item{perBase}{if TRUE, the per-base coverages i) only for targeted bases and
     ii) for all sequenced and/or targeted bases, are returned}
}
%\details{}

\value{If both \code{perTarget} and \code{perBase} are \code{FALSE}, only the average
      coverage over all target bases is returned.
      
      If at least one of \code{perTarget} and \code{perBase} are \code{TRUE} a list is
      returned with elements
  \item{avgTargetCoverage}{average coverage over all target bases}
  \item{targetCoverages}{input \code{RangedData} table \code{targets} with an additional 'values'
      column containing average coverages for each target; only returned if \code{perTarget} equals \code{TRUE}}
  \item{coverageAll}{\code{\link[IRanges:Rle-class]{RleList}} containing a \code{\link[IRanges:Rle-class]{Rle}}
        vector for each chromosome with coverages for all bases that are sequenced and/or within a targeted;
        only returned if \code{perBase} equals \code{TRUE}}
  \item{coverageTarget}{\code{\link[IRanges:Rle-class]{RleList}} containing a \code{\link[IRanges:Rle-class]{Rle}}
        vector for each chromosome with coverages for target bases only;
        only returned if \code{perBase} equals \code{TRUE}}
}

%\references{}

\author{Manuela Hummel \email{manuela.hummel@crg.es}}
%\note{}

\seealso{\code{\link{covered.k}}, \code{\link{coverage.hist}}, \code{\link{coverage.uniformity}}, \code{\link{coverage.plot}},
    \code{\link{coverage.targetlength.plot}}}

\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- paste(exptPath, "ExampleSet_Reads.bed", sep="/")
reads <- get.reads(readsfile, idcol=4, skip=0)
targetsfile <- paste(exptPath, "ExampleSet_Targets.bed", sep="/")
targets <- get.targets(targetsfile, skip=0)

## total average, per-base and per-target coverages
Coverage <- coverage.target(reads, targets)
}

\keyword{ math }
