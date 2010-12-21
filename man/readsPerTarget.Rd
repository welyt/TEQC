\name{readsPerTarget}
\alias{readsPerTarget}
\title{Numbers of reads per target}
\description{Counts the numbers of reads overlapping each target region}

\usage{
readsPerTarget(reads, targets, Offset = 0)
}

\arguments{
  \item{reads}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of sequenced reads, i.e. output from \code{\link{get.reads}}}
  \item{targets}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of target regions, i.e. output from \code{\link{get.targets}}}
  \item{Offset}{integer; add \code{Offset} bases on both sides to targeted regions and potentially
      collapse resulting overlapping target regions}
}
%\details{}

\value{The input \code{RangedData} table \code{targets} with an additional 'values'
      column containing numbers of reads overlapping each target}

%\references{}
\author{Manuela Hummel \email{manuela.hummel@crg.es}}

\note{As \code{reads} input also the \code{mappingReads} output of function \code{\link{fraction.reads.target}}
    can be used to speed up calculation. In this case, make sure that \code{targets} and \code{Offset} parameters were the
    same in \code{fraction.reads.target} as then specified in \code{readsPerTarget}.}

\seealso{\code{\link{coverage.target}}, \code{\link{fraction.reads.target}}, \code{\link{covered.k}}, \code{\link{coverage.hist}},
    \code{\link{coverage.uniformity}}, \code{\link{coverage.plot}}, \code{\link{coverage.targetlength.plot}}}

\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- paste(exptPath, "ExampleSet_Reads.bed", sep="/")
reads <- get.reads(readsfile, idcol=4, skip=0)
targetsfile <- paste(exptPath, "ExampleSet_Targets.bed", sep="/")
targets <- get.targets(targetsfile, skip=0)

## number of reads per target
readsPerTarget(reads, targets)
}

\keyword{ math }
