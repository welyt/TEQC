\name{make.wigfiles}
\alias{make.wigfiles}
\title{Creates wiggle files with per-base coverages}
\description{Prepares wiggle files with (non-zero) per-base coverages for the upload
    and visualization with genome browsers}

\usage{
make.wigfiles(coverageAll, chroms, trackname = "Coverage", filename = "Coverage")
}

\arguments{
  \item{coverageAll}{\code{\link[IRanges:Rle-class]{RleList}} containing \code{\link[IRanges:Rle-class]{Rle}}
        vectors of per-base coverages for each chromosome, i.e. \code{coverageAll} output of \code{\link{coverage.target}}}
  \item{chroms}{vector of chromosome names for which to produce wiggle files; if missing
      wiggle files will be produced for all chromosomes on which there are reads}
  \item{trackname}{trackname for wiggle file header}
  \item{filename}{part of output wiggle file name. Respective chromosome number and '.wig' will be added}
}
\details{Only non-zero coverages will be listed}

\value{One or more wiggle files listing per-base (non-zero) read coverages}

%\references{}
\author{Manuela Hummel \email{manuela.hummel@crg.es}}
%\note{}

\seealso{\code{\link{coverage.target}}, \code{\link{coverage.plot}}, \code{\link{covered.k}},
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

## create wiggle files for read coverages on chromsomes 13 and 17
make.wigfiles(Coverage$coverageAll, chroms=c("chr13", "chr17"))
}

\keyword{ file }
