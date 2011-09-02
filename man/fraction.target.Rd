\name{fraction.target}
\alias{fraction.target}
\title{Fraction of the target within the genome}
\description{Calculates the fraction of the reference genome that is targeted}

\usage{fraction.target(targets, Offset = 0, genome = c(NA, "hg19", "hg18"), genomesize)}

\arguments{
  \item{targets}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of target regions, i.e. output from \code{\link{get.targets}}}
  \item{Offset}{integer; add \code{Offset} bases on both sides to targeted regions and potentially
      collapse resulting overlapping target regions}
  \item{genome}{genome version targets were designed and reads aligned to.
      For the given options the total genome size is set automatically. For other genomes or versions,
      leave this option empty ('NA') and specify the genome size with option 'genomesize'}
  \item{genomesize}{integer: specify the total genome size manually. If 'genomesize' is given,
      option 'genome' will be ignored.}
}
%\details{}

\value{Returns the fraction of nucleotides within the genome that were targeted.}

%\references{}
\author{Manuela Hummel \email{manuela.hummel@crg.es}}

\note{With the output from \code{\link{fraction.target}} and \code{\link{fraction.reads.target}}
    the 'enrichment' of the target capture experiment can be calculated as
    'fraction of on-target reads / fraction of target within genome'}

\seealso{\code{\link{fraction.reads.target}}, \code{\link{get.targets}}}

\examples{
exptPath <- system.file("extdata", package="TEQC")
targetsfile <- file.path(exptPath, "ExampleSet_Targets.bed")
targets <- get.targets(targetsfile, skip=0)
fraction.target(targets, genome="hg19")
}

\keyword{ math }
