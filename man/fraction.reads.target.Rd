\name{fraction.reads.target}
\alias{fraction.reads.target}
\title{Target capture specificity}
\description{Calculates the fraction of reads that align to target regions.
    Can also be used to retrieve those reads mapping to targets.}

\usage{
fraction.reads.target(reads, targets, Offset = 0, mappingReads = FALSE)
}

\arguments{
  \item{reads}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of sequenced reads, i.e. output of \code{\link{get.reads}}. Alternatively, for paired-end data, it can be the output
      of \code{\link{reads2pairs}} when fraction of on-target read \emph{pairs} shall be calculated instead
      of fraction of single on-target reads.}
  \item{targets}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of target regions, i.e. output from \code{\link{get.targets}}}
  \item{Offset}{integer; add \code{Offset} bases on both sides to targeted regions and potentially
      collapse resulting overlapping target regions}
  \item{mappingReads}{if \code{TRUE}, reduced \code{RangedData} table with only those reads
      mapping to target regions is returned. When \code{reads} is output of \code{\link{reads2pairs}},
      \code{mappingReads} will be the corresponding subset of on-target read pairs.}
}
%\details{}

\value{If \code{mappingReads} equals \code{FALSE}, just the fraction of reads / read pairs mapping to targets is returned.
    When \code{reads} contains all single reads (i.e. is output of \code{\link{get.reads}}), this is the number of target-overlapping reads,
    divided by the number of all single reads. When \code{reads} contains read pairs (i.e. is output of \code{\link{reads2pairs}}),
    it is the number of read pairs with at least one target-overlapping read, divided by the
    number of read pairs (= half the number of reads). In case of small targets and large insert sizes
    the two reads of a pair could be located on both sides of the target without overlap, respectively.
    Still, the read pair will be counted as on-target, since the corresponding DNA molecule was covering the target.

    If \code{mappingReads} equals \code{TRUE}, a list is returned with elements
  \item{onTargetFraction}{fraction of reads / read pairs mapping to targets}
  \item{mappingReads}{\code{RangedData} table containing positions of the reads / read pairs mapping to target regions}
}
%\references{}
\author{Manuela Hummel \email{manuela.hummel@crg.es}}

\note{With the output from \code{\link{fraction.target}} and \code{\link{fraction.reads.target}}
    the 'enrichment' of the target capture experiment can be calculated as
    'fraction of on-target reads / fraction of target within genome'}

\seealso{\code{\link{fraction.target}}, \code{\link{get.reads}}, \code{\link{reads2pairs}}, \code{\link{get.targets}}}

\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- paste(exptPath, "ExampleSet_Reads.bed", sep="/")
reads <- get.reads(readsfile, idcol=4, skip=0)
targetsfile <- paste(exptPath, "ExampleSet_Targets.bed", sep="/")
targets <- get.targets(targetsfile, skip=0)

## fraction of on-target reads
fraction.reads.target(reads, targets)
}

\keyword{ math }
