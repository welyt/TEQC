\name{chrom.barplot}
\alias{chrom.barplot}
\title{Reads per chromosome barplot}
\description{Barplot of numbers (or fractions) of reads (and targets) falling on each chromosome}

\usage{chrom.barplot(reads, targets, col = c("darkgreen", "orange"), ylab, legendpos = "topright", ...)}

\arguments{
  \item{reads}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing read positions,
      i.e. output from \code{\link{get.reads}}. To ensure a useful ordering of the bars, the
      chromosome information ('spaces' of \code{reads}) should be given as "chr" plus a
      number/letter [plus further specification], e.g. "chr1", "chrX", "chr17_ctg5_hap1", "chrUn_gl000211".}
  \item{targets}{Optional \code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of target regions, i.e. output from \code{\link{get.targets}}. The chromosome information should match
      the one of \code{reads}. If \code{targets} is missing, only numbers of reads will be displayed.}
  \item{col}{color(s) of the bars}
  \item{ylab}{y-axis label}
  \item{legendpos}{Position of the legend. String from the list "bottomright", "bottom", "bottomleft",
     "left", "topleft", "top", "topright", "right" and "center". Ignored if \code{targets} is missing.}
  \item{\dots}{graphical parameters passed to \code{barplot}}
}
\details{If \code{targets} is not specified, absolute read counts per chromosome are shown in the barplot.
  If \code{targets} is provided, fractions of reads and targets are shown. For reads, this is the
  fraction within the total \emph{number} of reads (since reads are expected to have all the same length).
  In contrast, for the targets, the fraction of targeted bases on each chromosome is calculated.
  Since targets might vary in length it is reasonable to account for the actual target \emph{sizes} instead
  of considering merely numbers of targets per chromosome.}

\value{Barplot of reads and optionally targets per chromosome.}

%\references{}
\author{Manuela Hummel \email{manuela.hummel@crg.es}}
%\note{}

\seealso{\code{\link{get.reads}}}

\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- paste(exptPath, "ExampleSet_Reads.bed", sep="/")
reads <- get.reads(readsfile, idcol=4, skip=0)
targetsfile <- paste(exptPath, "ExampleSet_Targets.bed", sep="/")
targets <- get.targets(targetsfile, skip=0)

chrom.barplot(reads, targets)
}

\keyword{ hplot }
