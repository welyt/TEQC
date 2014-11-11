\name{coverage.GC}
\alias{coverage.GC}
\title{Bait coverage versus GC content plot}
\description{Calculates and plots average normalized coverage per hybridization probe versus GC content
    of the respective probe. A smoothing spline is added to the scatter plot.}

\usage{coverage.GC(coverageAll, baits, returnBaitValues = FALSE, linecol = "darkred", lwd, xlab, ylab, pch, col, cex, ...)}

\arguments{
  \item{coverageAll}{\code{\link[IRanges:AtomicList-class]{RleList}} containing \code{\link[S4Vectors:Rle-class]{Rle}}
        vectors of per-base coverages for each chromosome, i.e. \code{coverageAll} output of \code{\link{coverage.target}}}
  \item{baits}{A \code{\link[IRanges:RangedData-class]{RangedData}} table holding the hybridization probe
  ("bait") positions and sequences, i.e. output of\code{\link{get.baits}}}
  \item{returnBaitValues}{if \code{TRUE}, average coverage, average normalized coverage and GC content per bait are returned}
  \item{linecol, lwd}{color and width of spline curve}
  \item{xlab, ylab}{x- and y-axis labels}
  \item{pch}{plotting character}
  \item{col, cex}{color and size of plotting character}
  \item{\dots}{further graphical parameters passed to \code{plot}}
}
\details{The function calculates average normalized coverages for each bait: the average coverage
    over all bases within a bait is divided by the
    average coverage over all bait-covered bases. Normalized coverages are not dependent
    on the absolute quantity of reads and are hence better comparable between different samples
    or even different experiments.}

\value{A scatterplot with normalized per-bait coverages on the y-axis and
    GC content of respective baits on the x-axis. A smoothing spline is added to the plot.

    If \code{returnBaitValues = TRUE} average coverage, average normalized coverage and GC content per bait are returned
    as 'values' columns of the \code{baits} input \code{RangedData} table}

\references{Tewhey R, Nakano M, Wang X, Pabon-Pena C, Novak B, Giuffre A, Lin E, Happe S, Roberts DN, LeProust EM, Topol EJ, Harismendy O, Frazer KA.
Enrichment of sequencing targets from the human genome by solution hybridization.
Genome Biol. 2009; 10(10): R116.}

\author{Manuela Hummel \email{m.hummel@dkfz.de}}
%\note{}

\seealso{\code{\link{coverage.target}}, \code{\link{covered.k}}, \code{\link{coverage.hist}}, \code{\link{coverage.plot}},
    \code{\link{coverage.uniformity}}, \code{\link{coverage.targetlength.plot}}}

\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- file.path(exptPath, "ExampleSet_Reads.bed")
reads <- get.reads(readsfile, idcol=4, skip=0)
targetsfile <- file.path(exptPath, "ExampleSet_Targets.bed")
targets <- get.targets(targetsfile, skip=0)

## calculate per-base coverages
Coverage <- coverage.target(reads, targets, perBase=TRUE)

## get bait positions and sequences
baitsfile <- file.path(exptPath, "ExampleSet_Baits.txt")
baits <- get.baits(baitsfile, chrcol=3, startcol=4, endcol=5, seqcol=2)

## do coverage vs GC plot
coverage.GC(Coverage$coverageAll, baits)
}

\keyword{ hplot }
