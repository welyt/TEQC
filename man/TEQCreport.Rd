\name{TEQCreport}
\alias{TEQCreport}
\title{Creates an html report}
\description{Creates an automated html report for the complete TEQC analysis of one sample}
\usage{
TEQCreport(sampleName = "", targetsName = "", referenceName = "", destDir = "TEQCreport",
            reads = get.reads(), targets = get.targets(), Offset = 0, pairedend = FALSE, genome = c(NA, "hg19", "hg18"),
            genomesize, k=c(1, 2, 3, 5, 10, 20), CovUniformityPlot = FALSE, CovTargetLengthPlot = FALSE, CovGCPlot = FALSE,
            duplicatesPlot = FALSE, baits = get.baits(), WigFiles = FALSE, saveWorkspace = FALSE)
}

\arguments{
  \item{sampleName}{descriptive sample name; will be written on top of the html report}
  \item{targetsName}{descriptive name of the captured target; will be written on top of the html report}
  \item{referenceName}{descriptive name of the reference genome the reads were aligned against; will be written on top of the html report}
  \item{destDir}{directory where results and html documents shall be saved}
  \item{reads}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of sequenced reads, or call to \code{\link{get.reads}} to read in positions from a bed or BAM file}
  \item{targets}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing positions
      of target regions, or call to \code{\link{get.targets}} to read in positions from a bed file}
  \item{Offset}{integer; add \code{Offset} bases on both sides to targeted regions and potentially
      collapse resulting overlapping target regions}
  \item{pairedend}{if \code{TRUE}, data will be considered to be paired-end data, i.e. reads will be "merged" to read
      pairs, and chromosome bar plot, specificity, enrichment and duplicate analysis (if selected) will be
      based on read pairs rather than on single reads}
  \item{genome}{genome version targets were designed and reads aligned to.
      For the given options the total genome size is set automatically. For other genomes or versions,
      leave this option empty ('NA') and specify the genome size with option 'genomesize'}
  \item{genomesize}{integer: specify the total genome size manually. If 'genomesize' is given,
      option 'genome' will be ignored.}
  \item{k}{integer vector of \code{k}-values for which to show fraction
        of target bases with coverage >= \code{k} in 'Sensitivity' table. Passed to \code{\link{covered.k}}}
  \item{CovUniformityPlot}{if \code{TRUE}, a coverage uniformity plot is created,
        see \code{\link{coverage.uniformity}}}
  \item{CovTargetLengthPlot}{if \code{TRUE}, coverage vs target length plots are created,
        see \code{\link{coverage.targetlength.plot}} }
  \item{CovGCPlot}{if \code{TRUE}, a coverage vs GC content plot is created,
        see \code{\link{coverage.GC}}}
  \item{duplicatesPlot}{if \code{TRUE}, a duplicates barplot is created,
        see \code{\link{duplicates.barplot}}}
  \item{baits}{A \code{\link[IRanges:RangedData-class]{RangedData}} table holding the hybridization probe
        ("bait") positions and sequences, or call to \code{\link{get.baits}} to read in
        positions from a bed file. Only needed if \code{CovGCPlot = TRUE}.}
  \item{WigFiles}{if \code{TRUE}, wiggle files with per-base coverage are created for each chromosome}
  \item{saveWorkspace}{if \code{TRUE}, an R workspace with objects \code{reads}, \code{targets} and
        output of \code{\link{coverage.target}} and \code{\link{reads2pairs}}
        (in case \code{pairedend = TRUE}) are saved in \code{destDir} to be available for further analyses}
}


\details{TEQC analysis is preformed and files for an html report are created in \code{destDir}.
        The report can be viewed by opening \code{destDir}/index.html in a web browser. Images are saved in
        \code{destDir}/image. Wiggle files (in case \code{WigFiles = TRUE}) are saved in
        \code{destDir}/wiggle. The table with coverage values per target and the R workspace
        containing R objects for potential further analysis (in case \code{saveWorkspace = TRUE})
        are saved in \code{destDir}.}

\value{The function is invoked for its side effect}

\references{Hummel M, Bonnin S, Lowy E, Roma G. TEQC: an R-package for quality control in
    target capture experiments. Bioinformatics 2011; doi: 10.1093/bioinformatics/btr122
}

\author{Manuela Hummel \email{manuela.hummel@crg.es}}
%\note{}

\seealso{\code{\link{get.reads}}, \code{\link{get.targets}}, \code{\link{fraction.target}}, \code{\link{fraction.reads.target}},
    \code{\link{coverage.target}}, \code{\link{readsPerTarget}}, \code{\link{reads2pairs}},
    \code{\link{covered.k}}, \code{\link{coverage.hist}}, \code{\link{coverage.uniformity}},
    \code{\link{coverage.targetlength.plot}}, \code{\link{coverage.GC}}, \code{\link{get.baits}},
    \code{\link{make.wigfiles}}}

\examples{
## get reads and targets files
exptPath <- system.file("extdata", package="TEQC")
readsfile <- file.path(exptPath, "ExampleSet_Reads.bed")
targetsfile <- file.path(exptPath, "ExampleSet_Targets.bed")

## create report
\dontrun{
TEQCreport(sampleName="Test Sample", targetsName="Human Exome", referenceName="Human Genome",
          destDir="report", reads=get.reads(readsfile, skip=0, idcol=4),
          targets=get.targets(targetsfile, skip=0), genome="hg19")}
}

\keyword{ misc }
