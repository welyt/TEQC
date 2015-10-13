\name{multiTEQCreport}
\alias{multiTEQCreport}
\title{Creates an html report for multiple samples}
\description{Creates an automated html report comparing TEQC analysis results of several samples}
\usage{
multiTEQCreport(singleReportDirs, samplenames, projectName = "", targetsName = "",
          referenceName = "", destDir = "multiTEQCreport", k = c(1, 2, 3, 5, 10, 20), figureFormat = c("jpeg", "png", "tiff"))
}

\arguments{
  \item{singleReportDirs}{string of directory names; output directories of function TEQCreport(), launched beforehand for each single sample}
  \item{samplenames}{names of the samples that will be used in tables and figures}
  \item{projectName}{descriptive name for the project / collection of samples; will be written on top of the html report}
  \item{targetsName}{descriptive name of the captured target; will be written on top of the html report}
  \item{referenceName}{descriptive name of the reference genome the reads were aligned against; will be written on top of the html report}
  \item{destDir}{directory where results and html documents shall be saved}
  \item{k}{integer vector of \code{k}-values for which to show fraction
        of target bases with coverage >= \code{k} in section 'Sensitivity'}
  \item{figureFormat}{format of the figures produced for the html report (besides pdf graphs)}
}

\details{Before creating the html report for multiple samples, \code{TEQCreport} has to be run for each of the samples separately. 
  The output directories of those analyses are the input for \code{multiTEQCreport}. While the creation of single-sample reports is time and memory intensive, \code{multiTEQCreport} finishes quickly, since it just collects and summarizes the results from the single analyses.  
  The files for the multiple sample html report are created in \code{destDir}.
        The report can be viewed by opening \code{destDir}/index.html in a web browser. Images are saved in
        \code{destDir}/image.}

\value{The function is invoked for its side effect}

\references{Hummel M, Bonnin S, Lowy E, Roma G. TEQC: an R-package for quality control in
    target capture experiments. Bioinformatics 2011; 27(9):1316-7.
}

\author{Manuela Hummel \email{m.hummel@dkfz.de}}
%\note{}

\seealso{\code{\link{TEQCreport}}}

\examples{
## get reads and targets
exptPath <- system.file("extdata", package="TEQC")
readsfile <- file.path(exptPath, "ExampleSet_Reads.bed")
reads <- get.reads(readsfile, skip=0, idcol=4)
targetsfile <- file.path(exptPath, "ExampleSet_Targets.bed")
targets <- get.targets(targetsfile, skip=0) 

## simulated 2nd sample
r <- sample(nrow(reads), 0.1 * nrow(reads))
reads2 <- reads[-r,,drop=TRUE]

\dontrun{
## create single-sample reports
TEQCreport(sampleName="Test Sample A", targetsName="Human Exome", referenceName="Human Genome",
          destDir="./reportA", reads=reads, targets=targets, genome="hg19")
TEQCreport(sampleName="Test Sample B", targetsName="Human Exome", referenceName="Human Genome",
          destDir="./reportB", reads=reads2, targets=targets, genome="hg19")
          
## create multi-sample report
multiTEQCreport(singleReportDirs=c("./reportA", "./reportB"), samplenames=c("Sample A","Sample B"), 
          projectName="Test Project", targetsName="Human Exome", referenceName="Human Genome", 
          destDir="./multiTEQCreport")
          }
}

\keyword{ misc }
