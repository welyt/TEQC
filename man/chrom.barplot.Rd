\name{chrom.barplot}
\alias{chrom.barplot}
\title{Reads per chromosome barplot}
\description{Barplot of the absolute numbers of reads aligning to each chromosome}

\usage{chrom.barplot(reads, ylab, ...)}

\arguments{
  \item{reads}{\code{\link[IRanges:RangedData-class]{RangedData}} table containing read positions,
      i.e. output from \code{\link{get.reads}}. To ensure a useful ordering of the bars, the
      chromosome information ('spaces' of \code{reads}) should be given as "chr" plus a
      number/letter [plus further specification], e.g. "chr1", "chrX", "chr17_ctg5_hap1", "chrUn_gl000211".}
  \item{ylab}{y-axis label}
  \item{\dots}{graphical parameters passed to \code{barplot}}
}
%\details{}

\value{Barplot of reads per chromosome.}

%\references{}
\author{Manuela Hummel \email{manuela.hummel@crg.es}}
%\note{}

\seealso{\code{\link{get.reads}}}

\examples{
exptPath <- system.file("extdata", package="TEQC")
readsfile <- paste(exptPath, "ExampleSet_Reads.bed", sep="/")
reads <- get.reads(readsfile, idcol=4, skip=0)
chrom.barplot(reads)
}

\keyword{ hplot }
