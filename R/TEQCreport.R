
# html report generation

# get results values and plots
TEQCreport <- function(sampleName="", targetsName="", referenceName="", destDir="TEQCreport",
            reads=get.reads(), targets=get.targets(), Offset=0, pairedend=FALSE, genome=c(NA, "hg19", "hg18"),
            genomesize, k=c(1, 2, 3, 5, 10, 20), covthreshold=8, CovUniformityPlot=FALSE, CovTargetLengthPlot=FALSE, CovGCPlot=FALSE,
            duplicatesPlot=FALSE, baits=get.baits(), WigFiles=FALSE, saveWorkspace=FALSE, figureFormat=c("jpeg","png","tiff")){
# sampleName, targetsName, referenceName: names that can be chosen by user and will be placed on top of html report
# destDir: output directory
# reads, targets: reads/targets RangedData tables or commands how to read them
# Offset: option used in various functions to add ?Offset? bases on both sides of targeted regions
# pairedend: are the data paired-end reads?
# genome, genomesize: options as needed for 'fraction.target()'
# k: parameter for 'covered.k()'
# !! covthreshold: coverage threshold for 'coverage.hist()'
# CovUniformityPlot, CovTargetLengthPlot, CovGCPlot, duplicatesPlot: shall corresponding plots be created?
# baits: baits table or command how to read it
# WigFiles: shall wiggle files be created
# saveWorkspace: should R workspace with 'reads', 'targets' and output of 'coverage.target()' and 'reads2pairs()'
#   (if pairedend=T) bet saved for further analyses?
# figureFormat: format of the figures for the html report (besides pdf graphs)
  
  figureFormat <- match.arg(figureFormat)

  # output directory
  if (!file.exists(destDir))
    dir.create(destDir, recursive=TRUE)
  wd <- getwd()
  setwd(destDir)
  Path <- getwd()
  setwd(wd)
  print(paste("results and report are saved in folder", Path))

  # image directory
  imgDir <- file.path(destDir, "image")
  if (!file.exists(imgDir))
    dir.create(imgDir)

  # wiggle file directory
  if(WigFiles){
    wigDir <- file.path(destDir, "wiggle")
    if (!file.exists(wigDir))
      dir.create(wigDir)
  }

  # some checks
  genome <- match.arg(genome)
  if (missing(genomesize) & is.na(genome))
    stop("either 'genome' or 'genomesize' has to be specified")

  if(CovGCPlot){
    if(missing(baits))
      stop("if 'CovGCPlot = TRUE', a 'baits' table has to be specified")
    else if(data.class(baits) != "RangedData")
      stop("the 'baits' table has to be of class 'RangedData'")
  }

  print("reading data...")
  if(missing(reads))
    stop("'reads' have to be specified")
    
  if(data.class(reads) != "RangedData")
    stop("the 'reads' table has to be of class 'RangedData'")

  if(missing(targets))
    stop("'targets' have to be specified")

  if(data.class(targets) != "RangedData")
    stop("the 'targets' table has to be of class 'RangedData'")

  n.reads <- nrow(reads)
  n.targets <- nrow(targets)

  if(pairedend){
    print("collapsing reads to pairs...")
    readpairs <- reads2pairs(reads)
    if(is.list(readpairs)){
      n.pairs <- nrow(readpairs$readpairs)
      n.singles <- nrow(readpairs$singleReads)
    }
    else{
      n.pairs <- nrow(readpairs)
      n.singles <- 0
    }
  }

  # specificity and enrichment
  ft <- fraction.target(targets, Offset=Offset, genome=genome, genomesize=genomesize)
  if(pairedend){
    print("calculating fraction of on-target read pairs")
    fr <- fraction.reads.target(readpairs, targets, Offset=Offset)
  }
  else {
    print("calculating fraction of on-target reads")
    fr <- fraction.reads.target(reads, targets, Offset=Offset)
  }
  enr <- as.character(round(fr / ft))
  
  # coverage
  print("calculating coverage...")
  Coverage <- coverage.target(reads, targets, Offset=Offset)
  stats <- c(fr, Coverage$avgTargetCoverage, Coverage$targetCoverageSD, Coverage$targetCoverageQuantiles)
  names(stats) <- c("fractionReadsOnTarget", "avgCoverage", "coverageSD", "coverageMin", "coverageQuartile1", "medianCoverage", "coverageQuartile3", "coverageMax")
  write.table(data.frame(stats), file=file.path(destDir, "onTargetStatistics.txt"), sep="\t", col.names=FALSE, quote=FALSE)
  avgcov <- data.frame(round(Coverage$avgTargetCoverage, 2), round(Coverage$targetCoverageSD, 2), matrix(Coverage$targetCoverageQuantiles, ncol=5))
  names(avgcov) <- c("avgTargetCoverage", "targetCoverageSD", paste(names(Coverage$targetCoverageQuantiles), "quantile"))

  # coverage per target
  print("counting reads per target...")
  targetcov0 <- Coverage$targetCoverages
  targetcov0 <- readsPerTarget(reads, targetcov0, Offset=Offset)
  targetcov <- as.data.frame(targetcov0)
  write.table(targetcov, file=file.path(destDir, "target_coverage.txt"),
              sep="\t", row.names=F, quote=F)
  if(nrow(targetcov) > 20)
    targetcov <- rbind(apply(targetcov[1:20,], 2, as.character), "...")
    
  # sensitivity
  sensi0 <- covered.k(Coverage$coverageTarget, k=k)
  sensi <- round(sensi0 * 100, 2)
  N <- paste(">=", names(sensi), "X", sep="")
  sensi <- paste(sensi, "%", sep="")
  names(sensi) <- N

  # values for make.report
  print("generating figures...")
  values <-
        list(SAMPLE=sampleName,
             NREADS=as.character(nrow(reads)),
             TARGETS=targetsName,
             NTARGETS=as.character(nrow(targets)),
             REFERENCE=referenceName,
             OFFSET=as.character(Offset),
             SPECIFICITY=hwrite(paste(round(fr*100, 2), "%", sep="")),
             ENRICHMENT=hwrite(enr),
             CHROM_BARPLOT=htmlChromBarplot(destDir, reads, targets, figureFormat),
             AVGCOV=hwrite(avgcov),
             COVTARG=hwrite(targetcov),
             SENSITIVITY=hwrite(sensi),
             COV_HIST=htmlCoverageHist(destDir, Coverage$coverageTarget, figureFormat, covthreshold=covthreshold))

  # special output for paired-end data
  if(pairedend)
    values <- c(values, list(
             NPAIRS=as.character(n.pairs),
             NSINGLES=as.character(n.singles),
             ISIZEHIST=htmlInsertSizeHist(destDir, readpairs, figureFormat)))

  # optional additional plots
  if(CovUniformityPlot)
    values <- c(values, list(COV_UNIFORM=htmlCovUniformity(destDir, Coverage, figureFormat)))
  if(CovTargetLengthPlot)
    values <- c(values, list(COV_TARGLEN=htmlCovTargetLength(destDir, targetcov0, figureFormat)))
  if(CovGCPlot)
    values <- c(values, list(COV_GC=htmlCovGC(destDir, Coverage$coverageAll, baits, figureFormat)))
  if(duplicatesPlot){
    print("duplicates analysis...")
    if(pairedend)
      values <- c(values, list(DUPLICATES=htmlDuplicatesBarplot(destDir, readpairs, targets, ylab="Fraction of read pairs", figureFormat)))
    else
      values <- c(values, list(DUPLICATES=htmlDuplicatesBarplot(destDir, reads, targets, figureFormat)))
  }
  
  # create wiggle files
  if(WigFiles){
    make.wigfiles(Coverage$coverageAll, filename=file.path(destDir, "wiggle", "Coverage"))
    chroms <- names(Coverage$coverageAll)
    wignames <- paste("Coverage_", chroms, ".wig", sep="")
    values <- c(values, list(WIG=hwrite(matrix(wignames), link=file.path(".", "wiggle", wignames))))
  }

  # make html report
  print("generating html report...")
  make.report(destDir=destDir, values=values, pairedend=pairedend, CovUniformityPlot=CovUniformityPlot,
              CovTargetLengthPlot=CovTargetLengthPlot, CovGCPlot=CovGCPlot, duplicatesPlot=duplicatesPlot,
              WigFiles=WigFiles)
  
  # save R objects for further usage
  if(saveWorkspace){
    print("saving workspace...")
    if(pairedend)
      save(reads, targets, Coverage, readpairs, file=file.path(destDir, "results.RData"))
    else
      save(reads, targets, Coverage, file=file.path(destDir, "results.RData"))
  }
}



# write results to html report
make.report <- function(destDir, values, pairedend, CovUniformityPlot, CovTargetLengthPlot,
                        CovGCPlot, duplicatesPlot, WigFiles, ...){

  # different sections of html page, templates are in inst/template
  if(pairedend)
    fls <- c("0000-Header.html", "1000-Overview.html", "1001-Pairedend.html", "2001-SpeciEnrichment.html",
           "3000-Coverage.html")
  else
    fls <- c("0000-Header.html", "1000-Overview.html", "2000-SpeciEnrichment.html",
           "3000-Coverage.html")
  if(CovUniformityPlot)
    fls <- c(fls, "4000-CoverageUniformity.html")
  if(CovTargetLengthPlot)
    fls <- c(fls, "5000-CoverageLength.html")
  if(CovGCPlot)
    fls <- c(fls, "6000-CoverageGC.html")
  if(WigFiles)
    fls <- c(fls, "7000-Wiggle.html")
  if(duplicatesPlot & pairedend)
    fls <- c(fls, "8001-Duplicates.html")
  else if(duplicatesPlot)
    fls <- c(fls, "8000-Duplicates.html")
  fls <- c(fls, "9999-Footer.html")

  sections <- system.file("template", fls, package="TEQC")
   
  # cssFile: html settings template, QA.css template file taken from ShortRead package in inst/template
  cssFile <- c(QA.css=system.file("template", "QA.css", package="TEQC"))
  htmlFile <- file.path(destDir, "index.html")
  biocFile <- "bioclogo-small.jpg"
  values <- c(list(CSS=names(cssFile), DATE=date(), VERSION=packageDescription("TEQC")$Version), values)

  # open index.html, copy sections from templates and fill in values
  toConn <- file(htmlFile, "w")
  for (sec in sections) {
    fromConn <- file(sec, open="r")
    copySubstitute(sec, toConn, values)       # function from Biobase
    close(fromConn)
  }
  close(toConn)
  
  # copy QA.css and BioC image to destDir
  file.copy(cssFile, file.path(destDir, names(cssFile)))
  file.copy(system.file("template", "image", biocFile, package="TEQC"), file.path(destDir, "image", biocFile))
  htmlFile
}


# create jpeg and pdf figures and put them into report
html_img <- function(dir, file, fig, figureFormat, ...){
    figFile <- paste(file, figureFormat, sep=".")
    pdfFile <- paste(file, "pdf", sep=".")
    imgDir <- file.path(dir, "image")

    if(figureFormat == "jpeg")
      jpeg(file.path(imgDir, figFile), ...)
    else if(figureFormat == "png")
      png(file.path(imgDir, figFile), ...)
    if(figureFormat == "tiff")
      tiff(file.path(imgDir, figFile), ...)
    fig
    dev.off()

    pdf(file.path(imgDir, pdfFile), ...)    # pdf file is always empty!?!
    fig                                     # seems that only one device can be printed usefully like that...
    dev.off()

    # show jpeg in the report but link to pdf
    hwriteImage(file.path(".", "image", figFile), link=file.path(".", "image", pdfFile))  # function from hwriter
}


htmlChromBarplot <- function(dir, reads, targets, figureFormat, ...){
    figFile <- paste("chrom_barplot", figureFormat, sep=".")
    pdfFile <- "chrom_barplot.pdf"
    imgDir <- file.path(dir, "image")

    if(figureFormat == "jpeg")
      jpeg(file.path(imgDir, figFile), width=800, ...)
    else if(figureFormat == "png")
      png(file.path(imgDir, figFile), width=800, ...)
    if(figureFormat == "tiff")
      tiff(file.path(imgDir, figFile), width=800, ...)
    chrom.barplot(reads, targets)
    dev.off()

    pdf(file.path(imgDir, pdfFile), width=12, ...)
    chrom.barplot(reads, targets)
    dev.off()

    # show jpeg in the report but link to pdf
    hwriteImage(file.path(".", "image", figFile), link=file.path(".", "image", pdfFile))  # function from hwriter
}


htmlCoverageHist <- function(dir, coverageTarget, figureFormat, ...){
  figFile <- paste("coverage_histogram", figureFormat, sep=".")
    pdfFile <- "coverage_histogram.pdf"
    imgDir <- file.path(dir, "image")

  if(figureFormat == "jpeg")
    jpeg(file.path(imgDir, figFile))
  else if(figureFormat == "png")
    png(file.path(imgDir, figFile))
  if(figureFormat == "tiff")
    tiff(file.path(imgDir, figFile))
  sensi <- .coverage.hist(coverageTarget, ...)
    dev.off()
    write.table(sensi, file=file.path(dir, "sensitivity.txt"), sep="\t", row.names=FALSE, quote=FALSE)

    pdf(file.path(imgDir, pdfFile))
    coverage.hist(coverageTarget, ...)
    dev.off()

    hwriteImage(file.path(".", "image", figFile), link=file.path(".", "image", pdfFile))
}


htmlInsertSizeHist <- function(dir, readpairs, figureFormat, ...){
  figFile <- paste("insert_size_histogram", figureFormat, sep=".")
    pdfFile <- "insert_size_histogram.pdf"
    imgDir <- file.path(dir, "image")

  if(figureFormat == "jpeg")
    jpeg(file.path(imgDir, figFile), ...)
  else if(figureFormat == "png")
    png(file.path(imgDir, figFile), ...)
  if(figureFormat == "tiff")
    tiff(file.path(imgDir, figFile), ...)
  insert.size.hist(readpairs)
    dev.off()

    pdf(file.path(imgDir, pdfFile), ...)
    insert.size.hist(readpairs)
    dev.off()

    hwriteImage(file.path(".", "image", figFile), link=file.path(".", "image", pdfFile))
}

htmlCovUniformity <- function(dir, Coverage, figureFormat, ...){
  figFile <- paste("coverage_uniformity", figureFormat, sep=".")
    pdfFile <- "coverage_uniformity.pdf"
    imgDir <- file.path(dir, "image")

  if(figureFormat == "jpeg")
    jpeg(file.path(imgDir, figFile), ...)
  else if(figureFormat == "png")
    png(file.path(imgDir, figFile), ...)
  if(figureFormat == "tiff")
    tiff(file.path(imgDir, figFile), ...)
  coverage.uniformity(Coverage)
    dev.off()

    pdf(file.path(imgDir, pdfFile), ...)
    coverage.uniformity(Coverage)
    dev.off()

    hwriteImage(file.path(".", "image", figFile), link=file.path(".", "image", pdfFile))
}

htmlCovTargetLength <- function(dir, targetcov, figureFormat, ...){
  figFile <- paste("coverage_targetlength", figureFormat, sep=".")
    pdfFile <- "coverage_targetlength.pdf"
    imgDir <- file.path(dir, "image")

  if(figureFormat == "jpeg")
    jpeg(file.path(imgDir, figFile), width=1000, ...)
  else if(figureFormat == "png")
    png(file.path(imgDir, figFile), width=1000, ...)
  if(figureFormat == "tiff")
    tiff(file.path(imgDir, figFile), width=1000, ...)
    par(mfrow=c(1,2))
    coverage.targetlength.plot(targetcov, plotcolumn="nReads")
    coverage.targetlength.plot(targetcov, plotcolumn="avgCoverage")
    dev.off()

    pdf(file.path(imgDir, pdfFile), width=14, ...)
    par(mfrow=c(1,2))
    coverage.targetlength.plot(targetcov, plotcolumn="nReads")
    coverage.targetlength.plot(targetcov, plotcolumn="avgCoverage")
    dev.off()

    hwriteImage(file.path(".", "image", figFile), link=file.path(".", "image", pdfFile))
}

htmlCovGC <- function(dir, coverageAll, baits, figureFormat, ...){
  figFile <- paste("coverage_GC", figureFormat, sep=".")
    pdfFile <- "coverage_GC.pdf"
    imgDir <- file.path(dir, "image")

  if(figureFormat == "jpeg")
    jpeg(file.path(imgDir, figFile), ...)
  else if(figureFormat == "png")
    png(file.path(imgDir, figFile), ...)
  if(figureFormat == "tiff")
    tiff(file.path(imgDir, figFile), ...)
  coverage.GC(coverageAll, baits)
    dev.off()

    pdf(file.path(imgDir, pdfFile), ...)
    coverage.GC(coverageAll, baits)
    dev.off()

    hwriteImage(file.path(".", "image", figFile), link=file.path(".", "image", pdfFile))
}

htmlDuplicatesBarplot <- function(dir, reads, targets, figureFormat, ...){
  figFile <- paste("duplicates_barplot", figureFormat, sep=".")
  pdfFile <- "duplicates_barplot.pdf"
    imgDir <- file.path(dir, "image")

  if(figureFormat == "jpeg")
    jpeg(file.path(imgDir, figFile), ...)
  else if(figureFormat == "png")
    png(file.path(imgDir, figFile), ...)
  if(figureFormat == "tiff")
    tiff(file.path(imgDir, figFile), ...)
  duplicates.barplot(reads, targets, ...)
    dev.off()

    pdf(file.path(imgDir, pdfFile))
    duplicates.barplot(reads, targets, ...)
    dev.off()

    hwriteImage(file.path(".", "image", figFile), link=file.path(".", "image", pdfFile))
}




