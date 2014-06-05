
# generation of a html report for several target capture samples
# based on output from TEQCreport() (produced with TEQC version >3.0.0)

multiTEQCreport <- function(singleReportDirs, samplenames, projectName="", targetsName="",
          referenceName="", destDir="multiTEQCreport", k=c(1, 2, 3, 5, 10, 20), figureFormat=c("jpeg","png","tiff")){
# singleReportDirs: string of directory names: output directories of function TEQCreport()
# samplenames: names of individual samples, used for plots
# projectName, targetsName, referenceName: names that can be chosen by user and will be placed on top of html report
# destDir: output directory
# k: parameter for 'covered.k()': ?k?-values for which to show fraction of target bases with coverage >= ?k?
# !!figureFormat: format of the figures for the html report (besides pdf graphs)
  
  # !!
  figureFormat <- match.arg(figureFormat)
  # !!

  # get/check sample directories and names
  n.samples <- length(singleReportDirs)

  if(n.samples < 2)
    stop("more than one directory name should be provided in 'singleReportDirs'")

  if(missing(samplenames))
    samplenames <- paste("sample", 1:n.samples)

  if(n.samples != length(samplenames))
    stop("'singleReportDirs' and 'samplenames' have to be of same length")

  if(!all(file.exists(singleReportDirs)))
    stop("invalid directory name(s) provided in 'singleReportDirs'")

  sampleInfo <- data.frame(Sample=samplenames, Report_Directory=singleReportDirs)
  
  # get target names
  tmp <- read.table(paste(singleReportDirs[1], "target_coverage.txt", sep="/"), header=TRUE, sep="\t", as.is=TRUE)
  targets <- paste(tmp$space, tmp$start, tmp$end, sep=":")

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


  ## go through single report directories and collect results
  targetstats <- NULL
  sensilist <- list()
#  sensi <- matrix(0, nrow=length(k), ncol=n.samples)
#  dimnames(sensi) <- list(k, samplenames)
  perTargCov <- NULL
  for(i in 1:n.samples){
    sr <- singleReportDirs[i]
    # on-target statistics
    tmp <- read.table(paste(sr, "onTargetStatistics.txt", sep="/"), row.names=1, header=FALSE)
    targetstats <- cbind(targetstats, as.matrix(tmp))

    # sensitivity
    tmp <- read.table(paste(sr, "sensitivity.txt", sep="/"), header=TRUE)
    sensilist <- c(sensilist, list(tmp))
#    k2 <- intersect(k, tmp$coverage)
#    sensi[as.character(k2), i] <- tmp$fractionTargetBases[tmp$coverage %in% k2]

    # per-target coverage
    tmp <- read.table(paste(sr, "target_coverage.txt", sep="/"), header=TRUE, sep="\t", as.is=TRUE)

    # check if targets are the same in all single reports
    targ <- paste(tmp$space, tmp$start, tmp$end, sep=":")
    if(!identical(targ, targets))
      stop("targets were not the same in all single reports")

    perTargCov <- cbind(perTargCov, tmp$avgCoverage)
  }
  colnames(targetstats) <- samplenames
  dimnames(perTargCov) <- list(targets, samplenames)

  ## save tables
  # specificity
  speci <- data.frame(fractionReadsOnTarget=targetstats["fractionReadsOnTarget",])
  write.table(speci, file=file.path(destDir, "fractionReadsOnTarget.txt"), sep="\t", quote=FALSE)

  # average and median target coverage
#  targcov <- data.frame(avgTargetCoverage=targetstats["avgCoverage",], medianTargetCoverage=targetstats["medianCoverage",])
  targcov <- targetstats[-1,]
  write.table(targcov, file=file.path(destDir, "targetCoverageStats.txt"), sep="\t", quote=FALSE)

  # sensitivity
  m <- max(sapply(sensilist, function(x) max(x$coverage)))
  sensi0 <- matrix(0, nrow=m+1, ncol=n.samples)  # all coverage values
  dimnames(sensi0) <- list(0:m, samplenames)
  for(i in 1:n.samples){
    tmp <- sensilist[[i]]
    k0 <- intersect(0:m, tmp$coverage)
    sensi0[as.character(k0), i] <- tmp$fractionTargetBases[tmp$coverage %in% k0]
  }
  k2 <- intersect(k, 0:m)        # just selected 'k' values
  sensi <- sensi0[as.character(k2),]
  sensi <- data.frame(coverage=rownames(sensi), sensi)
  write.table(sensi, file=file.path(destDir, "sensitivity.txt"), sep="\t", row.names=FALSE, quote=FALSE)

  # per-target coverage
  write.table(perTargCov, file=file.path(destDir, "targetCoverage.txt"), sep="\t", quote=FALSE)


  ## values for make.multi.report
  values <-
    list(PROJECT=projectName,
         NSAMPLES=as.character(n.samples),
         TARGETS=targetsName,
         REFERENCE=referenceName,
         SAMPLES=hwrite(sampleInfo),
         SPECI_BARPLOT=htmlSpeciBarplot(destDir, speci, figureFormat),
         #SPECI=hwrite(speci),
         COV_BOXPLOT=htmlCovBoxplot(destDir, targcov, figureFormat),
         SENSI_BARPLOT=htmlSensiBarplot(destDir, sensi, figureFormat),
         UNIF_PLOT=htmlUniformityPlot(destDir, sensi0, avgCov=targcov["avgCoverage",], figureFormat),
         COV_CORPLOT=htmlCovCorrelationPlot(destDir, perTargCov, figureFormat)
         )
         
  # create report
  make.multi.report(destDir=destDir, values=values)
}
                       

# write results to html report
make.multi.report <- function(destDir, values, ...){
  
  fls <- c("0000-Header.html", "1000-Overview.html", "2000-Specificity.html", "3000-TargetCoverage.html",
             "4000-Sensitivity.html","5000-Uniformity.html","6000-Correlation.html","9999-Footer.html")
  
  sections <- system.file("templateMulti", fls, package="TEQC")
  
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
  
  
  
## plot functions

speciBarplot <- function(speci){
  par(mar=c(8,4,4,3), las=2)
  barplot(t(100*speci), ylab="% on-target reads", ylim=c(0,100), col="cornflowerblue")
}

htmlSpeciBarplot <- function(dir, speci, figureFormat, ...){
  figFile <- paste("specificity_barplot", figureFormat, sep=".")
  pdfFile <- "specificity_barplot.pdf"
    imgDir <- file.path(dir, "image")

  if(figureFormat == "jpeg")
    jpeg(file.path(imgDir, figFile), ...)
  else if(figureFormat == "png")
    png(file.path(imgDir, figFile), ...)
  if(figureFormat == "tiff")
    tiff(file.path(imgDir, figFile), ...)
  speciBarplot(speci)
    dev.off()

    pdf(file.path(imgDir, pdfFile), ...)
    speciBarplot(speci)
    dev.off()

    # show jpeg in the report but link to pdf
    hwriteImage(file.path(".", "image", figFile), link=file.path(".", "image", pdfFile))
}

#[
covplot <- function(targcov){
  samplenames <- rownames(targcov)
  n.samples <- length(samplenames)
  par(xpd=T, mar=c(8,4,4,8), las=2)
  plot(targcov$avgTargetCoverage, ylim=c(0, max(targcov)), ylab="target coverage", xlab="",
        col="cornflowerblue", pch=16, cex=2, xaxt="n")
  points(targcov$medianTargetCoverage, col="darkorange", pch=17, cex=2)
  axis(side=1, at=1:n.samples, samplenames)
  legend(x=n.samples+n.samples*.05, y=mean(targcov$avgTargetCoverage), c("average", "median"), col=c("cornflowerblue","darkorange"), pch=16:17, pt.cex=2)
}
#]

covBoxplot <- function(targcov){
  par(mar=c(8,4,4,3), las=2)
  boxplot(targcov[-(1:2),], outline=FALSE, whiskcol="transparent", boxcol="cornflowerblue", medcol="cornflowerblue", staplecol="transparent", ylab="Coverage")
  points(1:ncol(targcov), targcov["avgCoverage",], pch="*", cex=2, col="darkred")
  legend("bottomleft", "average coverage", pch="*", pt.cex=2, col="darkred")
}

htmlCovBoxplot <- function(dir, targcov, figureFormat, ...){
  figFile <- paste("targetCoverage_boxplot", figureFormat, sep=".")
  pdfFile <- "targetCoverage_boxplot.pdf"
    imgDir <- file.path(dir, "image")

  if(figureFormat == "jpeg")
    jpeg(file.path(imgDir, figFile), ...)
  else if(figureFormat == "png")
    png(file.path(imgDir, figFile), ...)
  if(figureFormat == "tiff")
    tiff(file.path(imgDir, figFile), ...)
  covBoxplot(targcov)
    dev.off()

    pdf(file.path(imgDir, pdfFile), ...)
    covBoxplot(targcov)
    dev.off()

    # show jpeg in the report but link to pdf
    hwriteImage(file.path(".", "image", figFile), link=file.path(".", "image", pdfFile))
}


sensiBarplot <- function(sensi){
  k <- rownames(sensi)
  tmp <- as.matrix(sensi[,-1])
  n.samples <- ncol(tmp)

  # transform to "non-cumulative" fractions such that bars get correct heights (summing up to 1 and not more)
  tmp <- tmp - rbind(tmp[-1,], 0)
  tmp <- tmp[order(as.numeric(k), decreasing=TRUE),]

  col <- colorRampPalette(c("darkblue", "ghostwhite"))(nrow(sensi))
  par(mar=c(8,4,4,6), las=2, xpd=T)
  barplot(100*tmp, col=col, ylim=c(0,100), ylab="% target bases")
  p <- par("usr")[2]
  legend(x=p+p*0.025, y=80, legend=k, fill=rev(col), title="Coverage")
}

htmlSensiBarplot <- function(dir, sensi, figureFormat, ...){
  figFile <- paste("sensitivity_barplot", figureFormat, sep=".")
  pdfFile <- "sensitivity_barplot.pdf"
    imgDir <- file.path(dir, "image")

  if(figureFormat == "jpeg")
    jpeg(file.path(imgDir, figFile), ...)
  else if(figureFormat == "png")
    png(file.path(imgDir, figFile), ...)
  if(figureFormat == "tiff")
    tiff(file.path(imgDir, figFile), ...)
  sensiBarplot(sensi)
    dev.off()

    pdf(file.path(imgDir, pdfFile), ...)
    sensiBarplot(sensi)
    dev.off()

    # show jpeg in the report but link to pdf
    hwriteImage(file.path(".", "image", figFile), link=file.path(".", "image", pdfFile))
}


unifplot <- function(sensi, avgCov){
  n.samples <- ncol(sensi)
  par(mfrow=c(1,2))
  lty <- rep(1:3, length.out=n.samples)
  
  # raw coverage
  plot(NA, xlim=as.numeric(range(rownames(sensi))), ylim=c(0,1), main="Raw coverage", xlab="Coverage", ylab="Cumulative fraction of target bases")
  for(i in 1:n.samples)
    lines(x=rownames(sensi), y=sensi[,i], col=i, lty=lty[i], lwd=2)
  legend("topright", colnames(sensi), col=1:n.samples, lty=lty, lwd=2)
  
  # normalized coverage
  plot(NA, xlim=c(0,1), ylim=c(0,1), main="Normalized coverage", xlab="Normalized coverage", ylab="Cumulative fraction of target bases")
  for(i in 1:n.samples){
    normcov <- as.numeric(rownames(sensi)) / avgCov[i]
    lines(x=normcov, y=sensi[,i], col=i, lty=lty[i], lwd=2)
  }
  legend("topright", colnames(sensi), col=1:n.samples, lty=lty, lwd=2)
}

htmlUniformityPlot <- function(dir, sensi, avgCov, figureFormat, ...){
  figFile <- paste("coverageUniformity_plot", figureFormat, sep=".")
  pdfFile <- "coverageUniformity_plot.pdf"
    imgDir <- file.path(dir, "image")

  if(figureFormat == "jpeg")
    jpeg(file.path(imgDir, figFile), width=1000, ...)
  else if(figureFormat == "png")
    png(file.path(imgDir, figFile), width=1000, ...)
  if(figureFormat == "tiff")
    tiff(file.path(imgDir, figFile),width=1000, ...)
   unifplot(sensi, avgCov)
    dev.off()

    pdf(file.path(imgDir, pdfFile), width=14, ...)
    unifplot(sensi, avgCov)
    dev.off()

    # show jpeg in the report but link to pdf
    hwriteImage(file.path(".", "image", figFile), link=file.path(".", "image", pdfFile))
}


plot.i <- function(x, y, ...){
  pch <- "."
  cex <- 2
  if(length(x) < 1000){
    pch <- 20
    cex <- 1
  }
  points(x, y, cex=cex, pch=pch, ...)
  abline(0, 1, lty=2, col="grey")
}

cortext <- function(x, y, ...){
  np <- min(x, na.rm=T)
  xp <- max(x, na.rm=T)
  mp <- mean(c(np,xp))
  Co <- cor(x, y, use="complete.obs")
  text(x=mp, y=mp, round(Co, digits=2), cex=2*Co, ...)
}

covcorplot <- function(perTargCov){
  pairs(perTargCov, panel=plot.i, lower.panel=cortext)
}

htmlCovCorrelationPlot <- function(dir, perTargCov, figureFormat, ...){
  figFile <- paste("perTargetCoverageCorrelation_plot", figureFormat, sep=".")
  pdfFile <- "perTargetCoverageCorrelation_plot.pdf"
    imgDir <- file.path(dir, "image")

  if(figureFormat == "jpeg")
    jpeg(file.path(imgDir, figFile), width=800, height=800, ...)
  else if(figureFormat == "png")
    png(file.path(imgDir, figFile), width=800, height=800, ...)
  if(figureFormat == "tiff")
    tiff(file.path(imgDir, figFile), width=800, height=800, ...)
    covcorplot(perTargCov)
    dev.off()

    pdf(file.path(imgDir, pdfFile), width=12, height=12, ...)
    covcorplot(perTargCov)
    dev.off()

    # show jpeg in the report but link to pdf
    hwriteImage(file.path(".", "image", figFile), link=file.path(".", "image", pdfFile))
}

