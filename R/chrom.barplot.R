chrom.barplot <-
function(reads, ylab, ...){

  tab0 <- table(space(reads))

  # order chromosomes
  chr <- substr(names(tab0), 4, 4)
  g <- grep("chr\\d{2}", names(tab0), perl=TRUE)
  chr[g] <- substr(names(tab0)[g], 4, 5)
  chrn <- as.numeric(chr[!(chr %in% c("M", "U", "X", "Y"))])
  tab <- tab0[!(chr %in% c("M", "U", "X", "Y"))]
  tab <- tab[order(chrn)]
  tab <- c(tab, tab0[chr %in% c("M", "U", "X", "Y")])

  if(missing(ylab)) ylab <- "Number of reads"
  options(scipen=10)
  B <- barplot(tab, names.arg=rep("", length(tab)), ylab=ylab, ...)
  text(B, par("usr")[3], labels=names(tab), srt=45, adj=c(1,0.8), xpd=T, cex=1, font=2)
}

