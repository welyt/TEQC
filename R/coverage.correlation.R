
# 'pairs' plot function adapted to Rle input; and more specific

Rlepairs <-
function (x, plotfrac=0.001, seed=123, labels, main, cex.labels, cex.pch=1,
                      cex.main=1.2, cex.corr, font.labels=1, font.main=2, ...){

    # panel functions
    localAxis <- function(side, x, y, xpd, bg, col=NULL, main, oma, ...) {
        if(side %%2 == 1) Axis(x, side=side, xpd=NA, ...)
        else Axis(y, side=side, xpd=NA, ...)
    }

    localUpperPanel <- function(x, y, cex.corr){
      Corr <- cor(x, y, use="complete.obs")
      if(missing(cex.corr))
        cex.corr <- 3 * Corr
      text(x=0.5, y=0.5, round(Corr, digits=2), cex=cex.corr)
    }

    localLowerPanel <- function(x, y, ...){
      points(x, y, ...)
      abline(a=0, b=1, col=2)
    }

    # random selection of 'plotfrac' values for scatter plots
    n <- length(x[[1]])
    set.seed(seed)
    r <- sample(n, plotfrac*n)

    # sample labels
    nc <- length(x)
    if(missing(labels)){
      labels <- names(x)
      if (is.null(labels))
        labels <- paste("sample", 1L:nc)
    }

    # outer margins
    dots <- list(...); nmdots <- names(dots)
    oma <- if("oma" %in% nmdots) dots$oma else NULL
    if (is.null(oma)) {
        oma <- c(4, 4, 4, 4)
        if (!missing(main)) oma[3L] <- 6
    }
    opar <- par(mfrow = c(nc, nc), mar = rep.int(.5, 4), oma = oma)
    on.exit(par(opar))

    # pairwise plots
    for (i in 1L:nc)
        for (j in 1L:nc) {
            # initialize plot
            plot(x[[j]], x[[i]], xlab = "", ylab = "", axes = FALSE, type = "n", main="", ...)
                box()

                # axis ticks
                if(i == 1  && j > i)
                    localAxis(3, x[[j]], x[[i]], ...)
                if(j == nc && j > i)
                    localAxis(4, x[[j]], x[[i]], ...)

                # sample labels (diagonal panels)
                if(i == j) {
                  par(usr = c(0, 1, 0, 1))
                  if(missing(cex.labels)) {
                    l.wid <- strwidth(labels, "user")
                    cex.labels <- max(0.8, min(2, .9 / max(l.wid)))
                  }
                  text(x=0.5, y=0.5, labels[i], cex=cex.labels, font=font.labels)
                }

                # scatter plot (upper panels)
                else if(i < j)
                  localLowerPanel(x[[j]][r], x[[i]][r], cex=cex.pch, ...)

                # correlations (lower panels)
                else {
                    par(usr = c(0, 1, 0, 1))
                    localUpperPanel(x[[j]], x[[i]], cex.corr=cex.corr)
                }
       }

    # main title
    if (!missing(main))
        mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
}



# extract (normalized) coverages, scatterplot of randomly selected fraction of coverage values,
# coverage correlation (using all values)

coverage.correlation <-
function(coveragelist, normalized=TRUE, plotfrac=0.001, seed=123, labels, main, pch=".", cex.labels,
                  cex.pch=2, cex.main=1.2, cex.corr, font.labels=1, font.main=2, ...){

  # check input
  n.samples <- length(coveragelist)
  listnames <- lapply(coveragelist, names)
  if(!all(sapply(listnames, function(x) "coverageTarget" %in% x)))
    stop("elements of input list need element 'coverageTarget' (run 'coverage.target()' with option 'perBase=TRUE')")

  coverages <- lapply(coveragelist, function(x) x$coverageTarget)

  chroms <- names(coverages[[1]])
  L <- sapply(coverages[[1]], length)
  if(!all(sapply(coverages, function(x) identical(names(x), chroms))))
    stop("elements of 'covlist' do not seem to be coverages for the same targets (there are positions on different chromosomes for different samples)")
  if(!all(sapply(coverages, function(x) identical(sapply(x, length), L))))
    stop("elements of 'covlist' do not seem to be coverages for the same targets (there are different numbers of target positions for different samples)")

  coverages <- lapply(coverages, unlist, use.names=FALSE)

  # normalized coverages
  if(normalized){
    for(i in 1:n.samples)
      coverages[[i]] <- coverages[[i]] / coveragelist[[i]]$avgTargetCoverage
  }

  Rlepairs(coverages, plotfrac=plotfrac, seed=seed, labels=labels, main=main, cex.labels=cex.labels, cex.pch=cex.pch,
            cex.main=cex.main, cex.corr=cex.corr, font.labels=font.labels, font.main=font.main, pch=pch, ...)
}
