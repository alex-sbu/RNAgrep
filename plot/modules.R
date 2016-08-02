ylabcex <- 1
library(scales)
library(IRanges)


diagram <- function(input, mars, xlimit = NULL){
    data <- read.csv(input)
    data <- na.omit(data)
    rec <- which(data$type == "rectangle")
    tbone <- which(data$type == "tbone")
    lab <- which(data$name != "")
    marks <- lapply(1:nrow(data), function(x) data[x, c(1:ncol(data))])
    ylimit <- c(min(c(data$y1, data$y2)),
    (max(c(data$y1, data$y2))+(max(c(data$y1, data$y2))*0.35)))
    if (is.null(xlimit)) {
    xlimit <- c(min(c(data$x1, data$x2)),
                max(c(data$x1, data$x2)))
    }
    ## mars[3] <- 0.5
    par(mar = mars)
    plot(1000, 1000, xlim = xlimit, ylim = ylimit,
         bty='n', col="white", yaxt='n', xaxt='n', ann=FALSE)
    for(i in 1:length(rec)){
        drawrect(marks[[rec[i]]])
    }
    if(length(tbone)>0) {
    for(i in 1:length(tbone)){
        drawtbone(marks[[tbone[i]]])
    }
    }
    for(i in 1:length(lab)){
        drawlab(marks[[lab[i]]])
    }
}



synP <- function(input, mars, xlimit = NULL, modx = NULL, ...){
    load(input)
    window <- data$window
    nonzero <- data$nonzero
    threshold <- data$threshold
    Ay <- data$Ay
    By <- data$By
    Cy <- data$Cy
    x <- data$x
    y <- data$y    
    if(is.null(xlimit)) xlimit <- c(0, length(x))
    ylimit <- c(1, min((Cy), (threshold^2)))
    ylimit[2] <- 1e-20
    par(lend = 1, ljoin = 2)
    par(mar = mars)
    plot(x, y, col = "white",
         ylab = "n", xlab = "n",
         xaxt="n", yaxt="n",
         ann = FALSE, xlim = xlimit,
         ylim = ylimit, bty="n", log="y")
    lines(xlimit, rep(threshold, length(xlimit)), lty="dashed", col="#616161", lwd = 0.5)
    lines(x[1:(length(x)-1)], Cy, col="#000000")
    cc0 = 1.2
    cc1 = 1.0
    cc2 = 0.8
    ## Plot p-values
    ## mtext (paste0("p-value"), cex=ylabcex,
    ##        at=10^mean(log10(ylimit)), side=2, xpd=TRUE,
    ##        las=0, line = 1)
    ## text(0.9,0.5,labels="p-value",cex=cc1,adj=0.5,srt=90)
    ylimit[2] <- 10^(ceiling(log10(ylimit[2])))
    axisat <- axTicks(side = 4, log = TRUE)
    axislabels <- pretty10exp(axisat, drop.1 = TRUE)
    ## drawsideaxis(axisat = axisat, axislabels = axislabels, axislas = 2)
    ## drawaxistitle("p-value", 10^mean(log10(ylimit)))
    box()
}



synOE <- function(input, mars, xlimit = NULL, ...){
    load(input)
    window <- data$window
    nonzero <- data$nonzero
    threshold <- data$threshold
    Ay <- data$Ay
    By <- data$By
    Cy <- data$Cy
    x <- data$x
    y <- data$y
    if(is.null(xlimit)) xlimit <- c(0, length(x))
    ylimit <- c(0, 2)
    par(lend = 1, ljoin = 2)
    par(mar = mars)
    plot(xlimit, ylimit, col="white", ylab="n", xlab="n", xaxt="n",
         yaxt="n", ann=FALSE, xlim=xlimit, ylim=ylimit, bty="n")
    lines(xlimit, rep(1, length(xlimit)), lty="dashed", col="#616161", lwd = 0.5)
    lines(x, y, col="#555753")
    cc0 = 1.2
    cc1 = 1.0
    cc2 = 0.8
    ## Plot p-values
    ## text(0.9,0.5,labels="obs/exp",cex=cc1,adj=0.5,srt=90)
    ## mtext (paste0("obs/exp"), cex=ylabcex,
    ##        at=(mean(ylimit)), side=2, xpd=TRUE,
    ##        las = 0,  line = 1)
    ylimit <- round(ylimit, 2)
    axispos <- 0
    axiswidth <- 1
    labelcex <- 0.8
    labelpos <- 0
    ## drawsideaxis()
    ## drawaxistitle("obs/exp", mean(ylimit))
    box()
}




barcode <- function(input, mars, xlimit = NULL, labs = NULL, ...) {
    name <- gsub("^.+/(.+).csv$", "\\1", input, perl = TRUE)
    data <- read.csv(input)
    data <- na.omit(data)
    nlanes <- ncol(data) - 1
    nbars <- length(data$pos)
    labs <- names(data)
    labs <- labs[2:length(labs)]    
    gety <- function(height = 10, gap = 1, lanes = nlanes){
        fheight <- seq(0, (height * nlanes), height)
        sheight <- sapply(fheight[2:length(fheight)], function(x) x - 1)
        sheight <- c(0, sheight)
        yout <- list()
        for(i in 1:nlanes){
            k <- i + 1
            ylane <- numeric(length = 2)
            ylane <- c(fheight[i], sheight[k])
            yout[[i]] <- ylane
        }
        return(yout)
    }    
    getx <- function(width = 1, bars = nbars, pos = NULL, lanes = nlanes){
        ## print(paste(width, bars, lanes, sep=" "))
        if (is.null(pos)) {
            fwidth <- seq(1, (width * nbars), width)
        }else{
            fwidth <- seq((pos[1]*width), (width * pos[nbars]), width)            
        }
        swidth <- c(fwidth[2:length(fwidth)], (fwidth[length(fwidth)]+1))
        xout <- list()
        for(i in 1:nlanes){
            k <- i + 1
            xlist <- list()
            for(j in 1:nbars){
                l <- j - 1
                xbar <- numeric(length = 2)
                xbar <- c(fwidth[j], swidth[j])
                xlist[[j]] <- xbar
            }
            xout[[i]] <- xlist
        }
        return(xout)
    }    
    getmarkx <- function(x, y){
        counter <- 0
        markx <- list()
        marky <- list()
        for(i in 1:length(x)){
            ni <- length(x[[i]])
            xlist <- list()
            for(j in 1:ni){
                counter <- counter + 1
                markx[[counter]] <- x[[i]][[j]]
            }
        }
        return(markx)
    }    
    getmarky <- function(x, y){
        counter <- 0
        markx <- list()
        marky <- list()
        for(i in 1:length(x)){
            ni <- length(x[[i]])
            xlist <- list()
            for(j in 1:ni){
                counter <- counter + 1
                marky[[counter]] <- y[[i]]
            }
        }
        return(marky)
    }    
    getmarks <- function(markx, marky){
        marks <- list()
        marks[[1]] <- numeric(length = length(markx))
        marks[[2]] <- numeric(length = length(markx))
        marks[[3]] <- numeric(length = length(markx))
        marks[[4]] <- numeric(length = length(markx))
        for(i in 1:length(markx)){
            marks[[1]][i] <- markx[[i]][1]
            marks[[2]][i] <- markx[[i]][2]
            marks[[3]][i] <- marky[[i]][1]
            marks[[4]][i] <- marky[[i]][2]
        }
        return(marks)
    }
    x <- getx(width = 1,
              bars = nbars,
              pos = as.numeric(as.character(data$pos)),
              lanes = nlanes)
    y <- gety(height = 10, gap = 1, lanes = nlanes)    
    if (is.null(xlimit)) {
        xlimit <- c(min(unlist(x)), max(unlist(x)))
    }
    ylimit <- c(min(unlist(y)), max(unlist(y)))
    markx <- getmarkx(x, y)
    marky <- getmarky(x, y)
    marks <- getmarks(markx, marky)
    names(marks) <- c("x1", "x2", "y1", "y2")
    ## marks <- lapply(1:nrow(data), function(x){
    ##     data[x, c(2:ncol(data))]})    
    ## ylimit <- c(min(c(data$y1, data$y2), na.rm = TRUE),
    ## (max(c(data$y1, data$y2), na.rm = TRUE)+(max(c(data$y1, data$y2), na.rm = TRUE)*0.35)))
    ## print(ylimit)
    ## xlimit <- c(min(c(data$x1, data$x2)),
    ##             max(c(data$x1, data$x2)))
    ## xlimit <- c(1, 3983)
    ## mars[1] <- 1
    par(mar = mars)
    plot(xlimit, ylimit,
         col = alpha("white", 0),
         ylab = "n", xlab = "n",
         xaxt = "n", yaxt = "n",
         xlim = xlimit, ylim = ylimit,
         bty = "n", ann = FALSE)
    counter <- 0
    for(i in 1:nlanes){
        k <- i + 1
        for(j in 1:nbars){
            counter <- counter + 1
            ## print(paste(k, j, data[k, j], sep = "--"))
            if(data[j, k] == 0){
                fixmark <- list((marks$x1[counter]),
                (marks$x2[counter]),
                marks$y1[counter],
                marks$y2[counter])
                names(fixmark) <- c("x1", "x2", "y1", "y2")
                drawpixels(fixmark)
            }
        }
    }
    ## print(marks)
    ## drawaxistitle("SNPs", mean(ylimit))
    if(is.null(labs)){
        labs <- names(data)
        labs <- labs[-1]
    }
    sideaxislabs <- numeric()
    sideaxisblanks <- character()
    for(i in 1:length(y)) {
        drawaxistitle(axtitle = labs[[i]],
                      axtitleat = mean(y[[i]]),
                      axtitleside = 4,
                      axtitlelas = 2,
                      axtitleline = 0.5,
                      axtitlecex = 0.8)
        sideaxislabs[i] <- mean(y[[i]])
        sideaxisblanks[i] <- " "
    }
    drawsideaxis(axisat = sideaxislabs, axislabels = sideaxisblanks)
    box()
    return(xlimit)
}



motifplot <- function(input, mars, xlimit = NULL, labs = NULL, col = "black",
                      sep = 0.5, colgradient = NULL,...) {
    load(input)
    x <- tout
    height <- 1
    if(is.null(colgradient)) colgradient <- rep("black", length(x))
    if(is.null(xlimit)){
        xlimit = x
        if (is(xlimit, "Ranges")){
            xlimit <- c(min(start(xlimit)), max(end(xlimit)))
        }
    }
    bins <- disjointBins(IRanges(start(x), end(x) + 1))
    ylimit <- c(0, max(bins)*(height + sep))
    ## plot.new()
    par(mar=mars)
    par(lend = 1, ljoin = 2, lwd = 0.75)
    par(oma = c(0,0,0,0))     
    plot(xlimit,
         ylimit,
         col = alpha("white", 0),
         ylab = "n", xlab = "n",
         xaxt = "n", yaxt = "n",
         bty = "n", ann = FALSE)
    ybottom <- bins * (sep + height) - height
    rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = colgradient, border = NA, ...)
    ylimit <- as.numeric(ylimit)
    ylimit <- round(ylimit, 2)
    axispos <- 0
    axiswidth <- 1
    labelcex <- 0.8
    labelpos <- 0
    ## drawsideaxis()
    ## drawaxistitle("count", mean(ylimit))
    if(!is.null(labs)) {
        drawaxistitle(axtitle = labs,
                      axtitleat = mean(ylimit),
                      axtitleside = 4,
                      axtitlelas = 2,
                      axtitleline = 0.5,
                      axtitlecex = 0.8)
}
    box()
}
