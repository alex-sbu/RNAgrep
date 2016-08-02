drawtbone <- function(marks){
    marks$x1 <- as.numeric(marks$x1)
    marks$x2 <- as.numeric(marks$x2)
    marks$y1 <- as.numeric(marks$y1)
    marks$y2 <- as.numeric(marks$y2)
    marks$size <- as.numeric(marks$size)
    marks$srt <- as.numeric(marks$srt)
    marks$name <- as.character(marks$name)
    marks$type <- as.character(marks$type)
    if(!is.null(marks$segcol)){
        segcol <- marks$segcol
    }else{
        segcol <- "black"
    }
    segments(marks$x1, marks$y1, marks$x1, marks$y2, lwd = 0.5, col = segcol)
    segments(marks$x2, marks$y1, marks$x2, marks$y2, lwd = 0.5, col = segcol)
    ymidline <- (marks$y1+marks$y2)/2
    xmidline <- (marks$x1+marks$x2)/2
    segments(marks$x1, ymidline, marks$x2, ymidline, lwd = 0.5)
}

drawrect <- function(marks){
    marks$x1 <- as.numeric(marks$x1)
    marks$x2 <- as.numeric(marks$x2)
    marks$y1 <- as.numeric(marks$y1)
    marks$y2 <- as.numeric(marks$y2)
    marks$size <- as.numeric(marks$size)
    marks$srt <- as.numeric(marks$srt)
    marks$name <- as.character(marks$name)
    marks$type <- as.character(marks$type)
    if(!is.null(marks$segcol)){
        segcol <- as.numeric(marks$segcol)
    }else{
        segcol <- "black"
    }
    segments(marks$x1, marks$y1, marks$x1, marks$y2, lwd = 0.5, col = segcol)
    segments(marks$x2, marks$y1, marks$x2, marks$y2, lwd = 0.5, col = segcol)
    segments(marks$x1, marks$y1, marks$x2, marks$y1, lwd = 0.5, col = segcol)
    segments(marks$x1, marks$y2, marks$x2, marks$y2, lwd = 0.5, col = segcol)
    ymidline <- (marks$y1+marks$y2)/2
    xmidline <- (marks$x1+marks$x2)/2
}

drawlab <- function(marks){
    marks$x1 <- as.numeric(marks$x1)
    marks$x2 <- as.numeric(marks$x2)
    marks$y1 <- as.numeric(marks$y1)
    marks$y2 <- as.numeric(marks$y2)
    marks$size <- as.numeric(marks$size)
    marks$srt <- as.numeric(marks$srt)
    marks$name <- as.character(marks$name)
    marks$type <- as.character(marks$type)
    marks$xlab <- as.numeric(marks$xlab)
    marks$ylab <- as.numeric(marks$ylab)
    marks$name <- as.character(marks$name)
    if(!is.null(marks$xlab)){
        xmidline <- as.numeric(marks$xlab)
    }else{
        xmidline <- (marks$x1+marks$x2)/2
    }
    if(!is.null(marks$xlab)){
        ymidline <- as.numeric(marks$ylab)
    }else{
        ymidline <- (marks$y1+marks$y2)/2
    }
    greekchar <- c("alpha", "beta")
    greek <- which(marks$name %in% greekchar)

    if(length(greek) > 0){
        extractgreek <- marks$name[greek]
        indexgreek <- which(greekchar %in% extractgreek)
        for(m in 1:length(indexgreek)){
            if(indexgreek[m] == 1){text(xmidline, ymidline, expression(alpha), srt = marks$srt, cex = marks$size)}
            if(indexgreek[m] == 2){text(xmidline, ymidline, expression(beta), srt = marks$srt, cex = marks$size)}
        }
    }else{
        text(xmidline, ymidline, marks$name,
             srt = marks$srt, cex = marks$size)        
    }
}

drawpixels <- function(marks){
    marks$x1 <- as.numeric(marks$x1)
    marks$x2 <- as.numeric(marks$x2)
    marks$y1 <- as.numeric(marks$y1)
    marks$y2 <- as.numeric(marks$y2)
    rect((marks$x1), marks$y1, (marks$x2), marks$y2,
         col = "black",
         lend = 1, border = NA)
}

drawpoints <- function(marks){
    marks$x1 <- as.numeric(marks$x1)
    marks$y1 <- as.numeric(marks$y1)
    points(marks$x1, marks$y1, pch = 15)
}
