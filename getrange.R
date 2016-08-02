suppressPackageStartupMessages(library(IRanges))

system(command = "./comp2lev.R -f compared-all.csv -h 743", intern = FALSE)
system(command = "./comp2lev.R -f compared-natural.csv -h 743", intern = FALSE)
system(command = "./comp2lev.R -f compared-synthetic.csv -h 743", intern = FALSE)
dirs <- c("compared-all", "compared-natural", "compared-synthetic")

plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, cn = NULL, ...){
    height <- 1
    if(is.null(cn)) cn <- rep("black", length(x))
    if (is(xlim, "Ranges"))
        xlim <- c(min(start(xlim)), max(end(xlim)))
    bins <- disjointBins(IRanges(start(x), end(x) + 1))
    plot.new()
    plot.window(xlim, c(0, max(bins)*(height + sep)))
    ybottom <- bins * (sep + height) - height
    print(bins)
    print(x)
    rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = cn, border = NA, ...)
    print(x)
    print(names(x))
    title(main)
    axis(1)
}


colsel <- c("#ff0000",
            "#ff0000",
            "#ff0000",
            "#000000",
            "#000000",
            "#859BC8",
            "#ABBAE3")
cgrad <- list()
for(x in 1:length(ir)){
    irnames <- character()
    irnames <- names(ir[[x]])
    csgrad[[x]] <- character(length=length(irnames))
    counter <- 0
    for(n in 1:length(irnames)){
        counter <- counter+1
        temp <- character()
        temp <- irnames[n]
        temp <- as.numeric((as.character(temp)))
        if(temp>=0 && temp<=2) csgrad[[x]][counter] <- colsel[1]
        if(temp>=3 && temp<=5) csgrad[[x]][counter] <- colsel[2]
        if(temp>=6 && temp<=8) csgrad[[x]][counter] <- colsel[3]
        if(temp>=9 && temp<=11) csgrad[[x]][counter] <- colsel[4]
        if(temp>=12) csgrad[[x]][counter] <- colsel[5]
        print(n)
    }
}


ir <- list()
for (n in 1:length(dirs)) {
    system(command = paste0("th ./levdistance.lua -d ", dirs[n]))
    s <- read.csv(file.path(getwd(), dirs[n], "start.csv"), header = FALSE)
    s <- as.numeric(unlist(s))
    e <- read.csv(file.path(getwd(), dirs[n], "end.csv"), header = FALSE)
    e <- as.numeric(unlist(e))
    r <- read.csv(file.path(getwd(), dirs[n], "structures.csv"), header = FALSE)
    r <- as.character(unlist(r))
    m <- read.csv(file.path(getwd(),dirs[1],"lev.csv"), header = FALSE)
    w <- sapply(1:length(s), function(x) length(seq(s[x],e[x], 1)))
    for(i in 1:ncol(m)){
        i <- 1
        ir[[n]] <- IRanges(start = s, end = e, width = w, names = m[,i])    
    }
    
}

for(i in 1:length(ir)) {
    tout <- ir[[i]]
    fileout = paste0(dirs[i],".RData")
    save(tout, file = fileout)
}
plotRanges(ir[[3]])
### ++++++ TODO integrate plotRanges into genomeplot, get syn and nonsyn data
### ++++++ TODO colored diagrams with best lengthwise scan (most hits across
### ++++++      bins).
