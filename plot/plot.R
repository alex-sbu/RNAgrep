### ++++++ TODO: integrate plotRanges into genedplot modules
### ++++++ TODO: save iranges output, move synplot output
### ++++++ TODO: plotit
### ++++++ TODO: figure legend
### ++++++ TODO: Results blurb
### ++++++ TODO: Supplementary (Pairwise comp, with colored min per range bin)
### ++++++ TODO: get max no. of lengthwise motifs (at n mismatches)
### ++++++ TODO: discussion
### ++++++ TODO: publish


################################################################################
### Gene plotter v2.0
################################################################################
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(sfsmisc))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(caTools))
suppressPackageStartupMessages(library(IRanges))

################################################################################
### Load plot functions
################################################################################
source("drawings.R")
source("modules.R")
source("axis.R")
source("axistitle.R")
plotmods <- function(module, ...) do.call(module, ...)





################################################################################
### Sequence features & data input
################################################################################
genestart <- 3386
geneend <- 7440
allseqs <- read.csv("polyprotein_na_aligned.fasta")
positions <- "diagram-P2P3.csv"
nonsynonmuts <- "NONSYN-MUTS_P2P3_na_aligned_wt_nogap_REF-pvm_wt_polyprotein_na.csv"
synonmuts <- "SYN-MUTS_P2P3_na_aligned_wt_nogap_REF-pvm_wt_polyprotein_na.csv"
synplotPVAL <- "synplot-pval-plot.RData"
synplotOE <- "synplot-obs-exp-plot.RData"
naturalmotifs <- "compared-natural.RData"
syntheticmotifs <- "compared-synthetic.RData"
mars <- c(0,0,0,1.35)
gint = c(genestart, geneend)

modlist <- list(diagram, synP, synOE, motifplot, motifplot, barcode, barcode)
inputlist <- list(
    list(input = positions, mars = mars, xlimit = gint),
    list(input = synplotPVAL, mars = mars, xlimit = gint),
    list(input = synplotOE, mars = mars, xlimit = gint),
    list(input = naturalmotifs, mars = mars, xlimit = gint),
    list(input = syntheticmotifs, mars = mars, xlimit = gint),
    list(input = synonmuts, mars = mars, xlimit = gint),
    list(input = nonsynonmuts, mars = mars, xlimit = gint))




################################################################################
### Plotting area config
## if splitOverride is TRUE then the split matrix for split.screen is saved
## as splitscreen.csv, and if mod_splitscreen.csv exists, it overwrites default
## split.screen positions.
################################################################################
axiswidth <- 1
labelcex <- 1
labelpos <- 1
mtickby <- 100
btickby <- 500
axispos <- 1
xlabels <- c(genestart, geneend)
count <- length(inputlist)
## count <- 5
alphalab <- toupper(rep(letters, count))
splitOverride <- TRUE
drawfig <- "pdf"

## dev.off()



################################################################################
### Print figure
## drawfig must be either pdf, png, or svg, or else nothing is drawn
################################################################################
width <- 9.68504
height <- 6.85039
if(!is.null(drawfig)) {
    if(drawfig == "pdf"){
        pdf(paste0("genomeplot.pdf"),
            width = width, height = height,
            pointsize = 9, family = "sans", colormodel = "rgb")
    }else if(drawfig == "png"){
        png(filename = paste0("genomeplot.png"),
            width = width, height = height,
            units = "in", pointsize = 9,
            bg = "white",  res = 800)
    }else if(drawfig == "svg"){
        svg(filename = paste0("genomeplot.svg"),
            width = width, height = height, pointsize = 8,
            bg = "white")
    }
}

par(mar = c(0,0,0,0))
par(oma = c(0,0,0,0))




################################################################################
### Generate split.screen plot
################################################################################
if(file.exists("mod_splitscreen.csv")) {
    f <- readLines("mod_splitscreen.csv")
    f <- sapply(f, strsplit, ",")
    f <- unlist(f)
    f <- as.numeric(f)
    f <- matrix(f, ncol = 4, nrow = count, byrow = TRUE)
} else {
    f <- matrix(0, ncol = 4, nrow = count)
    f[,1] <- 0.1                            # left
    f[,2] <- 0.9                            # right
    f[nrow(f), 3] <- 0.1                    # bottom
    f[1, 4] <- 0.95                         # top
    divf <- 0.85/count
    seqf <- seq(f[nrow(f), 3], f[1, 4], divf)
    seqf <- seqf[-1]
    seqf <- seqf[-length(seqf)]
    seqf <- rev(seqf)
    counter <- 0
    for(i in 2:nrow(f)){
        counter <- counter+1
        f[i, 4] <- seqf[counter]
        f[(i-1), 3] <- seqf[counter]
    }
}
if(splitOverride) {
    write.table(f, file = "splitscreen.csv", sep = ",", append = FALSE, col.names = FALSE, row.names = FALSE)
}

print("\n\n")
print("PANEL SPLIT:\n")
print(f)
print("\n")













################################################################################
### Plot area
################################################################################
split.screen(f)
for(i in 1:count) {
    screen(i)
    xlimit <- plotmods(modlist[[i]], inputlist[[i]])

    
    mtext(side = 2,                     # margin labels
          line = 1,
          padj = 0,
          text = alphalab[i],
          las = 2,
          cex = 1.8)
    
    if(i == count) {
        minitick <- c(seq(xlimit[1], xlimit[2], by = mtickby), xlimit[2])
        bigtick <- c(seq(xlimit[1], xlimit[2], by = btickby), xlimit[2])
        bigticklabel <- c(seq(xlabels[1], xlabels[2], by = btickby), xlabels[2])
        minitick<- minitick[-which(minitick%in%bigtick)]
        axis(1,
             at = bigtick,
             labels = FALSE,
             tck = -.075,
             line = 0,
             lwd = 1)
        axis(1,
             at = minitick,
             labels = FALSE,
             tck = -0.05,
             line = 0,
             lwd = 1)        
        axis(1,
             at = bigtick,
             labels = bigticklabel,
             tick = FALSE,
             line= -0.5,
             cex.axis=labelcex,
             )
        mtext("Nucleotide position", cex = 1, side = 1, padj = 4, line = -0.8)
    }

    ## if(i == 1) {
    ##     input <- inputlist[[1]]
    ##     data <- read.csv(input)
    ##     data <- na.omit(data)
    ##     xlimit <- c(min(c(data$x1, data$x2)),
    ##                 max(c(data$x1, data$x2)))            
    ##     minitick <- c(seq(xlimit[1], xlimit[2], by = mtickby), xlimit[2])
    ##     bigtick <- c(seq(xlimit[1], xlimit[2], by = btickby), xlimit[2])
    ##     minitick<- minitick[-which(minitick%in%bigtick)]
    ##     axis(3,
    ##          at = bigtick,
    ##          labels = FALSE,
    ##          tck = .075,
    ##          line = 0,
    ##          lwd = 1)        
    ##     axis(3,
    ##          at = minitick,
    ##          labels = FALSE,
    ##          tck = 0.05,
    ##          line = 0,
    ##          lwd = 1)
    ## }    
}







if(!is.null(drawfig)) {
    dev.off()
}
