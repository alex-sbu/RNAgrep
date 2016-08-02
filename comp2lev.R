#!/usr/bin/env r

### ++++++ options
suppressMessages(library(docopt))
doc <- "Usage: comp2lev.R [-s select] [-f filein] [-h shift]
-s --select STRING Comma seperated column headers which to select [default: all]
-f --filein FILE  Local path to files [default: compared-all.csv]
-h --shift NUMERIC  shift all positions by a constant [default: 0]"
opt <- docopt(doc)
filein <- opt$filein
sel <- opt$select
shift <- opt$shift
shift <- as.numeric(shift)


### ++++++ make output directory
dir <- gsub("\\.csv","", filein)
if (grepl("/", dir)) {
    dir <- gsub("^.+/([[:alpha:]]+$)", "\\1", dir)
}
dir <- file.path(getwd(), dir)
if(sel!="all") dir <- paste0(as.character(dir),"-", gsub(",","-",sel))
if( dir.exists(dir) ) {
    allf <- paste0( dir, "/", list.files(dir) )
} else {
        dir.create(dir) 
}

group <- as.numeric(c(1,2,4))

### ++++++ remove non-conserved structures
d <- read.csv(filein)
if(sel != "all"){
    n <- names(d)
    n <- sapply(n, function(x) gsub("\\..+$","",x))
    n <- as.character(n)
    cat("\n\nCOLNAMES:\n")
    cat(n, sep = ",")
    if(any(grepl(",", sel))) {
        selv <- unlist(strsplit(sel, split = ","))
        selv <- as.character(selv)
    }
    cat("\n\nSELECTED:\n")
    cat(selv, sep = ",")
    cat("\n\n",paste0("OUTPUT DIRECTORY:\n", dir, "\n"))
    selv <- as.character(selv)
    if(any(n%in%selv)){
        d <- d[,as.numeric(which(n%in%selv))]
    }
}


### ++++++ omit non-conserved
d <- na.omit(d)
### ++++++ order by start pos
d <- d[order(d[,2]),]

if (shift!=0) {
    d[,group[1]] <- as.numeric(as.character(d[,group[1]]))+shift
    d[,group[2]] <- as.numeric(as.character(d[,group[2]]))+shift
    print(d)
}

### ++++++ output conserved structures
rp <- file.path( dir, "structures.csv")
if(file.exists(rp)) file.remove(rp)
hp <- as.character(d[,group[3]])
write.table(as.character(hp), file = rp, sep = ",",
            col.names = FALSE, row.names = FALSE)


### ++++++ output start positions of structure intervals
sp <- file.path( dir, "start.csv" )
if(file.exists(sp)) file.remove(sp)
start <- as.numeric(as.character(d[,group[1]]))
write.table(as.numeric(as.character(start)), file = sp, sep = ",",
            col.names = FALSE, row.names = FALSE)


### ++++++ output end positions of structure intervals
ep <- file.path( dir, "end.csv" )
if(file.exists(ep)) file.remove(ep)
fin <- as.numeric(as.character(d[,group[2]]))
write.table(as.numeric(as.character(fin)), file = ep, sep = ",",
            col.names = FALSE, row.names = FALSE)
cat("\n\n SUCCESS\n\n")

### ++++++ TODO: write.table is inserting trailing newlines
