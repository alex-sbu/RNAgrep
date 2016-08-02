### ++++++ CUSTOM VARS
shift <- 740
dfnames <- c("pos", "SD23","MAX23","P3L37","23127","w2","lansing","CAV20-IH35")
filein <- "/files/projects/polio/assembly_of_the_poliovirion/supplemental/polyprotein_na_aligned_wt_nogap.fasta"
refname <- "pvm_wt_polyprotein_na"
fileout <- paste0("-MUTS_", gsub(".+/(.+)\\.fasta$", "\\1", filein),
                  "_REF-", refname, ".csv", collapse = "")

suppressPackageStartupMessages(library("seqinr"))





splitbycodon <- function(seq) {
    seq <- c2s(seq)
    seq <- paste0(unlist(seq), collapse = "")
    size <- 3
    cuts <- nchar(seq)%/%size
    seq <- substring(seq, 1, (cuts*size))
    splitseq <- unlist(strsplit(seq, split=""))
    seqmatrix <- matrix(splitseq, nrow=cuts, ncol=size, byrow=TRUE)
    cutseq <- apply(seqmatrix, 1, paste0, collapse="")
    return(cutseq)
}

### ++++++ TODO: manually manage gapping/alignment
seqs <- read.fasta(filein)
ref <- seqs[paste0(refname)]
ref <- unlist(ref)
seqnames <- names(seqs)
seqs <- seqs[-which(seqnames == refname)]
seqnames <- seqnames[-which(seqnames == refname)]
refaa <- seqinr::translate(ref)
aaseqs <- list()
for(i in 1:length(seqs)) aaseqs[[i]] <- seqinr::translate(seqs[[i]])
names(aaseqs) <- seqnames
codonseqs <- lapply(seqs, splitbycodon)
codonref <- splitbycodon(ref)

mutL <- list()
for(i in 1:length(codonseqs)) {
    mutL[[i]] <- character()
    counter <- 0
    for(j in 1:length(codonseqs[[i]])) {
        aamut <- logical()
        aamut <- aaseqs[[i]][j] == refaa[j]
        naseq <- character()
        naseq <- as.character(unlist(strsplit(codonseqs[[i]][j], split="")))
        naref <- character()
        naref <- as.character(unlist(strsplit(codonref[j], split="")))
        type <- character()
        type <- ifelse(aamut == TRUE, "S", "N")
        for(k in 1:length(naseq)) {
            counter <- counter+1
            if(naseq[k] == naref[k]){
                if(aamut == TRUE){
                    mutL[[i]][counter] <- "NOMUTATION"
                }else{
                    mutL[[i]][counter] <- type
                }
            }else{
                if(aamut == TRUE){
                    mutL[[i]][counter] <- type
                }else{
                    mutL[[i]][counter] <- type
                }
            }
            ## cat("\n\n",naseq[k],"   ",
            ## naref[k],"    ",
            ## aamut , ">>>>>>>>         ",
            ## mutL[[i]][counter], "\n\n")
        }
    }
}
names(mutL) <- seqnames

m <- matrix(unlist(mutL), ncol = length(mutL), nrow = length(mutL[[1]]))
pos <- 1:length(mutL[[1]])
if(!is.null(shift)){
    pos <- pos+shift
}




### ++++++ Synonymous output
df <- data.frame(pos, m, stringsAsFactors = FALSE)
names(df) <- seqnames
df[which(df == "NOMUTATION", arr.ind = TRUE)] <- 1
df[which(df == "N", arr.ind = TRUE)] <- 1
df[which(df == "S", arr.ind = TRUE)] <- 0
df <- apply(df,2,as.numeric)
if(!is.null(dfnames)){
    colnames(df) <- dfnames
}
synfileout <- paste0("SYN", fileout)
write.table(df,
            file = synfileout,
            row.names = FALSE,
            append = FALSE,
            sep = ",",
            eol = "\n",
            na = "NA",
            col.names = TRUE)





### ++++++ nonsynonymous output
df <- data.frame(pos, m, stringsAsFactors = FALSE)
names(df) <- seqnames
df[which(df == "NOMUTATION", arr.ind = TRUE)] <- 1
df[which(df == "S", arr.ind = TRUE)] <- 1
df[which(df == "N", arr.ind = TRUE)] <- 0
df <- apply(df,2,as.numeric)
if(!is.null(dfnames)){
    colnames(df) <- dfnames
}
nonsynfileout <- paste0("NONSYN", fileout)
write.table(df,
            file = nonsynfileout,
            row.names = FALSE,
            append = FALSE,
            sep = ",",
            eol = "\n",
            na = "NA",
            col.names = TRUE)
