library(scales)

args=commandArgs(trailing=TRUE)
assoc <- args[[1]]
link<- args[[2]]
genes <- args[[3]]
snps <- args[[4]]
out <- args [[5]]

assoc <- read.table(assoc, stringsAsFactors = FALSE, header = TRUE)
link <- read.table(link, stringsAsFactors = FALSE, header = TRUE)
genes <- read.delim(genes, stringsAsFactors = FALSE, header = TRUE)
snps <- read.table(snps, stringsAsFactors = FALSE, header = FALSE, col.names="SNP")
compl<-assoc[assoc$SNP %in% snps$SNP,]
sign <- compl[compl$P<0.00000005,]

setwd(out)

source("/home/nardone/software/LocusZooms/functions/locus_zoom.R")

for (line in dim(sign)[1]){
chr=sign[line,1]
snp=sign[line,2]
start=sign[line,3]
end=start + 1
locus.zoom(data = assoc,
           region = c(chr, start, end),
           offset_bp = 250000,
           ld.file = link,
           genes.data = genes,
           plot.title = "",
           file.name = paste(snp,".jpg", sep=""),
           secondary.snp = snp,
           secondary.label = TRUE)
}
