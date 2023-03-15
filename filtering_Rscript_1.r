
genes_strand=Sys.getenv("genes_strand_location") 
genes=Sys.getenv("genes_location") 
a<-read.table(genes_strand, header=F)
b<-read.table(genes, header=F)
c<-merge(a,b, by="V1")
write.table(c, file=genes, col.names=F, row.names=F , quote=F, sep="\t")

