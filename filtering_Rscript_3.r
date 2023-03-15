
bidirTSS=Sys.getenv("bidirTSS_loc") 
transcripts_5=Sys.getenv("transcripts_loc") 
transcripts_bidir=Sys.getenv("transcripts_bidir_loc")
a<-read.table(bidirTSS, header=F)
b<-read.table(transcripts_5, header=F)
c<-b[b$V4 %in% a$V4,]
write.table(c, file=transcripts_bidir, col.names=F, row.names=F , quote=F, sep="\t")

