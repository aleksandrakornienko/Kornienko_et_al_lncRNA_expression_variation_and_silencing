
transcripts_2=Sys.getenv("transcripts_loc") 
transcripts_3=Sys.getenv("transcripts_tofiltert_loc")
transcripts_4=Sys.getenv("transcripts_filtered")
a<-read.table(transcripts_2, header=F)
b<-read.table(transcripts_3, header=F)
c<-a[!(a$V4 %in% b$V1),]
write.table(c, file=transcripts_4, col.names=F, row.names=F , quote=F, sep="\t")

