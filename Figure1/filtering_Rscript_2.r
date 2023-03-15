
CPC=Sys.getenv("CPC_loc") 
CPC_coding=Sys.getenv("CPC_coding_loc") 

transcripts_2=Sys.getenv("transcripts_loc") 
transcripts_3=Sys.getenv("transcripts_filt_loc")
transcripts_4=Sys.getenv("transcripts_coding_loc")

a<-read.table(CPC, header=F)
b<-read.table(transcripts_2, header=F)
c<-b[b$V4 %in% a$V1,]
write.table(c, file=transcripts_3, col.names=F, row.names=F , quote=F, sep="\t")

a<-read.table(CPC_coding, header=F)
b<-read.table(transcripts_2, header=F)
c<-b[b$V4 %in% a$V1,]
write.table(c, file=transcripts_4, col.names=F, row.names=F , quote=F, sep="\t")
