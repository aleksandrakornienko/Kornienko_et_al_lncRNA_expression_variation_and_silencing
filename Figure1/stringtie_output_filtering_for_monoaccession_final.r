
bedloc<-Sys.getenv("bedloc")
abund<-Sys.getenv("abundance")
namepool<-Sys.getenv("name_pool")
workdir <-Sys.getenv("working_folder")
gtfloc <- Sys.getenv("gtf_loc")
artef<-Sys.getenv("artefacts")
gtf_filtered_loc<-Sys.getenv("gtf_loc_filter")
bed_filtered_loc<-Sys.getenv("bedloc_filter")


bed<-read.table (bedloc,header=F,  sep="\t",quote="")
abundance<-read.table (abund,header=F,  sep="\t",quote="")
gtf<-read.table (gtfloc,header=F,  sep="\t",quote="")
try({artefacts<-read.table (artef,header=F,  sep="\t",quote="")})
#names(bed)<-c("chr","start","end","name","score","strand",)

if(exists("artefacts")) {
	bed_with_TPM<-merge (bed,abundance[,c("V4","V9")],by="V4")
	#filter annotation for lowly expressed transcripts (TPM<0.5) and for short single-exon transcripts(with very high TPM cut-off) 
	bed_prefiltered<-bed_with_TPM[bed_with_TPM[,13]>=0.5 & !(bed_with_TPM$V10==1&(bed_with_TPM$V3-bed_with_TPM$V2)<400&bed_with_TPM[,13]<5)& !(bed_with_TPM$V10==1&(bed_with_TPM$V3-bed_with_TPM$V2)<200), ] 
	bed_filtered <- bed_prefiltered[!(bed_prefiltered$V4 %in% artefacts$V1),c(2,3,4,1,5,6,7,8,9,10,11,12)]

	#filter gtf 
	gtf_filtered<-gtf[gtf$V1 %in% bed_filtered$V4,2:10]
	write.table(gtf_filtered, gtf_filtered_loc,sep="\t", row.names=F, col.names=F,quote=F, eol="\n")
	write.table(bed_filtered, bed_filtered_loc,sep="\t", row.names=F, col.names=F,quote=F, eol="\n")
} else {
	bed_with_TPM<-merge (bed,abundance[,c("V4","V9")],by="V4")
	#filter annotation for lowly expressed transcripts (TPM<0.5) and for short single-exon transcripts(with very high TPM cut-off) 
	bed_prefiltered<-bed_with_TPM[bed_with_TPM[,13]>=0.5 & !(bed_with_TPM$V10==1&(bed_with_TPM$V3-bed_with_TPM$V2)<400&bed_with_TPM[,13]<5)& !(bed_with_TPM$V10==1&(bed_with_TPM$V3-bed_with_TPM$V2)<200), ] 
	bed_filtered <- bed_prefiltered[,c(2,3,4,1,5,6,7,8,9,10,11,12)]

	#filter gtf 
	gtf_filtered<-gtf[gtf$V1 %in% bed_filtered$V4,2:10]
	write.table(gtf_filtered, gtf_filtered_loc,sep="\t", row.names=F, col.names=F,quote=F, eol="\n")
	write.table(bed_filtered, bed_filtered_loc,sep="\t", row.names=F, col.names=F,quote=F, eol="\n")
}



