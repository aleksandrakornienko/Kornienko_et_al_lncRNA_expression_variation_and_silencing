lincRNAs_TE_coverage.TE_opposite_strand.TAIR10 <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/lincRNAs_TE_coverage.TE_opposite_strand.TAIR10.bed", header=FALSE)

lincRNAs_TE_coverage.TE_same_strand.TAIR10 <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/lincRNAs_TE_coverage.TE_same_strand.TAIR10.bed", header=FALSE)

lincRNAs_TE_coverage.TAIR10 <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/lincRNAs_TE_coverage.TAIR10.bed", header=FALSE)

names(lincRNAs_TE_coverage.TE_same_strand.TAIR10)<-c("gene","start_of_TE_piece","end_of_TE_piece","TE_type","direction_re_gene","proportion_of_gene","TE_length","start_of_gene","end_of_gene")
names(lincRNAs_TE_coverage.TE_opposite_strand.TAIR10)<-c("gene","start_of_TE_piece","end_of_TE_piece","TE_type","direction_re_gene","proportion_of_gene","TE_length","start_of_gene","end_of_gene")
names(lincRNAs_TE_coverage.TAIR10)<-c("gene","start_of_TE_piece","end_of_TE_piece","TE_type","direction_re_gene","proportion_of_gene","TE_length","start_of_gene","end_of_gene")

AS_TE_coverage.TAIR10 <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/AS_TE_coverage.TAIR10.bed", header=FALSE)
names(AS_TE_coverage.TAIR10)<-c("gene","start_of_TE_piece","end_of_TE_piece","TE_type","direction_re_gene","proportion_of_gene","TE_length","start_of_gene","end_of_gene")
PC_TE_coverage.TAIR10 <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/PC_TE_coverage.TAIR10.bed", header=FALSE)
names(PC_TE_coverage.TAIR10)<-c("gene","start_of_TE_piece","end_of_TE_piece","TE_type","direction_re_gene","proportion_of_gene","TE_length","start_of_gene","end_of_gene")
Ar11PC_TE_coverage.TAIR10 <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/Ar11PC_TE_coverage.TAIR10.bed", header=FALSE)
names(Ar11PC_TE_coverage.TAIR10)<-c("gene","start_of_TE_piece","end_of_TE_piece","TE_type","direction_re_gene","proportion_of_gene","TE_length","start_of_gene","end_of_gene")

TEgene_TE_coverage.TAIR10 <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/te_TE_coverage.TAIR10.bed", header=FALSE)
names(TEgene_TE_coverage.TAIR10)<-c("gene","start_of_TE_piece","end_of_TE_piece","TE_type","direction_re_gene","proportion_of_gene","TE_length","start_of_gene","end_of_gene")



linc_TE_cov_loci<-read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/lincRNAs_TE_coverage.TAIR10.summarized_per_locus.bed", header=FALSE)
names(linc_TE_cov_loci)<-c("gene","start","end","N_of_non_overlapping_TE_patches","total_length_of_TE_patches","length_of_gene","coverage")



lincTSS_TE_cov_loci<-read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/lincTSS_TE_coverage.TAIR10.summarized_per_locus.bed", header=FALSE)
names(lincTSS_TE_cov_loci)<-c("gene","start","end","N_of_non_overlapping_TE_patches","total_length_of_TE_patches","length_of_gene","coverage")




lincTES_TE_cov_loci<-read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/lincTES_TE_coverage.TAIR10.summarized_per_locus.bed", header=FALSE)
names(lincTES_TE_cov_loci)<-c("gene","start","end","N_of_non_overlapping_TE_patches","total_length_of_TE_patches","length_of_gene","coverage")




linc_TE_cov_all_loci_2cols<-linc_TE_cov_loci[,c(1,7)]
addition<-as.data.frame(cbind(as.vector(lncRNAs.intergenic.loci$gene[!(lncRNAs.intergenic.loci$gene %in% linc_TE_cov_all_loci_2cols$gene)]),0))
names(addition)<-c("gene","coverage")
linc_TE_cov_all_loci_2cols<-rbind(as.data.frame(linc_TE_cov_loci[,c(1,7)]),as.data.frame(addition))
linc_TE_cov_all_loci_2cols$coverage<-as.numeric(linc_TE_cov_all_loci_2cols$coverage)

as_TE_cov_loci<-read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/AS_TE_coverage.TAIR10.summarized_per_locus.bed", header=FALSE)
names(as_TE_cov_loci)<-c("gene","start","end","N_of_non_overlapping_TE_patches","total_length_of_TE_patches","length_of_gene","coverage")
as_TE_cov_all_loci_2cols<-as_TE_cov_loci[,c(1,7)]
addition<-as.data.frame(cbind(as.vector(lncRNAs.antisense.loci$gene[!(lncRNAs.antisense.loci$gene %in% as_TE_cov_all_loci_2cols$gene)]),0))
names(addition)<-c("gene","coverage")
as_TE_cov_all_loci_2cols<-rbind(as.data.frame(as_TE_cov_loci[,c(1,7)]),as.data.frame(addition))
as_TE_cov_all_loci_2cols$coverage<-as.numeric(as_TE_cov_all_loci_2cols$coverage)



pc_TE_cov_loci<-read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/PC_TE_coverage.TAIR10.summarized_per_locus.bed", header=FALSE)
names(pc_TE_cov_loci)<-c("gene","start","end","N_of_non_overlapping_TE_patches","total_length_of_TE_patches","length_of_gene","coverage")

Ar11pc_TE_cov_loci<-read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/Ar11PC_TE_coverage.TAIR10.summarized_per_locus.bed", header=FALSE)
names(Ar11pc_TE_cov_loci)<-c("gene","start","end","N_of_non_overlapping_TE_patches","total_length_of_TE_patches","length_of_gene","coverage")

te_TE_cov_loci<-read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/te_TE_coverage.TAIR10.summarized_per_locus.bed", header=FALSE)
names(te_TE_cov_loci)<-c("gene","start","end","N_of_non_overlapping_TE_patches","total_length_of_TE_patches","length_of_gene","coverage")


# nr of patches

nr_of_TEpatches <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/TE_content_20211013_annotation/nr_of_TEpatches.txt")





### tables with TAIR10 and 27 genomes


lincRNAs_TE_coverage_forward <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/lincRNAs_TE_coverage_forward.bed") 
lincRNAs_TE_coverage_reverse <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/lincRNAs_TE_coverage_reverse.bed") 
rownames(lincRNAs_TE_coverage_forward)<-lincRNAs_TE_coverage_forward$gene
lincRNAs_TE_coverage_forward<-as.data.frame(t(lincRNAs_TE_coverage_forward[,2:2247]))
rownames(lincRNAs_TE_coverage_reverse)<-lincRNAs_TE_coverage_reverse$gene
lincRNAs_TE_coverage_reverse<-as.data.frame(t(lincRNAs_TE_coverage_reverse[,2:2247]))


lincRNAs_TE_coverage_forward$sd<-apply(lincRNAs_TE_coverage_forward[2:28],1,sd,na.rm=T)
lincRNAs_TE_coverage_forward$mean<-apply(lincRNAs_TE_coverage_forward[2:28],1,mean,na.rm=T)
lincRNAs_TE_coverage_forward$min<-apply(lincRNAs_TE_coverage_forward[2:28],1,min,na.rm=T)
lincRNAs_TE_coverage_forward$max<-apply(lincRNAs_TE_coverage_forward[2:28],1,max,na.rm=T)
#lincRNAs_TE_coverage_forward$TE_presense_freq<-apply(lincRNAs_TE_coverage_forward[2:28],1,function(i) sum(i>0,na.rm=TRUE))
lincRNAs_TE_coverage_forward$TE_presense_nrlines<-apply(lincRNAs_TE_coverage_forward[2:28],1,function(i) sum(i>0,na.rm=TRUE))
lincRNAs_TE_coverage_forward$gene_presence_nrlines<-apply(lincRNAs_TE_coverage_forward[2:28],1,function(i) sum(!is.na(i)))
lincRNAs_TE_coverage_forward$gene<-rownames(lincRNAs_TE_coverage_forward)



lincRNAs_TE_coverage_reverse$sd<-apply(lincRNAs_TE_coverage_reverse[2:28],1,sd,na.rm=T)
lincRNAs_TE_coverage_reverse$mean<-apply(lincRNAs_TE_coverage_reverse[2:28],1,mean,na.rm=T)
lincRNAs_TE_coverage_reverse$min<-apply(lincRNAs_TE_coverage_reverse[2:28],1,min,na.rm=T)
lincRNAs_TE_coverage_reverse$max<-apply(lincRNAs_TE_coverage_reverse[2:28],1,max,na.rm=T)
#lincRNAs_TE_coverage_reverse$TE_presense_freq<-apply(lincRNAs_TE_coverage_reverse[2:28],1,function(i) sum(i>0,na.rm=TRUE))
lincRNAs_TE_coverage_reverse$TE_presense_nrlines<-apply(lincRNAs_TE_coverage_reverse[2:28],1,function(i) sum(i>0,na.rm=TRUE))
lincRNAs_TE_coverage_reverse$gene_presence_nrlines<-apply(lincRNAs_TE_coverage_reverse[2:28],1,function(i) sum(!is.na(i)))
lincRNAs_TE_coverage_reverse$gene<-rownames(lincRNAs_TE_coverage_reverse)




lincRNAs_TE_coverage_anystrand <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/lincRNAs_TE_coverage_anystrand.bed",header = T,na.strings = "") 
rownames(lincRNAs_TE_coverage_anystrand)<-lincRNAs_TE_coverage_anystrand$gene
lincRNAs_TE_coverage_anystrand<-as.data.frame(t(lincRNAs_TE_coverage_anystrand[,2:2247]))


AS_TE_coverage_anystrand <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/ASlncRNAs_TE_coverage_anystrand.bed") 
#remove genes with mistakes: 
a<-lncRNAs.antisense.loci[grep("-",lncRNAs.antisense.loci$gene),] 
a$g<-lapply(strsplit(as.character(a$gene),"-", fixed = T), "[",1)
AS_TE_coverage_anystrand<-AS_TE_coverage_anystrand[,!(names(AS_TE_coverage_anystrand) %in% a$g) ]

#AS_TE_coverage_anystrand<-AS_TE_coverage_anystrand1[,1:100]
rownames(AS_TE_coverage_anystrand)<-AS_TE_coverage_anystrand$gene
AS_TE_coverage_anystrand<-AS_TE_coverage_anystrand[,2:8192]
AS_TE_coverage_anystrand<-as.data.frame(t(AS_TE_coverage_anystrand))


PC_TE_coverage_anystrand <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/PC_TE_coverage_anystrand.bed") 
#remove genes with mistakes: 
a<-denovoPC.loci[grep("-",denovoPC.loci$gene),] 
a$g<-lapply(strsplit(as.character(a$gene),"-", fixed = T), "[",1)
PC_TE_coverage_anystrand<-PC_TE_coverage_anystrand[,!(names(PC_TE_coverage_anystrand) %in% a$g) ]

rownames(PC_TE_coverage_anystrand)<-PC_TE_coverage_anystrand$gene
PC_TE_coverage_anystrand<-PC_TE_coverage_anystrand[,2:23140]
PC_TE_coverage_anystrand<-as.data.frame(t(PC_TE_coverage_anystrand))


TEgene_TE_coverage_anystrand <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/TEgenes_TE_coverage_anystrand.bed") 
rownames(TEgene_TE_coverage_anystrand)<-TEgene_TE_coverage_anystrand$gene
TEgene_TE_coverage_anystrand<-as.data.frame(t(TEgene_TE_coverage_anystrand[,2:2131]))


lincRNAs_TE_coverage_anystrand$sd<-apply(lincRNAs_TE_coverage_anystrand[2:28],1,sd,na.rm=T)
lincRNAs_TE_coverage_anystrand$mean<-apply(lincRNAs_TE_coverage_anystrand[2:28],1,mean,na.rm=T)
lincRNAs_TE_coverage_anystrand$min<-apply(lincRNAs_TE_coverage_anystrand[2:28],1,min,na.rm=T)
lincRNAs_TE_coverage_anystrand$max<-apply(lincRNAs_TE_coverage_anystrand[2:28],1,max,na.rm=T)
#lincRNAs_TE_coverage_anystrand$TE_presense_freq<-apply(lincRNAs_TE_coverage_anystrand[2:28],1,function(i) sum(i>0,na.rm=TRUE))
lincRNAs_TE_coverage_anystrand$TE_presense_nrlines<-apply(lincRNAs_TE_coverage_anystrand[2:28],1,function(i) sum(i>0,na.rm=TRUE))
lincRNAs_TE_coverage_anystrand$gene_presence_nrlines<-apply(lincRNAs_TE_coverage_anystrand[2:28],1,function(i) sum(!is.na(i)))
lincRNAs_TE_coverage_anystrand$gene<-rownames(lincRNAs_TE_coverage_anystrand)

AS_TE_coverage_anystrand$sd<-apply(AS_TE_coverage_anystrand[2:28],1,sd,na.rm=T)
AS_TE_coverage_anystrand$mean<-apply(AS_TE_coverage_anystrand[2:28],1,mean,na.rm=T)
AS_TE_coverage_anystrand$min<-apply(AS_TE_coverage_anystrand[2:28],1,min,na.rm=T)
AS_TE_coverage_anystrand$max<-apply(AS_TE_coverage_anystrand[2:28],1,max,na.rm=T)
#AS_TE_coverage_anystrand$TE_presense_freq<-apply(AS_TE_coverage_anystrand[2:28],1,function(i) sum(i>0,na.rm=TRUE))
AS_TE_coverage_anystrand$TE_presense_nrlines<-apply(AS_TE_coverage_anystrand[2:28],1,function(i) sum(i>0,na.rm=TRUE))
AS_TE_coverage_anystrand$gene_presence_nrlines<-apply(AS_TE_coverage_anystrand[2:28],1,function(i) sum(!is.na(i)))
AS_TE_coverage_anystrand$gene<-rownames(AS_TE_coverage_anystrand)

PC_TE_coverage_anystrand$sd<-apply(PC_TE_coverage_anystrand[2:28],1,sd,na.rm=T)
PC_TE_coverage_anystrand$mean<-apply(PC_TE_coverage_anystrand[2:28],1,mean,na.rm=T)
PC_TE_coverage_anystrand$min<-apply(PC_TE_coverage_anystrand[2:28],1,min,na.rm=T)
PC_TE_coverage_anystrand$max<-apply(PC_TE_coverage_anystrand[2:28],1,max,na.rm=T)
#PC_TE_coverage_anystrand$TE_presense_freq<-apply(PC_TE_coverage_anystrand[2:28],1,function(i) sum(i>0,na.rm=TRUE))
PC_TE_coverage_anystrand$TE_presense_nrlines<-apply(PC_TE_coverage_anystrand[2:28],1,function(i) sum(i>0,na.rm=TRUE))

PC_TE_coverage_anystrand$gene_presence_nrlines<-apply(PC_TE_coverage_anystrand[2:28],1,function(i) sum(!is.na(i)))
PC_TE_coverage_anystrand$gene<-rownames(PC_TE_coverage_anystrand)


TEgene_TE_coverage_anystrand$sd<-apply(TEgene_TE_coverage_anystrand[2:28],1,sd,na.rm=T)
TEgene_TE_coverage_anystrand$mean<-apply(TEgene_TE_coverage_anystrand[2:28],1,mean,na.rm=T)
TEgene_TE_coverage_anystrand$min<-apply(TEgene_TE_coverage_anystrand[2:28],1,min,na.rm=T)
TEgene_TE_coverage_anystrand$max<-apply(TEgene_TE_coverage_anystrand[2:28],1,max,na.rm=T)
TEgene_TE_coverage_anystrand$gene_presence_nrlines<-apply(TEgene_TE_coverage_anystrand[2:28],1,function(i) sum(!is.na(i)))
TEgene_TE_coverage_anystrand$TE_presense_nrlines<-apply(TEgene_TE_coverage_anystrand[2:28],1,function(i) sum(i>0,na.rm=TRUE))
#TEgene_TE_coverage_anystrand$TE_presense_freq<-apply(TEgene_TE_coverage_anystrand[2:28],1,function(i) sum(i>0,na.rm=TRUE))

TEgene_TE_coverage_anystrand$gene<-rownames(TEgene_TE_coverage_anystrand)








hist(PC_TE_coverage_anystrand$TE_presense_freq,breaks=40)
hist(lincRNAs_TE_coverage_anystrand$TE_presense_freq,breaks=40)
boxplot(lincRNAs_TE_coverage_anystrand$sd,PC_TE_coverage_anystrand$sd)
vioplot(lincRNAs_TE_coverage_anystrand$TE_presense_freq[lincRNAs_TE_coverage_anystrand$TE_presense_freq>0],
        PC_TE_coverage_anystrand$TE_presense_freq[PC_TE_coverage_anystrand$TE_presense_freq>0])
boxplot(lincRNAs_TE_coverage_anystrand$TE_presense_freq[lincRNAs_TE_coverage_anystrand$TE_presense_freq>0],
        PC_TE_coverage_anystrand$TE_presense_freq[PC_TE_coverage_anystrand$TE_presense_freq>0])


library(pheatmap)

pheatmap(lincRNAs_TE_coverage_anystrand[complete.cases(lincRNAs_TE_coverage_anystrand),3:10],na_col = "gray")
pheatmap(lincRNAs_TE_coverage_anystrand[complete.cases(lincRNAs_TE_coverage_anystrand[,3:10])&apply(lincRNAs_TE_coverage_anystrand[,3:10],1,max,na.rm=TRUE)>0.4,3:10],na_col = "gray",scale = "row")

a<-is.na(lincRNAs_TE_coverage_anystrand[lincRNAs_TE_coverage_anystrand$max>0,3:10])
