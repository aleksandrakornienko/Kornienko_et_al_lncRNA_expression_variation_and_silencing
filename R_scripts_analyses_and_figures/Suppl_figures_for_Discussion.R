
setwd("Z:/01_POSTDOC/")

# AS-PC anticorrelation candidates 

# upload gene pairs

pairs_AS_Ar11PC<- read.delim("03_Projects/2018_lncRNA_variation_paper/positional_analysis/20211013_annotation/Ar11_PC_loci_closestAS_to_AS_RNAs.bed", header=FALSE)
names(pairs_AS_Ar11PC)<-c("chr_NC","start_NC", "end_NC","gene_NC","score_NC","strand_NC","chr_PC","start_PC", "end_PC","gene_PC","score_PC","strand_PC","distance")


a<-pairs_AS_Ar11PC[,c("gene_NC","gene_PC")]
a<-a[a$gene_PC %in% Araport11.TPMs.genes.1001G$gene,]
length(pairs_AS_Ar11PC$chr_NC)
#9061


a$corr_1001G<-0
a$NC_variance_1001G<-0
a$NC_max_1001G<-0
a$NC_mean_1001G<-0
a$NC_nrlines_1001G<-0
a$PC_variance_1001G<-0
a$PC_max_1001G<-0
a$PC_mean_1001G<-0
a$PC_nrlines_1001G<-0
a$corr_Cortijo<-0
a$corr_1001Gnew<-0
a$NC_max_1001Gnew<-0
a$PC_max_1001Gnew<-0
a$corr_EC<-0
a$NC_max_EC<-0
a$PC_max_EC<-0
a$corr_EC<-0
a$NC_max_EC<-0
a$PC_max_EC<-0
a$corr_EC_pollen<-0
a$NC_max_pollen<-0
a$PC_max_pollen<-0
a$corr_EC_flower<-0
a$NC_max_flower<-0
a$PC_max_flower<-0
a$corr_EC_rosette<-0
a$NC_max_ros<-0
a$PC_max_ros<-0
a$corr_EC_seedling<-0
a$NC_max_seedl<-0
a$PC_max_seedl<-0
a$corr_Vu<-0
a$corr_Vu_WT<-0
a$corr_Vu_stem<-0
a$Vu_mean_ddm1<-0
a$Vu_mean_pol4<-0
a$Vu_mean_WT<-0


for (i in 1:9026) {
 # print(AS)
  AS=as.character(a$gene_NC[i])
  print(AS)
  PC=as.character(a$gene_PC[i])
  a$corr_1001G[i]<-cor(t(denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene==AS,2:462]),t(Araport11.TPMs.genes.1001G[Araport11.TPMs.genes.1001G$gene==PC,2:462]))
  a$NC_variance_1001G[i]<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene==AS]
  a$NC_max_1001G[i]<-denovo2021.TPMs.genes.1001G$max[denovo2021.TPMs.genes.1001G$gene==AS]
  a$NC_nrlines_1001G[i]<-denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene==AS]
  a$PC_variance_1001G[i]<-Araport11.TPMs.genes.1001G$variance[Araport11.TPMs.genes.1001G$gene==PC]
  a$PC_max_1001G[i]<-Araport11.TPMs.genes.1001G$max[Araport11.TPMs.genes.1001G$gene==PC]
  a$PC_nrlines_1001G[i]<-Araport11.TPMs.genes.1001G$Nacc_where_expressed05[Araport11.TPMs.genes.1001G$gene==PC]
  
  
  a$corr_1001Gnew[i]<-cor(t(denovo2021.TPMs.genes.1001Gnew[denovo2021.TPMs.genes.1001Gnew$gene==AS,2:88]),t(Araport11.TPMs.genes.1001Gnew[Araport11.TPMs.genes.1001Gnew$gene==PC,2:88]))
  a$corr_Cortijo[i]<-cor(t(denovo2021.TPMs.genes.Cortijo[denovo2021.TPMs.genes.Cortijo$gene==AS,2:169]),t(Araport11.TPMs.genes.Cortijo[Araport11.TPMs.genes.Cortijo$gene==PC,2:169]))
  
  
  a$corr_EC[i]<-cor(t(denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene==AS,2:97]),t(Araport11.TPMs.genes.ERACAPS[Araport11.TPMs.genes.ERACAPS$gene==PC,2:97]))
  a$NC_max_EC[i]<-denovo2021.TPMs.genes.ERACAPS$max[denovo2021.TPMs.genes.ERACAPS$gene==AS]
  a$PC_max_EC[i]<-Araport11.TPMs.genes.ERACAPS$max[Araport11.TPMs.genes.ERACAPS$gene==PC]
  
  a$corr_EC_pollen[i]<-cor(t(denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene==AS,grep("P.",names(denovo2021.TPMs.genes.ERACAPS))]),t(Araport11.TPMs.genes.ERACAPS[Araport11.TPMs.genes.ERACAPS$gene==PC,grep("P.",names(Araport11.TPMs.genes.ERACAPS))]))
  a$NC_max_pollen[i]<-denovo2021.TPMs.genes.ERACAPS$max.pollen[denovo2021.TPMs.genes.ERACAPS$gene==AS]
  a$PC_max_pollen[i]<-Araport11.TPMs.genes.ERACAPS$max.pollen[Araport11.TPMs.genes.ERACAPS$gene==PC]
  
  a$corr_EC_flower[i]<-cor(t(denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene==AS,grep("F.",names(denovo2021.TPMs.genes.ERACAPS))]),t(Araport11.TPMs.genes.ERACAPS[Araport11.TPMs.genes.ERACAPS$gene==PC,grep("F.",names(Araport11.TPMs.genes.ERACAPS))]))
  a$NC_max_flower[i]<-denovo2021.TPMs.genes.ERACAPS$max.flower[denovo2021.TPMs.genes.ERACAPS$gene==AS]
  a$PC_max_flower[i]<-Araport11.TPMs.genes.ERACAPS$max.flower[Araport11.TPMs.genes.ERACAPS$gene==PC]
  
  a$corr_EC_rosette[i]<-cor(t(denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene==AS,grep("R",names(denovo2021.TPMs.genes.ERACAPS))]),t(Araport11.TPMs.genes.ERACAPS[Araport11.TPMs.genes.ERACAPS$gene==PC,grep("R",names(Araport11.TPMs.genes.ERACAPS))]))
  a$NC_max_ros[i]<-denovo2021.TPMs.genes.ERACAPS$max.rosette[denovo2021.TPMs.genes.ERACAPS$gene==AS]
  a$PC_max_ros[i]<-Araport11.TPMs.genes.ERACAPS$max.rosette[Araport11.TPMs.genes.ERACAPS$gene==PC]
  
  
  a$corr_EC_seedling[i]<-cor(t(denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene==AS,grep("S",names(denovo2021.TPMs.genes.ERACAPS))]),t(Araport11.TPMs.genes.ERACAPS[Araport11.TPMs.genes.ERACAPS$gene==PC,grep("S",names(Araport11.TPMs.genes.ERACAPS))]))
  a$NC_max_seedl[i]<-denovo2021.TPMs.genes.ERACAPS$max.seedl[denovo2021.TPMs.genes.ERACAPS$gene==AS]
  a$PC_max_seedl[i]<-Araport11.TPMs.genes.ERACAPS$max.seedl[Araport11.TPMs.genes.ERACAPS$gene==PC]
  
  a$corr_Vu[i]<-cor(t(Vu_denovo_TPM[Vu_denovo_TPM$gene==AS,2:27]),t(Vu_Araport11_TPM[Vu_Araport11_TPM$gene==PC,2:27]))
  a$corr_Vu_WT[i]<-cor(t(Vu_denovo_TPM[Vu_denovo_TPM$gene==AS,grep("WT",names(Vu_denovo_TPM))]),t(Vu_Araport11_TPM[Vu_Araport11_TPM$gene==PC,grep("WT",names(Vu_denovo_TPM))]))
  a$corr_Vu_stem[i]<-cor(t(Vu_denovo_TPM[Vu_denovo_TPM$gene==AS,grep("stem",names(Vu_denovo_TPM))]),t(Vu_Araport11_TPM[Vu_Araport11_TPM$gene==PC,grep("stem",names(Vu_denovo_TPM))]))
  a$Vu_mean_ddm1[i]<-cor(t(Vu_denovo_TPM[Vu_denovo_TPM$gene==AS,grep("ddm1",names(Vu_denovo_TPM))]),t(Vu_Araport11_TPM[Vu_Araport11_TPM$gene==PC,grep("ddm1",names(Vu_denovo_TPM))]))
  a$Vu_mean_pol4[i]<-cor(t(Vu_denovo_TPM[Vu_denovo_TPM$gene==AS,grep("pol4",names(Vu_denovo_TPM))]),t(Vu_Araport11_TPM[Vu_Araport11_TPM$gene==PC,grep("pol4",names(Vu_denovo_TPM))]))
  a$Vu_mean_WT[i]<-cor(t(Vu_denovo_TPM[Vu_denovo_TPM$gene==AS,grep("WT",names(Vu_denovo_TPM))]),t(Vu_Araport11_TPM[Vu_Araport11_TPM$gene==PC,grep("WT",names(Vu_denovo_TPM))]))
  
} 


for (i in 1:9026) {
  # print(AS)
  AS=as.character(a$gene_NC[i])
  print(AS)
  PC=as.character(a$gene_PC[i])
  a$NC_max_1001Gnew[i]<-denovo2021.TPMs.genes.1001Gnew$ma_x[denovo2021.TPMs.genes.1001Gnew$gene==AS]
  a$PC_max_1001Gnew[i]<-Araport11.TPMs.genes.1001Gnew$ma_x[Araport11.TPMs.genes.1001Gnew$gene==PC]
  
}

tmp<-denovo2021.TPMs.genes.ERACAPS[,names(denovo2021.TPMs.genes.ERACAPS)!="F.22004"]
tmpAr<-Araport11.TPMs.genes.ERACAPS[,names(Araport11.TPMs.genes.ERACAPS)!="F.22004"]
# F.22004 - outlier!! exclude! 
for (i in 1:9026) {
  # print(AS)
  AS=as.character(a$gene_NC[i])
  print(AS)
  PC=as.character(a$gene_PC[i])
  a$corr_EC_flower[i]<-cor(t(tmp[tmp$gene==AS,grep("F.",names(tmp))]),t(tmpAr[tmpAr$gene==PC,grep("F.",names(tmpAr))]))
  a$NC_max_flower[i]<-tmp$max.flower[tmp$gene==AS]
  a$PC_max_flower[i]<-tmpAr$max.flower[tmpAr$gene==PC]
}
  
  
hist (a$corr_1001G)
hist (a$corr_1001Gnew)
hist (a$corr_EC_seedling)
hist (a$corr_EC_rosette)
hist (a$corr_EC_flower)
hist (a$corr_EC_pollen)



AS="CUFF_NC.4367"
PC="AT2G35810"
  CUFF_NC.4367
AT2G35810


CUFF_NC.8035
AT4G18390

AS="CUFF_NC.8035"
PC="AT4G18390"



AS="CUFF_NC.417"
PC="AT1G12050"
AT1G12050
plot(t(denovo2021.TPMs.genes.1001Gnew[denovo2021.TPMs.genes.1001Gnew$gene==AS,2:88]),t(Araport11.TPMs.genes.1001Gnew[Araport11.TPMs.genes.1001Gnew$gene==PC,2:88]))



a[a$corr_1001G<(-0.4)&a$NC_max_1001G>1 & a$PC_max_1001G>1,c("gene_NC","gene_PC","corr_1001G","NC_max_1001G","PC_max_1001G"),]
a[a$corr_1001Gnew<(-0.4)&!is.na(a$corr_1001Gnew)&a$NC_max_1001Gnew>1 & a$PC_max_1001Gnew>1,c("gene_NC","gene_PC","corr_1001G","corr_1001Gnew","PC_max_1001Gnew","NC_max_1001Gnew"),]


a[a$corr_EC_seedling<(-0.4)&a$NC_max_seedl>1 & a$PC_max_seedl>1,c("gene_NC","gene_PC","corr_EC_seedling","NC_max_seedl","PC_max_seedl"),]

a[a$corr_EC_rosette<(-0.4)&a$NC_max_ros>1 & a$PC_max_ros>1,c("gene_NC","gene_PC","corr_EC_rosette","NC_max_ros","PC_max_ros"),]

a[a$corr_EC_flower<(-0.4)&a$NC_max_flower>1 & a$PC_max_flower>1,c("gene_NC","gene_PC","corr_EC_flower","NC_max_flower","PC_max_flower"),]

a[a$corr_EC_pollen<(-0.4)&a$NC_max_pollen>1 & a$PC_max_pollen>1,c("gene_NC","gene_PC","corr_EC_pollen","NC_max_pollen","PC_max_pollen"),]

cand_1001<-a$gene_NC[a$corr_1001G<(-0.4)&a$NC_max_1001G>1 & a$PC_max_1001G>1]
cand_1001new<-a$gene_NC[a$corr_1001Gnew<(-0.4)&!is.na(a$corr_1001Gnew)&a$NC_max_1001Gnew>1 & a$PC_max_1001Gnew>1]
cand_seedl<-a$gene_NC[a$corr_EC_seedling<(-0.4)&a$NC_max_seedl>1 & a$PC_max_seedl>1]
cand_ros<-a$gene_NC[a$corr_EC_rosette<(-0.4)&a$NC_max_ros>1 & a$PC_max_ros>1]
cand_fl<-a$gene_NC[a$corr_EC_flower<(-0.4)&a$NC_max_flower>1 & a$PC_max_flower>1]
cand_pol<-a$gene_NC[a$corr_EC_pollen<(-0.4)&a$NC_max_pollen>1 & a$PC_max_pollen>1]

intersect(cand_1001,cand_1001new)
intersect(cand_1001,cand_ros)
intersect(cand_1001,cand_seedl)
intersect(cand_1001,cand_fl)
intersect(cand_1001,cand_pol)
intersect(cand_ros,cand_pol)
intersect(cand_ros,cand_seedl)

library(scales)
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Discussion/Supplem_CUFF_NC.6822_vs_partnergene.1001G.pdf",height = 3.5,width = 3.5)
AS="CUFF_NC.6822"
PC="AT3G54363"
plot(t(denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene==AS,2:462]),t(Araport11.TPMs.genes.1001G[Araport11.TPMs.genes.1001G$gene==PC,2:462]), main="Expression in 461 accessions, rosette", ylab="AT3G54363 expression, TPM", xlab="CUFF_NC.6822 expression, TPM",pch=19, col=alpha("black",alpha=0.2) )
dev.off()


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Discussion/Supplem_CUFF_NC.6822_vs_partnergene.EC_pollen.pdf",height = 3.5,width = 3.5)
AS="CUFF_NC.6822"
PC="AT3G54363"
plot(t(denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene==AS,grep("P.",names(denovo2021.TPMs.genes.ERACAPS))]),t(Araport11.TPMs.genes.ERACAPS[Araport11.TPMs.genes.ERACAPS$gene==PC,grep("P.",names(Araport11.TPMs.genes.ERACAPS))]), main="Expression in 23 accessions, pollen", ylab="AT3G54363 expression, TPM", xlab="CUFF_NC.6822 expression, TPM",pch=19, col=alpha("black",alpha=0.4) )
dev.off()





pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Discussion/Supplem_CUFF_NC.6822_vs_partnergene.EC_pollen.pdf",height = 3.5,width = 3.5)
AS="CUFF_NC.6822"
PC="AT3G54363"
plot(t(denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene==AS,grep("P.",names(denovo2021.TPMs.genes.ERACAPS))]),t(Araport11.TPMs.genes.ERACAPS[Araport11.TPMs.genes.ERACAPS$gene==PC,grep("P.",names(Araport11.TPMs.genes.ERACAPS))]), main="Expression in 23 accessions, pollen", ylab="AT3G54363 expression, TPM", xlab="CUFF_NC.6822 expression, TPM",pch=19, col=alpha("black",alpha=0.4) )
dev.off()



library(scales)
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Discussion/Supplem_CUFF_NC.4367_vs_partnergene.1001G.pdf",height = 3.5,width = 3.5)
AS="CUFF_NC.4367"
PC="AT2G35810"
plot(t(denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene==AS,2:462]),t(Araport11.TPMs.genes.1001G[Araport11.TPMs.genes.1001G$gene==PC,2:462]), main="Expression in 461 accessions, rosette", ylab="AT2G35810 expression, TPM", xlab="CUFF_NC.4367 expression, TPM",pch=19, col=alpha("black",alpha=0.2) )
dev.off()

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Discussion/Supplem_CUFF_NC.11117_vs_partnergene.1001G.pdf",height = 3.5,width = 3.5)
AS="CUFF_NC.11117"
PC="AT5G59660"
plot(t(denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene==AS,2:462]),t(Araport11.TPMs.genes.1001G[Araport11.TPMs.genes.1001G$gene==PC,2:462]), main="Expression in 461 accessions, rosette", ylab="AT5G59660 expression, TPM", xlab="CUFF_NC.11117 expression, TPM",pch=19, col=alpha("black",alpha=0.2) )
dev.off()


CUFF_NC.7318	AT4G03100


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Discussion/Supplem_CUFF_NC.7318_vs_partnergene.EC_seedl.pdf",height = 3.5,width = 3.5)
AS="CUFF_NC.7318"
PC="AT4G03100"
plot(t(denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene==AS,grep("S.",names(denovo2021.TPMs.genes.ERACAPS))]),t(Araport11.TPMs.genes.ERACAPS[Araport11.TPMs.genes.ERACAPS$gene==PC,grep("S.",names(Araport11.TPMs.genes.ERACAPS))]), main="Expression in 25 accessions, 9-leaf rosette", ylab="AT4G03100 expression, TPM", xlab="CUFF_NC.7318 expression, TPM",pch=19, col=alpha("black",alpha=0.4) )
dev.off()


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Discussion/Supplem_CUFF_NC.1636_vs_partnergene.EC_rosette.pdf",height = 3.5,width = 3.5)
AS="CUFF_NC.1636"
PC="AT1G48970"
plot(t(denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene==AS,grep("R.",names(denovo2021.TPMs.genes.ERACAPS))]),t(Araport11.TPMs.genes.ERACAPS[Araport11.TPMs.genes.ERACAPS$gene==PC,grep("R.",names(Araport11.TPMs.genes.ERACAPS))]), main="Expression in 25 accessions, 9-leaf rosette", ylab="AT1G48970 expression, TPM", xlab="CUFF_NC.1636 expression, TPM",pch=19, col=alpha("black",alpha=0.4) )
dev.off()


CUFF_NC.2380	AT1G67790


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Discussion/Supplem_CUFF_NC.2380_vs_partnergene.EC_flowers.pdf",height = 3.5,width = 3.5)
AS="CUFF_NC.2380"
PC="AT1G67790"
plot(t(denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene==AS,grep("F.",names(denovo2021.TPMs.genes.ERACAPS))]),t(Araport11.TPMs.genes.ERACAPS[Araport11.TPMs.genes.ERACAPS$gene==PC,grep("F.",names(Araport11.TPMs.genes.ERACAPS))]), main="Expression in 22 accessions, flowers", ylab="AT1G67790 expression, TPM", xlab="CUFF_NC.2380 expression, TPM",pch=19, col=alpha("black",alpha=0.4) )
dev.off()




# other genes and TE pieces 



linc_noTE<-lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10==0] 
linc_50TE<-lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10>0 & lincRNAs_TE_coverage_anystrand$TAIR10<=0.5] 
linc_50_80TE<-lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10>0.5 & lincRNAs_TE_coverage_anystrand$TAIR10<=0.8] 
linc_80TE<-lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10>0.8] 

length(linc_noTE) #1070
length(linc_50TE)#703
length(linc_50_80TE)#187
length(linc_80TE)#286

as_noTE<-AS_TE_coverage_anystrand$gene[AS_TE_coverage_anystrand$TAIR10==0] 
as_50TE<-AS_TE_coverage_anystrand$gene[AS_TE_coverage_anystrand$TAIR10>0 & AS_TE_coverage_anystrand$TAIR10<=0.5] 
as_50_80TE<-AS_TE_coverage_anystrand$gene[AS_TE_coverage_anystrand$TAIR10>0.5 & AS_TE_coverage_anystrand$TAIR10<=0.8]
as_80TE<-AS_TE_coverage_anystrand$gene[AS_TE_coverage_anystrand$TAIR10>0.8] 

length(as_noTE) #6834
length(as_50TE)#1236
length(as_50_80TE)#80
length(as_80TE)#41

pc_noTE<-PC_TE_coverage_anystrand$gene[PC_TE_coverage_anystrand$TAIR10==0] 
pc_50TE<-PC_TE_coverage_anystrand$gene[PC_TE_coverage_anystrand$TAIR10>0 & PC_TE_coverage_anystrand$TAIR10<=0.5] 
pc_50_80TE<-PC_TE_coverage_anystrand$gene[PC_TE_coverage_anystrand$TAIR10>0.5 & PC_TE_coverage_anystrand$TAIR10<=0.8]
pc_80TE<-PC_TE_coverage_anystrand$gene[PC_TE_coverage_anystrand$TAIR10>0.8] 

length(pc_noTE) #17919
length(pc_50TE)#4872
length(pc_50_80TE)#214
length(pc_80TE)#134


#1001G
length(denovo2021.TPMs.genes.1001G$X6909[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% as_noTE]) #591
length(denovo2021.TPMs.genes.1001G$X6909[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% as_50TE])#151
length(denovo2021.TPMs.genes.1001G$X6909[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% as_50_80TE])#3
length(denovo2021.TPMs.genes.1001G$X6909[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% as_80TE])#0

#1001Gnew
length(denovo2021.TPMs.genes.1001Gnew$mean.6909[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% as_noTE]) #392
length(denovo2021.TPMs.genes.1001Gnew$mean.6909[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% as_50TE])#108
length(denovo2021.TPMs.genes.1001Gnew$mean.6909[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% as_50_80TE])#7
length(denovo2021.TPMs.genes.1001Gnew$mean.6909[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% as_80TE])#0
#pollen
length(denovo2021.TPMs.genes.ERACAPS$P.6909[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% as_noTE]) #622
length(denovo2021.TPMs.genes.ERACAPS$P.6909[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% as_50TE]) #127
length(denovo2021.TPMs.genes.ERACAPS$P.6909[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% as_50_80TE]) #7
length(denovo2021.TPMs.genes.ERACAPS$P.6909[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% as_80TE]) #6



#1001G
length(denovo2021.TPMs.genes.1001G$X6909[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% pc_noTE]) #12234
length(denovo2021.TPMs.genes.1001G$X6909[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% pc_50TE])#3082
length(denovo2021.TPMs.genes.1001G$X6909[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% pc_50_80TE])#56
length(denovo2021.TPMs.genes.1001G$X6909[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% pc_80TE])#11

#1001Gnew
length(denovo2021.TPMs.genes.1001Gnew$mean.6909[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% pc_noTE]) #12074
length(denovo2021.TPMs.genes.1001Gnew$mean.6909[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% pc_50TE])#2981
length(denovo2021.TPMs.genes.1001Gnew$mean.6909[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% pc_50_80TE])#56
length(denovo2021.TPMs.genes.1001Gnew$mean.6909[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% pc_80TE])#8

#pollen
length(denovo2021.TPMs.genes.ERACAPS$P.6909[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% pc_noTE]) #5579
length(denovo2021.TPMs.genes.ERACAPS$P.6909[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% pc_50TE]) #1332
length(denovo2021.TPMs.genes.ERACAPS$P.6909[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% pc_50_80TE]) #34
length(denovo2021.TPMs.genes.ERACAPS$P.6909[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% pc_80TE]) #13



########################################################
# do the TE pieces affect epigenetic state of lincRNAs? 
########################################################




#chipseq 

# chipseq coverage in lincRNA with different TE content 

#K9 boxplot PC + AS + lincs + TE genes
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Discussion/K9_cov_AS_PC_TEcontent.pdf",height = 4,width = 5)
###########################################################################

par(mar=c(6,4,4,2))
a1<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% pc_noTE]
a2<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% pc_50TE]
a3<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% pc_50_80TE]
a4<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% pc_80TE]

a12<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% as_noTE]
a22<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% as_50TE]
a32<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% as_50_80TE]
a42<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% as_80TE]

a13<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% linc_noTE]
a23<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% linc_50TE]
a33<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% linc_50_80TE]
a43<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% linc_80TE]
a5<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene]
boxplot(   a1,a2,a3,a4,-10, a12,a22,a32,a42,-10, a13,a23,a33,a43,-10,a5 ,main="H3K9me2 level\n(Col-0 rosette)",cex.main=1.2,ylim=c(-2,2.2),               col=c("#486EB4",  "#395891", "#395891",  "#395891",  "#395891",     "#90C473",  "#698f54",  "#698f54", "#698f54", "#698f54",   "#F2AB54","#C55B16","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","","noTE","TE<50%","TE50-80%","TE>80%","","noTE","TE<50%","TE50-80%","TE>80%","","TE genes"),las=2, notch = T, outline = F, ylab="log2(CHIP/INPUT)")
mtext (text = "lincRNAs",side = 3,at = 13,line = 0)
mtext (text = "AS lncRNAs",side = 3,at = 7,line = 0)
mtext (text = "PC genes",side = 3,at = 2,line = 0)

#################
#add p values   #
#################
d<-wilcox.test(a1,a2)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=2)
d<-wilcox.test(a2,a3)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=2)
d<-wilcox.test(a3,a4)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=2)

d<-wilcox.test(a12,a22)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=2)

d<-wilcox.test(a32,a22)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=2)
d<-wilcox.test(a32,a42)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=2)

d<-wilcox.test(a13,a23)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=11.5,y=2)
d<-wilcox.test(a23,a33)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=12.5,y=2)
d<-wilcox.test(a33,a43)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=13.5,y=2)

#################
dev.off()



#H1 boxplot PC + AS + lincs + TE genes
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Discussion/H1_cov_AS_PC_TEcontent.pdf",height = 4,width = 5)
###########################################################################

par(mar=c(6,4,4,2))
a1<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% pc_noTE]
a2<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% pc_50TE]
a3<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% pc_50_80TE]
a4<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% pc_80TE]

a12<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% as_noTE]
a22<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% as_50TE]
a32<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% as_50_80TE]
a42<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% as_80TE]

a13<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% linc_noTE]
a23<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% linc_50TE]
a33<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% linc_50_80TE]
a43<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% linc_80TE]
a5<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene]
boxplot(   a1,a2,a3,a4,-10, a12,a22,a32,a42,-10, a13,a23,a33,a43,-10,a5 ,main="H1 level\n(Col-0 rosette)",cex.main=1.2,ylim=c(-2,2.2),               col=c("#486EB4",  "#395891", "#395891",  "#395891",  "#395891",     "#90C473",  "#698f54",  "#698f54", "#698f54", "#698f54",   "#F2AB54","#C55B16","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","","noTE","TE<50%","TE50-80%","TE>80%","","noTE","TE<50%","TE50-80%","TE>80%","","TE genes"),las=2, notch = T, outline = F, ylab="log2(CHIP/INPUT)")
mtext (text = "lincRNAs",side = 3,at = 13,line = 0)
mtext (text = "AS lncRNAs",side = 3,at = 7,line = 0)
mtext (text = "PC genes",side = 3,at = 2,line = 0)

#################
#add p values   #
#################
d<-wilcox.test(a1,a2)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=2)
d<-wilcox.test(a2,a3)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=2)
d<-wilcox.test(a3,a4)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=2)

d<-wilcox.test(a12,a22)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=2)

d<-wilcox.test(a32,a22)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=2)
d<-wilcox.test(a32,a42)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=2)

d<-wilcox.test(a13,a23)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=11.5,y=2)
d<-wilcox.test(a23,a33)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=12.5,y=2)
d<-wilcox.test(a33,a43)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=13.5,y=2)

#################
dev.off()




#CG boxplot PC + AS + lincs + TE genes
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Discussion/CG_AS_PC_TEcontent.pdf",height = 4,width = 5)
###########################################################################

par(mar=c(6,4,4,2))
a1<-CG.1001new.denovo$mean.6909[ CG.1001new.denovo$transcript %in% pc_noTE]
a2<-CG.1001new.denovo$mean.6909[ CG.1001new.denovo$transcript %in% pc_50TE]
a3<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% pc_50_80TE]
a4<-CG.1001new.denovo$mean.6909[ CG.1001new.denovo$transcript %in% pc_80TE]

a12<-CG.1001new.denovo$mean.6909[ CG.1001new.denovo$transcript %in% as_noTE]
a22<-CG.1001new.denovo$mean.6909[ CG.1001new.denovo$transcript %in% as_50TE]
a32<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% as_50_80TE]
a42<-CG.1001new.denovo$mean.6909[ CG.1001new.denovo$transcript %in% as_80TE]

a13<-CG.1001new.denovo$mean.6909[ CG.1001new.denovo$transcript %in% linc_noTE]
a23<-CG.1001new.denovo$mean.6909[ CG.1001new.denovo$transcript %in% linc_50TE]
a33<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% linc_50_80TE]
a43<-CG.1001new.denovo$mean.6909[ CG.1001new.denovo$transcript %in% linc_80TE]
a5<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]
boxplot(   a1,a2,a3,a4,-10, a12,a22,a32,a42,-10, a13,a23,a33,a43,-10,a5 ,main="CG methylation level\n(Col-0 rosette)",cex.main=1.2,ylim=c(0,1.1),               col=c("#486EB4",  "#395891", "#395891",  "#395891",  "#395891",     "#90C473",  "#698f54",  "#698f54", "#698f54", "#698f54",   "#F2AB54","#C55B16","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","","noTE","TE<50%","TE50-80%","TE>80%","","noTE","TE<50%","TE50-80%","TE>80%","","TE genes"),las=2, notch = T, outline = F, ylab="methylation level")
mtext (text = "lincRNAs",side = 3,at = 13,line = 0)
mtext (text = "AS lncRNAs",side = 3,at = 7,line = 0)
mtext (text = "PC genes",side = 3,at = 2,line = 0)

#################
#add p values   #
#################
d<-wilcox.test(a1,a2)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
d<-wilcox.test(a2,a3)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
d<-wilcox.test(a3,a4)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)

d<-wilcox.test(a12,a22)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=1)

d<-wilcox.test(a32,a22)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
d<-wilcox.test(a32,a42)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=1)

d<-wilcox.test(a13,a23)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=11.5,y=1)
d<-wilcox.test(a23,a33)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=12.5,y=1)
d<-wilcox.test(a33,a43)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=13.5,y=1)

#################
dev.off()


#CHH boxplot PC + AS + lincs + TE genes
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Discussion/CHH_AS_PC_TEcontent.pdf",height = 4,width = 5)
###########################################################################

par(mar=c(6,4,4,2))
a1<-CHH.1001new.denovo$mean.6909[ CHH.1001new.denovo$transcript %in% pc_noTE]
a2<-CHH.1001new.denovo$mean.6909[ CHH.1001new.denovo$transcript %in% pc_50TE]
a3<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% pc_50_80TE]
a4<-CHH.1001new.denovo$mean.6909[ CHH.1001new.denovo$transcript %in% pc_80TE]

a12<-CHH.1001new.denovo$mean.6909[ CHH.1001new.denovo$transcript %in% as_noTE]
a22<-CHH.1001new.denovo$mean.6909[ CHH.1001new.denovo$transcript %in% as_50TE]
a32<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% as_50_80TE]
a42<-CHH.1001new.denovo$mean.6909[ CHH.1001new.denovo$transcript %in% as_80TE]

a13<-CHH.1001new.denovo$mean.6909[ CHH.1001new.denovo$transcript %in% linc_noTE]
a23<-CHH.1001new.denovo$mean.6909[ CHH.1001new.denovo$transcript %in% linc_50TE]
a33<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% linc_50_80TE]
a43<-CHH.1001new.denovo$mean.6909[ CHH.1001new.denovo$transcript %in% linc_80TE]
a5<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% TE_genes.loci$gene]
boxplot(   a1,a2,a3,a4,-10, a12,a22,a32,a42,-10, a13,a23,a33,a43,-10,a5 ,main="CHH methylation level\n(Col-0 rosette)",cex.main=1.2,ylim=c(0,0.32),               col=c("#486EB4",  "#395891", "#395891",  "#395891",  "#395891",     "#90C473",  "#698f54",  "#698f54", "#698f54", "#698f54",   "#F2AB54","#C55B16","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","","noTE","TE<50%","TE50-80%","TE>80%","","noTE","TE<50%","TE50-80%","TE>80%","","TE genes"),las=2, notch = T, outline = F, ylab="methylation level")
mtext (text = "lincRNAs",side = 3,at = 13,line = 0)
mtext (text = "AS lncRNAs",side = 3,at = 7,line = 0)
mtext (text = "PC genes",side = 3,at = 2,line = 0)

#################
#add p values   #
#################
d<-wilcox.test(a1,a2)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.2)
d<-wilcox.test(a2,a3)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.2)
d<-wilcox.test(a3,a4)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.2)

d<-wilcox.test(a12,a22)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=0.2)

d<-wilcox.test(a32,a22)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.2)
d<-wilcox.test(a32,a42)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=0.2)

d<-wilcox.test(a13,a23)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=11.5,y=0.2)
d<-wilcox.test(a23,a33)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=12.5,y=0.2)
d<-wilcox.test(a33,a43)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=13.5,y=0.2)

#################
dev.off()


#24nt boxplot PC + AS + lincs + TE genes
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Discussion/24nt_AS_PC_TEcontent.pdf",height = 4,width = 5)
###########################################################################

par(mar=c(6,4,4,2))
a1<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% pc_noTE]
a2<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% pc_50TE]
a3<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% pc_50_80TE]
a4<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% pc_80TE]

a12<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% as_noTE]
a22<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% as_50TE]
a32<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% as_50_80TE]
a42<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% as_80TE]

a13<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% linc_noTE]
a23<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% linc_50TE]
a33<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% linc_50_80TE]
a43<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% linc_80TE]
a5<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% TE_genes.loci$gene]
boxplot(   a1,a2,a3,a4,-10, a12,a22,a32,a42,-10, a13,a23,a33,a43,-10,a5 ,main="24nt sRNA coverage\n(Col-0 flowers)",cex.main=1.2,ylim=c(0,2.5),               col=c("#486EB4",  "#395891", "#395891",  "#395891",  "#395891",     "#90C473",  "#698f54",  "#698f54", "#698f54", "#698f54",   "#F2AB54","#C55B16","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","","noTE","TE<50%","TE50-80%","TE>80%","","noTE","TE<50%","TE50-80%","TE>80%","","TE genes"),las=2, notch = T, outline = F, ylab="24nt coverage, RPM")
mtext (text = "lincRNAs",side = 3,at = 13,line = 0)
mtext (text = "AS lncRNAs",side = 3,at = 7,line = 0)
mtext (text = "PC genes",side = 3,at = 2,line = 0)

#################
#add p values   #
#################
d<-wilcox.test(a1,a2)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.2)
d<-wilcox.test(a2,a3)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.2)
d<-wilcox.test(a3,a4)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.2)

d<-wilcox.test(a12,a22)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=0.2)

d<-wilcox.test(a32,a22)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.2)
d<-wilcox.test(a32,a42)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=0.2)

d<-wilcox.test(a13,a23)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=11.5,y=0.2)
d<-wilcox.test(a23,a33)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=12.5,y=0.2)
d<-wilcox.test(a33,a43)
d<-d$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=13.5,y=0.2)

#################
dev.off()











