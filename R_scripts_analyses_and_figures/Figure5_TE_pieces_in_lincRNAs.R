############# 
# Figure 5 TE content of lincRNAs 
###############

TE_types_sorted=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/TAIR10_TEs_name_7superfam_length.sorted.txt


#how many lincRNAs have a TE piece? 
#no TE piece 
length(lncRNAs.intergenic.loci$gene[!(lncRNAs.intergenic.loci$gene %in% unique(lincRNAs_TE_coverage.TAIR10$gene))])
# 1070 lincRNA loci 

#any strand
length(unique(lincRNAs_TE_coverage.TAIR10$gene))
# 1176 lincRNA loci 

# same strand 
length(unique(lincRNAs_TE_coverage.TE_same_strand.TAIR10$gene))
# 866
# opposite strand 
length(unique(lincRNAs_TE_coverage.TE_opposite_strand.TAIR10$gene))
# 938

#only same strand
length(unique(lincRNAs_TE_coverage.TE_same_strand.TAIR10$gene[!(lincRNAs_TE_coverage.TE_same_strand.TAIR10$gene %in% lincRNAs_TE_coverage.TE_opposite_strand.TAIR10$gene)]))
#238
#both same and opposite strands
length(unique(lincRNAs_TE_coverage.TE_same_strand.TAIR10$gene[lincRNAs_TE_coverage.TE_same_strand.TAIR10$gene %in% lincRNAs_TE_coverage.TE_opposite_strand.TAIR10$gene]))
#628
#only opposite strand
length(unique(lincRNAs_TE_coverage.TE_opposite_strand.TAIR10$gene[!(lincRNAs_TE_coverage.TE_opposite_strand.TAIR10$gene %in% lincRNAs_TE_coverage.TE_same_strand.TAIR10$gene)]))
#310

# other genes 



# TE content of genes 

hist (lincRNAs_TE_coverage_anystrand$TAIR10[lincRNAs_TE_coverage_anystrand$TAIR10>0]*100, breaks=10)
hist (AS_TE_coverage_anystrand$TAIR10*100, breaks=10)
hist (pc_TE_cov_loci$coverage*100, breaks=10)
hist (Ar11pc_TE_cov_loci$coverage*100, breaks=10)
hist (te_TE_cov_loci$coverage*100, breaks=10)



#for the main figure - lincRNAs 

#histogram TE content in lincRNAs - any direction (sense and antisense to lincRNA direction)
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_lincs_TEcontent_histogram.pdf",height = 3,width = 3)
############################################################
par(mar=c(4,4,4,2)) 
hist (lincRNAs_TE_coverage_anystrand$TAIR10[lincRNAs_TE_coverage_anystrand$TAIR10>0]*100, breaks=10, main="TE content\n of lincRNA loci", ylab="number of loci", xlab='% locus length\ncovered by TE patches',xaxt='n', col="#F2AB54",las=2,panel.first=grid(lty = 1, col = "darkgray"),ylim=c(0,250))
axis(side=1, at=seq(0,100, 10), labels=seq(0,100,10),las=2)
############################################################
dev.off()
#histogram TE content in lincRNAs - forward direction (sense to lincRNA direction)
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_lincs_TEcontent.FORW.histogram.pdf",height = 3,width = 3)
############################################################
par(mar=c(4,4,4,2)) 
hist (lincRNAs_TE_coverage_forward$TAIR10[lincRNAs_TE_coverage_forward$TAIR10>0]*100, breaks=10, main="TE content of lincRNA loci\nsense direction", ylab="number of loci", xlab='% locus length\ncovered by TE patches',xaxt='n', col="#F2AB54",las=2,panel.first=grid(lty = 1, col = "darkgray",lwd=0.6),ylim=c(0,250),lwd=0.6)
axis(side=1, at=seq(0,100, 10), labels=seq(0,100,10),las=2)
############################################################
dev.off()

#histogram TE content in lincRNAs - reverse direction (antisense to lincRNA direction)
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_lincs_TEcontent.Rev.histogram.pdf",height = 3,width = 3)
############################################################
par(mar=c(4,4,4,2)) 
hist (lincRNAs_TE_coverage_reverse$TAIR10[lincRNAs_TE_coverage_reverse$TAIR10>0]*100, breaks=10, main="TE content of lincRNA loci\nantisense direction", ylab="number of loci", xlab='% locus length\ncovered by TE patches',xaxt='n', col="#F2AB54",las=2,panel.first=grid(lty = 1, col = "darkgray",lwd=0.6),ylim=c(0,250),lwd=0.6)
axis(side=1, at=seq(0,100, 10), labels=seq(0,100,10),las=2)
############################################################
dev.off()


#



#What is the length of the TE piece inside lincRNAs? 
# TE patch length in lincRNAs (and AS and PC and TE genes)
############################################################
#histogram TE patch length in lincRNAs
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_lincs_TE_patch_length.pdf",height = 3,width = 3.2)
############################################################
par(mar=c(4,4,4,2)) 

hist (ifelse(lincRNAs_TE_coverage.TAIR10$TE_length>500,600,lincRNAs_TE_coverage.TAIR10$TE_length), breaks=10, main="length of TE patches\n within lincRNA loci", ylab="number of TE patches", xlab='TE patch length,bp',xaxt='n', col="#F2AB54",las=2,ylim=c(0,2000),panel.first=grid(lty = 1, col = "darkgray",lwd=0.6), xlim=c(0,600))
axis(side=1, at=seq(0,500, 50), labels=seq(0,500,50),las=2,lwd=0.6)
axis(side=1, at=575, labels=">500",las=2)
text(paste("median =", median(lincRNAs_TE_coverage.TAIR10$TE_length),"bp"),x=400,y=1500)
text(paste("min =", min(lincRNAs_TE_coverage.TAIR10$TE_length),"bp"),x=400,y=1300)
############################################################
dev.off()

#histogram TE patch length in AS lncRNAs
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_AS_TE_patch_length.pdf",height = 3,width = 3.2)
############################################################
par(mar=c(4,4,4,2)) 

hist (ifelse(AS_TE_coverage.TAIR10$TE_length>500,600,AS_TE_coverage.TAIR10$TE_length), breaks=10, main="length of TE patches\n within AS lncRNA loci", ylab="number of TE patches", xlab='TE patch length,bp',xaxt='n', col="#90C473",las=2,panel.first=grid(lty = 1, col = "darkgray",lwd=0.6), xlim=c(0,600))
axis(side=1, at=seq(0,500, 50), labels=seq(0,500,50),las=2,lwd=0.6)
axis(side=1, at=575, labels=">500",las=2)
text(paste("median =", median(AS_TE_coverage.TAIR10$TE_length,na.rm = T),"bp"),x=400,y=600)
text(paste("min =", min(AS_TE_coverage.TAIR10$TE_length,na.rm = T),"bp"),x=400,y=500)
############################################################
dev.off()

#histogram TE patch length in PC genes
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_PC_TE_patch_length.pdf",height = 3,width = 3.2)
############################################################
par(mar=c(4,4,4,2)) 

hist (ifelse(PC_TE_coverage.TAIR10$TE_length>500,600,PC_TE_coverage.TAIR10$TE_length), breaks=10, main="length of TE patches\n within PC loci", ylab="number of TE patches", xlab='TE patch length,bp',xaxt='n', col="#486EB4",las=2,panel.first=grid(lty = 1, col = "darkgray",lwd=0.6), xlim=c(0,600))
axis(side=1, at=seq(0,500, 50), labels=seq(0,500,50),las=2,lwd=0.6)
axis(side=1, at=575, labels=">500",las=2)
text(paste("median =", median(PC_TE_coverage.TAIR10$TE_length,na.rm = T),"bp"),x=400,y=2500)
text(paste("min =", min(PC_TE_coverage.TAIR10$TE_length,na.rm = T),"bp"),x=400,y=2300)
############################################################
dev.off()


#histogram TE patch length in TE genes
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_TEgene_TE_patch_length.pdf",height = 3,width = 3.2)
############################################################
par(mar=c(4,4,4,2)) 
hist (ifelse(TEgene_TE_coverage.TAIR10$TE_length>500,600,TEgene_TE_coverage.TAIR10$TE_length), breaks=10, main="length of TE patches\n within TE genes", ylab="number of TE patches", xlab='TE patch length,bp',xaxt='n', col="#673A8E",las=2,panel.first=grid(lty = 1, col = "darkgray",lwd=0.6), xlim=c(0,600))
axis(side=1, at=seq(0,500, 50), labels=seq(0,500,50),las=2,lwd=0.6)
axis(side=1, at=575, labels=">500",las=2)
text(paste("median =", median(TEgene_TE_coverage.TAIR10$TE_length,na.rm = T),"bp"),x=400,y=1500)
text(paste("min =", min(TEgene_TE_coverage.TAIR10$TE_length,na.rm = T),"bp"),x=400,y=1300)
############################################################
dev.off()





col=c("#486EB4","#90C473","#F2AB54","#673A8E")



pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_lincs_TE_patch_number.pdf",height = 3,width = 3.2)
par(mar=c(4,4,4,2)) 
barplot(table(ifelse(nr_of_TEpatches$nr.of.patches>5,6,nr_of_TEpatches$nr.of.patches)),space = 0, main="TE-patch number\n within lincRNA loci", ylab="number of loci", xlab='Number of TE patchs in locus', col="#F2AB54",las=1,ylim=c(0,500))
dev.off()


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_TE_genes_TE_patch_number.pdf",height = 3,width = 3.2)
par(mar=c(4,4,4,2)) 
barplot(table(ifelse(nr_of_TEpatches$nr.of.patches.1>5,6,nr_of_TEpatches$nr.of.patches.1)),space = 0, main="TE-patch number\n within TE genes", ylab="number of loci", xlab='Number of TE patchs in locus', col="#673A8E",las=1,ylim=c(0,1500))
dev.off()




# how many TE seq bp per 1 kb of lincRNA?
mean(lincRNAs_TE_coverage_anystrand$TAIR10[lincRNAs_TE_coverage_anystrand$TAIR10>0])
#0.4671137
mean(lincRNAs_TE_coverage_forward$TAIR10[lincRNAs_TE_coverage_forward$TAIR10>0])
#0.3087895
mean(lincRNAs_TE_coverage_reverse$TAIR10[lincRNAs_TE_coverage_reverse$TAIR10>0])
#0.4360636

mean(TEgene_TE_coverage_anystrand$TAIR10[TEgene_TE_coverage_anystrand$TAIR10>0])
#0.9100271




#SUPPLEMENT 
#what is the normal length of TE fragments? so what we observe inside lincRNAs and other genes - is that just normal TE fragments or really just small pieces of TEs?
################################################################################
# length of annotated TE fragments in TAIR10 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_TAIR10_TEfrags_length.pdf",height = 3,width = 3.2)
par(mar=c(4,4,4,2)) 

hist (ifelse((Ar11_TE_fragments$end -Ar11_TE_fragments$start) >500,600,(Ar11_TE_fragments$end -Ar11_TE_fragments$start)), breaks=10, main="length of TE fragment\nsannotated in TAIR10", ylab="number of TE patches", xlab='TE patch length,bp',xaxt='n', col="mediumpurple1",las=2,panel.first=grid(lty = 1, col = "darkgray",lwd=0.6), xlim=c(0,600))
axis(side=1, at=seq(0,500, 50), labels=seq(0,500,50),las=2,lwd=0.6)
axis(side=1, at=575, labels=">500",las=2)
text(paste("median =", median(Ar11_TE_fragments$end -Ar11_TE_fragments$start),"bp"),x=400,y=8000)
text(paste("min =", min(Ar11_TE_fragments$end -Ar11_TE_fragments$start),"bp"),x=400,y=7000)
dev.off()
################################################################################





#SUPPLEMENT 
# where are TE pieces located within the gene? 
############################################################
#normalize the length of the locus 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4_linc_vs_TEs_TEpatch_position_within_linc_all_dir.pdf",height = 5,width = 4)
############################################################
par(mar=c(4,4,4,2)) 

a<-lincRNAs_TE_coverage.TAIR10
a$start_norm<-100*a$start_of_TE_piece/a$end_of_gene
a$end_norm<-100*a$end_of_TE_piece/a$end_of_gene
a$TE_type_first<-sapply(strsplit(as.character(a$TE_type), "\\,"), `[`, 1)
a$col<-"black"
a$col[a$TE_type_first=="RC_Helitron"]<-"darkblue"
a$col[a$TE_type_first=="SINE_LINE"]<-"green"
a$col[a$TE_type_first=="DNA_MuDR"]<-"blue"
a$col[a$TE_type_first=="LTR_Copia"]<-"brown"
a$col[a$TE_type_first=="LTR_Gypsy"]<-"orange"
a$col[a$TE_type_first=="DNA_other"]<-"lightblue"
a<-a[order(a$col,a$start_norm,a$end_norm,decreasing = F),]

plot(1:3971,1:3971, ylim=c(0,3971), xlim=c(0,100), col="white" ,main="TE patch location within lincRNA gene\n TE patch in the any direction to the gene",ylab="TE patches", xlab="normalized locus length 0-100%")
for (i in 1:3971) {
  segments(a$start_norm[i],i,a$end_norm[i],i,col = as.character(a$col[i]))
}
############################################################

dev.off()



# plot the position of TE pieces (antisense to lincRNA direction) inside the lincRNA loci

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_linc_vs_TEs_TEpatch_position_within_linc_same_dir.pdf",height = 5,width = 4)
############################################################
par(mar=c(4,4,4,2)) 

a<-lincRNAs_TE_coverage.TE_same_strand.TAIR10
a$start_norm<-100*a$start_of_TE_piece/a$end_of_gene
a$end_norm<-100*a$end_of_TE_piece/a$end_of_gene
a$TE_type_first<-sapply(strsplit(as.character(a$TE_type), "\\,"), `[`, 1)
a$col<-"black"
a$col[a$TE_type_first=="RC_Helitron"]<-"darkblue"
a$col[a$TE_type_first=="SINE_LINE"]<-"green"
a$col[a$TE_type_first=="DNA_MuDR"]<-"blue"
a$col[a$TE_type_first=="LTR_Copia"]<-"brown"
a$col[a$TE_type_first=="LTR_Gypsy"]<-"orange"
a$col[a$TE_type_first=="DNA_other"]<-"lightblue"
a<-a[order(a$col,a$start_norm,a$end_norm,decreasing = F),]

plot(1:1901,1:1901, ylim=c(0,1901), xlim=c(0,100), col="white",main="TE patch location within lincRNA gene\n TE patch in the same direction as the gene",ylab="TE patches", xlab="normalized locus length 0-100%")
for (i in 1:1901) {
  segments(a$start_norm[i],i,a$end_norm[i],i,col = as.character(a$col[i]))
}
############################################################
dev.off()



# plot the position of TE pieces (antisense to lincRNA direction) inside the lincRNA loci
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4_linc_vs_TEs_TEpatch_position_within_linc_opp_dir.pdf",height = 5,width = 4)
###########################################################################
par(mar=c(4,4,4,2)) 

a<-lincRNAs_TE_coverage.TE_opposite_strand.TAIR10
a$start_norm<-100*a$start_of_TE_piece/a$end_of_gene
a$end_norm<-100*a$end_of_TE_piece/a$end_of_gene
a$TE_type_first<-sapply(strsplit(as.character(a$TE_type), "\\,"), `[`, 1)
a$col<-"black"
a$col[a$TE_type_first=="RC_Helitron"]<-"darkblue"
a$col[a$TE_type_first=="SINE_LINE"]<-"green"
a$col[a$TE_type_first=="DNA_MuDR"]<-"blue"
a$col[a$TE_type_first=="LTR_Copia"]<-"brown"
a$col[a$TE_type_first=="LTR_Gypsy"]<-"orange"
a$col[a$TE_type_first=="DNA_other"]<-"lightblue"
a<-a[order(a$col,a$start_norm,a$end_norm,decreasing = F),]

plot(1:2070,1:2070, ylim=c(0,2070), xlim=c(0,100), col="white",main="TE patch location within lincRNA gene\n TE patch in the opposite direction to the gene",ylab="TE patches", xlab="normalized locus length 0-100%")
for (i in 1:2070) {
  segments(a$start_norm[i],i,a$end_norm[i],i,col = as.character(a$col[i]))
}

############################################################
dev.off()


# TE content vs distance to the centromere 
# boxplot CG new 1001G distant to centromere close to centromere
###########################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPLEM.TEcontent.distant.nondist.pdf",height = 3,width = 5)
par(mar=c(6,6,2,2)) 
pc1<-PC_TE_coverage_anystrand$TAIR10[PC_TE_coverage_anystrand$gene %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere<1000000] ]
pc2<-PC_TE_coverage_anystrand$TAIR10[PC_TE_coverage_anystrand$gene %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere>1000000] ]

as1<-AS_TE_coverage_anystrand$TAIR10[AS_TE_coverage_anystrand$gene %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere<1000000]]
as2<-AS_TE_coverage_anystrand$TAIR10[AS_TE_coverage_anystrand$gene %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere>1000000]]

linc1<-lincRNAs_TE_coverage_anystrand$TAIR10[lincRNAs_TE_coverage_anystrand$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
linc2<-lincRNAs_TE_coverage_anystrand$TAIR10[lincRNAs_TE_coverage_anystrand$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

te1<-TEgene_TE_coverage_anystrand$TAIR10[TEgene_TE_coverage_anystrand$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere<1000000]]
te2<-TEgene_TE_coverage_anystrand$TAIR10[TEgene_TE_coverage_anystrand$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere>1000000]]

boxplot( pc1,as1,linc1,te1,pc2,as2,linc2,te2,
         col=c("#486EB4","#90C473","#F2AB54","#673A8E","#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE","PC","AS","linc","TE"),las=2, main="TE sequence content in TAIR10",ylab="portion of the gene coverd by TE pieces", ylim=c(0,1.1),outline = F, notch = T,rectCol="black", lineCol="black")

#################
#add p values   #
#################
len=min(length(pc1),length(as1),length(linc1),length(te1),length(pc2),length(as2),length(linc2),length(te2))
a<-wilcox.test(pc1,as1)
b<-wilcox.test(pc1,as1)
c<-wilcox.test(pc1,as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(linc1,as1)
b<-wilcox.test(linc1,as1)
c<-wilcox.test(linc1,as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(linc1,te1)
b<-wilcox.test(linc1,te1)
c<-wilcox.test(linc1,te1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
a<-wilcox.test(sample(pc2,2000),sample(as2,2000))
b<-wilcox.test(sample(pc2,2000),sample(as2,2000))
c<-wilcox.test(sample(pc2,2000),sample(as2,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.6)
a<-wilcox.test(sample(linc2,2000),sample(as2,2000))
b<-wilcox.test(sample(linc2,2000),sample(as2,2000))
c<-wilcox.test(sample(linc2,2000),sample(as2,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=0.6)
a<-wilcox.test(sample(linc2,2000),te2)
b<-wilcox.test(sample(linc2,2000),te2)
c<-wilcox.test(sample(linc2,2000),te2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.6)

a<-wilcox.test(pc1,sample(pc2,2000))
b<-wilcox.test(pc1,sample(pc2,2000))
c<-wilcox.test(pc1,sample(pc2,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=1.05)
a<-wilcox.test(sample(as2,2000),as1)
b<-wilcox.test(sample(as2,2000),as1)
c<-wilcox.test(sample(as2,2000),as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=1.05)
a<-wilcox.test(linc1,linc2)
b<-wilcox.test(linc1,linc2)
c<-wilcox.test(linc1,linc2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=1.05)
a<-wilcox.test(te1,te2)
b<-wilcox.test(te1,te2)
c<-wilcox.test(te1,te2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=1.05)
#################
dev.off()



#correlation between linc TE content and distance to centromere? 

a<-merge(lincRNAs_TE_coverage_anystrand[,c("gene","TAIR10")],lncRNAs.intergenic.loci)
plot(a$dist_from_centromere,a$TAIR10,pch=19,col=alpha("black",alpha=0.2),ylab=)
cor(a$dist_from_centromere,a$TAIR10)






col=c("#486EB4","#90C473","#F2AB54","#673A8E")
median (lincRNAs_TE_coverage.TAIR10$TE_length)#91







linc_noTE<-lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10==0] 
linc_50TE<-lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10>0 & lincRNAs_TE_coverage_anystrand$TAIR10<=0.5] 
linc_50_80TE<-lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10>0.5 & lincRNAs_TE_coverage_anystrand$TAIR10<=0.8] 
linc_80TE<-lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10>0.8] 

length(linc_noTE) #1070
length(linc_50TE)#703
length(linc_50_80TE)#187
length(linc_80TE)#286


########################################################
# do the TE pieces affect expression of lincRNAs? 
########################################################

#1001G
length(denovo2021.TPMs.genes.1001G$X6909[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% linc_noTE]) #62
length(denovo2021.TPMs.genes.1001G$X6909[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% linc_50TE])#39
length(denovo2021.TPMs.genes.1001G$X6909[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% linc_50_80TE])#8
length(denovo2021.TPMs.genes.1001G$X6909[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% linc_80TE])#3

#1001Gnew
length(denovo2021.TPMs.genes.1001Gnew$mean.6909[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% linc_noTE]) #43
length(denovo2021.TPMs.genes.1001Gnew$mean.6909[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% linc_50TE])#36
length(denovo2021.TPMs.genes.1001Gnew$mean.6909[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% linc_50_80TE])#8
length(denovo2021.TPMs.genes.1001Gnew$mean.6909[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% linc_80TE])#2



#pollen
length(denovo2021.TPMs.genes.ERACAPS$P.6909[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% linc_noTE]) #91
length(denovo2021.TPMs.genes.ERACAPS$P.6909[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% linc_50TE]) #101
length(denovo2021.TPMs.genes.ERACAPS$P.6909[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% linc_50_80TE]) #19
length(denovo2021.TPMs.genes.ERACAPS$P.6909[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% linc_80TE]) #27




########################################################
# do the TE pieces affect epigenetic state of lincRNAs? 
########################################################


#chipseq 

# chipseq coverage in lincRNA with different TE content 

#K9 boxplot lincs + TE genes
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_linc_vs_TEs_boxplot_K9_cov_lincs_TEcontent.pdf",height = 4,width = 3)
par(mar=c(6,4,4,2))
a1<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% linc_noTE]
a2<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% linc_50TE]
a3<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% linc_50_80TE]
a4<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% linc_80TE]

a5<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10, main="H3K9me2 level\n(Col-0 rosette)",cex.main=1.2,ylim=c(-2,2.2),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="log2(CHIP/INPUT)")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)

abline(h=median(chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(-2,2.2),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=2)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=2)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=2)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=2)
#################
dev.off()

#K27 boxplot lincs + TE genes
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_linc_vs_TEs_boxplot_K27_cov_lincs_TEcontent.pdf",height = 4,width = 3)
par(mar=c(6,4,4,2))
a1<-chip.denovo.log2$K27.6909[ chip.denovo.log2$gene %in% linc_noTE]
a2<-chip.denovo.log2$K27.6909[ chip.denovo.log2$gene %in% linc_50TE]
a3<-chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% linc_50_80TE]
a4<-chip.denovo.log2$K27.6909[ chip.denovo.log2$gene %in% linc_80TE]

a5<-chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10, main="H3K27me3 level\n(Col-0 rosette)",cex.main=1.2,ylim=c(-2,2),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="log2(CHIP/INPUT)")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)

abline(h=median(chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(-2,2),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1)
#################
dev.off()


#H1 boxplot lincs + TE genes
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_linc_vs_TEs_boxplot_H1_cov_lincs_TEcontent.pdf",height = 4,width = 3)
par(mar=c(6,4,4,2))
a1<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% linc_noTE]
a2<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% linc_50TE]
a3<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% linc_50_80TE]
a4<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% linc_80TE]

a5<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10, main="H1 level\n(Col-0 rosette)",cex.main=1.2,ylim=c(-1.2,1.3),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="log2(CHIP/INPUT)")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)

abline(h=median(chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(-1.2,1.3),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1)
#################
dev.off()


# methylation in lincRNA with different TE content 
#CG boxplot lincs + TE genes
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_linc_vs_TEs_boxplot_CGmeth_cov_lincs_TEcontent.pdf",height = 3,width = 3)
par(mar=c(6,4,4,2))
a1<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% linc_noTE]
a2<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% linc_50TE]
a3<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% linc_50_80TE]
a4<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% linc_80TE]

a5<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10, main="CG methylation level\n(Col-0 rosette)",cex.main=1.2,ylim=c(0,1),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="CG methylation")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)

abline(h=median(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,1),  
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1)
#################
dev.off()



#CHH boxplot lincs + TE genes
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_linc_vs_TEs_boxplot_CHHmeth_cov_lincs_TEcontent.pdf",height = 3,width = 3)
par(mar=c(6,4,4,2))
a1<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% linc_noTE]
a2<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% linc_50TE]
a3<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% linc_50_80TE]
a4<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% linc_80TE]

a5<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10, main="CHH methylation level\n(Col-0 rosette)",cex.main=1.2,ylim=c(0,0.15),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="CHH methylation")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)

abline(h=median(CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,0.15),  
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.1)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.1)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.1)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.1)
#################
dev.off()


# small RNA coverage in lincRNA with different TE content 
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_linc_vs_TEs_boxplot_sRNA_24nt_cov_lincs_TEcontent.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% linc_noTE]
a2<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% linc_50TE]
a3<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% linc_50_80TE]
a4<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% linc_80TE]

a5<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% TE_genes.loci$gene]

boxplot(        -1,-1,-1,-1,-1, ylim=c(0,2.7),main="24nt sRNA coverage\n(Col-0 flowers)",cex.main=1.2,col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="24nt small RNA coverage, RPM")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)

abline(h=median(sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
#text(x=6,y=median(sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% denovoPC.loci$gene])+0.07,adj=c(0.5,0.5),"PC genes",col="darkblue",cex=0.8)
abline(h=median(sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
#text(x=6,y=median(sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.antisense.loci$gene])+0.07,adj=c(0.5,0.5),"AS lncRNAs",col="darkgreen",cex=0.8)
abline(h=median(sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)
#text(x=6,y=median(sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene])+0.07,adj=c(0.5,0.5),"all lincRNAs",col="orange",cex=0.8)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,2.5),
col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
len=min(length(pc1),length(as1),length(linc1),length(te1),length(pc2),length(as2),length(linc2),length(te2))
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=2)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=2)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=2)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=2)
#################
dev.off()




########################################################
# do the TE pieces affect expression and epigenetic variability of lincRNAs? 
########################################################


#Expression variation 1001G lincs + TE genes
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_Expression_variation_1001G_lincs_TEcontent.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% linc_noTE & denovo2021.TPMs.genes.1001G$max>1 ]
a2<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% linc_50TE& denovo2021.TPMs.genes.1001G$max>1 ]
a3<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% linc_50_80TE& denovo2021.TPMs.genes.1001G$max>1 ]
a4<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% linc_80TE& denovo2021.TPMs.genes.1001G$max>1 ]

a5<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.loci$gene& denovo2021.TPMs.genes.1001G$max>1 ]

boxplot(        -10,-10,-10,-10,-10, main="Expression variation\n(461 accessions, rosette)",cex.main=1.2,ylim=c(0,22),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="coefficient of variance")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)

abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene],na.rm = T),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene],na.rm = T),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene],na.rm = T),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,22),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=20)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=20)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=20)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=20)
#################
dev.off()
wilcox.test(a1,a3)

#expression level control: 
#Expression variation 1001G lincs + TE genes : 1<TPMmax<2
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_Expression_variation_1001G_lincs_TEcontent_tpmmax_1_2.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% linc_noTE & denovo2021.TPMs.genes.1001G$max>1 & denovo2021.TPMs.genes.1001G$max<2]
a2<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% linc_50TE& denovo2021.TPMs.genes.1001G$max>1 & denovo2021.TPMs.genes.1001G$max<2]
a3<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% linc_50_80TE& denovo2021.TPMs.genes.1001G$max>1 & denovo2021.TPMs.genes.1001G$max<2]
a4<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% linc_80TE& denovo2021.TPMs.genes.1001G$max>1 & denovo2021.TPMs.genes.1001G$max<2]
a5<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.loci$gene& denovo2021.TPMs.genes.1001G$max>1 & denovo2021.TPMs.genes.1001G$max<2]

boxplot(        -10,-10,-10,-10,-10, main="Expression variation\n(461 accessions, rosette)",cex.main=1.2,ylim=c(0,22),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="coefficient of variance")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)
boxplot(      a1,a2,a3,a4,a5, ylim=c(0,22),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)
mtext(paste("N=",length(a1),sep=""),at=1,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a2),sep=""),at=2,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a3),sep=""),at=3,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a4),sep=""),at=4,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a5),sep=""),at=5,line=-1,side=1,cex=0.6)
mtext("1<TPMmax<2",at=0,line=3,side=1,cex=0.8)

#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=20)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=20)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=20)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=20)
#################
dev.off()
wilcox.test(a1,a3)
wilcox.test(a1,a4)


#Expression variation ERACAPS flowers lincs + TE genes : 1<TPMmax<2
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_Expression_variation_Eracaps_flowers_lincs_TEcontent_tpmmax1-2.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-denovo2021.TPMs.genes.ERACAPS$var.flowers[ denovo2021.TPMs.genes.ERACAPS$gene %in% linc_noTE & denovo2021.TPMs.genes.ERACAPS$max.flowers>1& denovo2021.TPMs.genes.ERACAPS$max.flowers<2 ]
a2<-denovo2021.TPMs.genes.ERACAPS$var.flowers[ denovo2021.TPMs.genes.ERACAPS$gene %in% linc_50TE& denovo2021.TPMs.genes.ERACAPS$max.flowers>1 & denovo2021.TPMs.genes.ERACAPS$max.flowers<2 ]
a3<-denovo2021.TPMs.genes.ERACAPS$var.flowers[denovo2021.TPMs.genes.ERACAPS$gene %in% linc_50_80TE& denovo2021.TPMs.genes.ERACAPS$max.flowers>1 & denovo2021.TPMs.genes.ERACAPS$max.flowers<2 ]
a4<-denovo2021.TPMs.genes.ERACAPS$var.flowers[ denovo2021.TPMs.genes.ERACAPS$gene %in% linc_80TE& denovo2021.TPMs.genes.ERACAPS$max.flowers>1 & denovo2021.TPMs.genes.ERACAPS$max.flowers<2 ]

a5<-denovo2021.TPMs.genes.ERACAPS$var.flowers[denovo2021.TPMs.genes.ERACAPS$gene %in% TE_genes.loci$gene& denovo2021.TPMs.genes.ERACAPS$max.flowers>1 & denovo2021.TPMs.genes.ERACAPS$max.flowers<2 ]

boxplot(        -10,-10,-10,-10,-10, main="Expression variation\n(23 accessions, flowers)",cex.main=1.2,ylim=c(0,5),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="coefficient of variance")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)
boxplot(      a1,a2,a3,a4,a5, ylim=c(0,5),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)
mtext(paste("N=",length(a1),sep=""),at=1,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a2),sep=""),at=2,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a3),sep=""),at=3,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a4),sep=""),at=4,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a5),sep=""),at=5,line=-1,side=1,cex=0.6)
mtext("1<TPMmax<2",at=0,line=3,side=1,cex=0.8)

#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=5)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=5)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=5)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=5)
#################
dev.off()
wilcox.test(a1,a3)
wilcox.test(a1,a4)



#Expression variation 1001GNEW lincs + TE genes : 1<TPMmax<2
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_Expression_variation_1001GNEW_lincs_TEcontent_tpmmax1-2.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-denovo2021.TPMs.genes.1001Gnew$variance_of_means[ denovo2021.TPMs.genes.1001Gnew$gene %in% linc_noTE &denovo2021.TPMs.genes.1001Gnew$ma_x>1&denovo2021.TPMs.genes.1001Gnew$ma_x<2]
a2<-denovo2021.TPMs.genes.1001Gnew$variance_of_means[ denovo2021.TPMs.genes.1001Gnew$gene %in% linc_50TE&denovo2021.TPMs.genes.1001Gnew$ma_x>1&denovo2021.TPMs.genes.1001Gnew$ma_x<2]
a3<-denovo2021.TPMs.genes.1001Gnew$variance_of_means[denovo2021.TPMs.genes.1001Gnew$gene %in% linc_50_80TE&denovo2021.TPMs.genes.1001Gnew$ma_x>1&denovo2021.TPMs.genes.1001Gnew$ma_x<2]
a4<-denovo2021.TPMs.genes.1001Gnew$variance_of_means[ denovo2021.TPMs.genes.1001Gnew$gene %in% linc_80TE&denovo2021.TPMs.genes.1001Gnew$ma_x>1&denovo2021.TPMs.genes.1001Gnew$ma_x<2]

a5<-denovo2021.TPMs.genes.1001Gnew$variance_of_means[denovo2021.TPMs.genes.1001Gnew$gene %in% TE_genes.loci$gene&denovo2021.TPMs.genes.1001Gnew$ma_x>1&denovo2021.TPMs.genes.1001Gnew$ma_x<2]

boxplot(        -10,-10,-10,-10,-10, main="Expression variation\n(28 accessions, rosette)",cex.main=1.2,ylim=c(0,6),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="coefficient of variance")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)
boxplot(      a1,a2,a3,a4,a5, ylim=c(0,6),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

mtext(paste("N=",length(a1),sep=""),at=1,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a2),sep=""),at=2,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a3),sep=""),at=3,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a4),sep=""),at=4,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a5),sep=""),at=5,line=-1,side=1,cex=0.6)
mtext("1<TPMmax<2",at=0,line=3,side=1,cex=0.8)

#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=6)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=6)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=6)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=6)
#################
dev.off()
wilcox.test(a1,a3)
wilcox.test(a1,a4)


#Intravariance expression 1001GNEW lincs + TE genes : 1<TPMmax<2
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_Expression_intravariance_1001GNEW_lincs_TEcontent_tpmmax1-2.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-denovo2021.TPMs.genes.1001Gnew$mean_intravariance[ denovo2021.TPMs.genes.1001Gnew$gene %in% linc_noTE &denovo2021.TPMs.genes.1001Gnew$ma_x>1&denovo2021.TPMs.genes.1001Gnew$ma_x<2]
a2<-denovo2021.TPMs.genes.1001Gnew$mean_intravariance[ denovo2021.TPMs.genes.1001Gnew$gene %in% linc_50TE&denovo2021.TPMs.genes.1001Gnew$ma_x>1&denovo2021.TPMs.genes.1001Gnew$ma_x<2]
a3<-denovo2021.TPMs.genes.1001Gnew$mean_intravariance[denovo2021.TPMs.genes.1001Gnew$gene %in% linc_50_80TE&denovo2021.TPMs.genes.1001Gnew$ma_x>1&denovo2021.TPMs.genes.1001Gnew$ma_x<2]
a4<-denovo2021.TPMs.genes.1001Gnew$mean_intravariance[ denovo2021.TPMs.genes.1001Gnew$gene %in% linc_80TE&denovo2021.TPMs.genes.1001Gnew$ma_x>1&denovo2021.TPMs.genes.1001Gnew$ma_x<2]

a5<-denovo2021.TPMs.genes.1001Gnew$mean_intravariance[denovo2021.TPMs.genes.1001Gnew$gene %in% TE_genes.loci$gene&denovo2021.TPMs.genes.1001Gnew$ma_x>1&denovo2021.TPMs.genes.1001Gnew$ma_x<2]

boxplot(        -10,-10,-10,-10,-10, main="Expression variation\n(2-4 replicates, rosette)",cex.main=1.2,ylim=c(0,6),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="coefficient of variance")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)

abline(h=median(denovo2021.TPMs.genes.1001Gnew$mean_intravariance[denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene&denovo2021.TPMs.genes.1001Gnew$ma_x>1],na.rm = T),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(denovo2021.TPMs.genes.1001Gnew$mean_intravariance[denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene&denovo2021.TPMs.genes.1001Gnew$ma_x>1],na.rm = T),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(denovo2021.TPMs.genes.1001Gnew$mean_intravariance[denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene&denovo2021.TPMs.genes.1001Gnew$ma_x>1],na.rm = T),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,6),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

mtext(paste("N=",length(a1),sep=""),at=1,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a2),sep=""),at=2,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a3),sep=""),at=3,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a4),sep=""),at=4,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a5),sep=""),at=5,line=-1,side=1,cex=0.6)
mtext("1<TPMmax<2",at=0,line=3,side=1,cex=0.8)
#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=6)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=6)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=6)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=6)
#################
dev.off()


#Noise lincs + TE genes  : 1<TPMmax<2
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_Noiise_lincs_TEcontent_tpmmax_1-2.pdf",height = 4,width = 3)
###########################################################################

par(mar=c(6,4,4,2))
a1<-denovo2021.TPMs.genes.Cortijo$noise_average[ denovo2021.TPMs.genes.Cortijo$gene %in% linc_noTE &denovo2021.TPMs.genes.Cortijo$max>1&denovo2021.TPMs.genes.Cortijo$max<2]
a2<-denovo2021.TPMs.genes.Cortijo$noise_average[ denovo2021.TPMs.genes.Cortijo$gene %in% linc_50TE &denovo2021.TPMs.genes.Cortijo$max>1&denovo2021.TPMs.genes.Cortijo$max<2]
a3<-denovo2021.TPMs.genes.Cortijo$noise_average[denovo2021.TPMs.genes.Cortijo$gene %in% linc_50_80TE &denovo2021.TPMs.genes.Cortijo$max>1&denovo2021.TPMs.genes.Cortijo$max<2]
a4<-denovo2021.TPMs.genes.Cortijo$noise_average[ denovo2021.TPMs.genes.Cortijo$gene %in% linc_80TE &denovo2021.TPMs.genes.Cortijo$max>1&denovo2021.TPMs.genes.Cortijo$max<2]

a5<-denovo2021.TPMs.genes.Cortijo$noise_average[denovo2021.TPMs.genes.Cortijo$gene %in% TE_genes.loci$gene &denovo2021.TPMs.genes.Cortijo$max>1&denovo2021.TPMs.genes.Cortijo$max<2]

boxplot(        -10,-10,-10,-10,-10, main="Expression noise\n(Col-0, seedling)",cex.main=1.2,ylim=c(0,5),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="coefficient of variance")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,5),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

mtext(paste("N=",length(a1),sep=""),at=1,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a2),sep=""),at=2,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a3),sep=""),at=3,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a4),sep=""),at=4,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a5),sep=""),at=5,line=-1,side=1,cex=0.6)
mtext("1<TPMmax<2",at=0,line=3,side=1,cex=0.8)
#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=5)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=5)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=5)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=5)
#################
dev.off()

#Noise lincs + TE genes  : TPMmax>1
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_Noiise_lincs_TEcontent.pdf",height = 4,width = 3)
###########################################################################

par(mar=c(6,4,4,2))
a1<-denovo2021.TPMs.genes.Cortijo$noise_average[ denovo2021.TPMs.genes.Cortijo$gene %in% linc_noTE &denovo2021.TPMs.genes.Cortijo$max>1]
a2<-denovo2021.TPMs.genes.Cortijo$noise_average[ denovo2021.TPMs.genes.Cortijo$gene %in% linc_50TE &denovo2021.TPMs.genes.Cortijo$max>1]
a3<-denovo2021.TPMs.genes.Cortijo$noise_average[denovo2021.TPMs.genes.Cortijo$gene %in% linc_50_80TE &denovo2021.TPMs.genes.Cortijo$max>1]
a4<-denovo2021.TPMs.genes.Cortijo$noise_average[ denovo2021.TPMs.genes.Cortijo$gene %in% linc_80TE &denovo2021.TPMs.genes.Cortijo$max>1]

a5<-denovo2021.TPMs.genes.Cortijo$noise_average[denovo2021.TPMs.genes.Cortijo$gene %in% TE_genes.loci$gene &denovo2021.TPMs.genes.Cortijo$max>1]

boxplot(        -10,-10,-10,-10,-10, main="Expression noise\n(Col-0, seedling)",cex.main=1.2,ylim=c(0,5),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="coefficient of variance")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,5),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

mtext(paste("N=",length(a1),sep=""),at=1,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a2),sep=""),at=2,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a3),sep=""),at=3,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a4),sep=""),at=4,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a5),sep=""),at=5,line=-1,side=1,cex=0.6)
mtext("1<TPMmax",at=0,line=3,side=1,cex=0.8)
#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=5)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=5)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=5)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=5)
#################
dev.off()


#Circadian variance lincs + TE genes : 1<TPMmax<2
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_circadian_variance_lincs_TEcontent_tpmmax1-2.pdf",height = 4,width = 3)
###########################################################################

par(mar=c(6,4,4,2))
a1<-denovo2021.TPMs.genes.Cortijo$variance_circadian[ denovo2021.TPMs.genes.Cortijo$gene %in% linc_noTE &denovo2021.TPMs.genes.Cortijo$max>1&denovo2021.TPMs.genes.Cortijo$max<2]
a2<-denovo2021.TPMs.genes.Cortijo$variance_circadian[ denovo2021.TPMs.genes.Cortijo$gene %in% linc_50TE &denovo2021.TPMs.genes.Cortijo$max>1&denovo2021.TPMs.genes.Cortijo$max<2]
a3<-denovo2021.TPMs.genes.Cortijo$variance_circadian[denovo2021.TPMs.genes.Cortijo$gene %in% linc_50_80TE &denovo2021.TPMs.genes.Cortijo$max>1&denovo2021.TPMs.genes.Cortijo$max<2]
a4<-denovo2021.TPMs.genes.Cortijo$variance_circadian[ denovo2021.TPMs.genes.Cortijo$gene %in% linc_80TE &denovo2021.TPMs.genes.Cortijo$max>1&denovo2021.TPMs.genes.Cortijo$max<2]

a5<-denovo2021.TPMs.genes.Cortijo$variance_circadian[denovo2021.TPMs.genes.Cortijo$gene %in% TE_genes.loci$gene &denovo2021.TPMs.genes.Cortijo$max>1&denovo2021.TPMs.genes.Cortijo$max<2]

boxplot(        -10,-10,-10,-10,-10, main="Circadian expression variation\n(Col-0, seedling)",cex.main=1.2,ylim=c(0,2.2),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="coefficient of variance")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,2.2),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

mtext(paste("N=",length(a1),sep=""),at=1,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a2),sep=""),at=2,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a3),sep=""),at=3,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a4),sep=""),at=4,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(a5),sep=""),at=5,line=-1,side=1,cex=0.6)
mtext("1<TPMmax<2",at=0,line=3,side=1,cex=0.8)
#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=2)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=2)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=2)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=2)
#################
dev.off()






###########################################################################
#Epigenetic variation  ChiPseq
###########################################################################

#K9
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_K9_variation_lincs_TEcontent.pdf",height = 4,width = 3)
par(mar=c(6,4,4,2))
###########################################################################
a1<-chip.denovo.quantstan$sd.key9[ chip.denovo.quantstan$gene %in% linc_noTE]
a2<-chip.denovo.quantstan$sd.key9[ chip.denovo.quantstan$gene %in% linc_50TE]
a3<-chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% linc_50_80TE]
a4<-chip.denovo.quantstan$sd.key9[ chip.denovo.quantstan$gene %in% linc_80TE]

a5<-chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10, main="H3K9me2 variation\n(13 accessions, rosette)",cex.main=1.2,ylim=c(0,1.2),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="standard deviation\n of normalized ChIP signal")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)

abline(h=median(chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% denovoPC.loci$gene],na.rm = T),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene],na.rm = T),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene],na.rm = T),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,1.2),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1)
#################
dev.off()

#K27
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_K27_variation_lincs_TEcontent.pdf",height = 4,width = 3)
par(mar=c(6,4,4,2))
###########################################################################
a1<-chip.denovo.quantstan$sd.key27[ chip.denovo.quantstan$gene %in% linc_noTE]
a2<-chip.denovo.quantstan$sd.key27[ chip.denovo.quantstan$gene %in% linc_50TE]
a3<-chip.denovo.quantstan$sd.key27[chip.denovo.quantstan$gene %in% linc_50_80TE]
a4<-chip.denovo.quantstan$sd.key27[ chip.denovo.quantstan$gene %in% linc_80TE]

a5<-chip.denovo.quantstan$sd.key27[chip.denovo.quantstan$gene %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10, main="H3K27me3 variation\n(12 accessions, rosette)",cex.main=1.2,ylim=c(0,1.2),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="standard deviation\n of normalized ChIP signal")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)

abline(h=median(chip.denovo.quantstan$sd.key27[chip.denovo.quantstan$gene %in% denovoPC.loci$gene],na.rm = T),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.quantstan$sd.key27[chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene],na.rm = T),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.quantstan$sd.key27[chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene],na.rm = T),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,1.2),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1)
#################
dev.off()


#H1
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_H1_variation_lincs_TEcontent.pdf",height = 4,width = 3)
par(mar=c(6,4,4,2))
###########################################################################
a1<-chip.denovo.quantstan$sd.hist1[ chip.denovo.quantstan$gene %in% linc_noTE]
a2<-chip.denovo.quantstan$sd.hist1[ chip.denovo.quantstan$gene %in% linc_50TE]
a3<-chip.denovo.quantstan$sd.hist1[chip.denovo.quantstan$gene %in% linc_50_80TE]
a4<-chip.denovo.quantstan$sd.hist1[ chip.denovo.quantstan$gene %in% linc_80TE]

a5<-chip.denovo.quantstan$sd.hist1[chip.denovo.quantstan$gene %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10, main="H1 variation\n(14 accessions, rosette)",cex.main=1.2,ylim=c(0,1.2),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="standard deviation\n of normalized ChIP signal")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)

abline(h=median(chip.denovo.quantstan$sd.hist1[chip.denovo.quantstan$gene %in% denovoPC.loci$gene],na.rm = T),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.quantstan$sd.hist1[chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene],na.rm = T),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.quantstan$sd.hist1[chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene],na.rm = T),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,1.2),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1)
#################
dev.off()




###########################################################################
#epigenetic variability methylation 
###########################################################################

#CG 1001G
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_CG_variation_1001G_vs_lincs_TEcontent.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,5,4,2))

a1<-CG.1001.denovo$sd[CG.1001.denovo$transcript %in% linc_noTE]
a2<-CG.1001.denovo$sd[CG.1001.denovo$transcript %in% linc_50TE]
a3<-CG.1001.denovo$sd[CG.1001.denovo$transcript %in% linc_50_80TE]
a4<-CG.1001.denovo$sd[CG.1001.denovo$transcript %in% linc_80TE]
a5<-CG.1001.denovo$sd[CG.1001.denovo$transcript %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10,main="CG methylation variation\n444 accessions",ylab="standard deviation of CG methylation level",cex.main=1.2,ylim=c(0,0.5),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F)

#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,0.5),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E"), las=2, notch = T, outline = F, add=T)
#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(a4,a3)
b<-wilcox.test(a4,a3)
c<-wilcox.test(a4,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.4)
#################
dev.off()

#CHH 1001G
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_CHH_variation_1001G_vs_lincs_TEcontent.pdf",height =4,width = 3)
###########################################################################
par(mar=c(6,5,4,2))

a1<-CHH.1001.denovo$sd[CHH.1001.denovo$transcript %in% linc_noTE]
a2<-CHH.1001.denovo$sd[CHH.1001.denovo$transcript %in% linc_50TE]
a3<-CHH.1001.denovo$sd[CHH.1001.denovo$transcript %in% linc_50_80TE]
a4<-CHH.1001.denovo$sd[CHH.1001.denovo$transcript %in% linc_80TE]
a5<-CHH.1001.denovo$sd[CHH.1001.denovo$transcript %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10,main="CHH methylation variation\n 450 accessions",ylab="standard deviation of CHH methylation level",cex.main=1.2,ylim=c(0,0.5),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F)

#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,0.5),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E"), las=2, notch = T, outline = F, add=T)
#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(a4,a3)
b<-wilcox.test(a4,a3)
c<-wilcox.test(a4,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.4)
#################
dev.off()


#CG 1001Gnew
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_CG_variation_1001GNEW_vs_lincs_TEcontent.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,5,4,2))

a1<-CG.1001new.denovo$sd_of_means[CG.1001.denovo$transcript %in% linc_noTE]
a2<-CG.1001new.denovo$sd_of_means[CG.1001.denovo$transcript %in% linc_50TE]
a3<-CG.1001new.denovo$sd_of_means[CG.1001.denovo$transcript %in% linc_50_80TE]
a4<-CG.1001new.denovo$sd_of_means[CG.1001.denovo$transcript %in% linc_80TE]
a5<-CG.1001new.denovo$sd_of_means[CG.1001.denovo$transcript %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10,main="CG methylation variation\n28 accessions",ylab="standard deviation of CG methylation level",cex.main=1.2,ylim=c(0,0.5),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F)

#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,0.5),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E"), las=2, notch = T, outline = F, add=T)
#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(a4,a3)
b<-wilcox.test(a4,a3)
c<-wilcox.test(a4,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.4)
#################
dev.off()

#CHH 1001Gnew
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_CHH_variation_1001GNEW_vs_lincs_TEcontent.pdf",height =4,width = 3)
###########################################################################
par(mar=c(6,5,4,2))

a1<-CHH.1001new.denovo$sd_of_means[CHH.1001new.denovo$transcript %in% linc_noTE]
a2<-CHH.1001new.denovo$sd_of_means[CHH.1001new.denovo$transcript %in% linc_50TE]
a3<-CHH.1001new.denovo$sd_of_means[CHH.1001new.denovo$transcript %in% linc_50_80TE]
a4<-CHH.1001new.denovo$sd_of_means[CHH.1001new.denovo$transcript %in% linc_80TE]
a5<-CHH.1001new.denovo$sd_of_means[CHH.1001new.denovo$transcript %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10,main="CHH methylation variation\n28 accessions",ylab="standard deviation of CHH methylation level",cex.main=1.2,ylim=c(0,0.5),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F)

#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,0.5),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E"), las=2, notch = T, outline = F, add=T)
#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(a4,a3)
b<-wilcox.test(a4,a3)
c<-wilcox.test(a4,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.4)
#################
dev.off()


#CG intravariation 1001Gnew
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_CG_intravariation_1001Gnew_lincs_TEcontent.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,5,4,2))

a1<-CG.1001new.denovo$mean_intra_sd[CG.1001new.denovo$transcript %in% linc_noTE]
a2<-CG.1001new.denovo$mean_intra_sd[CG.1001new.denovo$transcript %in% linc_50TE]
a3<-CG.1001new.denovo$mean_intra_sd[CG.1001new.denovo$transcript %in% linc_50_80TE]
a4<-CG.1001new.denovo$mean_intra_sd[CG.1001new.denovo$transcript %in% linc_80TE]
a5<-CG.1001new.denovo$mean_intra_sd[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10,main="CG methylation intravariation\n2-4 reps",ylab="standard deviation of CG methylation level",cex.main=1.2,ylim=c(0,0.25),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F)

#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,0.25),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E"), las=2, notch = T, outline = F, add=T)
#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.2)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.2)
a<-wilcox.test(a4,a3)
b<-wilcox.test(a4,a3)
c<-wilcox.test(a4,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.2)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.2)
#################
dev.off()

#CHH - intravariation 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_CHH_intravariation_1001Gnew_lincs_TEcontent.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,5,4,2))

a1<-CHH.1001new.denovo$mean_intra_sd[CHH.1001new.denovo$transcript %in% linc_noTE]
a2<-CHH.1001new.denovo$mean_intra_sd[CHH.1001new.denovo$transcript %in% linc_50TE]
a3<-CHH.1001new.denovo$mean_intra_sd[CHH.1001new.denovo$transcript %in% linc_50_80TE]
a4<-CHH.1001new.denovo$mean_intra_sd[CHH.1001new.denovo$transcript %in% linc_80TE]
a5<-CHH.1001new.denovo$mean_intra_sd[CHH.1001new.denovo$transcript %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10,main="CHH methylation intravariation\n2-4 reps",ylab="standard deviation of CHH methylation level",cex.main=1.2,ylim=c(0,0.4),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F)

#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a1,a2,a3,a4,a5, ylim=c(0,0.4),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E"), las=2, notch = T, outline = F, add=T)
#################
#add p values   #
#################
a<-wilcox.test(a1,a2)
b<-wilcox.test(a1,a2)
c<-wilcox.test(a1,a2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.2)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.2)
a<-wilcox.test(a4,a3)
b<-wilcox.test(a4,a3)
c<-wilcox.test(a4,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.2)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.2)
#################
dev.off()




###########################################################################
#### Control - close and far from centromeres 
###########################################################################


#K9 boxplot lincs + TE genes :  far and close to centromere 
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_linc_vs_TEs_boxplot_K9_vs_lincs_TEcontent.distant_nondist.pdf",height = 4,width = 5)
par(mar=c(6,5,4,2))
a11<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% linc_noTE & chip.denovo.log2$gene %in%lncRNAs.intergenic.loci$gene[ lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a12<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% linc_noTE & chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a21<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% linc_50TE& chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a22<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% linc_50TE& chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a31<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% linc_50_80TE& chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a32<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% linc_50_80TE& chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a41<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% linc_80TE & chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a42<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% linc_80TE & chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a51<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere<1000000]]
a52<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere>1000000]]

boxplot(        -10,-10,-10,-10,-10,-10,-10,-10,-10,-10, main="H3K9me2 level\n(Col-0 rosette)",cex.main=1.2,ylim=c(-2,2.5),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes","noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="log2(CHIP/INPUT)")
mtext (text = "distance to\n centromere",side = 3,at = -1,line = 0)

mtext (text = "close, <1Mb",side = 3,at = 3,line = 0)
mtext (text = "far, >1Mb",side = 3,at = 8,line = 0)

abline(h=median(chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a11,a21,a31,a41,a51,a12,a22,a32,a42,a52, ylim=c(-2,2.5),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a11,a21)
b<-wilcox.test(a11,a21)
c<-wilcox.test(a11,a21)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1.5)
a<-wilcox.test(a21,a31)
b<-wilcox.test(a21,a31)
c<-wilcox.test(a21,a31)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1.5)
a<-wilcox.test(a31,a41)
b<-wilcox.test(a31,a41)
c<-wilcox.test(a31,a41)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1.5)
a<-wilcox.test(a41,a51)
b<-wilcox.test(a41,a51)
c<-wilcox.test(a41,a51)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1.5)

a<-wilcox.test(a12,a22)
b<-wilcox.test(a12,a22)
c<-wilcox.test(a12,a22)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=1)
a<-wilcox.test(a22,a32)
b<-wilcox.test(a22,a32)
c<-wilcox.test(a22,a32)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
a<-wilcox.test(a32,a42)
b<-wilcox.test(a32,a42)
c<-wilcox.test(a32,a42)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=1)
a<-wilcox.test(a42,a52)
b<-wilcox.test(a42,a52)
c<-wilcox.test(a42,a52)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=9.5,y=1.5)

a<-wilcox.test(a11,a12)
b<-wilcox.test(a11,a12)
c<-wilcox.test(a11,a12)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=2.5)
a<-wilcox.test(a21,a22)
b<-wilcox.test(a21,a22)
c<-wilcox.test(a21,a22)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=2.5)
a<-wilcox.test(a31,a32)
b<-wilcox.test(a31,a32)
c<-wilcox.test(a31,a32)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=2.5)
a<-wilcox.test(a41,a42)
b<-wilcox.test(a41,a42)
c<-wilcox.test(a41,a42)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=2.5)
a<-wilcox.test(a51,a52)
b<-wilcox.test(a51,a52)
c<-wilcox.test(a51,a52)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5,y=2.5)
#################
dev.off()

#H1 boxplot lincs + TE genes : far and close to centromere 
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_linc_vs_TEs_boxplot_H1_vs_lincs_TEcontent.distant_nondist.pdf",height = 4,width = 5)
par(mar=c(6,5,4,2))
a11<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% linc_noTE & chip.denovo.log2$gene %in%lncRNAs.intergenic.loci$gene[ lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a12<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% linc_noTE & chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a21<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% linc_50TE& chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a22<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% linc_50TE& chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a31<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% linc_50_80TE& chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a32<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% linc_50_80TE& chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a41<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% linc_80TE & chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a42<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% linc_80TE & chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a51<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere<1000000]]
a52<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere>1000000]]

boxplot(        -10,-10,-10,-10,-10,-10,-10,-10,-10,-10, main="H1 level\n(Col-0 rosette)",cex.main=1.2,ylim=c(-1,1.5),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes","noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="log2(CHIP/INPUT)")
mtext (text = "distance to\n centromere",side = 3,at = -1,line = 0)

mtext (text = "close, <1Mb",side = 3,at = 3,line = 0)
mtext (text = "far, >1Mb",side = 3,at = 8,line = 0)

abline(h=median(chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a11,a21,a31,a41,a51,a12,a22,a32,a42,a52, ylim=c(-1,1.5),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a11,a21)
b<-wilcox.test(a11,a21)
c<-wilcox.test(a11,a21)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
a<-wilcox.test(a21,a31)
b<-wilcox.test(a21,a31)
c<-wilcox.test(a21,a31)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
a<-wilcox.test(a31,a41)
b<-wilcox.test(a31,a41)
c<-wilcox.test(a31,a41)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(a41,a51)
b<-wilcox.test(a41,a51)
c<-wilcox.test(a41,a51)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1)

a<-wilcox.test(a12,a22)
b<-wilcox.test(a12,a22)
c<-wilcox.test(a12,a22)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=1)
a<-wilcox.test(a22,a32)
b<-wilcox.test(a22,a32)
c<-wilcox.test(a22,a32)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
a<-wilcox.test(a32,a42)
b<-wilcox.test(a32,a42)
c<-wilcox.test(a32,a42)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=1)
a<-wilcox.test(a42,a52)
b<-wilcox.test(a42,a52)
c<-wilcox.test(a42,a52)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=9.5,y=1)

a<-wilcox.test(a11,a12)
b<-wilcox.test(a11,a12)
c<-wilcox.test(a11,a12)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=1.5)
a<-wilcox.test(a21,a22)
b<-wilcox.test(a21,a22)
c<-wilcox.test(a21,a22)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=1.5)
a<-wilcox.test(a31,a32)
b<-wilcox.test(a31,a32)
c<-wilcox.test(a31,a32)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=1.5)
a<-wilcox.test(a41,a42)
b<-wilcox.test(a41,a42)
c<-wilcox.test(a41,a42)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=1.5)
a<-wilcox.test(a51,a52)
b<-wilcox.test(a51,a52)
c<-wilcox.test(a51,a52)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5,y=1.5)
#################
dev.off()



# boxplot 24nt small RNA coverage lincs + TE genes: far and close to centromere 
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_24ntRPM_vs_lincs_TEcontent.distant_nondist.pdf",height = 4,width = 5)
###########################################################################
par(mar=c(6,5,4,2))
a11<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% linc_noTE & sRNA.24nt.denovo2021.RPM$gene %in%lncRNAs.intergenic.loci$gene[ lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a12<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% linc_noTE & sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a21<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% linc_50TE& sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a22<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% linc_50TE& sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a31<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% linc_50_80TE& sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a32<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% linc_50_80TE& sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a41<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% linc_80TE & sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a42<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% linc_80TE & sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a51<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere<1000000]]
a52<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere>1000000]]

boxplot(        -10,-10,-10,-10,-10,-10,-10,-10,-10,-10, main="24nt sRNA coverage\n(Col-0 flowers)",cex.main=1.2,ylim=c(0,3.5),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes","noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="24nt coverage, RPM")
mtext (text = "distance to\n centromere",side = 3,at = -1,line = 0)

mtext (text = "close, <1Mb",side = 3,at = 3,line = 0)
mtext (text = "far, >1Mb",side = 3,at = 8,line = 0)

abline(h=median(sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a11,a21,a31,a41,a51,a12,a22,a32,a42,a52, ylim=c(0,3.5),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a11,a21)
b<-wilcox.test(a11,a21)
c<-wilcox.test(a11,a21)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
a<-wilcox.test(a21,a31)
b<-wilcox.test(a21,a31)
c<-wilcox.test(a21,a31)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
a<-wilcox.test(a31,a41)
b<-wilcox.test(a31,a41)
c<-wilcox.test(a31,a41)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(a41,a51)
b<-wilcox.test(a41,a51)
c<-wilcox.test(a41,a51)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1)

a<-wilcox.test(a12,a22)
b<-wilcox.test(a12,a22)
c<-wilcox.test(a12,a22)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=2)
a<-wilcox.test(a22,a32)
b<-wilcox.test(a22,a32)
c<-wilcox.test(a22,a32)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=2)
a<-wilcox.test(a32,a42)
b<-wilcox.test(a32,a42)
c<-wilcox.test(a32,a42)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=2)
a<-wilcox.test(a42,a52)
b<-wilcox.test(a42,a52)
c<-wilcox.test(a42,a52)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=9.5,y=2)

a<-wilcox.test(a11,a12)
b<-wilcox.test(a11,a12)
c<-wilcox.test(a11,a12)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=3.5)
a<-wilcox.test(a21,a22)
b<-wilcox.test(a21,a22)
c<-wilcox.test(a21,a22)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=3.5)
a<-wilcox.test(a31,a32)
b<-wilcox.test(a31,a32)
c<-wilcox.test(a31,a32)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=3.5)
a<-wilcox.test(a41,a42)
b<-wilcox.test(a41,a42)
c<-wilcox.test(a41,a42)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=3.5)
a<-wilcox.test(a51,a52)
b<-wilcox.test(a51,a52)
c<-wilcox.test(a51,a52)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5,y=3.5)
#################
dev.off()


linc_centromere<-lncRNAs.intergenic.loci$gene[ lncRNAs.intergenic.loci$dist_from_centromere<1000000]
linc_arm<-lncRNAs.intergenic.loci$gene[ lncRNAs.intergenic.loci$dist_from_centromere>1000000]


# boxplot CHH lincs + TE genes: far and close to centromere 
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_CHH_vs_lincs_TEcontent.distant_nondist.pdf",height = 4,width = 5)
par(mar=c(6,5,4,2))
a11<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% linc_noTE & CHH.1001.denovo$transcript %in%lncRNAs.intergenic.loci$gene[ lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a12<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% linc_noTE & CHH.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a21<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% linc_50TE& CHH.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a22<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% linc_50TE& CHH.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a31<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% linc_50_80TE& CHH.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a32<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% linc_50_80TE& CHH.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a41<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% linc_80TE & CHH.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a42<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% linc_80TE & CHH.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a51<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere<1000000]]
a52<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere>1000000]]

boxplot(        -10,-10,-10,-10,-10,-10,-10,-10,-10,-10, main="CHH methylation level\n(Col-0 rosette)",cex.main=1.2,ylim=c(0,0.16),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes","noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="CHH methylation")
mtext (text = "distance to\n centromere",side = 3,at = -1,line = 0)

mtext (text = "close, <1Mb",side = 3,at = 3,line = 0)
mtext (text = "far, >1Mb",side = 3,at = 8,line = 0)

abline(h=median(CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a11,a21,a31,a41,a51,a12,a22,a32,a42,a52, ylim=c(0,0.16),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a11,a21)
b<-wilcox.test(a11,a21)
c<-wilcox.test(a11,a21)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.1)
a<-wilcox.test(a21,a31)
b<-wilcox.test(a21,a31)
c<-wilcox.test(a21,a31)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.1)
a<-wilcox.test(a31,a41)
b<-wilcox.test(a31,a41)
c<-wilcox.test(a31,a41)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.1)
a<-wilcox.test(a41,a51)
b<-wilcox.test(a41,a51)
c<-wilcox.test(a41,a51)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.1)

a<-wilcox.test(a12,a22)
b<-wilcox.test(a12,a22)
c<-wilcox.test(a12,a22)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=0.1)
a<-wilcox.test(a22,a32)
b<-wilcox.test(a22,a32)
c<-wilcox.test(a22,a32)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.1)
a<-wilcox.test(a32,a42)
b<-wilcox.test(a32,a42)
c<-wilcox.test(a32,a42)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=0.1)
a<-wilcox.test(a42,a52)
b<-wilcox.test(a42,a52)
c<-wilcox.test(a42,a52)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=9.5,y=0.1)

a<-wilcox.test(a11,a12)
b<-wilcox.test(a11,a12)
c<-wilcox.test(a11,a12)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=0.16)
a<-wilcox.test(a21,a22)
b<-wilcox.test(a21,a22)
c<-wilcox.test(a21,a22)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=0.16)
a<-wilcox.test(a31,a32)
b<-wilcox.test(a31,a32)
c<-wilcox.test(a31,a32)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=0.16)
a<-wilcox.test(a41,a42)
b<-wilcox.test(a41,a42)
c<-wilcox.test(a41,a42)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=0.16)
a<-wilcox.test(a51,a52)
b<-wilcox.test(a51,a52)
c<-wilcox.test(a51,a52)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5,y=0.16)
#################
dev.off()

# boxplot CG lincs + TE genes: far and close to centromere 
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_CG_vs_lincs_TEcontent.distant_nondist.pdf",height = 4,width = 5)
par(mar=c(6,5,4,2))
a11<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% linc_noTE & CG.1001.denovo$transcript %in%lncRNAs.intergenic.loci$gene[ lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a12<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% linc_noTE & CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a21<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% linc_50TE& CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a22<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% linc_50TE& CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a31<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% linc_50_80TE& CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a32<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% linc_50_80TE& CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a41<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% linc_80TE & CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
a42<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% linc_80TE & CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

a51<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere<1000000]]
a52<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere>1000000]]

boxplot(        -10,-10,-10,-10,-10,-10,-10,-10,-10,-10, main="CG methylation level\n(Col-0 rosette)",cex.main=1.2,ylim=c(0,1),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes","noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="CG methylation")
mtext (text = "distance to\n centromere",side = 3,at = -1,line = 0)

mtext (text = "close, <1Mb",side = 3,at = 3,line = 0)
mtext (text = "far, >1Mb",side = 3,at = 8,line = 0)

abline(h=median(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a11,a21,a31,a41,a51,a12,a22,a32,a42,a52, ylim=c(0,1),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a11,a21)
b<-wilcox.test(a11,a21)
c<-wilcox.test(a11,a21)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.9)
a<-wilcox.test(a21,a31)
b<-wilcox.test(a21,a31)
c<-wilcox.test(a21,a31)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.9)
a<-wilcox.test(a31,a41)
b<-wilcox.test(a31,a41)
c<-wilcox.test(a31,a41)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.9)
a<-wilcox.test(a41,a51)
b<-wilcox.test(a41,a51)
c<-wilcox.test(a41,a51)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.9)

a<-wilcox.test(a12,a22)
b<-wilcox.test(a12,a22)
c<-wilcox.test(a12,a22)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=0.9)
a<-wilcox.test(a22,a32)
b<-wilcox.test(a22,a32)
c<-wilcox.test(a22,a32)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.9)
a<-wilcox.test(a32,a42)
b<-wilcox.test(a32,a42)
c<-wilcox.test(a32,a42)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=0.9)
a<-wilcox.test(a42,a52)
b<-wilcox.test(a42,a52)
c<-wilcox.test(a42,a52)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=9.5,y=0.9)

a<-wilcox.test(a11,a12)
b<-wilcox.test(a11,a12)
c<-wilcox.test(a11,a12)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=1)
a<-wilcox.test(a21,a22)
b<-wilcox.test(a21,a22)
c<-wilcox.test(a21,a22)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=1)
a<-wilcox.test(a31,a32)
b<-wilcox.test(a31,a32)
c<-wilcox.test(a31,a32)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=1)
a<-wilcox.test(a41,a42)
b<-wilcox.test(a41,a42)
c<-wilcox.test(a41,a42)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=1)
a<-wilcox.test(a51,a52)
b<-wilcox.test(a51,a52)
c<-wilcox.test(a51,a52)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5,y=1)
#################
dev.off()


length(a11)#66

length(a21)#71

length(a31)#38

length(a41)#69

length(a51)#655

length(a12)#1004

length(a22)#632

length(a32)#149

length(a42)#217

length(a52)#1475




#EXpression variation  boxplot lincs + TE genes :  far and close to centromere 
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_ExpVariability_vs_lincs_TEcontent.distant_nondist.pdf",height = 4,width = 5)
###########################################################################
par(mar=c(6,5,4,2))
a11<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% linc_noTE & denovo2021.TPMs.genes.1001G$gene %in% linc_centromere]
a12<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% linc_noTE & denovo2021.TPMs.genes.1001G$gene %in% linc_arm]
a21<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% linc_50TE& denovo2021.TPMs.genes.1001G$gene %in% linc_centromere]
a22<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% linc_50TE& denovo2021.TPMs.genes.1001G$gene %in% linc_arm]
a31<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% linc_50_80TE& denovo2021.TPMs.genes.1001G$gene %in% linc_centromere]
a32<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% linc_50_80TE& denovo2021.TPMs.genes.1001G$gene %in% linc_arm]
a41<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% linc_80TE & denovo2021.TPMs.genes.1001G$gene %in% linc_centromere]
a42<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% linc_80TE & denovo2021.TPMs.genes.1001G$gene %in% linc_arm]
a51<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere<1000000]]
a52<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere>1000000]]

boxplot(        -10,-10,-10,-10,-10,-10,-10,-10,-10,-10, main="Expression variabiltiy\n(461 acc. rosette)",cex.main=1.2,ylim=c(0,22),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes","noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="coeficient of variance")
mtext (text = "distance to\n centromere",side = 3,at = -1,line = 0)

mtext (text = "close, <1Mb",side = 3,at = 3,line = 0)
mtext (text = "far, >1Mb",side = 3,at = 8,line = 0)

abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a11,a21,a31,a41,a51,a12,a22,a32,a42,a52, ylim=c(0,22),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a11,a21)
b<-wilcox.test(a11,a21)
c<-wilcox.test(a11,a21)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=20)
a<-wilcox.test(a21,a31)
b<-wilcox.test(a21,a31)
c<-wilcox.test(a21,a31)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=20)
a<-wilcox.test(a31,a41)
b<-wilcox.test(a31,a41)
c<-wilcox.test(a31,a41)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=20)
a<-wilcox.test(a41,a51)
b<-wilcox.test(a41,a51)
c<-wilcox.test(a41,a51)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=20)

a<-wilcox.test(a12,a22)
b<-wilcox.test(a12,a22)
c<-wilcox.test(a12,a22)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=20)
a<-wilcox.test(a22,a32)
b<-wilcox.test(a22,a32)
c<-wilcox.test(a22,a32)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=20)
a<-wilcox.test(a32,a42)
b<-wilcox.test(a32,a42)
c<-wilcox.test(a32,a42)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=20)
a<-wilcox.test(a42,a52)
b<-wilcox.test(a42,a52)
c<-wilcox.test(a42,a52)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=9.5,y=20)

a<-wilcox.test(a11,a12)
b<-wilcox.test(a11,a12)
c<-wilcox.test(a11,a12)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=22)
a<-wilcox.test(a21,a22)
b<-wilcox.test(a21,a22)
c<-wilcox.test(a21,a22)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=22)
a<-wilcox.test(a31,a32)
b<-wilcox.test(a31,a32)
c<-wilcox.test(a31,a32)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=22)
a<-wilcox.test(a41,a42)
b<-wilcox.test(a41,a42)
c<-wilcox.test(a41,a42)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=22)
a<-wilcox.test(a51,a52)
b<-wilcox.test(a51,a52)
c<-wilcox.test(a51,a52)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5,y=22)
#################
dev.off()



#inter-rep expression variation  boxplot lincs + TE genes 
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_INTRA_ExpVariability_vs_lincs_TEcontent.pdf",height = 4,width = 5)
###########################################################################
par(mar=c(6,5,4,2))
a11<-denovo2021.TPMs.genes.1001Gnew$mean_intravariance[denovo2021.TPMs.genes.1001Gnew$ma_x>1 &denovo2021.TPMs.genes.1001Gnew$gene %in% linc_noTE ]

a21<-denovo2021.TPMs.genes.1001Gnew$mean_intravariance[denovo2021.TPMs.genes.1001Gnew$ma_x>1 & denovo2021.TPMs.genes.1001Gnew$gene %in% linc_50TE]

a31<-denovo2021.TPMs.genes.1001Gnew$mean_intravariance[denovo2021.TPMs.genes.1001Gnew$ma_x>1 &denovo2021.TPMs.genes.1001Gnew$gene %in% linc_50_80TE]
a41<-denovo2021.TPMs.genes.1001Gnew$mean_intravariance[denovo2021.TPMs.genes.1001Gnew$ma_x>1 & denovo2021.TPMs.genes.1001Gnew$gene %in% linc_80TE ]

a51<-denovo2021.TPMs.genes.1001Gnew$mean_intravariance[denovo2021.TPMs.genes.1001Gnew$ma_x>1 &denovo2021.TPMs.genes.1001Gnew$gene %in% TE_genes.loci$gene]

boxplot(        -10,-10,-10,-10,-10,main="Inter-replicate expression variability",cex.main=1.2,ylim=c(0,2.3),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes"),las=2, notch = T, outline = F, ylab="coeficient of variance")

#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a11,a21,a31,a41,a51, ylim=c(0,2.3),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a11,a21)
b<-wilcox.test(a11,a21)
c<-wilcox.test(a11,a21)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=2.2)
a<-wilcox.test(a21,a31)
b<-wilcox.test(a21,a31)
c<-wilcox.test(a21,a31)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=2.2)
a<-wilcox.test(a31,a41)
b<-wilcox.test(a31,a41)
c<-wilcox.test(a31,a41)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=2.2)
a<-wilcox.test(a41,a51)
b<-wilcox.test(a41,a51)
c<-wilcox.test(a41,a51)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=2.2)

a<-wilcox.test(a12,a22)
b<-wilcox.test(a12,a22)
c<-wilcox.test(a12,a22)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=20)
a<-wilcox.test(a22,a32)
b<-wilcox.test(a22,a32)
c<-wilcox.test(a22,a32)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=20)
a<-wilcox.test(a32,a42)
b<-wilcox.test(a32,a42)
c<-wilcox.test(a32,a42)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=20)
a<-wilcox.test(a42,a52)
b<-wilcox.test(a42,a52)
c<-wilcox.test(a42,a52)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=9.5,y=20)

a<-wilcox.test(a11,a12)
b<-wilcox.test(a11,a12)
c<-wilcox.test(a11,a12)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=22)
a<-wilcox.test(a21,a22)
b<-wilcox.test(a21,a22)
c<-wilcox.test(a21,a22)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=22)
a<-wilcox.test(a31,a32)
b<-wilcox.test(a31,a32)
c<-wilcox.test(a31,a32)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=22)
a<-wilcox.test(a41,a42)
b<-wilcox.test(a41,a42)
c<-wilcox.test(a41,a42)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=22)
a<-wilcox.test(a51,a52)
b<-wilcox.test(a51,a52)
c<-wilcox.test(a51,a52)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5,y=22)
#################
dev.off()

#noise  boxplot lincs + TE genes 
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_Noise_vs_lincs_TEcontent.pdf",height = 4,width = 5)
par(mar=c(6,5,4,2))
a01<-denovo2021.TPMs.genes.Cortijo$noise_average[denovo2021.TPMs.genes.Cortijo$max>1 &denovo2021.TPMs.genes.Cortijo$gene %in% lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TE_presense_nrlines==0] ]
  
a11<-denovo2021.TPMs.genes.Cortijo$noise_average[denovo2021.TPMs.genes.Cortijo$max>1 &denovo2021.TPMs.genes.Cortijo$gene %in% linc_noTE ]

a21<-denovo2021.TPMs.genes.Cortijo$noise_average[denovo2021.TPMs.genes.Cortijo$max>1 & denovo2021.TPMs.genes.Cortijo$gene %in% linc_50TE]

a31<-denovo2021.TPMs.genes.Cortijo$noise_average[denovo2021.TPMs.genes.Cortijo$max>1 &denovo2021.TPMs.genes.Cortijo$gene %in% linc_50_80TE]
a41<-denovo2021.TPMs.genes.Cortijo$noise_average[denovo2021.TPMs.genes.Cortijo$max>1 & denovo2021.TPMs.genes.Cortijo$gene %in% linc_80TE ]

a51<-denovo2021.TPMs.genes.Cortijo$noise_average[denovo2021.TPMs.genes.Cortijo$max>1 &denovo2021.TPMs.genes.Cortijo$gene %in% TE_genes.loci$gene]
a61<-denovo2021.TPMs.genes.Cortijo$noise_average[denovo2021.TPMs.genes.Cortijo$max>1 &denovo2021.TPMs.genes.Cortijo$gene %in% TE_frags.transcripts$gene]

boxplot(        -10,-10,-10,-10,-10,-10,main="Expression noise (Cortijo, seedling)",cex.main=1.2,ylim=c(0,5),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E"), names=c("noTE","TE<50%","TE50-80%","TE>80%","TE genes","TE frags"),las=2, notch = T, outline = F, ylab="coeficient of variance")

#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]),col="darkblue",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]),col="darkgreen",lty=2,lwd = 1.6)
#abline(h=median(denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]),col="orange",lty=2,lwd = 1.6)

boxplot(      a11,a21,a31,a41,a51,a61, ylim=c(0,5),
              col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E"), las=2, notch = T, outline = F, add=T)

#################
#add p values   #
#################
a<-wilcox.test(a11,a21)
b<-wilcox.test(a11,a21)
c<-wilcox.test(a11,a21)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=2.2)
a<-wilcox.test(a21,a31)
b<-wilcox.test(a21,a31)
c<-wilcox.test(a21,a31)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=2.2)
a<-wilcox.test(a31,a41)
b<-wilcox.test(a31,a41)
c<-wilcox.test(a31,a41)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=2.2)
a<-wilcox.test(a41,a51)
b<-wilcox.test(a41,a51)
c<-wilcox.test(a41,a51)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=2.2)
a<-wilcox.test(a61,a51)
b<-wilcox.test(a61,a51)
c<-wilcox.test(a61,a51)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=2.2)

#################
dev.off()


