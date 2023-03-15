######################################################
####### Figure 3: epigenetic patterns of lncRNAs 
######################################################



#install.packages("rlang")
#install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_0.4.10.tar.gz")
install.packages("Z:/programs/Rpackages/rlang_1.0.6.tar.gz", repos = NULL,lib ="Z:/programs/Rpackages" )
install.packages("Z:/programs/Rpackages/vctrs_0.3.8.tar.gz", repos = NULL,lib ="Z:/programs/Rpackages" )
install.packages ("magrittr",lib ="Z:/programs/Rpackages")
install.packages ("ggthemes",lib ="Z:/programs/Rpackages")
install.packages("pillar",lib ="Z:/programs/Rpackages")
install.packages("ggplot2",lib ="Z:/programs/Rpackages")
install.packages("ggfortify",lib ="Z:/programs/Rpackages")
install.packages("ggbiplot",lib ="Z:/programs/Rpackages")
install.packages("plotly",lib ="Z:/programs/Rpackages")
install.packages("Z:/programs/Rpackages/vctrs_0.4.1.tar.gz" , repos=NULL, lib ="Z:/programs/Rpackages")
#install.packages("devtools")
library (vctrs,lib.loc = "Z:/programs/Rpackages")
library (rlang,lib.loc = "Z:/programs/Rpackages")
library (ggthemes,lib ="Z:/programs/Rpackages")
library (pillar,lib.loc = "Z:/programs/Rpackages")
library (ggplot2,lib.loc = "Z:/programs/Rpackages")
library (vioplot)
library(scales)
library(RColorBrewer)
#library(ggfortify,lib.loc = "Z:/programs/Rpackages")
#library(ggbiplot)
library(vioplot)

library(devtools)
install_github("ggbiplot", "vqv")
library(ggbiplot)

########################
#CHIPSEQ quality check
#########################

# PCA on all samples 
#transpose the chipseq coverage 
a<-as.data.frame(t(chip.denovo.allsamples.quantstan[,names(chip.denovo.allsamples.quantstan) %in% names_chip_r$V1]))

#remove columns without variation
a<-a[,apply(a,2,sd)>0]
#do PCA
pca_res <- prcomp(a, scale. = TRUE)
autoplot(pca_res)
PCA_table_chipseq<-as.data.frame(pca_res$x)
PCA_table_chipseq$accession<- as.factor(unlist(lapply(strsplit(as.character(rownames(PCA_table_chipseq)),".", fixed = T), "[",1))) 
PCA_table_chipseq$mark<-as.factor(unlist(lapply(strsplit(as.character(rownames(PCA_table_chipseq)),".", fixed = T), "[",3)))

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_SUPPL_PCA_CHIPSEQ.pdf",height = 4,width = 4)
plot(-200,-200,ylim=c(min(PCA_table_chipseq$PC2),max(PCA_table_chipseq$PC2)),xlim=c(min(PCA_table_chipseq$PC1),max(PCA_table_chipseq$PC1)),ylab="PC2(16.9%)",xlab="PC1 (39.08%)",main="PCA of ChIP-seq samples \nPCA on coverage on 36,877 loci \n from this study's annotation",cex.main =1 )
points(PCA_table_chipseq[PCA_table_chipseq$mark=="H3K9me2",1],PCA_table_chipseq[PCA_table_chipseq$mark=="H3K9me2",2],col="darkmagenta",pch=20)
points(PCA_table_chipseq[PCA_table_chipseq$mark=="H1",1],PCA_table_chipseq[PCA_table_chipseq$mark=="H1",2],col="aquamarine4",pch=20)
points(PCA_table_chipseq[PCA_table_chipseq$mark=="H3K27me3",1],PCA_table_chipseq[PCA_table_chipseq$mark=="H3K27me3",2],col="red",pch=20)
points(PCA_table_chipseq[PCA_table_chipseq$mark=="H3K36me3",1],PCA_table_chipseq[PCA_table_chipseq$mark=="H3K36me3",2],col="goldenrod4",pch=20)
points(PCA_table_chipseq[PCA_table_chipseq$mark=="H3K4me3",1],PCA_table_chipseq[PCA_table_chipseq$mark=="H3K4me3",2],col="darkgreen",pch=20)
dev.off()

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_SUPPL_PCAggplot_CHIPSEQ.pdf",height = 3,width = 3.5)
ggplot(data = PCA_table_chipseq,  mapping = aes(x = PC1, y = PC2)) +
  geom_point(aes(color=mark))
dev.off()


###########################################################
#METHYLATION###############################################
###########################################################
library(vioplot)
#1001G 
pc<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% denovoPC.loci$gene]
as<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% TE_genes.loci$gene]

# methylation levels on different genes 

#CG methylation 1001G vioplot 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_vioplot_CGmeth_6909_from1001Gdata.pdf",height = 3,width = 3.5)
###########################################################
par(mar=c(6,6,2,2)) 
vioplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, main="CG methylation in Col-0",ylab="CG methylation in the locus", rectCol="black", lineCol="black")

#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()



# CG methylation 1001Gnew data (this study) vioplot
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_vioplot_CGmeth_6909_from1001GNEWdata.pdf",height = 3,width = 3.5)
###########################################################
par(mar=c(6,6,2,2)) 
pc<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]
as<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]

vioplot( pc,as,linc,te,
         col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, main="CG methylation in Col-0",ylab="CG methylation in the locus", rectCol="black", lineCol="black")
###########################################################
dev.off()


#CG 1001G new data boxplot
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CGmeth_6909_1001GNEW.pdf",height = 3,width = 4)
###########################################################
par(mar=c(6,6,2,2)) 
pc<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]
as<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        notch = T,outline = F,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, main="CG methylation in Col-0",ylab="CG methylation in the locus")
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
#################
dev.off()


#CHG 1001G new data boxplot
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CHGmeth_6909_1001GNEW.pdf",height = 3,width = 4)
###########################################################
par(mar=c(6,6,2,2)) 
pc<-CHG.1001new.denovo$mean.6909[CHG.1001new.denovo$transcript %in% denovoPC.loci$gene]
as<-CHG.1001new.denovo$mean.6909[CHG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CHG.1001new.denovo$mean.6909[CHG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CHG.1001new.denovo$mean.6909[CHG.1001new.denovo$transcript %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        notch = T,outline = F,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, main="CHG methylation in Col-0",ylab="CHG methylation in the locus")
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.6)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.6)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.6)
#################
dev.off()


#CHH 1001G new data boxplot
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CHHmeth_6909_1001GNEW.pdf",height = 3,width = 4)
###########################################################
par(mar=c(6,6,2,2)) 
pc<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% denovoPC.loci$gene]
as<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        notch = T,outline = F,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, main="CHH methylation in Col-0",ylab="CHH methylation in the locus")
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.1)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.1)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.1)
#################
dev.off()


#CG 1001G boxplot
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CGmeth_6909_from1001Gdata.pdf",height = 3,width = 4)
###########################################################
par(mar=c(6,6,2,2)) 
pc<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% denovoPC.loci$gene]
as<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        notch = T,outline = F,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, main="CG methylation in Col-0",ylab="CG methylation in the locus")

#stripchart(pc,method = "jitter",pch=1,col = alpha("black", 0.4), bg=c("orange","black"),jitter=0.1,vertical = TRUE,add = TRUE,cex = 0.6,at=1) 

#stripchart(as,method = "jitter",pch=1,col = alpha("black", 0.4), bg=c("orange","black"),jitter=0.1,vertical = TRUE,add = TRUE,cex = 0.6,at=2) 

#stripchart(linc,method = "jitter",pch=1,col = alpha("black", 0.4), bg=c("orange","black"),jitter=0.1,vertical = TRUE,add = TRUE,cex = 0.6,at=3) 

#stripchart(te,method = "jitter",pch=1,col = alpha("black", 0.4), bg=c("orange","black"),jitter=0.1,vertical = TRUE,add = TRUE,cex = 0.6,at=4) 

#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
#################
dev.off()

#CHG 1001G boxplot
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CHGmeth_6909_from1001Gdata.pdf",height = 3,width = 4)
###########################################################
par(mar=c(6,6,2,2)) 
pc<-CHG.1001.denovo$X6909[CHG.1001.denovo$transcript %in% denovoPC.loci$gene]
as<-CHG.1001.denovo$X6909[CHG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CHG.1001.denovo$X6909[CHG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CHG.1001.denovo$X6909[CHG.1001.denovo$transcript %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        notch = T,outline = F,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, main="CHG methylation in Col-0",ylab="CHG methylation in the locus")

#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()


#CHH1001G boxplot
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CHHmeth_6909_from1001Gdata.pdf",height = 3,width = 4)
###########################################################
par(mar=c(6,6,2,2)) 
pc<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% denovoPC.loci$gene]
as<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        notch = T,outline = F,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, main="CHH methylation in Col-0",ylab="CHH methylation in the locus")

#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.1)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.1)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.1)
#################
dev.off()






median(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% denovoPC.loci$gene],na.rm = T)#0.04732935
median(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene])# 0.00558659
median(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene])#0.1640385
median(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% TE_genes.loci$gene])#0.8577525

length(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% denovoPC.loci$gene])#23676
length(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene])# 8195
length(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene])#2246
length(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% TE_genes.loci$gene])#2130


col=c("#486EB4","#90C473","#F2AB54","#673A8E")

# bimodal distribution of lincRNAs 

#1001G data 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/densityplot_CG_from1001Gdata.pdf",height = 3.5,width = 3.5)
###########################################################
par(mar=c(3,3,2,2),mfrow=c(2,2)) 
plot(density(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% denovoPC.loci$gene]),xlab="CG methylation level",main="",lwd=2,las=2, col="#486EB4")
mtext("PC genes",side=3,at=0.6,line = -1)
plot(density(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]),xlab="CG methylation level",main="",lwd=2,las=2, col="#90C473")
mtext("AS lncRNAs",side=3,at=0.6,line = -1)
plot(density(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),xlab="CG methylation level",main="",lwd=2,las=2, col="#F2AB54")
mtext("lincRNAs",side=3,at=0.6,line = -1)
plot(density(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% TE_genes.loci$gene]),xlab="CG methylation level",lwd=2,las=2, main="",col="#673A8E")
mtext("TE genes",side=3,at=0.2,line = -1)
###########################################################
dev.off()

#1001G new data  CG
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/densityplot_CG_from1001GNEW.pdf",height = 3.5,width = 3.5)
###########################################################
par(mar=c(3,3,2,2),mfrow=c(2,2)) 
plot(density(CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),xlab="CG methylation level",main="",lwd=2,las=2, col="#486EB4")
mtext("PC genes",side=3,at=0.6,line = -1)
plot(density(CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),xlab="CG methylation level",main="",lwd=2,las=2, col="#90C473")
mtext("AS lncRNAs",side=3,at=0.6,line = -1)
plot(density(CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),xlab="CG methylation level",main="",lwd=2,las=2, col="#F2AB54")
mtext("lincRNAs",side=3,at=0.6,line = -1)
plot(density(CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),xlab="CG methylation level",lwd=2,las=2, main="",col="#673A8E")
mtext("TE genes",side=3,at=0.2,line = -1)
###########################################################
dev.off()

#1001G new data CHG
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/densityplot_CHG_from1001GNEW.pdf",height = 3.5,width = 3.5)
###########################################################
par(mar=c(3,3,2,2),mfrow=c(2,2)) 
plot(density(CHG.1001new.denovo$mean.6909[CHG.1001new.denovo$transcript %in% denovoPC.loci$gene]),xlab="CG methylation level",main="",lwd=2,las=2, col="#486EB4")
mtext("PC genes",side=3,at=0.6,line = -1)
plot(density(CHG.1001new.denovo$mean.6909[CHG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),xlab="CG methylation level",main="",lwd=2,las=2, col="#90C473")
mtext("AS lncRNAs",side=3,at=0.6,line = -1)
plot(density(CHG.1001new.denovo$mean.6909[CHG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),xlab="CG methylation level",main="",lwd=2,las=2, col="#F2AB54")
mtext("lincRNAs",side=3,at=0.6,line = -1)
plot(density(CHG.1001new.denovo$mean.6909[CHG.1001new.denovo$transcript %in% TE_genes.loci$gene]),xlab="CG methylation level",lwd=2,las=2, main="",col="#673A8E")
mtext("TE genes",side=3,at=0.2,line = -1)
###########################################################
dev.off()


#1001G new data CHH
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/densityplot_CHH_from1001GNEW.pdf",height = 3.5,width = 3.5)
###########################################################
par(mar=c(3,3,2,2),mfrow=c(2,2)) 
plot(density(CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% denovoPC.loci$gene]),xlab="CG methylation level",main="",lwd=2,las=2, col="#486EB4")
mtext("PC genes",side=3,at=0.6,line = -1)
plot(density(CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),xlab="CG methylation level",main="",lwd=2,las=2, col="#90C473")
mtext("AS lncRNAs",side=3,at=0.6,line = -1)
plot(density(CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),xlab="CG methylation level",main="",lwd=2,las=2, col="#F2AB54")
mtext("lincRNAs",side=3,at=0.6,line = -1)
plot(density(CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% TE_genes.loci$gene]),xlab="CG methylation level",lwd=2,las=2, main="",col="#673A8E")
mtext("TE genes",side=3,at=0.2,line = -1)
###########################################################
dev.off()

# bimodal distribution 
plot(density(CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]))
plot(density(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]))
plot(density(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% denovoPC.loci$gene]))
plot(density(CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% TE_genes.loci$gene]))




#CHG  1001G data vioplot
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_vioplot_CHGmeth_6909_from1001Gdata.pdf",height = 3,width = 3.5)
#########################################
par(mar=c(6,6,2,2)) 
pc<-CHG.1001.denovo$X6909[CHG.1001.denovo$transcript %in% denovoPC.loci$gene]
as<-CHG.1001.denovo$X6909[CHG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CHG.1001.denovo$X6909[CHG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CHG.1001.denovo$X6909[CHG.1001.denovo$transcript %in% TE_genes.loci$gene]

vioplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, main="CHG methylation in Col-0",ylab="CHG methylation in the locus", rectCol="black", lineCol="black")

#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()

# CHG 1001Gnew  data
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_vioplot_CHGmeth_6909_from1001GNEWdata.pdf",height = 3,width = 3.5)
#########################################
par(mar=c(6,6,2,2)) 
pc<-CHG.1001new.denovo$mean.6909[CHG.1001new.denovo$transcript %in% denovoPC.loci$gene]
as<-CHG.1001new.denovo$mean.6909[CHG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CHG.1001new.denovo$mean.6909[CHG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CHG.1001new.denovo$mean.6909[CHG.1001new.denovo$transcript %in% TE_genes.loci$gene]

vioplot( pc,as,linc,te,
         col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, main="CHG methylation in Col-0",ylab="CHG methylation in the locus", rectCol="black", lineCol="black")

#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()


# CHH 1001G data vioplot
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_vioplot_CHHmeth_6909_from1001Gdata.pdf",height = 3,width = 3.5)
#########################################
par(mar=c(6,6,2,2)) 
pc<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% denovoPC.loci$gene]
as<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% TE_genes.loci$gene]
vioplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, main="CHH methylation in Col-0",ylab="CHH methylation in the locus", rectCol="black", lineCol="black")

#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()


# CHH 1001G NEW data vioplot
pc<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% denovoPC.loci$gene]
as<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% TE_genes.loci$gene]

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_vioplot_CHHmeth_6909_from1001GNEWdata.pdf",height = 3,width = 3.5)
#########################################
par(mar=c(6,6,2,2)) 
vioplot( pc,as,linc,te,
         col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, main="CHH methylation in Col-0",ylab="CHH methylation in the locus", rectCol="black", lineCol="black")
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()




#CG 1001G vioplot 
pc<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% denovoPC.loci$gene]
as<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% TE_genes.loci$gene]

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_vioplot_CGmeth_6909_from1001Gdata.pdf",height = 3,width = 3.5)
#########################################
par(mar=c(6,6,2,2)) 
vioplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, main="CG methylation in Col-0",ylab="CG methylation in the locus", rectCol="black", lineCol="black")

#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################

dev.off()


###################################################################
########### plot CG/CHG/CHH methylation along chromosomes in 6909 1001Gnew 
###################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/CG_CHG_CHH_alongchr.Chrposition.1001Gnew.pdf",height = 3,width =8)
###################################################################
par(mar=c(3,3,3,0) + 0.1)
par(mgp=c(1.5,1,0))
par(mfrow=c(1,4))

#PC   
a<-merge(denovoPC.loci,CG.1001new.denovo[,c("transcript","mean.6909")],by.x="gene",by.y="transcript")
a<-merge(a,CHG.1001new.denovo[,c("transcript","mean.6909")],by.x="gene",by.y="transcript")
a<-merge(a,CHH.1001new.denovo[,c("transcript","mean.6909")],by.x="gene",by.y="transcript")
a$start[a$chr=="Chr1"]<-a$start[a$chr=="Chr1"]
a$end[a$chr=="Chr1"]<-a$end[a$chr=="Chr1"]
a$start[a$chr=="Chr2"]<-a$start[a$chr=="Chr2"]+30427671
a$end[a$chr=="Chr2"]<-a$end[a$chr=="Chr2"]+30427671
a$start[a$chr=="Chr3"]<-a$start[a$chr=="Chr3"]+30427671+19698289
a$end[a$chr=="Chr3"]<-a$end[a$chr=="Chr3"]+30427671+19698289
a$start[a$chr=="Chr4"]<-a$start[a$chr=="Chr4"]+30427671+19698289+23459830
a$end[a$chr=="Chr4"]<-a$end[a$chr=="Chr4"]+30427671+19698289+23459830
a$start[a$chr=="Chr5"]<-a$start[a$chr=="Chr5"]+30427671+19698289+23459830+18585056
a$end[a$chr=="Chr5"]<-a$end[a$chr=="Chr5"]+30427671+19698289+23459830+18585056

plot(a$start,a$mean.6909.x+2,pch=20,cex=0.4,col=alpha("tomato4",alpha=0.4),ylim=c(-2,3),ylab=" methylation level",xlab="gene start position",main="PC genes",yaxt = "n",xaxt = "n")
axis(1,at = c(30427671,30427671+19698289,30427671+19698289+23459830,30427671+19698289+23459830+18585056,30427671+19698289+23459830+18585056+26975502 ), labels =c("","","","","") ,las=2)
axis(1,at = c(30427671/2,30427671+19698289/2,30427671+19698289+23459830/2,30427671+19698289+23459830+18585056/2,30427671+19698289+23459830+18585056+26975502/2 ), labels =c("Chr1","Chr2","Chr3","Chr4","Chr5") ,tick = F,las=1,mgp=c(0,0,0),cex.axis=0.8)
axis(2,at = -2:3, labels =c(0,1,0,1,0,1) ,las=2)
axis(2,at = c(-1.5,0.5,2.5), labels =c("CHH","CHG","CG") ,las=3,mgp=c(0.3,0,0),tick = F)
points(a$start,a$mean.6909.y,pch=20,cex=0.4,col=alpha("tomato3",alpha=0.4))
points(a$start,a$mean.6909-2,pch=20,cex=0.4,col=alpha("tomato1",alpha=0.4))
#add centromeres
#add (approximate) centromere center positions
abline(v=15000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+4700000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+13000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+3800000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+18585056+11900000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
#AS   
a<-merge(lncRNAs.antisense.loci,CG.1001new.denovo[,c("transcript","mean.6909")],by.x="gene",by.y="transcript")
a<-merge(a,CHG.1001new.denovo[,c("transcript","mean.6909")],by.x="gene",by.y="transcript")
a<-merge(a,CHH.1001new.denovo[,c("transcript","mean.6909")],by.x="gene",by.y="transcript")
a$start[a$chr=="Chr1"]<-a$start[a$chr=="Chr1"]
a$end[a$chr=="Chr1"]<-a$end[a$chr=="Chr1"]
a$start[a$chr=="Chr2"]<-a$start[a$chr=="Chr2"]+30427671
a$end[a$chr=="Chr2"]<-a$end[a$chr=="Chr2"]+30427671
a$start[a$chr=="Chr3"]<-a$start[a$chr=="Chr3"]+30427671+19698289
a$end[a$chr=="Chr3"]<-a$end[a$chr=="Chr3"]+30427671+19698289
a$start[a$chr=="Chr4"]<-a$start[a$chr=="Chr4"]+30427671+19698289+23459830
a$end[a$chr=="Chr4"]<-a$end[a$chr=="Chr4"]+30427671+19698289+23459830
a$start[a$chr=="Chr5"]<-a$start[a$chr=="Chr5"]+30427671+19698289+23459830+18585056
a$end[a$chr=="Chr5"]<-a$end[a$chr=="Chr5"]+30427671+19698289+23459830+18585056

plot(a$start,a$mean.6909.x+2,pch=20,cex=0.4,col=alpha("tomato4",alpha=0.4),ylim=c(-2,3),ylab=" methylation level",xlab="gene start position",main="AS lncRNAs",yaxt = "n",xaxt = "n")
axis(1,at = c(30427671,30427671+19698289,30427671+19698289+23459830,30427671+19698289+23459830+18585056,30427671+19698289+23459830+18585056+26975502 ), labels =c("","","","","") ,las=2)
axis(1,at = c(30427671/2,30427671+19698289/2,30427671+19698289+23459830/2,30427671+19698289+23459830+18585056/2,30427671+19698289+23459830+18585056+26975502/2 ), labels =c("Chr1","Chr2","Chr3","Chr4","Chr5") ,tick = F,las=1,mgp=c(0,0,0),cex.axis=0.8)
axis(2,at = -2:3, labels =c(0,1,0,1,0,1) ,las=2)
axis(2,at = c(-1.5,0.5,2.5), labels =c("CHH","CHG","CG") ,las=3,mgp=c(0.3,0,0),tick = F)
points(a$start,a$mean.6909.y,pch=20,cex=0.4,col=alpha("tomato3",alpha=0.4))
points(a$start,a$mean.6909-2,pch=20,cex=0.4,col=alpha("tomato1",alpha=0.4))
#add (approximate) centromere center positions
abline(v=15000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+4700000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+13000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+3800000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+18585056+11900000,col=alpha("black",alpha=0.5), lty=2, lwd=1)

#linc  
a<-merge(lncRNAs.intergenic.loci,CG.1001new.denovo[,c("transcript","mean.6909")],by.x="gene",by.y="transcript")
a<-merge(a,CHG.1001new.denovo[,c("transcript","mean.6909")],by.x="gene",by.y="transcript")
a<-merge(a,CHH.1001new.denovo[,c("transcript","mean.6909")],by.x="gene",by.y="transcript")
a$start[a$chr=="Chr1"]<-a$start[a$chr=="Chr1"]
a$end[a$chr=="Chr1"]<-a$end[a$chr=="Chr1"]
a$start[a$chr=="Chr2"]<-a$start[a$chr=="Chr2"]+30427671
a$end[a$chr=="Chr2"]<-a$end[a$chr=="Chr2"]+30427671
a$start[a$chr=="Chr3"]<-a$start[a$chr=="Chr3"]+30427671+19698289
a$end[a$chr=="Chr3"]<-a$end[a$chr=="Chr3"]+30427671+19698289
a$start[a$chr=="Chr4"]<-a$start[a$chr=="Chr4"]+30427671+19698289+23459830
a$end[a$chr=="Chr4"]<-a$end[a$chr=="Chr4"]+30427671+19698289+23459830
a$start[a$chr=="Chr5"]<-a$start[a$chr=="Chr5"]+30427671+19698289+23459830+18585056
a$end[a$chr=="Chr5"]<-a$end[a$chr=="Chr5"]+30427671+19698289+23459830+18585056

plot(a$start,a$mean.6909.x+2,pch=20,cex=0.4,col=alpha("tomato4",alpha=0.4),ylim=c(-2,3),ylab=" methylation level",xlab="gene start position",main="lincRNAs",yaxt = "n",xaxt = "n")
axis(1,at = c(30427671,30427671+19698289,30427671+19698289+23459830,30427671+19698289+23459830+18585056,30427671+19698289+23459830+18585056+26975502 ), labels =c("","","","","") ,las=2)
axis(1,at = c(30427671/2,30427671+19698289/2,30427671+19698289+23459830/2,30427671+19698289+23459830+18585056/2,30427671+19698289+23459830+18585056+26975502/2 ), labels =c("Chr1","Chr2","Chr3","Chr4","Chr5") ,tick = F,las=1,mgp=c(0,0,0),cex.axis=0.8)
axis(2,at = -2:3, labels =c(0,1,0,1,0,1) ,las=2)
axis(2,at = c(-1.5,0.5,2.5), labels =c("CHH","CHG","CG") ,las=3,mgp=c(0.3,0,0),tick = F)
points(a$start,a$mean.6909.y,pch=20,cex=0.4,col=alpha("tomato3",alpha=0.4))
points(a$start,a$mean.6909-2,pch=20,cex=0.4,col=alpha("tomato1",alpha=0.4))
#add (approximate) centromere center positions
abline(v=15000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+4700000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+13000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+3800000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+18585056+11900000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
#TE genes   
a<-merge(TE_genes.loci,CG.1001new.denovo[,c("transcript","mean.6909")],by.x="gene",by.y="transcript")
a<-merge(a,CHG.1001new.denovo[,c("transcript","mean.6909")],by.x="gene",by.y="transcript")
a<-merge(a,CHH.1001new.denovo[,c("transcript","mean.6909")],by.x="gene",by.y="transcript")
a$start[a$chr=="Chr1"]<-a$start[a$chr=="Chr1"]
a$end[a$chr=="Chr1"]<-a$end[a$chr=="Chr1"]
a$start[a$chr=="Chr2"]<-a$start[a$chr=="Chr2"]+30427671
a$end[a$chr=="Chr2"]<-a$end[a$chr=="Chr2"]+30427671
a$start[a$chr=="Chr3"]<-a$start[a$chr=="Chr3"]+30427671+19698289
a$end[a$chr=="Chr3"]<-a$end[a$chr=="Chr3"]+30427671+19698289
a$start[a$chr=="Chr4"]<-a$start[a$chr=="Chr4"]+30427671+19698289+23459830
a$end[a$chr=="Chr4"]<-a$end[a$chr=="Chr4"]+30427671+19698289+23459830
a$start[a$chr=="Chr5"]<-a$start[a$chr=="Chr5"]+30427671+19698289+23459830+18585056
a$end[a$chr=="Chr5"]<-a$end[a$chr=="Chr5"]+30427671+19698289+23459830+18585056

plot(a$start,a$mean.6909.x+2,pch=20,cex=0.4,col=alpha("tomato4",alpha=0.4),ylim=c(-2,3),ylab=" methylation level",xlab="gene start position",main="TE genes",yaxt = "n",xaxt = "n")
axis(1,at = c(30427671,30427671+19698289,30427671+19698289+23459830,30427671+19698289+23459830+18585056,30427671+19698289+23459830+18585056+26975502 ), labels =c("","","","","") ,las=2)
axis(1,at = c(30427671/2,30427671+19698289/2,30427671+19698289+23459830/2,30427671+19698289+23459830+18585056/2,30427671+19698289+23459830+18585056+26975502/2 ), labels =c("Chr1","Chr2","Chr3","Chr4","Chr5") ,tick = F,las=1,mgp=c(0,0,0),cex.axis=0.8)
axis(2,at = -2:3, labels =c(0,1,0,1,0,1) ,las=2)
axis(2,at = c(-1.5,0.5,2.5), labels =c("CHH","CHG","CG") ,las=3,mgp=c(0.3,0,0),tick = F)
points(a$start,a$mean.6909.y,pch=20,cex=0.4,col=alpha("tomato3",alpha=0.4))
points(a$start,a$mean.6909-2,pch=20,cex=0.4,col=alpha("tomato1",alpha=0.4))
#add (approximate) centromere center positions
abline(v=15000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+4700000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+13000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+3800000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+18585056+11900000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
##############################
dev.off()



# boxplot CG  1001Gnew distant to centromere close to centromere
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CGmeth_6909_from1001GNEWdata.distant.nondist.pdf",height = 3,width = 5)
###########################################
par(mar=c(6,6,2,2)) 
pc1<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere<2000000] ]
pc2<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere>2000000] ]
as1<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere<2000000]]
as2<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere>2000000]]
linc1<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<2000000]]
linc2<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>2000000]]
te1<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere<2000000]]
te2<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere>2000000]]
boxplot( pc1,as1,linc1,te1,pc2,as2,linc2,te2,
         col=c("#486EB4","#90C473","#F2AB54","#673A8E","#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE","PC","AS","linc","TE"),las=2, main="CG methylation in Col-0",ylab="CG methylation in the locus", ylim=c(0,1.1),outline = F, notch = T)
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
text(b,x=1.5,y=0.7)
a<-wilcox.test(linc1,as1)
b<-wilcox.test(linc1,as1)
c<-wilcox.test(linc1,as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.7)
a<-wilcox.test(linc1,te1)
b<-wilcox.test(linc1,te1)
c<-wilcox.test(linc1,te1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.7)
a<-wilcox.test(sample(pc2,1612),sample(as2,1612))
b<-wilcox.test(sample(pc2,1612),sample(as2,1612))
c<-wilcox.test(sample(pc2,1612),sample(as2,1612))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.7)
a<-wilcox.test(sample(linc2,1612),sample(as2,1612))
b<-wilcox.test(sample(linc2,1612),sample(as2,1612))
c<-wilcox.test(sample(linc2,1612),sample(as2,1612))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=0.7)
a<-wilcox.test(sample(linc2,815),te2)
b<-wilcox.test(sample(linc2,815),te2)
c<-wilcox.test(sample(linc2,815),te2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.7)

a<-wilcox.test(pc1,sample(pc2,1012))
b<-wilcox.test(pc1,sample(pc2,1012))
c<-wilcox.test(pc1,sample(pc2,1012))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=1.1)
a<-wilcox.test(sample(as2,333),as1)
b<-wilcox.test(sample(as2,333),as1)
c<-wilcox.test(sample(as2,333),as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=1.1)
a<-wilcox.test(linc1,linc2)
b<-wilcox.test(linc1,linc2)
c<-wilcox.test(linc1,linc2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=1.1)
a<-wilcox.test(te1,te2)
b<-wilcox.test(te1,te2)
c<-wilcox.test(te1,te2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=1.1)
#################
dev.off()

# boxplot CHH  1001Gnew distant to centromere close to centromere
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CHHmeth_6909_from1001GNEWdata.distant.nondist.pdf",height = 3,width = 5)
###########################################
par(mar=c(6,6,2,2)) 
pc1<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere<2000000] ]
pc2<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere>2000000] ]
as1<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere<2000000]]
as2<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere>2000000]]
linc1<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<2000000]]
linc2<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>2000000]]
te1<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere<2000000]]
te2<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere>2000000]]
boxplot( pc1,as1,linc1,te1,pc2,as2,linc2,te2,
         col=c("#486EB4","#90C473","#F2AB54","#673A8E","#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE","PC","AS","linc","TE"),las=2, main="CG methylation in Col-0",ylab="CG methylation in the locus", ylim=c(0,0.25),outline = F, notch = T)
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
text(b,x=1.5,y=0.14)
a<-wilcox.test(linc1,as1)
b<-wilcox.test(linc1,as1)
c<-wilcox.test(linc1,as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.14)
a<-wilcox.test(linc1,te1)
b<-wilcox.test(linc1,te1)
c<-wilcox.test(linc1,te1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.14)
a<-wilcox.test(sample(pc2,1612),sample(as2,1612))
b<-wilcox.test(sample(pc2,1612),sample(as2,1612))
c<-wilcox.test(sample(pc2,1612),sample(as2,1612))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.14)
a<-wilcox.test(sample(linc2,1612),sample(as2,1612))
b<-wilcox.test(sample(linc2,1612),sample(as2,1612))
c<-wilcox.test(sample(linc2,1612),sample(as2,1612))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=0.14)
a<-wilcox.test(sample(linc2,815),te2)
b<-wilcox.test(sample(linc2,815),te2)
c<-wilcox.test(sample(linc2,815),te2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.14)

a<-wilcox.test(pc1,sample(pc2,1012))
b<-wilcox.test(pc1,sample(pc2,1012))
c<-wilcox.test(pc1,sample(pc2,1012))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=0.2)
a<-wilcox.test(sample(as2,333),as1)
b<-wilcox.test(sample(as2,333),as1)
c<-wilcox.test(sample(as2,333),as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=0.2)
a<-wilcox.test(linc1,linc2)
b<-wilcox.test(linc1,linc2)
c<-wilcox.test(linc1,linc2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=0.2)
a<-wilcox.test(te1,te2)
b<-wilcox.test(te1,te2)
c<-wilcox.test(te1,te2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=0.2)
#################
dev.off()




# vioplot CG 1001G distant to centromere close to centromere
############################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_vioplot_CGmeth_6909_from1001data.distant.nondist.pdf",height = 3,width = 5)
############################
par(mar=c(6,6,2,2)) 
pc1<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere<1000000] ]
pc2<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere>1000000] ]

as1<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere<1000000]]
as2<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere>1000000]]

linc1<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
linc2<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

te1<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere<1000000]]
te2<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere>1000000]]

vioplot( pc1,as1,linc1,te1,pc2,as2,linc2,te2,
         col=c("#486EB4","#90C473","#F2AB54","#673A8E","#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE","PC","AS","linc","TE"),las=2, main="CG methylation in Col-0",ylab="CG methylation in the locus", ylim=c(0,1.1),rectCol="black", lineCol="black",pchMed=18,colMed = "red")

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


# vioplot CHH 1001G distant to centromere close to centromere
############################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_vioplot_CHHmeth_6909_from1001data.distant.nondist.pdf",height = 3,width = 5)
############################
library(vioplot)
par(mar=c(6,6,2,2)) 
pc1<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere<1000000] ]
pc2<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere>1000000] ]

as1<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere<1000000]]
as2<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere>1000000]]

linc1<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<1000000]]
linc2<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>1000000]]

te1<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere<1000000]]
te2<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere>1000000]]

boxplot( pc1,as1,linc1,te1,pc2,as2,linc2,te2,outline = F,notch = T,
         col=c("#486EB4","#90C473","#F2AB54","#673A8E","#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE","PC","AS","linc","TE"),las=2, main="CHH methylation in Col-0",ylab="CHH methylation in the locus", ylim=c(0,0.12),rectCol="black", lineCol="black",pchMed=18,colMed = "red")

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
text(b,x=1.5,y=0.09)
a<-wilcox.test(linc1,as1)
b<-wilcox.test(linc1,as1)
c<-wilcox.test(linc1,as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.09)
a<-wilcox.test(linc1,te1)
b<-wilcox.test(linc1,te1)
c<-wilcox.test(linc1,te1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.09)
a<-wilcox.test(sample(pc2,2000),sample(as2,2000))
b<-wilcox.test(sample(pc2,2000),sample(as2,2000))
c<-wilcox.test(sample(pc2,2000),sample(as2,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.09)
a<-wilcox.test(sample(linc2,2000),sample(as2,2000))
b<-wilcox.test(sample(linc2,2000),sample(as2,2000))
c<-wilcox.test(sample(linc2,2000),sample(as2,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=0.09)
a<-wilcox.test(sample(linc2,2000),te2)
b<-wilcox.test(sample(linc2,2000),te2)
c<-wilcox.test(sample(linc2,2000),te2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.09)

a<-wilcox.test(pc1,sample(pc2,2000))
b<-wilcox.test(pc1,sample(pc2,2000))
c<-wilcox.test(pc1,sample(pc2,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=0.12)
a<-wilcox.test(sample(as2,2000),as1)
b<-wilcox.test(sample(as2,2000),as1)
c<-wilcox.test(sample(as2,2000),as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=0.12)
a<-wilcox.test(linc1,linc2)
b<-wilcox.test(linc1,linc2)
c<-wilcox.test(linc1,linc2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=0.12)
a<-wilcox.test(te1,te2)
b<-wilcox.test(te1,te2)
c<-wilcox.test(te1,te2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=0.12)
#################
dev.off()



############### 
#Methylation gene on gene off 
###############

#boxplot CG gene body  1001G data
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CG_expressed_silent_6909_allgenes.1001Gdata.pdf",height = 3,width =3)
###################################################
par(mar=c(6,3,3,2)) 
a1<-  CG.1001.denovo$X6909 [CG.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]]
a2<-     CG.1001.denovo$X6909 [CG.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5 & denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]]
a3<-      CG.1001.denovo$X6909 [CG.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]]
a4<-     CG.1001.denovo$X6909 [CG.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]]
a5<-  CG.1001.denovo$X6909 [CG.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]]
a6<-     CG.1001.denovo$X6909 [CG.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]]
a7<-     CG.1001.denovo$X6909 [CG.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]]
a8<-  CG.1001.denovo$X6909 [CG.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5& denovo2021.TPMs.genes.1001G$gene %in% TE_genes.loci$gene]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,
        col=c("#486EB4","#8aa7de","#90C473","#d4edc5","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="CG meth. level in the locus")
mtext('CG', side=3, line=-1, at=1,cex=0.9)
###################################################
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
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(a5,a6)
b<-wilcox.test(a5,a6)
c<-wilcox.test(a5,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
#################
dev.off()

#boxplot CG promoter 1001G data
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CG_TSS_expressed_silent_6909_allgenes.1001Gdata.pdf",height = 3,width =3)
###################################################
par(mar=c(6,3,3,2)) 
a1<-  CG.1001.denovo_TSS$X6909 [CG.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]]
a2<-     CG.1001.denovo_TSS$X6909 [CG.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5 & denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]]
a3<-      CG.1001.denovo_TSS$X6909 [CG.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]]
a4<-     CG.1001.denovo_TSS$X6909 [CG.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]]
a5<-  CG.1001.denovo_TSS$X6909 [CG.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]]
a6<-     CG.1001.denovo_TSS$X6909 [CG.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]]
a7<-     CG.1001.denovo_TSS$X6909 [CG.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]]
a8<-  CG.1001.denovo_TSS$X6909 [CG.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5& denovo2021.TPMs.genes.1001G$gene %in% TE_genes.loci$gene]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,
        col=c("#486EB4","#8aa7de","#90C473","#d4edc5","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="CG meth. level in the locus")
mtext('CG', side=3, line=-1, at=1,cex=0.9)
###################################################
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
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(a5,a6)
b<-wilcox.test(a5,a6)
c<-wilcox.test(a5,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
#################
dev.off()


#boxplot CHH gene body 1001G data
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CHH_expressed_silent_6909_allgenes.1001Gdata.pdf",height = 3,width =3)
###################################################
par(mar=c(6,3,3,2)) 
a1<-  CHH.1001.denovo$X6909 [CHH.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]]
a2<-     CHH.1001.denovo$X6909 [CHH.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5 & denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]]
a3<-      CHH.1001.denovo$X6909 [CHH.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]]
a4<-     CHH.1001.denovo$X6909 [CHH.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]]
a5<-  CHH.1001.denovo$X6909 [CHH.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]]
a6<-     CHH.1001.denovo$X6909 [CHH.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]]
a7<-     CHH.1001.denovo$X6909 [CHH.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]]
a8<-  CHH.1001.denovo$X6909 [CHH.1001.denovo$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5& denovo2021.TPMs.genes.1001G$gene %in% TE_genes.loci$gene]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,
        col=c("#486EB4","#8aa7de","#90C473","#d4edc5","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="CHH meth. level in the locus")
mtext('CHH', side=3, line=-1, at=1,cex=0.9)
###################################################
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
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.1)
a<-wilcox.test(a5,a6)
b<-wilcox.test(a5,a6)
c<-wilcox.test(a5,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.1)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.1)
#################
dev.off()

#boxplot CHH promoter 1001G data
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CHH_TSS_expressed_silent_6909_allgenes.1001Gdata.pdf",height = 3,width =3)
###################################################
par(mar=c(6,3,3,2)) 
a1<-  CHH.1001.denovo_TSS$X6909 [CHH.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]]
a2<-     CHH.1001.denovo_TSS$X6909 [CHH.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5 & denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]]
a3<-      CHH.1001.denovo_TSS$X6909 [CHH.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]]
a4<-     CHH.1001.denovo_TSS$X6909 [CHH.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]]
a5<-  CHH.1001.denovo_TSS$X6909 [CHH.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]]
a6<-     CHH.1001.denovo_TSS$X6909 [CHH.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]]
a7<-     CHH.1001.denovo_TSS$X6909 [CHH.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]]
a8<-  CHH.1001.denovo_TSS$X6909 [CHH.1001.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<0.5& denovo2021.TPMs.genes.1001G$gene %in% TE_genes.loci$gene]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,
        col=c("#486EB4","#8aa7de","#90C473","#d4edc5","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="CHH meth. level in the locus")
mtext('CHH', side=3, line=-1, at=1,cex=0.9)
###################################################
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
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.1)
a<-wilcox.test(a5,a6)
b<-wilcox.test(a5,a6)
c<-wilcox.test(a5,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.1)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.1)
#################
dev.off()




#boxplot CG gene body  1001GNEW data
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Suppl_Fig3_boxplot_CG_expressed_silent.1001GNEWdata.6909mean.pdf",height = 3,width =3)
###################################################
par(mar=c(6,3,3,2)) 
a1<-  CG.1001new.denovo$mean.6909  [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]]
a2<-     CG.1001new.denovo$mean.6909 [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]]
a3<-      CG.1001new.denovo$mean.6909 [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]]
a4<-     CG.1001new.denovo$mean.6909 [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]]
a5<-  CG.1001new.denovo$mean.6909 [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a6<-     CG.1001new.denovo$mean.6909 [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a7<-     CG.1001new.denovo$mean.6909 [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a8<-  CG.1001new.denovo$mean.6909 [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<0.5& denovo2021.TPMs.genes.1001Gnew$gene %in% TE_genes.loci$gene]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,
        col=c("#486EB4","#8aa7de","#90C473","#d4edc5","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="CG meth. level in the locus")
mtext("CG", side=3, line=-1, at=1,cex=0.9)
###################################################
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
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(a5,a6)
b<-wilcox.test(a5,a6)
c<-wilcox.test(a5,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
#################
dev.off()

#boxplot CG gene body  1001GNEW data
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Suppl_Fig3_boxplot_CG_expressed_silent.1001GNEWdata.6909exp1.pdf",height = 3,width =3)
###################################################
par(mar=c(6,3,3,2)) 
a1<-  CG.1001new.denovo$r14.Exp1.6909  [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]]
a2<-     CG.1001new.denovo$r14.Exp1.6909 [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]]
a3<-      CG.1001new.denovo$r14.Exp1.6909 [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]]
a4<-     CG.1001new.denovo$r14.Exp1.6909 [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]]
a5<-  CG.1001new.denovo$r14.Exp1.6909 [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a6<-     CG.1001new.denovo$r14.Exp1.6909 [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a7<-     CG.1001new.denovo$r14.Exp1.6909 [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a8<-  CG.1001new.denovo$r14.Exp1.6909 [CG.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5& denovo2021.TPMs.genes.1001Gnew$gene %in% TE_genes.loci$gene]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,
        col=c("#486EB4","#8aa7de","#90C473","#d4edc5","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="CG meth. level in the locus")
mtext("CG", side=3, line=-1, at=1,cex=0.9)
###################################################
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
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(a5,a6)
b<-wilcox.test(a5,a6)
c<-wilcox.test(a5,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
#################
dev.off()



#boxplot CG promoter  1001GNEW data
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Suppl_Fig3_boxplot_CG_expressed_silent.1001GNEWdata.6909mean.PROMOTER.pdf",height = 3,width =3)
###################################################
par(mar=c(6,3,3,2)) 
a1<-  CG.1001new.denovo_TSS$mean.6909  [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]]
a2<-     CG.1001new.denovo_TSS$mean.6909 [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]]
a3<-      CG.1001new.denovo_TSS$mean.6909 [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]]
a4<-     CG.1001new.denovo_TSS$mean.6909 [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]]
a5<-  CG.1001new.denovo_TSS$mean.6909 [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a6<-     CG.1001new.denovo_TSS$mean.6909 [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a7<-     CG.1001new.denovo_TSS$mean.6909 [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a8<-  CG.1001new.denovo_TSS$mean.6909 [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<0.5& denovo2021.TPMs.genes.1001Gnew$gene %in% TE_genes.loci$gene]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,
        col=c("#486EB4","#8aa7de","#90C473","#d4edc5","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="CG meth. level in the locus")
mtext("CG", side=3, line=-1, at=1,cex=0.9)
###################################################
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
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(a5,a6)
b<-wilcox.test(a5,a6)
c<-wilcox.test(a5,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
#################
dev.off()


#boxplot CG promoter 1001GNEW data
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Suppl_Fig3_boxplot_CG_expressed_silent.1001GNEWdata.6909exp1.pdf",height = 3,width =3)
###################################################
par(mar=c(6,3,3,2)) 
a1<-  CG.1001new.denovo_TSS$r14.Exp1.6909  [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]]
a2<-     CG.1001new.denovo_TSS$r14.Exp1.6909 [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]]
a3<-      CG.1001new.denovo_TSS$r14.Exp1.6909 [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]]
a4<-     CG.1001new.denovo_TSS$r14.Exp1.6909 [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]]
a5<-  CG.1001new.denovo_TSS$r14.Exp1.6909 [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a6<-     CG.1001new.denovo_TSS$r14.Exp1.6909 [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a7<-     CG.1001new.denovo_TSS$r14.Exp1.6909 [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a8<-  CG.1001new.denovo_TSS$r14.Exp1.6909 [CG.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5& denovo2021.TPMs.genes.1001Gnew$gene %in% TE_genes.loci$gene]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,
        col=c("#486EB4","#8aa7de","#90C473","#d4edc5","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="CG meth. level in the locus")
mtext("CG", side=3, line=-1, at=1,cex=0.9)
###################################################
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
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(a5,a6)
b<-wilcox.test(a5,a6)
c<-wilcox.test(a5,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
#################
dev.off()




#boxplot CHH gene body  1001GNEW data
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Suppl_Fig3_boxplot_CHH_expressed_silent.1001GNEWdata.6909mean.pdf",height = 3,width =3)
###################################################
par(mar=c(6,3,3,2)) 
a1<-  CHH.1001new.denovo$mean.6909  [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]]
a2<-     CHH.1001new.denovo$mean.6909 [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]]
a3<-      CHH.1001new.denovo$mean.6909 [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]]
a4<-     CHH.1001new.denovo$mean.6909 [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]]
a5<-  CHH.1001new.denovo$mean.6909 [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a6<-     CHH.1001new.denovo$mean.6909 [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a7<-     CHH.1001new.denovo$mean.6909 [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a8<-  CHH.1001new.denovo$mean.6909 [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<0.5& denovo2021.TPMs.genes.1001Gnew$gene %in% TE_genes.loci$gene]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,
        col=c("#486EB4","#8aa7de","#90C473","#d4edc5","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="CHH meth. level in the locus")
mtext("CHH", side=3, line=-1, at=1,cex=0.9)
###################################################
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
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.1)
a<-wilcox.test(a5,a6)
b<-wilcox.test(a5,a6)
c<-wilcox.test(a5,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.1)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.1)
#################
dev.off()

#boxplot CHH gene body  1001GNEW data
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Suppl_Fig3_boxplot_CHH_expressed_silent.1001GNEWdata.6909exp1.pdf",height = 3,width =3)
###################################################
par(mar=c(6,3,3,2)) 
a1<-  CHH.1001new.denovo$r14.Exp1.6909  [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]]
a2<-     CHH.1001new.denovo$r14.Exp1.6909 [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]]
a3<-      CHH.1001new.denovo$r14.Exp1.6909 [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]]
a4<-     CHH.1001new.denovo$r14.Exp1.6909 [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]]
a5<-  CHH.1001new.denovo$r14.Exp1.6909 [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a6<-     CHH.1001new.denovo$r14.Exp1.6909 [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a7<-     CHH.1001new.denovo$r14.Exp1.6909 [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a8<-  CHH.1001new.denovo$r14.Exp1.6909 [CHH.1001new.denovo$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5& denovo2021.TPMs.genes.1001Gnew$gene %in% TE_genes.loci$gene]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,
        col=c("#486EB4","#8aa7de","#90C473","#d4edc5","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="CHH meth. level in the locus")
mtext('CHH', side=3, line=-1, at=1,cex=0.9)
###################################################
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
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.1)
a<-wilcox.test(a5,a6)
b<-wilcox.test(a5,a6)
c<-wilcox.test(a5,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.1)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.1)
#################
dev.off()



#boxplot CHH gene body  1001GNEW data
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Suppl_Fig3_boxplot_CHH_expressed_silent.1001GNEWdata.6909exp1.PROMOTER.pdf",height = 3,width =3)
###################################################
par(mar=c(6,3,3,2)) 
a1<-  CHH.1001new.denovo_TSS$r14.Exp1.6909  [CHH.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]]
a2<-     CHH.1001new.denovo_TSS$r14.Exp1.6909 [CHH.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]]
a3<-      CHH.1001new.denovo_TSS$r14.Exp1.6909 [CHH.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]]
a4<-     CHH.1001new.denovo_TSS$r14.Exp1.6909 [CHH.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]]
a5<-  CHH.1001new.denovo_TSS$r14.Exp1.6909 [CHH.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a6<-     CHH.1001new.denovo_TSS$r14.Exp1.6909 [CHH.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a7<-     CHH.1001new.denovo_TSS$r14.Exp1.6909 [CHH.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]]
a8<-  CHH.1001new.denovo_TSS$r14.Exp1.6909 [CHH.1001new.denovo_TSS$transcript %in% denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$r14.Exp1.6909<0.5& denovo2021.TPMs.genes.1001Gnew$gene %in% TE_genes.loci$gene]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,
        col=c("#486EB4","#8aa7de","#90C473","#d4edc5","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="CHH meth. level in the locus")
mtext('CHH', side=3, line=-1, at=1,cex=0.9)
###################################################
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
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.1)
a<-wilcox.test(a5,a6)
b<-wilcox.test(a5,a6)
c<-wilcox.test(a5,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.1)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.1)
#################
dev.off()







#######################################################################
#######################################################################
#CHIP-SEQ
#######################################################################
#######################################################################
# chip-seq levels on different genes 

#K9 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K9_6909.pdf",height = 3.5,width = 2.5)
#################################################
pc<-chip.denovo.log2$K9.6909 [chip.denovo.log2$gene %in% denovoPC.loci$gene]
as<-chip.denovo.log2$K9.6909 [chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo.log2$K9.6909 [chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo.log2$K9.6909 [chip.denovo.log2$gene %in% TE_genes.loci$gene]

par(mar=c(8,4,3,2)) 
boxplot( pc,as,linc,te    ,   col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H3K9me2', side=3, line=-1, at=2,cex=0.9)
a<-wilcox.test(sample(pc,2000),sampleas)
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()

# K27 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K27_6909.pdf",height = 3.5,width = 2.5)
#################################################
pc<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% denovoPC.loci$gene]
as<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% TE_genes.loci$gene]

par(mar=c(8,4,3,2)) 
boxplot(pc,as,linc,te    , 
  col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H3K27me3', side=3, line=-1, at=2,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()

#K36
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K36_6909.pdf",height = 3.5,width = 2.5)
#################################################
par(mar=c(8,4,3,2)) 
pc<-chip.denovo.log2$K36.6909 [chip.denovo.log2$gene %in% denovoPC.loci$gene]
as<-chip.denovo.log2$K36.6909 [chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo.log2$K36.6909 [chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo.log2$K36.6909 [chip.denovo.log2$gene %in% TE_genes.loci$gene]

boxplot(pc,as,linc,te    , 
  col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H3K36me3', side=3, line=-1, at=4,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()

#K4
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K4_6909.pdf",height = 4,width = 2.5)
#################################################
par(mar=c(8,4,3,2)) 
pc<-chip.denovo.log2$K4.6909 [chip.denovo.log2$gene %in% denovoPC.loci$gene]
as<-chip.denovo.log2$K4.6909 [chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo.log2$K4.6909 [chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo.log2$K4.6909 [chip.denovo.log2$gene %in% TE_genes.loci$gene]

boxplot(pc,as,linc,te    , 
  col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H3K4me3', side=3, line=-1, at=2,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()

#H1
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_H1_6909.pdf",height = 3.5,width = 2.5)
#################################################
par(mar=c(8,4,3,2)) 
pc<-chip.denovo.log2$H1.6909 [chip.denovo.log2$gene %in% denovoPC.loci$gene]
as<-chip.denovo.log2$H1.6909 [chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo.log2$H1.6909 [chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo.log2$H1.6909 [chip.denovo.log2$gene %in% TE_genes.loci$gene]
boxplot(    pc,as,linc,te    ,  col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H1', side=3, line=-1, at=1,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()


#######################
#histone marks on promoters (TSS +/-200bp)
##################################

#K9 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K9_6909_TSS.pdf",height = 3.5,width = 2.5)
#################################################
pc<-chip.denovo_TSS.log2$K9.6909 [chip.denovo_TSS.log2$gene %in% denovoPC.loci$gene]
as<-chip.denovo_TSS.log2$K9.6909 [chip.denovo_TSS.log2$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo_TSS.log2$K9.6909 [chip.denovo_TSS.log2$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo_TSS.log2$K9.6909 [chip.denovo_TSS.log2$gene %in% TE_genes.loci$gene]

par(mar=c(8,4,3,2)) 
boxplot( pc,as,linc,te    ,   col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H3K9me2 TSS', side=3, line=-1, at=2,cex=0.9)

#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()

# K27 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K27_6909_TSS.pdf",height = 3.5,width = 2.5)
#################################################
pc<-chip.denovo_TSS.log2$K27.6909 [chip.denovo_TSS.log2$gene %in% denovoPC.loci$gene]
as<-chip.denovo_TSS.log2$K27.6909 [chip.denovo_TSS.log2$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo_TSS.log2$K27.6909 [chip.denovo_TSS.log2$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo_TSS.log2$K27.6909 [chip.denovo_TSS.log2$gene %in% TE_genes.loci$gene]

par(mar=c(8,4,3,2)) 
boxplot(pc,as,linc,te    , 
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H3K27me3 TSS', side=3, line=-1, at=2,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()

#K36
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K36_6909_TSS.pdf",height = 3.5,width = 2.5)
#################################################
par(mar=c(8,4,3,2)) 
pc<-chip.denovo_TSS.log2$K36.6909 [chip.denovo_TSS.log2$gene %in% denovoPC.loci$gene]
as<-chip.denovo_TSS.log2$K36.6909 [chip.denovo_TSS.log2$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo_TSS.log2$K36.6909 [chip.denovo_TSS.log2$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo_TSS.log2$K36.6909 [chip.denovo_TSS.log2$gene %in% TE_genes.loci$gene]

boxplot(pc,as,linc,te    , 
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H3K36me3 TSS', side=3, line=-1, at=4,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()

#K4
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K4_6909_TSS.pdf",height = 3.5,width = 2.5)
#################################################
par(mar=c(8,4,3,2)) 
pc<-chip.denovo_TSS.log2$K4.6909 [chip.denovo_TSS.log2$gene %in% denovoPC.loci$gene]
as<-chip.denovo_TSS.log2$K4.6909 [chip.denovo_TSS.log2$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo_TSS.log2$K4.6909 [chip.denovo_TSS.log2$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo_TSS.log2$K4.6909 [chip.denovo_TSS.log2$gene %in% TE_genes.loci$gene]

boxplot(pc,as,linc,te    , 
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H3K4me3 TSS', side=3, line=-1, at=2,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()


#H1
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_H1_6909_TSS.pdf",height = 3.5,width = 2.5)
#################################################
par(mar=c(8,4,3,2)) 
pc<-chip.denovo_TSS.log2$H1.6909 [chip.denovo_TSS.log2$gene %in% denovoPC.loci$gene]
as<-chip.denovo_TSS.log2$H1.6909 [chip.denovo_TSS.log2$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo_TSS.log2$H1.6909 [chip.denovo_TSS.log2$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo_TSS.log2$H1.6909 [chip.denovo_TSS.log2$gene %in% TE_genes.loci$gene]
boxplot(    pc,as,linc,te    ,  col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H1 TSS', side=3, line=-1, at=1,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()





######################################
## ChIP seq marks - genes off genes on#
#######################################

# H1 PC,AS, linc,TE
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_H1_expressed_silent_6909_allgenes.pdf",height = 3,width =3)
#################################################
par(mar=c(5,3,2,2)) 
a1<-chip.denovo.log2$H1.6909 [chip.denovo.log2$gene %in% pc_sasa$gene[pc_sasa$mean.6909>0.5]]
a2<-chip.denovo.log2$H1.6909 [chip.denovo.log2$gene %in% pc_sasa$gene[pc_sasa$mean.6909<0.5]]
a3<-chip.denovo.log2$H1.6909 [chip.denovo.log2$gene %in% as_sasa$gene[as_sasa$mean.6909>0.5]]
a4<-chip.denovo.log2$H1.6909 [chip.denovo.log2$gene %in% as_sasa$gene[as_sasa$mean.6909<0.5]]
a5<-chip.denovo.log2$H1.6909 [chip.denovo.log2$gene %in% linc_sasa$gene[linc_sasa$mean.6909>0.5]]
a6<-chip.denovo.log2$H1.6909 [chip.denovo.log2$gene %in% linc_sasa$gene[linc_sasa$mean.6909<0.5]]
a7<-chip.denovo.log2$H1.6909 [chip.denovo.log2$gene %in% te_sasa$gene[te_sasa$mean.6909>0.5]]
a8<-chip.denovo.log2$H1.6909 [chip.denovo.log2$gene %in% te_sasa$gene[te_sasa$mean.6909<0.5]]
boxplot(a1,a2,a3,a4,a5,a6,a7,a8,
        col=c("#486EB4","#bbcded","#90C473","#d0edc0","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H1', side=3, line=-1, at=1,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(a1,7602),sample(a2,7602))
b<-wilcox.test(sample(a1,7602),sample(a2,7602))
c<-wilcox.test(sample(a1,7602),sample(a2,7602))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(a3,238),sample(a4,238))
b<-wilcox.test(sample(a3,238),sample(a4,238))
c<-wilcox.test(sample(a3,238),sample(a4,238))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
a<-wilcox.test(sample(a5,61),sample(a6,61))
b<-wilcox.test(sample(a5,61),sample(a6,61))
c<-wilcox.test(sample(a5,61),sample(a6,61))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.4)
a<-wilcox.test(sample(a7,46),sample(a8,46))
b<-wilcox.test(sample(a7,46),sample(a8,46))
c<-wilcox.test(sample(a7,46),sample(a8,46))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.4)
#################
dev.off()

# H1 PC,AS, linc,TE
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_H1_expressed_silent_6909_allgenes_TSS.pdf",height = 3,width =3)
#################################################
par(mar=c(5,3,2,2)) 
a1<-chip.denovo_TSS.log2$H1.6909 [chip.denovo_TSS.log2$gene %in% pc_sasa$gene[pc_sasa$mean.6909>0.5]]
a2<-chip.denovo_TSS.log2$H1.6909 [chip.denovo_TSS.log2$gene %in% pc_sasa$gene[pc_sasa$mean.6909<0.5]]
a3<-chip.denovo_TSS.log2$H1.6909 [chip.denovo_TSS.log2$gene %in% as_sasa$gene[as_sasa$mean.6909>0.5]]
a4<-chip.denovo_TSS.log2$H1.6909 [chip.denovo_TSS.log2$gene %in% as_sasa$gene[as_sasa$mean.6909<0.5]]
a5<-chip.denovo_TSS.log2$H1.6909 [chip.denovo_TSS.log2$gene %in% linc_sasa$gene[linc_sasa$mean.6909>0.5]]
a6<-chip.denovo_TSS.log2$H1.6909 [chip.denovo_TSS.log2$gene %in% linc_sasa$gene[linc_sasa$mean.6909<0.5]]
a7<-chip.denovo_TSS.log2$H1.6909 [chip.denovo_TSS.log2$gene %in% te_sasa$gene[te_sasa$mean.6909>0.5]]
a8<-chip.denovo_TSS.log2$H1.6909 [chip.denovo_TSS.log2$gene %in% te_sasa$gene[te_sasa$mean.6909<0.5]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,
        col=c("#486EB4","#bbcded","#90C473","#d0edc0","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H1 TSS', side=3, line=-1, at=1,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(a1,7602),sample(a2,7602))
b<-wilcox.test(sample(a1,7602),sample(a2,7602))
c<-wilcox.test(sample(a1,7602),sample(a2,7602))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(a3,238),sample(a4,238))
b<-wilcox.test(sample(a3,238),sample(a4,238))
c<-wilcox.test(sample(a3,238),sample(a4,238))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
a<-wilcox.test(sample(a5,61),sample(a6,61))
b<-wilcox.test(sample(a5,61),sample(a6,61))
c<-wilcox.test(sample(a5,61),sample(a6,61))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.4)
a<-wilcox.test(sample(a7,46),sample(a8,46))
b<-wilcox.test(sample(a7,46),sample(a8,46))
c<-wilcox.test(sample(a7,46),sample(a8,46))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.4)
#################
dev.off()



# H3K27me3 gene body PC,AS, linc,TE
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K27_expressed_silent_6909_allgenes.pdf",height = 3,width =3)
#################################################
par(mar=c(5,3,2,2)) 
a1<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% pc_sasa$gene[pc_sasa$mean.6909>0.5]]
a2<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% pc_sasa$gene[pc_sasa$mean.6909<0.5]]

a3<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% as_sasa$gene[as_sasa$mean.6909>0.5]]
a4<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% as_sasa$gene[as_sasa$mean.6909<0.5]]
a5<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% linc_sasa$gene[linc_sasa$mean.6909>0.5]]
a6<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% linc_sasa$gene[linc_sasa$mean.6909<0.5]]
a7<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% te_sasa$gene[te_sasa$mean.6909>0.5]]
a8<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% te_sasa$gene[te_sasa$mean.6909<0.5]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,    col=c("#486EB4","#bbcded","#90C473","#d0edc0","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H3K27me3', side=3, line=-1, at=1,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(a1,7602),sample(a2,7602))
b<-wilcox.test(sample(a1,7602),sample(a2,7602))
c<-wilcox.test(sample(a1,7602),sample(a2,7602))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(a3,238),sample(a4,238))
b<-wilcox.test(sample(a3,238),sample(a4,238))
c<-wilcox.test(sample(a3,238),sample(a4,238))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
a<-wilcox.test(sample(a5,61),sample(a6,61))
b<-wilcox.test(sample(a5,61),sample(a6,61))
c<-wilcox.test(sample(a5,61),sample(a6,61))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.4)
a<-wilcox.test(sample(a7,46),sample(a8,46))
b<-wilcox.test(sample(a7,46),sample(a8,46))
c<-wilcox.test(sample(a7,46),sample(a8,46))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.4)
#################
dev.off()


# H3K27me3 promoter PC,AS, linc,TE
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K27_expressed_silent_6909_allgenes_TSS.pdf",height = 3,width =3)
#################################################
par(mar=c(5,3,2,2)) 
a1<-chip.denovo_TSS.log2$K27.6909 [chip.denovo_TSS.log2$gene %in% pc_sasa$gene[pc_sasa$mean.6909>0.5]]
a2<-chip.denovo_TSS.log2$K27.6909 [chip.denovo_TSS.log2$gene %in% pc_sasa$gene[pc_sasa$mean.6909<0.5]]
a3<-chip.denovo_TSS.log2$K27.6909 [chip.denovo_TSS.log2$gene %in% as_sasa$gene[as_sasa$mean.6909>0.5]]
a4<-chip.denovo_TSS.log2$K27.6909 [chip.denovo_TSS.log2$gene %in% as_sasa$gene[as_sasa$mean.6909<0.5]]
a5<-chip.denovo_TSS.log2$K27.6909 [chip.denovo_TSS.log2$gene %in% linc_sasa$gene[linc_sasa$mean.6909>0.5]]
a6<-chip.denovo_TSS.log2$K27.6909 [chip.denovo_TSS.log2$gene %in% linc_sasa$gene[linc_sasa$mean.6909<0.5]]
a7<-chip.denovo_TSS.log2$K27.6909 [chip.denovo_TSS.log2$gene %in% te_sasa$gene[te_sasa$mean.6909>0.5]]
a8<-chip.denovo_TSS.log2$K27.6909 [chip.denovo_TSS.log2$gene %in% te_sasa$gene[te_sasa$mean.6909<0.5]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,
        col=c("#486EB4","#bbcded","#90C473","#d0edc0","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H3K27me3_TSS', side=3, line=-1, at=1,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(a1,7602),sample(a2,7602))
b<-wilcox.test(sample(a1,7602),sample(a2,7602))
c<-wilcox.test(sample(a1,7602),sample(a2,7602))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(a3,238),sample(a4,238))
b<-wilcox.test(sample(a3,238),sample(a4,238))
c<-wilcox.test(sample(a3,238),sample(a4,238))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
a<-wilcox.test(sample(a5,61),sample(a6,61))
b<-wilcox.test(sample(a5,61),sample(a6,61))
c<-wilcox.test(sample(a5,61),sample(a6,61))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.4)
a<-wilcox.test(sample(a7,46),sample(a8,46))
b<-wilcox.test(sample(a7,46),sample(a8,46))
c<-wilcox.test(sample(a7,46),sample(a8,46))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.4)
#################
dev.off()



# H3K9me2 genebody PC,AS, linc,TE
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K9_expressed_silent_6909_allgenes.pdf",height = 3,width =3)
#################################################
par(mar=c(5,3,2,2)) 
a1<-chip.denovo.log2$K9.6909 [chip.denovo.log2$gene %in% pc_sasa$gene[pc_sasa$mean.6909>0.5]]
a2<-chip.denovo.log2$K9.6909 [chip.denovo.log2$gene %in% pc_sasa$gene[pc_sasa$mean.6909<0.5]]

a3<-chip.denovo.log2$K9.6909 [chip.denovo.log2$gene %in% as_sasa$gene[as_sasa$mean.6909>0.5]]
a4<-chip.denovo.log2$K9.6909 [chip.denovo.log2$gene %in% as_sasa$gene[as_sasa$mean.6909<0.5]]
a5<-chip.denovo.log2$K9.6909 [chip.denovo.log2$gene %in% linc_sasa$gene[linc_sasa$mean.6909>0.5]]
a6<-chip.denovo.log2$K9.6909 [chip.denovo.log2$gene %in% linc_sasa$gene[linc_sasa$mean.6909<0.5]]
a7<-chip.denovo.log2$K9.6909 [chip.denovo.log2$gene %in% te_sasa$gene[te_sasa$mean.6909>0.5]]
a8<-chip.denovo.log2$K9.6909 [chip.denovo.log2$gene %in% te_sasa$gene[te_sasa$mean.6909<0.5]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,  
        col=c("#486EB4","#bbcded","#90C473","#d0edc0","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H3K9me2', side=3, line=-1, at=1,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(a1,7602),sample(a2,7602))
b<-wilcox.test(sample(a1,7602),sample(a2,7602))
c<-wilcox.test(sample(a1,7602),sample(a2,7602))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(a3,238),sample(a4,238))
b<-wilcox.test(sample(a3,238),sample(a4,238))
c<-wilcox.test(sample(a3,238),sample(a4,238))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
a<-wilcox.test(sample(a5,61),sample(a6,61))
b<-wilcox.test(sample(a5,61),sample(a6,61))
c<-wilcox.test(sample(a5,61),sample(a6,61))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.4)
a<-wilcox.test(sample(a7,46),sample(a8,46))
b<-wilcox.test(sample(a7,46),sample(a8,46))
c<-wilcox.test(sample(a7,46),sample(a8,46))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.4)
#################
dev.off()


# H3K9me2 promoter PC,AS, linc,TE
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K9_expressed_silent_6909_allgenes_TSS.pdf",height = 3,width =3)
#################################################
par(mar=c(5,3,2,2)) 
a1<-chip.denovo_TSS.log2$K9.6909 [chip.denovo_TSS.log2$gene %in% pc_sasa$gene[pc_sasa$mean.6909>0.5]]
a2<-chip.denovo_TSS.log2$K9.6909 [chip.denovo_TSS.log2$gene %in% pc_sasa$gene[pc_sasa$mean.6909<0.5]]
a3<-chip.denovo_TSS.log2$K9.6909 [chip.denovo_TSS.log2$gene %in% as_sasa$gene[as_sasa$mean.6909>0.5]]
a4<-chip.denovo_TSS.log2$K9.6909 [chip.denovo_TSS.log2$gene %in% as_sasa$gene[as_sasa$mean.6909<0.5]]
a5<-chip.denovo_TSS.log2$K9.6909 [chip.denovo_TSS.log2$gene %in% linc_sasa$gene[linc_sasa$mean.6909>0.5]]
a6<-chip.denovo_TSS.log2$K9.6909 [chip.denovo_TSS.log2$gene %in% linc_sasa$gene[linc_sasa$mean.6909<0.5]]
a7<-chip.denovo_TSS.log2$K9.6909 [chip.denovo_TSS.log2$gene %in% te_sasa$gene[te_sasa$mean.6909>0.5]]
a8<-chip.denovo_TSS.log2$K9.6909 [chip.denovo_TSS.log2$gene %in% te_sasa$gene[te_sasa$mean.6909<0.5]]

boxplot(a1,a2,a3,a4,a5,a6,a7,a8,       col=c("#486EB4","#bbcded","#90C473","#d0edc0","#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("PC ON","PC OFF","AS ON","AS OFF","linc ON","linc OFF","TE ON","TE OFF"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H3K9me2 TSS', side=3, line=-1, at=1,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(a1,7602),sample(a2,7602))
b<-wilcox.test(sample(a1,7602),sample(a2,7602))
c<-wilcox.test(sample(a1,7602),sample(a2,7602))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(a3,238),sample(a4,238))
b<-wilcox.test(sample(a3,238),sample(a4,238))
c<-wilcox.test(sample(a3,238),sample(a4,238))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
a<-wilcox.test(sample(a5,61),sample(a6,61))
b<-wilcox.test(sample(a5,61),sample(a6,61))
c<-wilcox.test(sample(a5,61),sample(a6,61))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.4)
a<-wilcox.test(sample(a7,46),sample(a8,46))
b<-wilcox.test(sample(a7,46),sample(a8,46))
c<-wilcox.test(sample(a7,46),sample(a8,46))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.4)
#################
dev.off()






################################################
#### chromosomal position of chipseq marks 
################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/chipseq_H1_K9_K27_alongchr.pdf",height = 3,width =8)
par(mar=c(3,3,3,0) + 0.1)
par(mgp=c(2,1,0))
par(mfrow=c(1,4))
plot(chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% denovoPC.loci$gene]-5,pch=20,cex=0.4,col=alpha("red",alpha=0.4),ylim=c(-7,7),ylab="log2(ChIP/Input)",xlab="number of the gene (Chr1->5)",main="PC genes",yaxt = "n")
axis(2,at = -7:7, labels =c(-2,-1,0,1,2,3,-1,0,1,2,-2,-1,0,1,2) ,las=2)
points(chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% denovoPC.loci$gene]+5,pch=20,cex=0.4,col=alpha("aquamarine4",alpha=0.4))
points(chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% denovoPC.loci$gene]+1,pch=20,cex=0.4,col=alpha("darkmagenta",alpha=0.4))

plot(chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene]-5,pch=20,cex=0.4,col=alpha("red",alpha=0.4),ylim=c(-7,7),ylab="log2(ChIP/Input)",xlab="number of the gene (Chr1->5)",main="antisense lncRNAs",yaxt = "n")
axis(2,at = -7:7, labels =c(-2,-1,0,1,2,3,-1,0,1,2,-2,-1,0,1,2) ,las=2)
points(chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene]+5,pch=20,cex=0.4,col=alpha("aquamarine4",alpha=0.4))
points(chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene]+1,pch=20,cex=0.4,col=alpha("darkmagenta",alpha=0.4))

plot(chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene]-5,pch=20,cex=0.4,col=alpha("red",alpha=0.4),ylim=c(-7,7),ylab="log2(ChIP/Input)",xlab="number of the gene (Chr1->5)",main="lincRNAs",yaxt = "n")
axis(2,at = -7:7, labels =c(-2,-1,0,1,2,3,-1,0,1,2,-2,-1,0,1,2),las=2 )
points(chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene]+5,pch=20,cex=0.4,col=alpha("aquamarine4",alpha=0.4))
points(chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene]+1,pch=20,cex=0.4,col=alpha("darkmagenta",alpha=0.4))

plot(chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene]-5,pch=20,cex=0.4,col=alpha("red",alpha=0.4),ylim=c(-7,7),ylab="log2(ChIP/Input)",xlab="number of the gene (Chr1->5)",main="TE genes",yaxt = "n")
axis(2,at = -7:7, labels =c(-2,-1,0,1,2,3,-1,0,1,2,-2,-1,0,1,2),las=2 )
points(chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene]+5,pch=20,cex=0.4,col=alpha("aquamarine4",alpha=0.4))
points(chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene]+1,pch=20,cex=0.4,col=alpha("darkmagenta",alpha=0.4))

dev.off()


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/chipseq_H1_K9_K27_alongchr.Chrposition.pdf",height = 3,width =8)
par(mar=c(3,3,3,0) + 0.1)
par(mgp=c(2,1,0))
par(mfrow=c(1,4))

#PC   
a<-merge(denovoPC.loci,chip.denovo.log2[,c("gene","H1.6909","K9.6909","K27.6909")])
a$start[a$chr=="Chr1"]<-a$start[a$chr=="Chr1"]
a$end[a$chr=="Chr1"]<-a$end[a$chr=="Chr1"]
a$start[a$chr=="Chr2"]<-a$start[a$chr=="Chr2"]+30427671
a$end[a$chr=="Chr2"]<-a$end[a$chr=="Chr2"]+30427671
a$start[a$chr=="Chr3"]<-a$start[a$chr=="Chr3"]+30427671+19698289
a$end[a$chr=="Chr3"]<-a$end[a$chr=="Chr3"]+30427671+19698289
a$start[a$chr=="Chr4"]<-a$start[a$chr=="Chr4"]+30427671+19698289+23459830
a$end[a$chr=="Chr4"]<-a$end[a$chr=="Chr4"]+30427671+19698289+23459830
a$start[a$chr=="Chr5"]<-a$start[a$chr=="Chr5"]+30427671+19698289+23459830+18585056
a$end[a$chr=="Chr5"]<-a$end[a$chr=="Chr5"]+30427671+19698289+23459830+18585056

plot(a$start,a$K27.6909-5,pch=20,cex=0.4,col=alpha("red",alpha=0.4),ylim=c(-7,7),ylab="log2(ChIP/Input)",xlab="gene start position",main="PC genes",yaxt = "n",xaxt = "n")
axis(1,at = c(30427671,30427671+19698289,30427671+19698289+23459830,30427671+19698289+23459830+18585056,30427671+19698289+23459830+18585056+26975502 ), labels =c("","","","","") ,las=2)
axis(1,at = c(30427671/2,30427671+19698289/2,30427671+19698289+23459830/2,30427671+19698289+23459830+18585056/2,30427671+19698289+23459830+18585056+26975502/2 ), labels =c("Chr1","Chr2","Chr3","Chr4","Chr5") ,tick = F,las=1,mgp=c(0,0,0),cex.axis=0.8)
axis(2,at = -7:7, labels =c(-2,-1,0,1,2,3,-1,0,1,2,-2,-1,0,1,2) ,las=2)
points(a$start,a$H1.6909+5,pch=20,cex=0.4,col=alpha("aquamarine4",alpha=0.4))
points(a$start,a$K9.6909+1,pch=20,cex=0.4,col=alpha("darkmagenta",alpha=0.4))
#add (approximate) centromere center positions
abline(v=15000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+4700000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+13000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+3800000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+18585056+11900000,col=alpha("black",alpha=0.5), lty=2, lwd=1)

#AS   
a<-merge(lncRNAs.antisense.loci,chip.denovo.log2[,c("gene","H1.6909","K9.6909","K27.6909")])
a$start[a$chr=="Chr1"]<-a$start[a$chr=="Chr1"]
a$end[a$chr=="Chr1"]<-a$end[a$chr=="Chr1"]
a$start[a$chr=="Chr2"]<-a$start[a$chr=="Chr2"]+30427671
a$end[a$chr=="Chr2"]<-a$end[a$chr=="Chr2"]+30427671
a$start[a$chr=="Chr3"]<-a$start[a$chr=="Chr3"]+30427671+19698289
a$end[a$chr=="Chr3"]<-a$end[a$chr=="Chr3"]+30427671+19698289
a$start[a$chr=="Chr4"]<-a$start[a$chr=="Chr4"]+30427671+19698289+23459830
a$end[a$chr=="Chr4"]<-a$end[a$chr=="Chr4"]+30427671+19698289+23459830
a$start[a$chr=="Chr5"]<-a$start[a$chr=="Chr5"]+30427671+19698289+23459830+18585056
a$end[a$chr=="Chr5"]<-a$end[a$chr=="Chr5"]+30427671+19698289+23459830+18585056
plot(a$start,a$K27.6909-5,pch=20,cex=0.4,col=alpha("red",alpha=0.4),ylim=c(-7,7),ylab="log2(ChIP/Input)",xlab="gene start position",main="AS lncRNAs",yaxt = "n",xaxt = "n")
axis(1,at = c(30427671,30427671+19698289,30427671+19698289+23459830,30427671+19698289+23459830+18585056,30427671+19698289+23459830+18585056+26975502 ), labels =c("","","","","") ,las=2)
axis(1,at = c(30427671/2,30427671+19698289/2,30427671+19698289+23459830/2,30427671+19698289+23459830+18585056/2,30427671+19698289+23459830+18585056+26975502/2 ), labels =c("Chr1","Chr2","Chr3","Chr4","Chr5") ,tick = F,las=1,mgp=c(0,0,0),cex.axis=0.8)
axis(2,at = -7:7, labels =c(-2,-1,0,1,2,3,-1,0,1,2,-2,-1,0,1,2) ,las=2)
points(a$start,a$H1.6909+5,pch=20,cex=0.4,col=alpha("aquamarine4",alpha=0.4))
points(a$start,a$K9.6909+1,pch=20,cex=0.4,col=alpha("darkmagenta",alpha=0.4))
#add (approximate) centromere center positions
abline(v=15000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+4700000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+13000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+3800000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+18585056+11900000,col=alpha("black",alpha=0.5), lty=2, lwd=1)

#linc  
a<-merge(lncRNAs.intergenic.loci,chip.denovo.log2[,c("gene","H1.6909","K9.6909","K27.6909")])
a$start[a$chr=="Chr1"]<-a$start[a$chr=="Chr1"]
a$end[a$chr=="Chr1"]<-a$end[a$chr=="Chr1"]
a$start[a$chr=="Chr2"]<-a$start[a$chr=="Chr2"]+30427671
a$end[a$chr=="Chr2"]<-a$end[a$chr=="Chr2"]+30427671
a$start[a$chr=="Chr3"]<-a$start[a$chr=="Chr3"]+30427671+19698289
a$end[a$chr=="Chr3"]<-a$end[a$chr=="Chr3"]+30427671+19698289
a$start[a$chr=="Chr4"]<-a$start[a$chr=="Chr4"]+30427671+19698289+23459830
a$end[a$chr=="Chr4"]<-a$end[a$chr=="Chr4"]+30427671+19698289+23459830
a$start[a$chr=="Chr5"]<-a$start[a$chr=="Chr5"]+30427671+19698289+23459830+18585056
a$end[a$chr=="Chr5"]<-a$end[a$chr=="Chr5"]+30427671+19698289+23459830+18585056
plot(a$start,a$K27.6909-5,pch=20,cex=0.4,col=alpha("red",alpha=0.4),ylim=c(-7,7),ylab="log2(ChIP/Input)",xlab="gene start position",main="lincRNAs",yaxt = "n",xaxt = "n")
axis(1,at = c(30427671,30427671+19698289,30427671+19698289+23459830,30427671+19698289+23459830+18585056,30427671+19698289+23459830+18585056+26975502 ), labels =c("","","","","") ,las=2)
axis(1,at = c(30427671/2,30427671+19698289/2,30427671+19698289+23459830/2,30427671+19698289+23459830+18585056/2,30427671+19698289+23459830+18585056+26975502/2 ), labels =c("Chr1","Chr2","Chr3","Chr4","Chr5") ,tick = F,las=1,mgp=c(0,0,0),cex.axis=0.8)
axis(2,at = -7:7, labels =c(-2,-1,0,1,2,3,-1,0,1,2,-2,-1,0,1,2) ,las=2)
points(a$start,a$H1.6909+5,pch=20,cex=0.4,col=alpha("aquamarine4",alpha=0.4))
points(a$start,a$K9.6909+1,pch=20,cex=0.4,col=alpha("darkmagenta",alpha=0.4))
#add (approximate) centromere center positions
abline(v=15000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+4700000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+13000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+3800000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+18585056+11900000,col=alpha("black",alpha=0.5), lty=2, lwd=1)


#TE genes   
a<-merge(TE_genes.loci,chip.denovo.log2[,c("gene","H1.6909","K9.6909","K27.6909")])
a$start[a$chr=="Chr1"]<-a$start[a$chr=="Chr1"]
a$end[a$chr=="Chr1"]<-a$end[a$chr=="Chr1"]
a$start[a$chr=="Chr2"]<-a$start[a$chr=="Chr2"]+30427671
a$end[a$chr=="Chr2"]<-a$end[a$chr=="Chr2"]+30427671
a$start[a$chr=="Chr3"]<-a$start[a$chr=="Chr3"]+30427671+19698289
a$end[a$chr=="Chr3"]<-a$end[a$chr=="Chr3"]+30427671+19698289
a$start[a$chr=="Chr4"]<-a$start[a$chr=="Chr4"]+30427671+19698289+23459830
a$end[a$chr=="Chr4"]<-a$end[a$chr=="Chr4"]+30427671+19698289+23459830
a$start[a$chr=="Chr5"]<-a$start[a$chr=="Chr5"]+30427671+19698289+23459830+18585056
a$end[a$chr=="Chr5"]<-a$end[a$chr=="Chr5"]+30427671+19698289+23459830+18585056
plot(a$start,a$K27.6909-5,pch=20,cex=0.4,col=alpha("red",alpha=0.4),ylim=c(-7,7),ylab="log2(ChIP/Input)",xlab="gene start position",main="TE genes",yaxt = "n",xaxt = "n")
axis(1,at = c(30427671,30427671+19698289,30427671+19698289+23459830,30427671+19698289+23459830+18585056,30427671+19698289+23459830+18585056+26975502 ), labels =c("","","","","") ,las=2)
axis(1,at = c(30427671/2,30427671+19698289/2,30427671+19698289+23459830/2,30427671+19698289+23459830+18585056/2,30427671+19698289+23459830+18585056+26975502/2 ), labels =c("Chr1","Chr2","Chr3","Chr4","Chr5") ,tick = F,las=1,mgp=c(0,0,0),cex.axis=0.8)
axis(2,at = -7:7, labels =c(-2,-1,0,1,2,3,-1,0,1,2,-2,-1,0,1,2) ,las=2)
points(a$start,a$H1.6909+5,pch=20,cex=0.4,col=alpha("aquamarine4",alpha=0.4))
points(a$start,a$K9.6909+1,pch=20,cex=0.4,col=alpha("darkmagenta",alpha=0.4))
#add (approximate) centromere center positions
abline(v=15000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+4700000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+13000000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+3800000,col=alpha("black",alpha=0.5), lty=2, lwd=1)
abline(v=30427671+19698289+23459830+18585056+11900000,col=alpha("black",alpha=0.5), lty=2, lwd=1)

dev.off()


#centro_start <- c(14364752, 3602775, 12674550, 2919690, 11668616);
#centro_end   <- c(15750321, 3735247, 13674767, 4011692, 12082583);
#############################################################################

# chipseq dependent on distance from centromere 

# H3K9me2 distant to centromere close to centromere
###########################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K9.distant.nondist.pdf",height = 3,width = 5)
###############################################

par(mar=c(6,6,2,2)) 
pc1<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere<2000000] ]
pc2<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere>2000000] ]
as1<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere<2000000]]
as2<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere>2000000]]
linc1<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<2000000]]
linc2<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>2000000]]
te1<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere<2000000]]
te2<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere>2000000]]

boxplot( pc1,as1,linc1,te1,pc2,as2,linc2,te2,
         col=c("#486EB4","#90C473","#F2AB54","#673A8E","#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE","PC","AS","linc","TE"),las=2, main="H3K9me2 Col-0",ylab="log2(ChIP/Input)", notch = T,outline = F)

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
text(b,x=1.5,y=0.14)
a<-wilcox.test(linc1,as1)
b<-wilcox.test(linc1,as1)
c<-wilcox.test(linc1,as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.14)
a<-wilcox.test(linc1,te1)
b<-wilcox.test(linc1,te1)
c<-wilcox.test(linc1,te1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.14)
a<-wilcox.test(sample(pc2,1612),sample(as2,1612))
b<-wilcox.test(sample(pc2,1612),sample(as2,1612))
c<-wilcox.test(sample(pc2,1612),sample(as2,1612))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.14)
a<-wilcox.test(sample(linc2,1612),sample(as2,1612))
b<-wilcox.test(sample(linc2,1612),sample(as2,1612))
c<-wilcox.test(sample(linc2,1612),sample(as2,1612))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=0.14)
a<-wilcox.test(sample(linc2,815),te2)
b<-wilcox.test(sample(linc2,815),te2)
c<-wilcox.test(sample(linc2,815),te2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.14)

a<-wilcox.test(pc1,sample(pc2,1012))
b<-wilcox.test(pc1,sample(pc2,1012))
c<-wilcox.test(pc1,sample(pc2,1012))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=0.2)
a<-wilcox.test(sample(as2,333),as1)
b<-wilcox.test(sample(as2,333),as1)
c<-wilcox.test(sample(as2,333),as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=0.2)
a<-wilcox.test(linc1,linc2)
b<-wilcox.test(linc1,linc2)
c<-wilcox.test(linc1,linc2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=0.2)
a<-wilcox.test(te1,te2)
b<-wilcox.test(te1,te2)
c<-wilcox.test(te1,te2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=0.2)
#################
dev.off()


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K27.distant.nondist.pdf",height = 3,width = 5)
#######################################################
par(mar=c(6,6,2,2)) 
pc1<-chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere<2000000] ]
pc2<-chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere>2000000] ]

as1<-chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere<2000000]]
as2<-chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere>2000000]]

linc1<-chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<2000000]]
linc2<-chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>2000000]]

te1<-chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere<2000000]]
te2<-chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere>2000000]]

boxplot( pc1,as1,linc1,te1,pc2,as2,linc2,te2,
         col=c("#486EB4","#90C473","#F2AB54","#673A8E","#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE","PC","AS","linc","TE"),las=2, main="H3K27me3  Col-0",ylab="log2(ChIP/Input)", notch = T,outline = F)

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
text(b,x=1.5,y=0.14)
a<-wilcox.test(linc1,as1)
b<-wilcox.test(linc1,as1)
c<-wilcox.test(linc1,as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.14)
a<-wilcox.test(linc1,te1)
b<-wilcox.test(linc1,te1)
c<-wilcox.test(linc1,te1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.14)
a<-wilcox.test(sample(pc2,1612),sample(as2,1612))
b<-wilcox.test(sample(pc2,1612),sample(as2,1612))
c<-wilcox.test(sample(pc2,1612),sample(as2,1612))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.14)
a<-wilcox.test(sample(linc2,1612),sample(as2,1612))
b<-wilcox.test(sample(linc2,1612),sample(as2,1612))
c<-wilcox.test(sample(linc2,1612),sample(as2,1612))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=0.14)
a<-wilcox.test(sample(linc2,815),te2)
b<-wilcox.test(sample(linc2,815),te2)
c<-wilcox.test(sample(linc2,815),te2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.14)

a<-wilcox.test(pc1,sample(pc2,1012))
b<-wilcox.test(pc1,sample(pc2,1012))
c<-wilcox.test(pc1,sample(pc2,1012))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=0.2)
a<-wilcox.test(sample(as2,333),as1)
b<-wilcox.test(sample(as2,333),as1)
c<-wilcox.test(sample(as2,333),as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=0.2)
a<-wilcox.test(linc1,linc2)
b<-wilcox.test(linc1,linc2)
c<-wilcox.test(linc1,linc2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=0.2)
a<-wilcox.test(te1,te2)
b<-wilcox.test(te1,te2)
c<-wilcox.test(te1,te2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=0.2)
#################
dev.off()

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_H1.distant.nondist.pdf",height = 3,width = 5)
###############################################
par(mar=c(6,6,2,2)) 
pc1<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere<2000000] ]
pc2<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% denovoPC.loci$gene[denovoPC.loci$dist_from_centromere>2000000] ]

as1<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere<2000000]]
as2<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$dist_from_centromere>2000000]]

linc1<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere<2000000]]
linc2<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$dist_from_centromere>2000000]]

te1<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere<2000000]]
te2<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% TE_genes.loci$gene[TE_genes.loci$dist_from_centromere>2000000]]

boxplot( pc1,as1,linc1,te1,pc2,as2,linc2,te2,
         col=c("#486EB4","#90C473","#F2AB54","#673A8E","#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE","PC","AS","linc","TE"),las=2, main="H1 Col-0",ylab="log2(ChIP/Input)", notch = T,outline = F)

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
text(b,x=1.5,y=0.14)
a<-wilcox.test(linc1,as1)
b<-wilcox.test(linc1,as1)
c<-wilcox.test(linc1,as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.14)
a<-wilcox.test(linc1,te1)
b<-wilcox.test(linc1,te1)
c<-wilcox.test(linc1,te1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.14)
a<-wilcox.test(sample(pc2,1612),sample(as2,1612))
b<-wilcox.test(sample(pc2,1612),sample(as2,1612))
c<-wilcox.test(sample(pc2,1612),sample(as2,1612))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.14)
a<-wilcox.test(sample(linc2,1612),sample(as2,1612))
b<-wilcox.test(sample(linc2,1612),sample(as2,1612))
c<-wilcox.test(sample(linc2,1612),sample(as2,1612))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=0.14)
a<-wilcox.test(sample(linc2,815),te2)
b<-wilcox.test(sample(linc2,815),te2)
c<-wilcox.test(sample(linc2,815),te2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.14)

a<-wilcox.test(pc1,sample(pc2,1012))
b<-wilcox.test(pc1,sample(pc2,1012))
c<-wilcox.test(pc1,sample(pc2,1012))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=0.2)
a<-wilcox.test(sample(as2,333),as1)
b<-wilcox.test(sample(as2,333),as1)
c<-wilcox.test(sample(as2,333),as1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=0.2)
a<-wilcox.test(linc1,linc2)
b<-wilcox.test(linc1,linc2)
c<-wilcox.test(linc1,linc2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=0.2)
a<-wilcox.test(te1,te2)
b<-wilcox.test(te1,te2)
c<-wilcox.test(te1,te2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=0.2)
#################
dev.off()






############################################
#########################################
#small RNA levels on different genes 
#############################################
############################################

#################
#24 nt coverage boxplot Col0 flowers - 4 gene types (this study)
###################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_sRNA24nt_6909.pdf",height = 4,width = 2.5)
###################################
par(mar=c(8,4,3,2)) 
pc<-sRNA.24nt.denovo2021.RPM$X6909 [sRNA.24nt.denovo2021.RPM$gene %in% denovoPC.loci$gene]
as<-sRNA.24nt.denovo2021.RPM$X6909 [sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.antisense.loci$gene]
linc<-sRNA.24nt.denovo2021.RPM$X6909 [sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene]
te<-sRNA.24nt.denovo2021.RPM$X6909 [sRNA.24nt.denovo2021.RPM$gene %in% TE_genes.loci$gene]

boxplot(pc,as,linc,te,col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="RPM, 24nt")
mtext('24nt sRNA', side=3, line=-1, at=2,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()


#################
#21 nt coverage boxplot Col0 flowers - 4 gene types (this study)
###################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_sRNA21_22nt_6909.pdf",height = 4,width = 2.5)
############################
par(mar=c(8,4,3,2))
pc<-sRNA.21nt.denovo2021.RPM$X6909 [sRNA.21nt.denovo2021.RPM$gene %in% denovoPC.loci$gene]
as<-sRNA.21nt.denovo2021.RPM$X6909 [sRNA.21nt.denovo2021.RPM$gene %in% lncRNAs.antisense.loci$gene]
linc<-sRNA.21nt.denovo2021.RPM$X6909 [sRNA.21nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene]
te<-sRNA.21nt.denovo2021.RPM$X6909 [sRNA.21nt.denovo2021.RPM$gene %in% TE_genes.loci$gene]

boxplot( pc,as,linc,te,col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="RPM, 21-22nt")
mtext('21-22nt sRNA', side=3, line=-1, at=2,cex=0.9)

#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.08)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.08)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.08)
#################
dev.off()


#Ranj data 
#################
#24 nt coverage boxplot Col0 early heart- 4 gene types (Ranj)
###################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Suppl_Fig3_boxplot_sRNA24nt_6909_WTeh_RANJdata.pdf",height = 4,width = 2.5)
###################################
par(mar=c(8,4,3,2)) 
pc<-sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart [sRNA.24nt.denovo2021.RPM.Ranj$gene %in% denovoPC.loci$gene]
as<-sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart [sRNA.24nt.denovo2021.RPM.Ranj$gene %in% lncRNAs.antisense.loci$gene]
linc<-sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart [sRNA.24nt.denovo2021.RPM.Ranj$gene %in% lncRNAs.intergenic.loci$gene]
te<-sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart [sRNA.24nt.denovo2021.RPM.Ranj$gene %in% TE_genes.loci$gene]

boxplot(pc,as,linc,te,col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="RPM, 24nt",main="24nt sRNA level\n early heart, Col-0")
#mtext('24nt sRNA', side=3, line=-1, at=2,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1.4)
#################
dev.off()

#################
#21-22 nt coverage boxplot Col0 early heart- 4 gene types (Ranj)
###################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Suppl_Fig3_boxplot_sRNA21nt_6909_WTeh_RANJdata.pdf",height = 4,width = 2.5)
###################################
par(mar=c(8,4,3,2)) 
pc<-sRNA.2122nt.denovo2021.RPM.Ranj$WT.eheart [sRNA.2122nt.denovo2021.RPM.Ranj$gene %in% denovoPC.loci$gene]
as<-sRNA.2122nt.denovo2021.RPM.Ranj$WT.eheart [sRNA.2122nt.denovo2021.RPM.Ranj$gene %in% lncRNAs.antisense.loci$gene]
linc<-sRNA.2122nt.denovo2021.RPM.Ranj$WT.eheart [sRNA.2122nt.denovo2021.RPM.Ranj$gene %in% lncRNAs.intergenic.loci$gene]
te<-sRNA.2122nt.denovo2021.RPM.Ranj$WT.eheart [sRNA.2122nt.denovo2021.RPM.Ranj$gene %in% TE_genes.loci$gene]

boxplot(pc,as,linc,te,col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="RPM, 24nt",main="21-22nt sRNA level\n early heart, Col-0")
#mtext('21-22nt sRNA', side=3, line=-1, at=2,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()



#################
#24 nt coverage boxplot Col0 leaf- 4 gene types (Ranj)
###################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Suppl_Fig3_boxplot_sRNA24nt_6909_WTleaf_RANJdata.pdf",height = 4,width = 2.5)
###################################
par(mar=c(8,4,3,2)) 
pc<-sRNA.24nt.denovo2021.RPM.Ranj$WT.leaf [sRNA.24nt.denovo2021.RPM.Ranj$gene %in% denovoPC.loci$gene]
as<-sRNA.24nt.denovo2021.RPM.Ranj$WT.leaf [sRNA.24nt.denovo2021.RPM.Ranj$gene %in% lncRNAs.antisense.loci$gene]
linc<-sRNA.24nt.denovo2021.RPM.Ranj$WT.leaf [sRNA.24nt.denovo2021.RPM.Ranj$gene %in% lncRNAs.intergenic.loci$gene]
te<-sRNA.24nt.denovo2021.RPM.Ranj$WT.leaf [sRNA.24nt.denovo2021.RPM.Ranj$gene %in% TE_genes.loci$gene]

boxplot(pc,as,linc,te,col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="RPM, 24nt",main="24nt sRNA level\n leaf, Col-0")
#mtext('24nt sRNA', side=3, line=-1, at=2,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.4)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.4)
#################
dev.off()

#################
#21-22 nt coverage boxplot Col0 leaf - 4 gene types (Ranj)
###################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Suppl_Fig3_boxplot_sRNA21nt_6909_WTleaf_RANJdata.pdf",height = 4,width = 2.5)
###################################
par(mar=c(8,4,3,2)) 
pc<-sRNA.2122nt.denovo2021.RPM.Ranj$WT.leaf [sRNA.2122nt.denovo2021.RPM.Ranj$gene %in% denovoPC.loci$gene]
as<-sRNA.2122nt.denovo2021.RPM.Ranj$WT.leaf [sRNA.2122nt.denovo2021.RPM.Ranj$gene %in% lncRNAs.antisense.loci$gene]
linc<-sRNA.2122nt.denovo2021.RPM.Ranj$WT.leaf [sRNA.2122nt.denovo2021.RPM.Ranj$gene %in% lncRNAs.intergenic.loci$gene]
te<-sRNA.2122nt.denovo2021.RPM.Ranj$WT.leaf [sRNA.2122nt.denovo2021.RPM.Ranj$gene %in% TE_genes.loci$gene]

boxplot(pc,as,linc,te,col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="RPM, 24nt",main="21-22nt sRNA level\n leaves, Col-0")
#mtext('21-22nt sRNA', side=3, line=-1, at=2,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.2)
a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.2)
a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.2)
#################
dev.off()









##################################################################################
# smallRNA coverage expressed-silent 
###################################################################################

linc_EC <-
te_EC
##########################
#24 nt coverage boxplot - lincs and TEs on off 
##########################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_sRNA24nt_expressed_silent_6909.pdf",height = 4,width =3)
##########################
par(mar=c(8,4,3,2)) 
l1<-sRNA.24nt.denovo2021.RPM$X6909 [sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene &sRNA.24nt.denovo2021.RPM$gene %in% denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$R.6909>0.5]]
l2<-   sRNA.24nt.denovo2021.RPM$X6909 [sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene & sRNA.24nt.denovo2021.RPM$gene %in% denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$R.6909<0.5]]
l3<-     sRNA.24nt.denovo2021.RPM$X6909 [sRNA.24nt.denovo2021.RPM$gene %in% TE_genes.loci$gene &sRNA.24nt.denovo2021.RPM$gene %in% denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$R.6909>0.5]]
l4<-     sRNA.24nt.denovo2021.RPM$X6909 [sRNA.24nt.denovo2021.RPM$gene %in% TE_genes.loci$gene &sRNA.24nt.denovo2021.RPM$gene %in% denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$R.6909<0.5]]
        boxplot(l1,l2,l3,l4,
        col=c("#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("lincRNAs\nON","lincRNAs\nOFF","TEgenes\nON","TEgenes\nOFF"),las=2, notch = T, outline = F, ylab="RPM, 24nt")
mtext('24nt sRNA', side=3, line=-1, at=1,cex=0.9)
#################
#add p values   #
#################
a<-wilcox.test(l1,l2)
b<-wilcox.test(l1,l2)
c<-wilcox.test(l1,l2)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.2)
a<-wilcox.test(l3,l4)
b<-wilcox.test(l3,l4)
c<-wilcox.test(l3,l4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.2)

#################
dev.off()




##########################

##########################
#21 nt coverage boxplot - lincs and TEs on off 
##########################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_sRNA21nt_expressed_silent_6909.pdf",height = 4,width =3)
par(mar=c(8,4,3,2)) 
boxplot(sRNA.21nt.denovo2021.RPM$X6909 [sRNA.21nt.denovo2021.RPM$gene %in% linc_EC$gene[linc_EC$R.6909>0.5]],
        sRNA.21nt.denovo2021.RPM$X6909 [sRNA.21nt.denovo2021.RPM$gene %in% linc_EC$gene[linc_EC$R.6909<0.5]],
        sRNA.21nt.denovo2021.RPM$X6909 [sRNA.21nt.denovo2021.RPM$gene %in% te_EC$gene[te_EC$R.6909>0.5]],
        sRNA.21nt.denovo2021.RPM$X6909 [sRNA.21nt.denovo2021.RPM$gene %in% te_EC$gene[te_EC$R.6909<0.5]],
        
        col=c("#F2AB54","#f7d1a3","#673A8E","#d6b9f0"), names=c("lincRNAs\nON","lincRNAs\nOFF","TEgenes\nON","TEgenes\nOFF"),las=2, notch = T, outline = F, ylab="RPM, 24nt")
mtext('21-22nt sRNA', side=3, line=-1, at=2,cex=0.9)
dev.off()
##########################






#####################################################
# confirm epigenetic patterns in other accessions
#####################################################


#methylation 
#boxplots in accessions 1741, 5784, 9905 ,9888 (4 randomly picked accessions)

#CG 1001G new data boxplot 4 accessions
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/boxplot_CGmeth_4accessions_1001GNEW.pdf",height = 3.5,width = 5)
###########################################################
par(mar=c(8,4,3,2)) 
pc1<-CG.1001new.denovo$mean.1741[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]
as1<-CG.1001new.denovo$mean.1741[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc1<-CG.1001new.denovo$mean.1741[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te1<-CG.1001new.denovo$mean.1741[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]

pc2<-CG.1001new.denovo$mean.5784[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]
as2<-CG.1001new.denovo$mean.5784[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc2<-CG.1001new.denovo$mean.5784[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te2<-CG.1001new.denovo$mean.5784[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]

pc3<-CG.1001new.denovo$mean.9905[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]
as3<-CG.1001new.denovo$mean.9905[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc3<-CG.1001new.denovo$mean.9905[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te3<-CG.1001new.denovo$mean.9905[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]

pc4<-CG.1001new.denovo$mean.9888[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]
as4<-CG.1001new.denovo$mean.9888[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc4<-CG.1001new.denovo$mean.9888[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te4<-CG.1001new.denovo$mean.9888[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]

boxplot( pc1,pc2,pc3,pc4,as1,as2,as3,as4,linc1,linc2,linc3,linc4,te1 ,te2 ,te3 ,te4    ,   col=c("#486EB4","#486EB4","#486EB4","#486EB4","#90C473","#90C473","#90C473","#90C473","#F2AB54","#F2AB54","#F2AB54","#F2AB54","#673A8E","#673A8E","#673A8E","#673A8E"),las=2, notch = T, outline = F, ylab="CG methylation in the locus",names=c("1741","5784","9905","9888","1741","5784","9905","9888","1741","5784","9905","9888","1741","5784","9905","9888"))
mtext('CG', side=3, line=2, at=6,cex=0.9)
mtext('PC genes', side=3, line=0, at=2,cex=0.9)
mtext('AS lncRNAs', side=3, line=0, at=7,cex=0.9)
mtext('lincRNAs', side=3, line=0, at=12,cex=0.9)
mtext('TE genes', side=3, line=0, at=16,cex=0.9)

###############################################
dev.off()

#CHH 1001G new data boxplot 4 accessions
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/boxplot_CHHmeth_4accessions_1001GNEW.pdf",height = 3.5,width = 5)
###########################################################
par(mar=c(8,4,3,2)) 
pc1<-CHH.1001new.denovo$mean.1741[CHH.1001new.denovo$transcript %in% denovoPC.loci$gene]
as1<-CHH.1001new.denovo$mean.1741[CHH.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc1<-CHH.1001new.denovo$mean.1741[CHH.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te1<-CHH.1001new.denovo$mean.1741[CHH.1001new.denovo$transcript %in% TE_genes.loci$gene]

pc2<-CHH.1001new.denovo$mean.5784[CHH.1001new.denovo$transcript %in% denovoPC.loci$gene]
as2<-CHH.1001new.denovo$mean.5784[CHH.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc2<-CHH.1001new.denovo$mean.5784[CHH.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te2<-CHH.1001new.denovo$mean.5784[CHH.1001new.denovo$transcript %in% TE_genes.loci$gene]

pc3<-CHH.1001new.denovo$mean.9905[CHH.1001new.denovo$transcript %in% denovoPC.loci$gene]
as3<-CHH.1001new.denovo$mean.9905[CHH.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc3<-CHH.1001new.denovo$mean.9905[CHH.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te3<-CHH.1001new.denovo$mean.9905[CHH.1001new.denovo$transcript %in% TE_genes.loci$gene]

pc4<-CHH.1001new.denovo$mean.9888[CHH.1001new.denovo$transcript %in% denovoPC.loci$gene]
as4<-CHH.1001new.denovo$mean.9888[CHH.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc4<-CHH.1001new.denovo$mean.9888[CHH.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te4<-CHH.1001new.denovo$mean.9888[CHH.1001new.denovo$transcript %in% TE_genes.loci$gene]

boxplot( pc1,pc2,pc3,pc4,as1,as2,as3,as4,linc1,linc2,linc3,linc4,te1 ,te2 ,te3 ,te4    ,   col=c("#486EB4","#486EB4","#486EB4","#486EB4","#90C473","#90C473","#90C473","#90C473","#F2AB54","#F2AB54","#F2AB54","#F2AB54","#673A8E","#673A8E","#673A8E","#673A8E"),las=2, notch = T, outline = F, ylab="CG methylation in the locus",names=c("1741","5784","9905","9888","1741","5784","9905","9888","1741","5784","9905","9888","1741","5784","9905","9888"))
mtext('CHH', side=3, line=2, at=6,cex=0.9)
mtext('PC genes', side=3, line=0, at=2,cex=0.9)
mtext('AS lncRNAs', side=3, line=0, at=7,cex=0.9)
mtext('lincRNAs', side=3, line=0, at=12,cex=0.9)
mtext('TE genes', side=3, line=0, at=16,cex=0.9)

###############################################
dev.off()

#CHG 1001G new data boxplot 4 accessions
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/boxplot_CHGmeth_4accessions_1001GNEW.pdf",height = 3.5,width = 5)
###########################################################
par(mar=c(8,4,3,2)) 
pc1<-CHG.1001new.denovo$mean.1741[CHG.1001new.denovo$transcript %in% denovoPC.loci$gene]
as1<-CHG.1001new.denovo$mean.1741[CHG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc1<-CHG.1001new.denovo$mean.1741[CHG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te1<-CHG.1001new.denovo$mean.1741[CHG.1001new.denovo$transcript %in% TE_genes.loci$gene]

pc2<-CHG.1001new.denovo$mean.5784[CHG.1001new.denovo$transcript %in% denovoPC.loci$gene]
as2<-CHG.1001new.denovo$mean.5784[CHG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc2<-CHG.1001new.denovo$mean.5784[CHG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te2<-CHG.1001new.denovo$mean.5784[CHG.1001new.denovo$transcript %in% TE_genes.loci$gene]

pc3<-CHG.1001new.denovo$mean.9905[CHG.1001new.denovo$transcript %in% denovoPC.loci$gene]
as3<-CHG.1001new.denovo$mean.9905[CHG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc3<-CHG.1001new.denovo$mean.9905[CHG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te3<-CHG.1001new.denovo$mean.9905[CHG.1001new.denovo$transcript %in% TE_genes.loci$gene]

pc4<-CHG.1001new.denovo$mean.9888[CHG.1001new.denovo$transcript %in% denovoPC.loci$gene]
as4<-CHG.1001new.denovo$mean.9888[CHG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc4<-CHG.1001new.denovo$mean.9888[CHG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te4<-CHG.1001new.denovo$mean.9888[CHG.1001new.denovo$transcript %in% TE_genes.loci$gene]

boxplot( pc1,pc2,pc3,pc4,as1,as2,as3,as4,linc1,linc2,linc3,linc4,te1 ,te2 ,te3 ,te4    ,   col=c("#486EB4","#486EB4","#486EB4","#486EB4","#90C473","#90C473","#90C473","#90C473","#F2AB54","#F2AB54","#F2AB54","#F2AB54","#673A8E","#673A8E","#673A8E","#673A8E"),las=2, notch = T, outline = F, ylab="CHG methylation in the locus",names=c("1741","5784","9905","9888","1741","5784","9905","9888","1741","5784","9905","9888","1741","5784","9905","9888"))
mtext('CHG', side=3, line=2, at=6,cex=0.9)
mtext('PC genes', side=3, line=0, at=2,cex=0.9)
mtext('AS lncRNAs', side=3, line=0, at=7,cex=0.9)
mtext('lincRNAs', side=3, line=0, at=12,cex=0.9)
mtext('TE genes', side=3, line=0, at=16,cex=0.9)

###############################################
dev.off()





#1001G new data  CG ditribution density plot - 14 accessions 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/densityplot_CG_from1001GNEW.pdf",height =5,width =5)
###########################################################
par(mar=c(4,4,3,2),mfrow=c(2,2)) 
plot(density(CG.1001new.denovo$mean.5210[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),xlab="CG methylation level",main="",lwd=1,las=2, col=alpha(colour = "#486EB4",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.1741[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),lwd=1, col=alpha(colour = "#486EB4",alpha = 0.5))
mtext("PC genes",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.4807[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),lwd=1, col=alpha(colour = "#486EB4",alpha = 0.5))
mtext("PC genes",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.8366[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),lwd=1, col=alpha(colour = "#486EB4",alpha = 0.5))
mtext("PC genes",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.5772[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),lwd=1, col=alpha(colour = "#486EB4",alpha = 0.5))
mtext("PC genes",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.5856[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),lwd=1, col=alpha(colour = "#486EB4",alpha = 0.5))
mtext("PC genes",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.6021[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),lwd=1, col=alpha(colour = "#486EB4",alpha = 0.5))
mtext("PC genes",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.6220[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),lwd=1, col=alpha(colour = "#486EB4",alpha = 0.5))
mtext("PC genes",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.9057[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),lwd=1, col=alpha(colour = "#486EB4",alpha = 0.5))
mtext("PC genes",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.6966[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),lwd=1, col=alpha(colour = "#486EB4",alpha = 0.5))
mtext("PC genes",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.8244[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),lwd=1, col=alpha(colour = "#486EB4",alpha = 1))
mtext("PC genes",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.9412[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),lwd=1, col=alpha(colour = "#486EB4",alpha = 1))
mtext("PC genes",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),lwd=1, col=alpha(colour = "#486EB4",alpha = 1))
mtext("PC genes",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]),lwd=1, col=alpha(colour = "black",alpha = 1))
mtext("PC genes",side=3,at=0.6,line = -1)

plot(density(CG.1001new.denovo$mean.5210[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),xlab="CG methylation level",main="",lwd=1,las=2, col=alpha(colour = "#90C473",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.1741[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),lwd=1, col=alpha(colour = "#90C473",alpha = 0.5))
mtext("AS lncRNAs",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.4807[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),lwd=1, col=alpha(colour = "#90C473",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.8366[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),lwd=1, col=alpha(colour = "#90C473",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.5772[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),lwd=1, col=alpha(colour = "#90C473",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.5856[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),lwd=1, col=alpha(colour = "#90C473",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.6021[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),lwd=1, col=alpha(colour = "#90C473",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.6220[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),lwd=1, col=alpha(colour = "#90C473",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.9057[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),lwd=1, col=alpha(colour = "#90C473",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.6966[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),lwd=1, col=alpha(colour = "#90C473",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.8244[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),lwd=1, col=alpha(colour = "#90C473",alpha = 1))
lines(density(CG.1001new.denovo$mean.9412[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),lwd=1, col=alpha(colour = "#90C473",alpha = 1))
lines(density(CG.1001new.denovo$mean.[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),lwd=1, col=alpha(colour = "#90C473",alpha = 1))
lines(density(CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]),lwd=1, col=alpha(colour = "black",alpha = 1))

plot(density(CG.1001new.denovo$mean.5210[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),xlab="CG methylation level",main="",lwd=1,las=2, col=alpha(colour = "#F2AB54",alpha = 0.5), ylim=c(0,2.5))
lines(density(CG.1001new.denovo$mean.1741[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),lwd=1, col=alpha(colour = "#F2AB54",alpha = 0.5))
mtext("lincRNAs",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.4807[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),lwd=1, col=alpha(colour = "#F2AB54",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.8366[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),lwd=1, col=alpha(colour = "#F2AB54",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.5772[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),lwd=1, col=alpha(colour = "#F2AB54",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.5856[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),lwd=1, col=alpha(colour = "#F2AB54",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.6021[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),lwd=1, col=alpha(colour = "#F2AB54",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.6220[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),lwd=1, col=alpha(colour = "#F2AB54",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.9057[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),lwd=1, col=alpha(colour = "#F2AB54",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.6966[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),lwd=1, col=alpha(colour = "#F2AB54",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.8244[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),lwd=1, col=alpha(colour = "#F2AB54",alpha = 1))
lines(density(CG.1001new.denovo$mean.9412[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),lwd=1, col=alpha(colour = "#F2AB54",alpha = 1))
lines(density(CG.1001new.denovo$mean.[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),lwd=1, col=alpha(colour = "#F2AB54",alpha = 1))
lines(density(CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]),lwd=1, col=alpha(colour = "black",alpha = 1))

plot(density(CG.1001new.denovo$mean.5210[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),xlab="CG methylation level",main="",lwd=1,las=2, col=alpha(colour = "#673A8E",alpha = 0.5), ylim=c(0,7))
lines(density(CG.1001new.denovo$mean.1741[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),lwd=1, col=alpha(colour = "#673A8E",alpha = 0.5))
mtext("TE genes",side=3,at=0.6,line = -1)
lines(density(CG.1001new.denovo$mean.4807[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),lwd=1, col=alpha(colour = "#673A8E",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.8366[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),lwd=1, col=alpha(colour = "#673A8E",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.5772[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),lwd=1, col=alpha(colour = "#673A8E",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.5856[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),lwd=1, col=alpha(colour = "#673A8E",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.6021[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),lwd=1, col=alpha(colour = "#673A8E",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.6220[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),lwd=1, col=alpha(colour = "#673A8E",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.9057[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),lwd=1, col=alpha(colour = "#673A8E",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.6966[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),lwd=1, col=alpha(colour = "#673A8E",alpha = 0.5))
lines(density(CG.1001new.denovo$mean.8244[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),lwd=1, col=alpha(colour = "#673A8E",alpha = 1))
lines(density(CG.1001new.denovo$mean.9412[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),lwd=1, col=alpha(colour = "#673A8E",alpha = 1))
lines(density(CG.1001new.denovo$mean.[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),lwd=1, col=alpha(colour = "#673A8E",alpha = 1))
lines(density(CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]),lwd=1, col=alpha(colour = "black",alpha = 1))
###########################################################
dev.off()



#chipseq 



#K9 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/other_4_accessions_boxplot_K9_6909.pdf",height = 3.5,width = 5)
#################################################
pc1<-chip.denovo.quantstan$K9.1741 [chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as1<-chip.denovo.quantstan$K9.1741 [chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc1<-chip.denovo.quantstan$K9.1741 [chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te1<-chip.denovo.quantstan$K9.1741 [chip.denovo.quantstan$gene %in% TE_genes.loci$gene]
pc2<-chip.denovo.quantstan$K9.5784 [chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as2<-chip.denovo.quantstan$K9.5784 [chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc2<-chip.denovo.quantstan$K9.5784 [chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te2<-chip.denovo.quantstan$K9.5784 [chip.denovo.quantstan$gene %in% TE_genes.loci$gene]
pc3<-chip.denovo.quantstan$K9.9905 [chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as3<-chip.denovo.quantstan$K9.9905 [chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc3<-chip.denovo.quantstan$K9.9905 [chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te3<-chip.denovo.quantstan$K9.9905 [chip.denovo.quantstan$gene %in% TE_genes.loci$gene]
pc4<-chip.denovo.quantstan$K9.9888 [chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as4<-chip.denovo.quantstan$K9.9888 [chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc4<-chip.denovo.quantstan$K9.9888 [chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te4<-chip.denovo.quantstan$K9.9888 [chip.denovo.quantstan$gene %in% TE_genes.loci$gene]

par(mar=c(8,4,3,2)) 
names=c("PC genes","AS lncRNAs","lincRNAs","TE genes")
boxplot( pc1,pc2,pc3,pc4,as1,as2,as3,as4,linc1,linc2,linc3,linc4,te1 ,te2 ,te3 ,te4    ,   col=c("#486EB4","#486EB4","#486EB4","#486EB4","#90C473","#90C473","#90C473","#90C473","#F2AB54","#F2AB54","#F2AB54","#F2AB54","#673A8E","#673A8E","#673A8E","#673A8E"),las=2, notch = T, outline = F, ylab="quantile normalized log2(ChIP/Input)",names=c("1741","5784","9905","9888","1741","5784","9905","9888","1741","5784","9905","9888","1741","5784","9905","9888"))
mtext('H3K9me2', side=3, line=2, at=6,cex=0.9)
mtext('PC genes', side=3, line=0, at=2,cex=0.9)
mtext('AS lncRNAs', side=3, line=0, at=7,cex=0.9)
mtext('lincRNAs', side=3, line=0, at=12,cex=0.9)
mtext('TE genes', side=3, line=0, at=16,cex=0.9)

a<-wilcox.test(sample(pc,2000),sampleas)

dev.off()

# K27 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K27_6909.pdf",height = 3.5,width = 2.5)
#################################################
pc<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% denovoPC.loci$gene]
as<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo.log2$K27.6909 [chip.denovo.log2$gene %in% TE_genes.loci$gene]

par(mar=c(8,4,3,2)) 
boxplot(pc,as,linc,te    , 
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext('H3K27me3', side=3, line=-1, at=2,cex=0.9)
dev.off()


#H1 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/other_4_accessions_boxplot_H1_6909.pdf",height = 3.5,width = 5)
#################################################
pc1<-chip.denovo.quantstan$H1.1741 [chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as1<-chip.denovo.quantstan$H1.1741 [chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc1<-chip.denovo.quantstan$H1.1741 [chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te1<-chip.denovo.quantstan$H1.1741 [chip.denovo.quantstan$gene %in% TE_genes.loci$gene]
pc2<-chip.denovo.quantstan$H1.5784 [chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as2<-chip.denovo.quantstan$H1.5784 [chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc2<-chip.denovo.quantstan$H1.5784 [chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te2<-chip.denovo.quantstan$H1.5784 [chip.denovo.quantstan$gene %in% TE_genes.loci$gene]
pc3<-chip.denovo.quantstan$H1.9905 [chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as3<-chip.denovo.quantstan$H1.9905 [chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc3<-chip.denovo.quantstan$H1.9905 [chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te3<-chip.denovo.quantstan$H1.9905 [chip.denovo.quantstan$gene %in% TE_genes.loci$gene]
pc4<-chip.denovo.quantstan$H1.9888 [chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as4<-chip.denovo.quantstan$H1.9888 [chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc4<-chip.denovo.quantstan$H1.9888 [chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te4<-chip.denovo.quantstan$H1.9888 [chip.denovo.quantstan$gene %in% TE_genes.loci$gene]

par(mar=c(8,4,3,2)) 
names=c("PC genes","AS lncRNAs","lincRNAs","TE genes")
boxplot( pc1,pc2,pc3,pc4,as1,as2,as3,as4,linc1,linc2,linc3,linc4,te1 ,te2 ,te3 ,te4    ,   col=c("#486EB4","#486EB4","#486EB4","#486EB4","#90C473","#90C473","#90C473","#90C473","#F2AB54","#F2AB54","#F2AB54","#F2AB54","#673A8E","#673A8E","#673A8E","#673A8E"),las=2, notch = T, outline = F, ylab="quantile normalized log2(ChIP/Input)",names=c("1741","5784","9905","9888","1741","5784","9905","9888","1741","5784","9905","9888","1741","5784","9905","9888"))
mtext('H1', side=3, line=2, at=6,cex=0.9)
mtext('PC genes', side=3, line=0, at=2,cex=0.9)
mtext('AS lncRNAs', side=3, line=0, at=7,cex=0.9)
mtext('lincRNAs', side=3, line=0, at=12,cex=0.9)
mtext('TE genes', side=3, line=0, at=16,cex=0.9)

a<-wilcox.test(sample(pc,2000),sampleas)
####################################
dev.off()

#K36 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/other_4_accessions_boxplot_K36_6909.pdf",height = 3.5,width = 5)
#################################################
pc1<-chip.denovo.quantstan$K36.1741 [chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as1<-chip.denovo.quantstan$K36.1741 [chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc1<-chip.denovo.quantstan$K36.1741 [chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te1<-chip.denovo.quantstan$K36.1741 [chip.denovo.quantstan$gene %in% TE_genes.loci$gene]
pc2<-chip.denovo.quantstan$K36.5784 [chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as2<-chip.denovo.quantstan$K36.5784 [chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc2<-chip.denovo.quantstan$K36.5784 [chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te2<-chip.denovo.quantstan$K36.5784 [chip.denovo.quantstan$gene %in% TE_genes.loci$gene]
pc3<-chip.denovo.quantstan$K36.9905 [chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as3<-chip.denovo.quantstan$K36.9905 [chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc3<-chip.denovo.quantstan$K36.9905 [chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te3<-chip.denovo.quantstan$K36.9905 [chip.denovo.quantstan$gene %in% TE_genes.loci$gene]
pc4<-chip.denovo.quantstan$K36.9888 [chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as4<-chip.denovo.quantstan$K36.9888 [chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc4<-chip.denovo.quantstan$K36.9888 [chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te4<-chip.denovo.quantstan$K36.9888 [chip.denovo.quantstan$gene %in% TE_genes.loci$gene]

par(mar=c(8,4,3,2)) 
names=c("PC genes","AS lncRNAs","lincRNAs","TE genes")
boxplot( pc1,pc2,pc3,pc4,as1,as2,as3,as4,linc1,linc2,linc3,linc4,te1 ,te2 ,te3 ,te4    ,   col=c("#486EB4","#486EB4","#486EB4","#486EB4","#90C473","#90C473","#90C473","#90C473","#F2AB54","#F2AB54","#F2AB54","#F2AB54","#673A8E","#673A8E","#673A8E","#673A8E"),las=2, notch = T, outline = F, ylab="quantile normalized log2(ChIP/Input)",names=c("1741","5784","9905","9888","1741","5784","9905","9888","1741","5784","9905","9888","1741","5784","9905","9888"))
mtext('H3K36me3', side=3, line=2, at=6,cex=0.9)
mtext('PC genes', side=3, line=0, at=2,cex=0.9)
mtext('AS lncRNAs', side=3, line=0, at=7,cex=0.9)
mtext('lincRNAs', side=3, line=0, at=12,cex=0.9)
mtext('TE genes', side=3, line=0, at=16,cex=0.9)

a<-wilcox.test(sample(pc,2000),sampleas)
####################################
dev.off()



#24 nt 




#################
#24 nt coverage boxplot other accessions flowers - 4 gene types (this study)
###################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/otheraccessions_boxplot_sRNA24nt_6909.pdf",height = 3.5,width = 5)
###################################
par(mar=c(8,4,3,2)) 
pc1<-sRNA.24nt.denovo2021.RPM$X1741 [sRNA.24nt.denovo2021.RPM$gene %in% denovoPC.loci$gene]
as1<-sRNA.24nt.denovo2021.RPM$X1741 [sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.antisense.loci$gene]
linc1<-sRNA.24nt.denovo2021.RPM$X1741 [sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene]
te1<-sRNA.24nt.denovo2021.RPM$X1741 [sRNA.24nt.denovo2021.RPM$gene %in% TE_genes.loci$gene]

pc2<-sRNA.24nt.denovo2021.RPM$X10002 [sRNA.24nt.denovo2021.RPM$gene %in% denovoPC.loci$gene]
as2<-sRNA.24nt.denovo2021.RPM$X10002 [sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.antisense.loci$gene]
linc2<-sRNA.24nt.denovo2021.RPM$X10002 [sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene]
te2<-sRNA.24nt.denovo2021.RPM$X10002 [sRNA.24nt.denovo2021.RPM$gene %in% TE_genes.loci$gene]

pc3<-sRNA.24nt.denovo2021.RPM$X9905 [sRNA.24nt.denovo2021.RPM$gene %in% denovoPC.loci$gene]
as3<-sRNA.24nt.denovo2021.RPM$X9905 [sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.antisense.loci$gene]
linc3<-sRNA.24nt.denovo2021.RPM$X9905 [sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene]
te3<-sRNA.24nt.denovo2021.RPM$X9905 [sRNA.24nt.denovo2021.RPM$gene %in% TE_genes.loci$gene]

pc4<-sRNA.24nt.denovo2021.RPM$X9888 [sRNA.24nt.denovo2021.RPM$gene %in% denovoPC.loci$gene]
as4<-sRNA.24nt.denovo2021.RPM$X9888 [sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.antisense.loci$gene]
linc4<-sRNA.24nt.denovo2021.RPM$X9888 [sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene]
te4<-sRNA.24nt.denovo2021.RPM$X9888 [sRNA.24nt.denovo2021.RPM$gene %in% TE_genes.loci$gene]

boxplot( pc1,pc2,pc3,pc4,as1,as2,as3,as4,linc1,linc2,linc3,linc4,te1 ,te2 ,te3 ,te4    ,   col=c("#486EB4","#486EB4","#486EB4","#486EB4","#90C473","#90C473","#90C473","#90C473","#F2AB54","#F2AB54","#F2AB54","#F2AB54","#673A8E","#673A8E","#673A8E","#673A8E"),las=2, notch = T, outline = F, ylab="RPM, 24nt",names=c("1741","10002","9905","9888","1741","10002","9905","9888","1741","10002","9905","9888","1741","10002","9905","9888"))

mtext('PC genes', side=3, line=0, at=2,cex=0.9)
mtext('AS lncRNAs', side=3, line=0, at=7,cex=0.9)
mtext('lincRNAs', side=3, line=0, at=12,cex=0.9)
mtext('TE genes', side=3, line=0, at=16,cex=0.9)

#######################################
dev.off()





