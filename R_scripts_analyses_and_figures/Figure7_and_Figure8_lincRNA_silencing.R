##################################################################
######## Analyses for figure 7 and its supplements 
##################################################################
install.packages("scales")
library(scales)
library(pheatmap)
#k9 vs k27
pc_sasa<-denovo2021.TPMs.genes.1001Gnew[denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene,]
as_sasa<-denovo2021.TPMs.genes.1001Gnew[denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene,]
linc_sasa<-denovo2021.TPMs.genes.1001Gnew[denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene,]
te_sasa<-denovo2021.TPMs.genes.1001Gnew[denovo2021.TPMs.genes.1001Gnew$gene %in% TE_genes.loci$gene,]

#k27 vs k9

######################################################
#Scatter plot K9 vs K27 : PC, AS, linc and TE
###########################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/scatterplot_pc_as_linc_te_K27_vs_K9_6909.pdf",height = 2.5,width =10)
###########################################################
par(mar=c(4,4,2,1)) 
par(mfrow=c(1,4))
#PC
a<-chip.denovo.log2[chip.denovo.log2$gene %in% pc_sasa$gene[pc_sasa$mean.6909<0.5],]

plot (a$K9.6909,a$K27.6909, ylab="H3K27me3, log2(ChIP/Input)",xlab="H3K9me2, log2(ChIP/Input)",col=rgb(red=10, green=10, blue=10, alpha=50, maxColorValue=255),pch=20, main="PC genes",ylim=c(-2.2,2.5),xlim=c(-2.5,2.5))

points(a$K9.6909[a$K9.6909<0 & a$K27.6909>0],a$K27.6909[a$K9.6909<0 & a$K27.6909>0], ylab="H3K27me3, log2(ChIP/Input)",xlab="H3K9me2, log2(ChIP/Input)",col=rgb(red=180, green=10, blue=10, alpha=120, maxColorValue=255),pch=20)

points(a$K9.6909[a$K9.6909>0 & a$K27.6909<0],a$K27.6909[a$K9.6909>0 & a$K27.6909<0], ylab="H3K27me3, log2(ChIP/Input)",xlab="H3K9me2, log2(ChIP/Input)",col=rgb(red=10, green=10, blue=180, alpha=120, maxColorValue=255),pch=20)

text(x=2,y=-2,adj=c(0.5,0.5),col=rgb(red=10, green=10, blue=180,  maxColorValue=255),
     length(a$K27.6909[a$K9.6909>0 & a$K27.6909<0]))
text(x=-2,y=2,adj=c(0.5,0.5),col=rgb(red=180, green=10, blue=10,  maxColorValue=255),
     length(a$K9.6909[a$K9.6909<0 & a$K27.6909>0]))
length(a$K9.6909[a$K9.6909>0 & a$K27.6909<0])#579 K9
length(a$K27.6909[a$K9.6909<0 & a$K27.6909>0])#3586 K27
length(a$K27.6909) # 7602 all silent 
#AS
a<-chip.denovo.log2[chip.denovo.log2$gene %in% as_sasa$gene[as_sasa$mean.6909<0.5],]
plot (a$K9.6909,a$K27.6909, ylab="H3K27me3, log2(ChIP/Input)",xlab="H3K9me2, log2(ChIP/Input)",col=rgb(red=10, green=10, blue=10, alpha=50, maxColorValue=255),pch=20, main="AS lncRNAs",ylim=c(-2.2,2.5),xlim=c(-2.5,2.5))

points(a$K9.6909[a$K9.6909<0 & a$K27.6909>0],a$K27.6909[a$K9.6909<0 & a$K27.6909>0], ylab="H3K27me3, log2(ChIP/Input)",xlab="H3K9me2, log2(ChIP/Input)",col=rgb(red=180, green=10, blue=10, alpha=120, maxColorValue=255),pch=20)

points(a$K9.6909[a$K9.6909>0 & a$K27.6909<0],a$K27.6909[a$K9.6909>0 & a$K27.6909<0], ylab="H3K27me3, log2(ChIP/Input)",xlab="H3K9me2, log2(ChIP/Input)",col=rgb(red=10, green=10, blue=180, alpha=120, maxColorValue=255),pch=20)

text(x=2,y=-2,adj=c(0.5,0.5),col=rgb(red=10, green=10, blue=180,  maxColorValue=255),
     length(a$K27.6909[a$K9.6909>0 & a$K27.6909<0]))
text(x=-2,y=2,adj=c(0.5,0.5),col=rgb(red=180, green=10, blue=10,  maxColorValue=255),
     length(a$K9.6909[a$K9.6909<0 & a$K27.6909>0]))
length(a$K9.6909[a$K9.6909>0 & a$K27.6909<0])#258 K9
length(a$K27.6909[a$K9.6909<0 & a$K27.6909>0])#1522 K27
length(a$K27.6909) # 7688 all silent

#linc
a<-chip.denovo.log2[chip.denovo.log2$gene %in%  linc_sasa$gene[linc_sasa$mean.6909<0.5],]
plot (a$K9.6909,a$K27.6909, ylab="H3K27me3, log2(ChIP/Input)",xlab="H3K9me2, log2(ChIP/Input)",col=rgb(red=10, green=10, blue=10, alpha=50, maxColorValue=255),pch=20, main="lincRNAs",ylim=c(-2.2,2.5),xlim=c(-2.5,2.5))

points(a$K9.6909[a$K9.6909<0 & a$K27.6909>0],a$K27.6909[a$K9.6909<0 & a$K27.6909>0], ylab="H3K27me3, log2(ChIP/Input)",xlab="H3K9me2, log2(ChIP/Input)",col=rgb(red=180, green=10, blue=10, alpha=120, maxColorValue=255),pch=20)

points(a$K9.6909[a$K9.6909>0 & a$K27.6909<0],a$K27.6909[a$K9.6909>0 & a$K27.6909<0], ylab="H3K27me3, log2(ChIP/Input)",xlab="H3K9me2, log2(ChIP/Input)",col=rgb(red=10, green=10, blue=180, alpha=120, maxColorValue=255),pch=20)

text(x=2,y=-2,adj=c(0.5,0.5),col=rgb(red=10, green=10, blue=180,  maxColorValue=255),
     length(a$K27.6909[a$K9.6909>0 & a$K27.6909<0]))
text(x=-2,y=2,adj=c(0.5,0.5),col=rgb(red=180, green=10, blue=10,  maxColorValue=255),
     length(a$K9.6909[a$K9.6909<0 & a$K27.6909>0]))
length(a$K9.6909[a$K9.6909>0 & a$K27.6909<0])#656 K9
length(a$K27.6909[a$K9.6909<0 & a$K27.6909>0])#399 K27
length(a$K27.6909) # 2157 all silent


#te
a<-chip.denovo.log2[chip.denovo.log2$gene %in% te_sasa$gene[te_sasa$mean.6909<0.5],]
plot (a$K9.6909,a$K27.6909, ylab="H3K27me3, log2(ChIP/Input)",xlab="H3K9me2, log2(ChIP/Input)",col=rgb(red=10, green=10, blue=10, alpha=50, maxColorValue=255),pch=20, main="TE genes",ylim=c(-2.2,2.5),xlim=c(-2.5,2.5))

points(a$K9.6909[a$K9.6909<0 & a$K27.6909>0],a$K27.6909[a$K9.6909<0 & a$K27.6909>0], ylab="H3K27me3, log2(ChIP/Input)",xlab="H3K9me2, log2(ChIP/Input)",col=rgb(red=180, green=10, blue=10, alpha=120, maxColorValue=255),pch=20)

points(a$K9.6909[a$K9.6909>0 & a$K27.6909<0],a$K27.6909[a$K9.6909>0 & a$K27.6909<0], ylab="H3K27me3, log2(ChIP/Input)",xlab="H3K9me2, log2(ChIP/Input)",col=rgb(red=10, green=10, blue=180, alpha=120, maxColorValue=255),pch=20)

text(x=2,y=-2,adj=c(0.5,0.5),col=rgb(red=10, green=10, blue=180,  maxColorValue=255),
     length(a$K27.6909[a$K9.6909>0 & a$K27.6909<0]))
text(x=-2,y=2,adj=c(0.5,0.5),col=rgb(red=180, green=10, blue=10,  maxColorValue=255),
     length(a$K9.6909[a$K9.6909<0 & a$K27.6909>0]))

length(a$K27.6909[a$K9.6909>0 & a$K27.6909<0])#1865 K9
length(a$K9.6909[a$K9.6909<0 & a$K27.6909>0])#55 K27
length(a$K27.6909) # 2071 all silent


###########################################################
dev.off()

# make the list of k9 and k27 lincRNAs : only silent lincRNAs! 
a<-chip.denovo.log2[chip.denovo.log2$gene %in%  linc_sasa$gene[linc_sasa$mean.6909<0.5],]
length(a$K9.6909[a$K9.6909>0 & a$K27.6909<0])#656 K9
length(a$K27.6909[a$K9.6909<0 & a$K27.6909>0])#399 K27
################################################################
#define K9 and K27 lincRNAs
#################################################################
k27_lincs<-as.vector(a$gene[a$K9.6909<0 & a$K27.6909>0])
k9_lincs<-as.vector(a$gene[a$K9.6909>0 & a$K27.6909<0])
################################################################
length(k27_lincs)#399
length(k9_lincs)#656

# TE content
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Fig7_linc_boxplot_te_coverage_k27_k9.pdf",height = 3,width = 2.2)
####################################################
par(mar=c(5,4,3,2)) 
a1<-100*linc_TE_cov_all_loci_2cols$coverage[linc_TE_cov_all_loci_2cols$gene %in% k27_lincs  ]
a2<-100*linc_TE_cov_all_loci_2cols$coverage[linc_TE_cov_all_loci_2cols$gene %in% k9_lincs  ]
boxplot(a1,a2, main="TE content of the lincRNA locus",cex.main=1,
col=c("pink","lightblue","#673A8E"), names=c("K27","K9"),las=1, notch = T, outline = F,cex.lab=0.8, ylab="TE-sequence content, % of locus length ")
#add p values

a<-wilcox.test(sample(a2,399),a1)
b<-wilcox.test(sample(a2,399),a1)
c<-wilcox.test(sample(a2,399),a1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=90)
####################################################
dev.off()

# CHH 1001G
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Fig5_linc_boxplot_CHH_1001G_k27_k9.pdf",height = 3,width = 2.2)
####################################################
par(mar=c(5,4,3,2)) 
a1<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% k27_lincs  ]
a2<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% k9_lincs  ]
boxplot(a1,a2, main="CHH methylation of lincRNA loci\n (Col-0, rosette, 1001Genomes)",cex.main=1,
        col=c("pink","lightblue","#673A8E"), names=c("K27","K9"),las=1, notch = T, outline = F,cex.lab=0.8, ylab="CHH methylation level")
#add p values
a<-wilcox.test(sample(a2,399),a1)
b<-wilcox.test(sample(a2,399),a1)
c<-wilcox.test(sample(a2,399),a1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.10)
####################################################
dev.off()

# CG 1001G
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Fig5_linc_boxplot_CG_1001G_k27_k9.pdf",height = 3,width = 2.2)
####################################################
par(mar=c(5,4,3,2)) 
a1<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% k27_lincs  ]
a2<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% k9_lincs  ]
boxplot(a1,a2, main="CG methylation of lincRNA loci\n (Col-0, rosette, 1001Genomes)",cex.main=1,
        col=c("pink","lightblue","#673A8E"), names=c("K27","K9"),las=1, notch = T, outline = F,cex.lab=0.8, ylab="CG methylation level")
#add p values
a<-wilcox.test(sample(a2,399),a1)
b<-wilcox.test(sample(a2,399),a1)
c<-wilcox.test(sample(a2,399),a1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.98)
####################################################
dev.off()

# CHH 1001G NEW
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Fig5_linc_boxplot_CHH_1001GNEW_k27_k9.pdf",height = 3,width = 2.2)
####################################################
par(mar=c(5,4,3,2)) 
a1<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% k27_lincs]
a2<-CHH.1001new.denovo$mean.6909[CHH.1001new.denovo$transcript %in% k9_lincs]
boxplot(a1,a2, main="CHH methylation of lincRNA loci",cex.main=1, ylim=c(0,0.2),
        col=c("pink","lightblue","#673A8E"), names=c("K27","K9"),las=1, notch = T, outline = F,cex.lab=0.8, ylab="CHH methylation level")
#add p values
a<-wilcox.test(sample(a2,399),a1)
b<-wilcox.test(sample(a2,399),a1)
c<-wilcox.test(sample(a2,399),a1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.15)
####################################################
dev.off()

# CG 1001G NEW
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Fig5_linc_boxplot_CG_1001GNEW_k27_k9.pdf",height = 3,width = 2.2)
####################################################
par(mar=c(5,4,3,2)) 
a1<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% k27_lincs  ]
a2<-CG.1001new.denovo$mean.6909[CG.1001new.denovo$transcript %in% k9_lincs  ]
boxplot(a1,a2, main="CG methylation of lincRNA loci",cex.main=1,
        col=c("pink","lightblue","#673A8E"), names=c("K27","K9"),las=1, notch = T, outline = F,cex.lab=0.8, ylab="CG methylation level")
#add p values
a<-wilcox.test(sample(a2,399),a1)
b<-wilcox.test(sample(a2,399),a1)
c<-wilcox.test(sample(a2,399),a1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.98)
####################################################
dev.off()



# CHH 1001G variation
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Fig5_linc_boxplot_CHHvariation_1001G_k27_k9.pdf",height = 3,width = 2.2)
####################################################
par(mar=c(5,4,3,2)) 
a1<-CHH.1001.denovo$sd[CHH.1001.denovo$transcript %in% k27_lincs  ]
a2<-CHH.1001.denovo$sd[CHH.1001.denovo$transcript %in% k9_lincs  ]
boxplot(a1,a2, main="CHH methylation variation\n (444 accessions, 1001G)",cex.main=1,
        col=c("pink","lightblue","#673A8E"), names=c("K27","K9"),las=1, notch = T, outline = F,cex.lab=0.8, ylab="CHH methylation level")
#add p values
a<-wilcox.test(sample(a2,399),a1)
b<-wilcox.test(sample(a2,399),a1)
c<-wilcox.test(sample(a2,399),a1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.40)
####################################################
dev.off()

# CG 1001G variation
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Fig5_linc_boxplot_CGvariation_1001G_k27_k9.pdf",height = 3,width = 2.2)
####################################################
par(mar=c(5,4,3,2)) 
a1<-CG.1001.denovo$sd[CG.1001.denovo$transcript %in% k27_lincs  ]
a2<-CG.1001.denovo$sd[CG.1001.denovo$transcript %in% k9_lincs  ]
boxplot(a1,a2, main="CG methylation variation\n (444 accessions, 1001G)",cex.main=1,
        col=c("pink","lightblue","#673A8E"), names=c("K27","K9"),las=1, notch = T, outline = F,cex.lab=0.8, ylab="CG methylation level")
#add p values
a<-wilcox.test(sample(a2,399),a1)
b<-wilcox.test(sample(a2,399),a1)
c<-wilcox.test(sample(a2,399),a1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
####################################################
dev.off()



# CG 1001G NEW - variation
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Fig5_linc_boxplot_CGvariation_1001GNEW_k27_k9.pdf",height = 3,width = 2.2)
####################################################
par(mar=c(5,4,3,2)) 
a1<-CG.1001new.denovo$sd_of_means[CG.1001new.denovo$transcript %in% k27_lincs  ]
a2<-CG.1001new.denovo$sd_of_means[CG.1001new.denovo$transcript %in% k9_lincs  ]
boxplot(a1,a2, main="CG meth variation\n (28 accessions)",cex.main=1,
        col=c("pink","lightblue","#673A8E"), names=c("K27 linc","K9 linc"),las=1, notch = T, outline = F,cex.lab=0.8, ylab="CG methylation level")
#add p values
a<-wilcox.test(sample(a2,399),a1)
b<-wilcox.test(sample(a2,399),a1)
c<-wilcox.test(sample(a2,399),a1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
####################################################
dev.off()

# CHH 1001G NEW variation
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Fig5_linc_boxplot_CHHvariation_1001GNEW_k27_k9.pdf",height = 3,width = 2.2)
####################################################
par(mar=c(5,4,3,2)) 
a1<-CHH.1001new.denovo$sd_of_means[CHH.1001new.denovo$transcript %in% k27_lincs]
a2<-CHH.1001new.denovo$sd_of_means[CHH.1001new.denovo$transcript %in% k9_lincs]
boxplot(a1,a2, main="CHH methylation variation\n (28 accessions)",cex.main=1,
        col=c("pink","lightblue","#673A8E"), names=c("K27","K9"),las=1, notch = T, outline = F,cex.lab=0.8, ylab="CHH methylation level")
#add p values
a<-wilcox.test(sample(a2,399),a1)
b<-wilcox.test(sample(a2,399),a1)
c<-wilcox.test(sample(a2,399),a1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.4)
####################################################
dev.off()



# expression variation 1001G New
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Fig5_linc_boxplot_Expression_variation_1001GNEW_k27_k9.pdf",height = 3,width = 2.2)
####################################################
par(mar=c(5,4,3,2)) 
a1<-denovo2021.TPMs.genes.1001Gnew$variance_of_means[denovo2021.TPMs.genes.1001Gnew$gene %in% k27_lincs]
a2<-denovo2021.TPMs.genes.1001Gnew$variance_of_means[denovo2021.TPMs.genes.1001Gnew$gene %in% k9_lincs]
boxplot(a1,a2, main="Expression variation\n (28 accessions)",cex.main=1,
        col=c("pink","lightblue","#673A8E"), names=c("K27","K9"),las=1, notch = T, outline = F,cex.lab=0.8, ylab="coefficient of variance")
#add p values
a<-wilcox.test(sample(a2,399),a1)
b<-wilcox.test(sample(a2,399),a1)
c<-wilcox.test(sample(a2,399),a1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=5.4)
####################################################
dev.off()

# expression variation 1001G 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Fig5_linc_boxplot_Expression_variation_1001GNEWk27_k9.pdf",height = 3,width = 2.2)
####################################################
par(mar=c(5,4,3,2)) 
a1<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% k27_lincs]
a2<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% k9_lincs]
boxplot(a1,a2, main="Expression variation\n (461 accessions)",cex.main=1,
        col=c("pink","lightblue","#673A8E"), names=c("K27","K9"),las=1, notch = T, outline = F,cex.lab=0.8, ylab="coefficient of variance")
#add p values
a<-wilcox.test(sample(a2,399),a1)
b<-wilcox.test(sample(a2,399),a1)
c<-wilcox.test(sample(a2,399),a1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=18)
####################################################
dev.off()



# copy number in TAIR10 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Fig7_linc_boxplot_CNinTAIR10_k27_k9.pdf",height = 3,width = 2.2)
####################################################
par(mar=c(5,4,3,2)) 
a1<-CN_linc_27genomes$TAIR10[CN_linc_27genomes$gene %in% k27_lincs]
a2<-CN_linc_27genomes$TAIR10[CN_linc_27genomes$gene %in% k9_lincs]
boxplot(a1,a2, main="Copy number",cex.main=1,
        col=c("pink","lightblue","#673A8E"), names=c("K27","K9"),las=1, notch = T, outline = F,cex.lab=0.8, ylab="copy number in TAIR10")
#add p values
###################
a<-wilcox.test(sample(a2,399),a1)
b<-wilcox.test(sample(a2,399),a1)
c<-wilcox.test(sample(a2,399),a1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=10)
####################################################
dev.off()



# 24nt siRNAs flowers WT
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Fig7_linc_boxplot_24ntflowers_k27_k9.pdf",height = 3,width = 2.2)
####################################################
par(mar=c(5,4,3,2)) 
a1<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% k27_lincs]
a2<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% k9_lincs]
boxplot(a1,a2, main="24nt sRNA coverage\n (flowers,Col-0)",cex.main=1,
        col=c("pink","lightblue","#673A8E"), names=c("K27","K9"),las=1, notch = T, outline = F,cex.lab=0.8, ylab="RPM")
#add p values
###################
a<-wilcox.test(sample(a2,399),a1)
b<-wilcox.test(sample(a2,399),a1)
c<-wilcox.test(sample(a2,399),a1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1.2)
####################################################
dev.off()

# 24nt siRNAs e.heart WT
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Fig7_linc_boxplot_24nt_eHeart_Ranj_k27_k9.pdf",height = 3,width = 2.2)
####################################################
par(mar=c(5,4,3,2)) 
a1<-sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart[sRNA.24nt.denovo2021.RPM.Ranj$gene %in% k27_lincs]
a2<-sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart[sRNA.24nt.denovo2021.RPM.Ranj$gene %in% k9_lincs]
boxplot(a1,a2, main="24nt sRNA coverage\n (early heart,Col-0)",cex.main=1,
        col=c("pink","lightblue","#673A8E"), names=c("K27","K9"),las=1, notch = T, outline = F,cex.lab=0.8, ylab="RPM")
#add p values
###################
a<-wilcox.test(sample(a2,399),a1)
b<-wilcox.test(sample(a2,399),a1)
c<-wilcox.test(sample(a2,399),a1)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=2.2)
####################################################
dev.off()



# lincRNA sRNA targeting

# cutoff for sRNA targeting is > 90% quantile for PC gene level 

#PC highest in early heart 
quantile(sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart[sRNA.24nt.denovo2021.RPM.Ranj$gene %in% denovoPC.loci$gene],probs = seq(0.8, 1, 0.05))
80%          85%          90%          95%         100% 
0.002788510  0.004390543  0.012933982  0.109533447 50.003458937 

#PC highest in flowers
quantile(sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% denovoPC.loci$gene],probs = seq(0.8, 1, 0.05))
80%          85%          90%          95%         100% 
4.588468e-03 6.969082e-03 1.577103e-02 8.838389e-02 1.811259e+02 


#early heart 
quantile(sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart[sRNA.24nt.denovo2021.RPM.Ranj$gene %in% TE_genes.loci$gene],probs = seq(0,0.3, 0.05))
0%          5%         10%         15%         20%         25%         30% 
0.000000000 0.007265675 0.105471920 0.173522379 0.223007178 0.266412217 0.313032137 
#flowers
quantile(sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% TE_genes.loci$gene],probs = seq(0,0.3, 0.05))
0%          5%         10%         15%         20%         25%         30% 
0.000000000 0.003611994 0.017613828 0.027582618 0.038201127 0.046826086 0.056264461 


# choose cutoff of 0.03 

#lincRNAs 
#e.heart
length(sRNA.24nt.denovo2021.RPM.Ranj$gene[sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart>0.03 & sRNA.24nt.denovo2021.RPM.Ranj$gene %in% lncRNAs.intergenic.loci$gene])
#1207 lincRNAs are targeted by 24nt sRNAs in eHeart
1207*100/2246
#53.73998 % of all lincRNAs (all=2246)




length(sRNA.24nt.denovo2021.RPM.Ranj$gene[sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart>0.03 & sRNA.24nt.denovo2021.RPM.Ranj$gene %in% k9_lincs])
#644 
#98.2 % of k9 lincRNAs are targeted by 24nt in e heart
length(k9_lincs)
#656
length(sRNA.24nt.denovo2021.RPM.Ranj$gene[sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart>0.03 & sRNA.24nt.denovo2021.RPM.Ranj$gene %in% k27_lincs])
#92 
# 23.1% of k27 lincRNAs are targeted by 24nt in flowers
length(k27_lincs)
#399

#flowers
length(sRNA.24nt.denovo2021.RPM$gene[sRNA.24nt.denovo2021.RPM$X6909>0.03 & sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene])
#1131 lincRNAs are targeted by 24nt sRNAs in flowers
length(sRNA.24nt.denovo2021.RPM$gene[sRNA.24nt.denovo2021.RPM$X6909>0.03 & sRNA.24nt.denovo2021.RPM$gene %in% k9_lincs])
#595 
#90.7 % of k9 lincRNAs are targeted by 24nt in flowers
length(k9_lincs)
#656
length(sRNA.24nt.denovo2021.RPM$gene[sRNA.24nt.denovo2021.RPM$X6909>0.03 & sRNA.24nt.denovo2021.RPM$gene %in% k27_lincs])
#94 
# 23.6% of k27 lincRNAs are targeted by 24nt in flowers
length(k27_lincs)
#399







############################################
## rddm knockout from Ranj - Papareddy et al 
############################################


sRNAs24nt_disappearing_in_KOs_linc<-as.vector(
  sRNA.24nt.denovo2021.RPM.Ranj$gene [sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart>0.03 & sRNA.24nt.denovo2021.RPM.Ranj$nrpda3.eheart<0.03 & sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart>(3*sRNA.24nt.denovo2021.RPM.Ranj$nrpda3.eheart) & sRNA.24nt.denovo2021.RPM.Ranj$gene %in% lncRNAs.intergenic.loci$gene])

length(sRNAs24nt_disappearing_in_KOs_linc) #827
length(sRNAs24nt_disappearing_in_KOs_linc)*100/length(lncRNAs.intergenic.loci$gene) #36.82102% lincrnas 

a<-merge(sRNA.24nt.denovo2021.RPM.Ranj[,c("gene","WT.eheart" ,"nrpda3.eheart", "nrpda3.fb", "dcl234.fb" ,"WT.leaf", "nrpda3.leaf","dcl234.leaf")],sRNA.24nt.denovo2021.RPM [,c("gene","X6909" )],by="gene")
a$max<-apply(a[ ,c("WT.eheart" ,"nrpda3.eheart", "X6909","nrpda3.fb", "dcl234.fb" ,"WT.leaf", "nrpda3.leaf","dcl234.leaf")],1,max)
a$max1<-apply(a[ ,c("WT.eheart" ,"nrpda3.eheart", "nrpda3.fb", "dcl234.fb" )],1,max)


#24nt disapperance in pol4 and dicer mutants


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_Ranj_sRNAs_nrpdKO_eH_fb.pdf",height = 3,width =3)
par(mar=c(6,6,6,1)) 
#plot lincRNAs that are targeted in e heart by 24nt (off 0.03)
b<-a[a$gene %in% lncRNAs.intergenic.loci$gene & a$WT.eheart>0.03,c("WT.eheart" ,"nrpda3.eheart", "nrpda3.fb", "dcl234.fb" )]
pheatmap(b,main="24nt coverage of lincRNAs",cluster_cols = F,scale = "row",labels_row = F)
dev.off()

b<-a[a$gene %in% lncRNAs.intergenic.loci$gene & a$X6909>0.03,c("X6909","nrpda3.fb", "dcl234.fb","WT.leaf", "nrpda3.leaf" )]
pheatmap(b,main="24nt coverage of lincRNAs",cluster_cols = F,scale = "row",labels_row = F)


#Fig7H
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_Ranj_sRNAs_nrpdKO_eH_3reps.pdf",height = 3,width =3)
par(mar=c(6,6,6,1)) 
#plot lincRNAs that are targeted in e heart by 24nt (cutoff 0.03)
#plot lincRNAs that are targeted in the WT embryo  sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart>0.03
b<-sRNA.24nt.denovo2021.RPM.Ranj[sRNA.24nt.denovo2021.RPM.Ranj$gene %in% lncRNAs.intergenic.loci$gene & sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart>0.03,c("WT.eheart.r1" ,"WT.eheart.r2", "WT.eheart.r3", "nrpda3.eheart.r1", "nrpda3.eheart.r3", "nrpda3.eheart.r3" )]
length(b$WT.eheart.r1) # 1207 lincRNAs
pheatmap(b,main="24nt coverage of lincRNAs",cluster_cols = F,scale = "row",labels_row = F)
dev.off()

# number of lincRNAs targeted in e heart by 24nt (cutoff 0.03)
b<-sRNA.24nt.denovo2021.RPM.Ranj[sRNA.24nt.denovo2021.RPM.Ranj$gene %in% lncRNAs.intergenic.loci$gene & sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart.r1>0.03,c("WT.eheart.r1" , "nrpda3.eheart.r1", "nrpda3.eheart.r3", "nrpda3.eheart.r3" )]
length(b$WT.eheart.r1) # 1207 lincRNAs
# how many lose targeting? 


# loss of 24nt targeting in flowers 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_Ranj_sRNAs_nrpdKO_flowers_SUPPL.pdf",height = 4,width =3)
par(mar=c(6,6,6,1)) 
#plot lincRNAs that are targeted in flowers by 24nt (cutoff 0.03)
#plot lincRNAs that are targeted in the WT flowers  sRNA.24nt.denovo2021.RPM$X6909>0.03
a<-merge(sRNA.24nt.denovo2021.RPM.Ranj[,c("gene","WT.eheart" ,"nrpda3.eheart", "nrpda3.fb.r1","nrpda3.fb.r2","nrpda3.fb.r3","WT.leaf.r1" ,"WT.leaf.r2", "WT.leaf.r3", "nrpda3.leaf.r1" ,"nrpda3.leaf.r2", "nrpda3.leaf.r3")],sRNA.24nt.denovo2021.RPM [,c("gene","X6909" )],by="gene")
a$max<-apply(a[ ,c("X6909","nrpda3.fb.r1","nrpda3.fb.r2","nrpda3.fb.r3","WT.leaf.r1" ,"WT.leaf.r2", "WT.leaf.r3", "nrpda3.leaf.r1" ,"nrpda3.leaf.r2", "nrpda3.leaf.r3")],1,max)
b<-a[a$gene %in% lncRNAs.intergenic.loci$gene & a$X6909>0.03,c("X6909","nrpda3.fb.r1","nrpda3.fb.r2","nrpda3.fb.r3")]
pheatmap(b,main="24nt coverage of lincRNAs",cluster_cols = F,scale = "row",labels_row = F,labels_col = c("WT.flowers(1rep)","nrpd1a.flowers.rep1","nrpd1a.flowers.rep2","nrpd1a.flowers.rep3"))
dev.off()
length (b$X6909) # 1131 lincRNAs targeted in flowers (data from this study)
# how many lincRNAs  lose targeting? 
length (b$X6909[apply(b[c("nrpda3.fb.r1","nrpda3.fb.r2","nrpda3.fb.r3")],1,mean)<0.03]) # 1062 lincRNAs lose targeting in flowers (data from Papareddy et al)
#94% lose targeting


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_Ranj_sRNAs_nrpdKO_flowers_SUPPL.pdf",height = 4,width =3)
par(mar=c(6,6,6,1)) 
#plot lincRNAs that are targeted in flowers by 24nt (cutoff 0.03)
#plot lincRNAs that are targeted in the WT flowers  sRNA.24nt.denovo2021.RPM$X6909>0.03
a<-merge(sRNA.24nt.denovo2021.RPM.Ranj[,c("gene","WT.eheart" ,"nrpda3.eheart", "nrpda3.fb.r1","nrpda3.fb.r2","nrpda3.fb.r3","WT.leaf.r1" ,"WT.leaf.r2", "WT.leaf.r3", "nrpda3.leaf.r1" ,"nrpda3.leaf.r2", "nrpda3.leaf.r3")],sRNA.24nt.denovo2021.RPM [,c("gene","X6909" )],by="gene")
a$max<-apply(a[ ,c("X6909","nrpda3.fb.r1","nrpda3.fb.r2","nrpda3.fb.r3","WT.leaf.r1" ,"WT.leaf.r2", "WT.leaf.r3", "nrpda3.leaf.r1" ,"nrpda3.leaf.r2", "nrpda3.leaf.r3")],1,max)
b<-a[a$gene %in% lncRNAs.intergenic.loci$gene & a$X6909>0.03,c("X6909","nrpda3.fb.r1","nrpda3.fb.r2","nrpda3.fb.r3","WT.eheart","nrpda3.eheart")]
pheatmap(b,main="24nt coverage of lincRNAs",cluster_cols = F,scale = "row",labels_row = F,labels_col = c("WT.flowers(1rep)","nrpd1a.flowers.rep1","nrpd1a.flowers.rep2","nrpd1a.flowers.rep3","WT.embryo(3reps)","nrpd1a.embryo(3reps)"))
dev.off()

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_Ranj_sRNAs_nrpdKO.pdf",height = 3,width =3)
par(mar=c(6,6,6,1)) 
pheatmap(sRNA.24nt.denovo2021.RPM.Ranj[
  apply(sRNA.24nt.denovo2021.RPM.Ranj[,c("WT.eheart" ,"nrpda3.eheart", "nrpda3.fb", "dcl234.fb" ,"WT.leaf", "nrpda3.leaf","dcl234.leaf")],1,max)>0.03,
  c("WT.eheart" ,"nrpda3.eheart", "nrpda3.fb", "dcl234.fb" ,"WT.leaf", "nrpda3.leaf","dcl234.leaf")],main="lincRNAs\n ddm1 KO rosette",cluster_cols = F,scale = "row",labels_row = F)
dev.off()

a<-sRNA.24nt.denovo2021.RPM.Ranj[apply(sRNA.24nt.denovo2021.RPM.Ranj[,c("WT.eheart" ,"nrpda3.eheart", "nrpda3.fb", "dcl234.fb" ,"WT.leaf", "nrpda3.leaf","dcl234.leaf")],1,max)>0.05 & sRNA.24nt.denovo2021.RPM.Ranj$gene %in% lncRNAs.intergenic.loci$gene,c("WT.eheart" ,"nrpda3.eheart", "nrpda3.fb", "dcl234.fb" ,"WT.leaf", "nrpda3.leaf","dcl234.leaf")]

pheatmap(a,main="lincRNAs\n ddm1 KO rosette",cluster_cols = F,scale = "row",labels_row = F)







#################################################################################
### silencing by DDM 1 Bhagy's data Osakabe et al https://www.nature.com/articles/s41556-021-00658-1 
#################################################################################

#denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1
#rosette

#lincRNAs 

lincRNAs_reexpressed_in_ddm1_rosette<-as.vector(
  linc_Bhagy$gene [linc_Bhagy$ddm1>0.5 & linc_Bhagy$WT<0.5 & linc_Bhagy$ddm1>(3*linc_Bhagy$WT)])

length(lincRNAs_reexpressed_in_ddm1_rosette) #147
length(lincRNAs_reexpressed_in_ddm1_rosette)*100/length(lncRNAs.intergenic.loci$gene) #6.544969% lincrnas are reexpressed

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_Bhagy_ddm1_lincRNA_reexpressed.pdf",height = 3,width =3)
par(mar=c(6,6,6,1)) 
pheatmap(linc_Bhagy[lincRNAs_reexpressed_in_ddm1_rosette,c("Col_WT.rep1" ,"Col_WT.rep2", "Col_WT.rep3", "Col_ddm1.rep1" ,"Col_ddm1.rep2", "Col_ddm1.rep3")],main="lincRNAs\n ddm1 KO rosette",cluster_cols = F,scale = "row",labels_row = F)
dev.off()


#PC
PC_reexpressed_in_ddm1_rosette<-as.vector(
  pc_Bhagy$gene [    pc_Bhagy$ddm1>0.5    & pc_Bhagy$WT<0.5     & pc_Bhagy$ddm1>(3*pc_Bhagy$WT)])
length(PC_reexpressed_in_ddm1_rosette) #729
length(PC_reexpressed_in_ddm1_rosette)*100/length(denovoPC.loci$gene) #3.079067% PC genes are reexpressed
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_Bhagy_ddm1_PC_reexpressed.pdf",height = 3,width =3)
par(mar=c(6,6,6,1)) 
pheatmap(pc_Bhagy[PC_reexpressed_in_ddm1_rosette,c("Col_WT.rep1" ,"Col_WT.rep2", "Col_WT.rep3", "Col_ddm1.rep1" ,"Col_ddm1.rep2", "Col_ddm1.rep3")],main="PC genes\n ddm1 KO rosette",cluster_cols = F,scale = "row",labels_row = F)
dev.off()

#AS
AS_reexpressed_in_ddm1_rosette<-as.vector(
  as_Bhagy$gene [ as_Bhagy$ddm1>0.5  & as_Bhagy$WT<0.5    & as_Bhagy$ddm1>(3*as_Bhagy$WT)])
length(AS_reexpressed_in_ddm1_rosette) # 323
length(AS_reexpressed_in_ddm1_rosette)*100/length(lncRNAs.antisense.loci$gene) #3.941428% AS rnas are reexpressed
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_Bhagy_ddm1_AS_reexpressed.pdf",height = 3,width =3)
par(mar=c(6,6,6,1)) 
pheatmap(as_Bhagy[AS_reexpressed_in_ddm1_rosette,c("Col_WT.rep1" ,"Col_WT.rep2", "Col_WT.rep3", "Col_ddm1.rep1" ,"Col_ddm1.rep2", "Col_ddm1.rep3")],main="AS lncRNAs\n ddm1 KO rosette",cluster_cols = F,scale = "row",labels_row = F)
dev.off()

#TE genes 
TE_reexpressed_in_ddm1_rosette<-as.vector(
  te_Bhagy$gene [    te_Bhagy$ddm1>0.5  & te_Bhagy$WT<0.5   & te_Bhagy$ddm1>(3*te_Bhagy$WT)])
length(TE_reexpressed_in_ddm1_rosette) #  800
length(TE_reexpressed_in_ddm1_rosette)*100/length(TE_genes.loci$gene) #37.55869%  TE genes are reexpressed
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_Bhagy_ddm1_TE_reexpressed.pdf",height = 3,width =3)
par(mar=c(6,6,6,1)) 
pheatmap(te_Bhagy[TE_reexpressed_in_ddm1_rosette,c("Col_WT.rep1" ,"Col_WT.rep2", "Col_WT.rep3", "Col_ddm1.rep1" ,"Col_ddm1.rep2", "Col_ddm1.rep3")],main="TE genes\n ddm1 KO rosette",cluster_cols = F,scale = "row",labels_row = F)
dev.off()


lincRNAs_NOTreexpressed_in_ddm1_rosette<-as.vector(lncRNAs.intergenic.loci$gene[!(lncRNAs.intergenic.loci$gene %in% lincRNAs_reexpressed_in_ddm1_rosette)])
length(lincRNAs_NOTreexpressed_in_ddm1_rosette)#2099


#K9 or k 27 lincRNAs?
length(lincRNAs_reexpressed_in_ddm1_rosette) #147
length(intersect(lincRNAs_reexpressed_in_ddm1_rosette,k27_lincs)) #6
length(intersect(lincRNAs_reexpressed_in_ddm1_rosette,k9_lincs)) #90

length(lincRNAs_NOTreexpressed_in_ddm1_rosette) #2097
length(intersect(lincRNAs_NOTreexpressed_in_ddm1_rosette,k27_lincs)) #406
length(intersect(lincRNAs_NOTreexpressed_in_ddm1_rosette,k9_lincs)) # 575

#TE content 
length(intersect(lincRNAs_reexpressed_in_ddm1_rosette,lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10==0]))#47 
length(intersect(lincRNAs_reexpressed_in_ddm1_rosette,lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10>0])) #100

length(intersect(lincRNAs_NOTreexpressed_in_ddm1_rosette,lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10==0])) #1023
length(intersect(lincRNAs_NOTreexpressed_in_ddm1_rosette,lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10>0])) #1076


# copy number 
length(intersect(lincRNAs_reexpressed_in_ddm1_rosette,lincs_1copy))#96 
length(intersect(lincRNAs_reexpressed_in_ddm1_rosette,lincs_2copies))#20 
length(intersect(lincRNAs_reexpressed_in_ddm1_rosette,lincs3_10_copies))#27
length(intersect(lincRNAs_reexpressed_in_ddm1_rosette,lincs_more10copies))#4

length(intersect(lincRNAs_NOTreexpressed_in_ddm1_rosette,lincs_1copy))#1529 
length(intersect(lincRNAs_NOTreexpressed_in_ddm1_rosette,lincs_2copies))#192 
length(intersect(lincRNAs_NOTreexpressed_in_ddm1_rosette,lincs3_10_copies))#202
length(intersect(lincRNAs_NOTreexpressed_in_ddm1_rosette,lincs_more10copies))#176





####################################################
# data from Vu Nguyen et al - ddm1 knockout and heat stress
#######################################################

#lincRNAs 

lincRNAs_reexpressed_in_ddm1_stem_7d<-as.vector(
  linc_Vu$gene [( linc_Vu$stem_mock_7D_ddm1>0.5 & linc_Vu$stem_mock_7D_WT<0.5
                  & linc_Vu$stem_mock_7D_ddm1>3*linc_Vu$stem_mock_7D_WT) | ( linc_Vu$stem_heat_7D_ddm1>0.5 & linc_Vu$stem_heat_7D_WT<0.5 & linc_Vu$stem_heat_7D_ddm1>3*linc_Vu$stem_heat_7D_WT) ] )

length(lincRNAs_reexpressed_in_ddm1_stem_7d) #410
length(lincRNAs_reexpressed_in_ddm1_stem_7d)*100/length(lncRNAs.intergenic.loci$gene) #18.25467% lincrnas are reexpressed


lincRNAs_NOTreexpressed_in_ddm1_stem_7d<-as.vector(lncRNAs.intergenic.loci$gene[!(lncRNAs.intergenic.loci$gene %in% lincRNAs_reexpressed_in_ddm1_stem_7d)])
length(lincRNAs_NOTreexpressed_in_ddm1_stem_7d)#1836

#linc 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_Vu_ddm1_stem_7D_lincRNA_reexpressed.pdf",height = 4,width =3)
par(mar=c(6,6,6,1)) 
pheatmap(linc_Vu[lincRNAs_reexpressed_in_ddm1_stem_7d,c("stem_mock_7D_WT","stem_mock_7D_ddm1","stem_heat_7D_WT", "stem_heat_7D_ddm1" )],main="lincRNAs\n ddm1 KO stem cells",cluster_cols = F,scale = "row",labels_row = F)
dev.off()

#AS
as_reexpressed_in_ddm1_stem_7d<-as.vector(
  as_Vu$gene [
    apply(as_Vu[,c("stem_mock_7D_ddm1","stem_heat_7D_ddm1")],1,max)>0.5
    & apply(as_Vu[,c("stem_mock_7D_WT","stem_heat_7D_WT")],1,max)<0.5 
    & apply(as_Vu[,c("stem_mock_7D_ddm1","stem_heat_7D_ddm1")],1,max)>3* apply(as_Vu[,c("stem_mock_7D_WT","stem_heat_7D_WT")],1,max)])
length(as_reexpressed_in_ddm1_stem_7d) #216
length(as_reexpressed_in_ddm1_stem_7d)*100/length(lncRNAs.antisense.loci$gene) #2.635754% ASrnas are reexpressed

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_Vu_ddm1_stem_7D_AS_reexpressed.pdf",height = 4,width =3)
par(mar=c(6,6,6,1)) 
pheatmap(as_Vu[as_reexpressed_in_ddm1_stem_7d,c("stem_mock_7D_WT","stem_mock_7D_ddm1","stem_heat_7D_WT", "stem_heat_7D_ddm1")],main="AS lncRNAs\n ddm1 KO stem cells",cluster_cols = F,scale = "row",labels_row = F)
dev.off()

#PC
pc_reexpressed_in_ddm1_stem_7d<-as.vector(
  pc_Vu$gene [
    apply(pc_Vu[,c("stem_mock_7D_ddm1","stem_heat_7D_ddm1")],1,max)>0.5
    & apply(pc_Vu[,c("stem_mock_7D_WT","stem_heat_7D_WT")],1,max)<0.5 
    & apply(pc_Vu[,c("stem_mock_7D_ddm1","stem_heat_7D_ddm1")],1,max)>3* apply(pc_Vu[,c("stem_mock_7D_WT","stem_heat_7D_WT")],1,max)])
length(pc_reexpressed_in_ddm1_stem_7d) #594
length(pc_reexpressed_in_ddm1_stem_7d)*100/length(denovoPC.loci$gene) #2.50887% PC are reexpressed

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_Vu_ddm1_stem_7D_PC_reexpressed.pdf",height = 4,width =3)
par(mar=c(6,6,6,1)) 
pheatmap(pc_Vu[pc_reexpressed_in_ddm1_stem_7d,c("stem_mock_7D_WT","stem_mock_7D_ddm1","stem_heat_7D_WT", "stem_heat_7D_ddm1")],main="PC genes\n ddm1 KO stem cells",cluster_cols = F,scale = "row",labels_row = F)
dev.off()



#TE genes 
te_reexpressed_in_ddm1_stem_7d<-as.vector(
  te_Vu$gene [
    apply(te_Vu[,c("stem_mock_7D_ddm1","stem_heat_7D_ddm1")],1,max)>0.5
    & apply(te_Vu[,c("stem_mock_7D_WT","stem_heat_7D_WT")],1,max)<0.5 
    & apply(te_Vu[,c("stem_mock_7D_ddm1","stem_heat_7D_ddm1")],1,max)>3* apply(te_Vu[,c("stem_mock_7D_WT","stem_heat_7D_WT")],1,max)])
length(te_reexpressed_in_ddm1_stem_7d) #1160
length(te_reexpressed_in_ddm1_stem_7d)*100/length(TE_genes.loci$gene) #54.46009% TE genes  are reexpressed

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_Vu_ddm1_stem_7D_TE_reexpressed.pdf",height = 4,width =3)
par(mar=c(6,6,6,1)) 
pheatmap(te_Vu[te_reexpressed_in_ddm1_stem_7d,c("stem_mock_7D_WT","stem_mock_7D_ddm1","stem_heat_7D_WT", "stem_heat_7D_ddm1")],main="TE genes\n ddm1 KO stem cells",cluster_cols = F,scale = "row",labels_row = F)
dev.off()


#K9 or k 27 lincRNAs?
length(lincRNAs_reexpressed_in_ddm1_stem_7d) #410
length(intersect(lincRNAs_reexpressed_in_ddm1_stem_7d,k27_lincs)) #24
length(intersect(lincRNAs_reexpressed_in_ddm1_stem_7d,k9_lincs)) #269

length(lincRNAs_NOTreexpressed_in_ddm1_stem_7d) #1836
length(intersect(lincRNAs_NOTreexpressed_in_ddm1_stem_7d,k27_lincs)) #375
length(intersect(lincRNAs_NOTreexpressed_in_ddm1_stem_7d,k9_lincs)) #387

#TE content 
length(intersect(lincRNAs_reexpressed_in_ddm1_stem_7d,lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10==0]))#131 
length(intersect(lincRNAs_reexpressed_in_ddm1_stem_7d,lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10>0])) #279

length(intersect(lincRNAs_NOTreexpressed_in_ddm1_stem_7d,lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10==0])) #939
length(intersect(lincRNAs_NOTreexpressed_in_ddm1_stem_7d,lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10>0])) #897


# copy number 
length(intersect(lincRNAs_reexpressed_in_ddm1_stem_7d,lincs_1copy))#226 
length(intersect(lincRNAs_reexpressed_in_ddm1_stem_7d,lincs_2copies))#47 
length(intersect(lincRNAs_reexpressed_in_ddm1_stem_7d,lincs3_10_copies))#81
length(intersect(lincRNAs_reexpressed_in_ddm1_stem_7d,lincs_more10copies))#56

length(intersect(lincRNAs_NOTreexpressed_in_ddm1_stem_7d,lincs_1copy))#1399 
length(intersect(lincRNAs_NOTreexpressed_in_ddm1_stem_7d,lincs_2copies))#165 
length(intersect(lincRNAs_NOTreexpressed_in_ddm1_stem_7d,lincs3_10_copies))#148
length(intersect(lincRNAs_NOTreexpressed_in_ddm1_stem_7d,lincs_more10copies))#124




#################################################################################
### silencing by methylases He et al https://www.nature.com/articles/s41467-022-28940-2
#################################################################################

denovo_Oct2021.TPMs.genes.He_et_al <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/denovo_Oct2021.TPMs.genes.He_et_al.bed")
rownames(denovo_Oct2021.TPMs.genes.He_et_al)<-denovo_Oct2021.TPMs.genes.He_et_al$gene

denovo_Oct2021.TPMs.genes.He_et_al$ddcc<-apply(denovo_Oct2021.TPMs.genes.He_et_al[,2:4],1,mean)
denovo_Oct2021.TPMs.genes.He_et_al$mddcc<-apply(denovo_Oct2021.TPMs.genes.He_et_al[,5:7],1,mean)
denovo_Oct2021.TPMs.genes.He_et_al$met1<-apply(denovo_Oct2021.TPMs.genes.He_et_al[,8:10],1,mean)
denovo_Oct2021.TPMs.genes.He_et_al$wt<-apply(denovo_Oct2021.TPMs.genes.He_et_al[,11:13],1,mean)
denovo_Oct2021.TPMs.genes.He_et_al$wt_5w<-apply(denovo_Oct2021.TPMs.genes.He_et_al[,14:16],1,mean)

library (pheatmap)

#wt - 2 weeks old 
#ddcc - 2 weeks old 
#met1 - 2 weeks old 
#wt_5w - 5 weeks old 
#mddcc - 5 weeks old 


#AS RNAs
a<-denovo_Oct2021.TPMs.genes.He_et_al[denovo_Oct2021.TPMs.genes.He_et_al$gene %in% lncRNAs.antisense.loci$gene, ]

library (pheatmap)


############### linc lncRNAs reexpressed upon methylation deletion 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_He_et_al_lincRNA_reexpressed.pdf",height = 3,width =3)
###################
par(mar=c(6,6,6,1)) 
a<-denovo_Oct2021.TPMs.genes.He_et_al[denovo_Oct2021.TPMs.genes.He_et_al$gene %in% lncRNAs.intergenic.loci$gene, ]
pheatmap(a[
  (apply(a[,c("wt","ddcc","met1")],1,max)>0.5 & 
     apply(a[,c("ddcc","met1")],1,max)>3*a$wt & a$wt<0.5) |
    (apply(a[,c("wt_5w","mddcc")],1,max)>0.5 & a$mddcc>3*a$wt_5w & a$wt_5w<0.5) ,
  c("wt","ddcc","met1","wt_5w","mddcc")],scale="row",labels_row = F,cluster_cols = F,main="lincRNAs\n He et al data")
len=length(a[(apply(a[,c("wt","ddcc","met1")],1,max)>0.5 & apply(a[,c("ddcc","met1")],1,max)>3*a$wt & a$wt<0.5) |
               (apply(a[,c("wt_5w","mddcc")],1,max)>0.5 & a$mddcc>3*a$wt_5w& a$wt_5w<0.5),1])
#274
274*100/2246 #12.19947 lincRNAs reexpressed
####################
dev.off()

############### AS lncRNAs reexpressed upon methylation deletion 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_He_et_al_AS_RNA_reexpressed.pdf",height = 3,width =3)
##################################
par(mar=c(6,6,6,1)) 
a<-denovo_Oct2021.TPMs.genes.He_et_al[denovo_Oct2021.TPMs.genes.He_et_al$gene %in% lncRNAs.antisense.loci$gene, ]
pheatmap(a[
  (apply(a[,c("wt","ddcc","met1")],1,max)>0.5 & 
     apply(a[,c("ddcc","met1")],1,max)>3*a$wt& a$wt<0.5) |
    (apply(a[,c("wt_5w","mddcc")],1,max)>0.5 & a$mddcc>3*a$wt_5w& a$wt_5w<0.5) ,
  c("wt","ddcc","met1","wt_5w","mddcc")],scale="row",labels_row = F,cluster_cols = F,main="AS lncRNAs\n He et al data")
len=length(a[(apply(a[,c("wt","ddcc","met1")],1,max)>0.5 & apply(a[,c("ddcc","met1")],1,max)>3*a$wt& a$wt<0.5) |
               (apply(a[,c("wt_5w","mddcc")],1,max)>0.5 & a$mddcc>3*a$wt_5w& a$wt_5w<0.5),1])
#395
##################################
dev.off()

############### PC genes reexpressed upon methylation deletion 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_He_et_al_PC_reexpressed.pdf",height = 3,width =3)
##################################
par(mar=c(6,6,6,1)) 
a<-denovo_Oct2021.TPMs.genes.He_et_al[denovo_Oct2021.TPMs.genes.He_et_al$gene %in% sample(denovoPC.loci$gene,9000), ]
pheatmap(a[
  (apply(a[,c("wt","ddcc","met1")],1,max)>0.5 & 
     apply(a[,c("ddcc","met1")],1,max)>3*a$wt & a$wt<0.5)|
    (apply(a[,c("wt_5w","mddcc")],1,max)>0.5 & a$mddcc>3*a$wt_5w & a$wt_5w<0.5) ,
  c("wt","ddcc","met1","wt_5w","mddcc")],scale="row",labels_row = F,cluster_cols = F,main="PC genes\n He et al data")
len=length(a[(apply(a[,c("wt","ddcc","met1")],1,max)>0.5 & apply(a[,c("ddcc","met1")],1,max)>3*a$wt) |
               (apply(a[,c("wt_5w","mddcc")],1,max)>0.5 & a$mddcc>3*a$wt_5w),1])

length(a[  (apply(a[,c("wt","ddcc","met1")],1,max)>0.5 &  apply(a[,c("ddcc","met1")],1,max)>3*a$wt & a$wt<0.5)|    (apply(a[,c("wt_5w","mddcc")],1,max)>0.5 & a$mddcc>3*a$wt_5w & a$wt_5w<0.5) , 1])
#541
##################################
dev.off()


############### TE genes reexpressed upon methylation deletion 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_He_et_al_TEgenes_reexpressed.pdf",height = 3,width =3)
##################################
par(mar=c(6,6,6,1)) 
a<-denovo_Oct2021.TPMs.genes.He_et_al[denovo_Oct2021.TPMs.genes.He_et_al$gene %in% TE_genes.loci$gene, ]
pheatmap(a[
  (apply(a[,c("wt","ddcc","met1")],1,max)>0.5 & 
     apply(a[,c("ddcc","met1")],1,max)>3*a$wt& a$wt<0.5) |
    (apply(a[,c("wt_5w","mddcc")],1,max)>0.5 & a$mddcc>3*a$wt_5w& a$wt_5w<0.5) ,
  c("wt","ddcc","met1","wt_5w","mddcc")],scale="row",labels_row = F,cluster_cols = F,main="TE genes\n He et al data")
len=length(a[(apply(a[,c("wt","ddcc","met1")],1,max)>0.5 & apply(a[,c("ddcc","met1")],1,max)>3*a$wt& a$wt<0.5) |
               (apply(a[,c("wt_5w","mddcc")],1,max)>0.5 & a$mddcc>3*a$wt_5w& a$wt_5w<0.5),1])
#1132
##################################
dev.off()



# what is the K9/K27 distribution? TE content? copy number? 
a<-denovo_Oct2021.TPMs.genes.He_et_al[denovo_Oct2021.TPMs.genes.He_et_al$gene %in% lncRNAs.intergenic.loci$gene, ]
lincRNAs_reexpressed_in_methmutants<-as.vector(a$gene[(apply(a[,c("wt","ddcc","met1")],1,max)>0.5  & 
                                                         apply(a[,c("ddcc","met1")],1,max)>3*a$wt & a$wt<0.5) |
                                                        (apply(a[,c("wt_5w","mddcc")],1,max)>0.5 & a$mddcc>3*a$wt_5w & a$wt_5w<0.5)])
length(lincRNAs_reexpressed_in_methmutants)
#274
lincRNAs_NOTreexpressed_in_methmutants<-lncRNAs.intergenic.loci$gene[!(lncRNAs.intergenic.loci$gene %in% lincRNAs_reexpressed_in_methmutants) ]
length(lincRNAs_NOTreexpressed_in_methmutants)
#1972

#K9 or k 27 lincRNAs?
length(lincRNAs_reexpressed_in_methmutants) #274
length(intersect(lincRNAs_reexpressed_in_methmutants,k27_lincs)) #20
length(intersect(lincRNAs_reexpressed_in_methmutants,k9_lincs)) #173

length(lincRNAs_NOTreexpressed_in_methmutants) #1972
length(intersect(lincRNAs_NOTreexpressed_in_methmutants,k27_lincs)) #392
length(intersect(lincRNAs_NOTreexpressed_in_methmutants,k9_lincs)) # 492

#TE content 
length(intersect(lincRNAs_NOTreexpressed_in_methmutants,lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10==0]))#982 
length(intersect(lincRNAs_NOTreexpressed_in_methmutants,lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10>0])) #990

length(intersect(lincRNAs_reexpressed_in_methmutants,lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10==0])) #88
length(intersect(lincRNAs_reexpressed_in_methmutants,lincRNAs_TE_coverage_anystrand$gene[lincRNAs_TE_coverage_anystrand$TAIR10>0])) #186






###################################################
# number of lincRNAs that can be expressed in col-0 
######################################################

# all datasets
denovo_Oct2021.TPMs.genes.He_et_al 
denovo2021.TPMs.genes.ERACAPS
denovo2021.TPMs.genes.Cortijo
denovo2021.TPMs.genes.1001G
denovo2021.TPMs.genes.1001Gnew
denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1
denovo_Oct2021.TPMs.genes.Vu.PEsamples


a<- merge(denovo2021.TPMs.genes.ERACAPS[,c("gene","S.6909","R.6909","F.6909","P.6909")],denovo_Oct2021.TPMs.genes.He_et_al[,c("gene","ddcc","mddcc","met1","wt", "wt_5w")])
a<-merge(a,denovo2021.TPMs.genes.1001G[,c("gene","X6909")])
a<-merge(a,denovo2021.TPMs.genes.1001Gnew[,c("gene","mean.6909")])
a<-merge(a,denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1[,c("gene","WT","ddm1")])
a<-merge(a,denovo_Oct2021.TPMs.genes.Vu.PEsamples[,c("gene","stem_mock_7D_WT","stem_mock_7D_ddm1","stem_heat_7D_WT","stem_heat_7D_ddm1")])
a<-merge(a,denovo2021.TPMs.genes.Cortijo[,c("gene","mean")])

#,"nonstem_mock_7D_WT","nonstem_heat_7D_WT"
library(pheatmap)
b<-a[a$gene %in% lncRNAs.intergenic.loci$gene,]
b$max<-apply(b[,2:19],1,max)

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/heatmap_lincRNAs_expressed_across_all_samples_Col0.pdf",height = 5,width =5)
pheatmap (b[b$max>0.5,2:19],scale = "row",labels_row = F,clustering_method =  "median")
dev.off()


length(b[b$max>0.5,1])
#1194
#0.5316118





###################################################################################


# do lincRNAs with several pieces have the same TE type? 


a<-lincRNAs_TE_coverage.TAIR10[,c(1,4)]
head (a)
gene                        TE_type
1 CUFF_NC.13                    RC_Helitron
2 CUFF_NC.31                      SINE_LINE
3 CUFF_NC.32           DNA_MuDR,RC_Helitron
4 CUFF_NC.32 DNA_MuDR,DNA_other,RC_Helitron
5 CUFF_NC.32 DNA_MuDR,DNA_other,RC_Helitron


foo <- data.frame(do.call('rbind', strsplit(as.character(a$TE_type),',',fixed=TRUE))) 
a<-cbind(a,foo)
head (a)
gene                        TE_type          X1          X2          X3          X4          X5
1 CUFF_NC.13                    RC_Helitron RC_Helitron RC_Helitron RC_Helitron RC_Helitron RC_Helitron
2 CUFF_NC.31                      SINE_LINE   SINE_LINE   SINE_LINE   SINE_LINE   SINE_LINE   SINE_LINE
3 CUFF_NC.32           DNA_MuDR,RC_Helitron    DNA_MuDR RC_Helitron    DNA_MuDR RC_Helitron    DNA_MuDR
4 CUFF_NC.32 DNA_MuDR,DNA_other,RC_Helitron    DNA_MuDR   DNA_other RC_Helitron    DNA_MuDR   DNA_other

nr<-as.data.frame(t(table (a$gene)))
head (nr)
Var1          Var2 Freq
1    A CUFF_NC.10002    1
2    A CUFF_NC.10004    1
3    A CUFF_NC.10005    1
4    A CUFF_NC.10010    1
5    A CUFF_NC.10013    1
6    A CUFF_NC.10014    2

a<-a[a$gene %in% nr$Var2[nr$Freq>1],]
head (a)
gene                        TE_type          X1          X2          X3          X4          X5
3 CUFF_NC.32           DNA_MuDR,RC_Helitron    DNA_MuDR RC_Helitron    DNA_MuDR RC_Helitron    DNA_MuDR
4 CUFF_NC.32 DNA_MuDR,DNA_other,RC_Helitron    DNA_MuDR   DNA_other RC_Helitron    DNA_MuDR   DNA_other
5 CUFF_NC.32 DNA_MuDR,DNA_other,RC_Helitron    DNA_MuDR   DNA_other RC_Helitron    DNA_MuDR   DNA_other

loci<-as.data.frame(unique(a$gene))
loci$mix<-"tba"
a$type<-"a"
loci$mix_Ntypes<-0
loci$mix_DNA_RNA<-"a"
for (i in 1:740){
  linc=as.character(loci$`unique(a$gene)`[i])
  b<-a[a$gene==linc,]
  df2 <- as.vector(as.matrix(b[,3:7]))
  loci$mix_Ntypes[i]<-length(unique(df2))
  b <- data.frame(lapply(b, function(x) {   gsub("DNA_MuDR", "DNA", x)}))
  b <- data.frame(lapply(b, function(x) {   gsub("RC_Helitron", "DNA", x)}))
  b <- data.frame(lapply(b, function(x) {   gsub("DNA_other", "DNA", x)}))
  b <- data.frame(lapply(b, function(x) {   gsub("LTR_Gypsy", "RNA", x)}))
  b <- data.frame(lapply(b, function(x) {   gsub("LTR_Copia", "RNA", x)}))
  b <- data.frame(lapply(b, function(x) {   gsub("SINE_LINE ", "RNA", x)}))
  loci$mix_DNA_RNA[i]<-any(b=="DNA") & any(b=="RNA") 
}


length(loci$mix_DNA_RNA[loci$mix_DNA_RNA=="TRUE"])
#198  - 27% 
length(loci$mix_DNA_RNA) #740 
write.table(a,"Z:/01_POSTDOC/te_pieces.txt",col.names = F,row.names = F,sep = "\t",quote = F)
length(loci$mix_Ntypes[loci$mix_Ntypes>1])
#548    - 74%


################################################################################
# Figure 8   ###################################################################
################################################################################


################################################################################
# are lincRNAs silenced by the same mechanisms as the TE pieces inside them ? 
###############################################################################

linc_noTE<-as.vector(unique(lncRNAs.intergenic.loci$gene [!(lncRNAs.intergenic.loci$gene %in%lincRNAs_TE_coverage.TAIR10$gene) ] ) )
length(linc_noTE)#1070 
linc_Helitron<-as.vector(unique(lincRNAs_TE_coverage.TAIR10$gene [grep("RC_Helitron",lincRNAs_TE_coverage.TAIR10$TE_type)] ) )
length(linc_Helitron) #763

a<-lincRNAs_TE_coverage.TAIR10 [grep("RC_Helitron",lincRNAs_TE_coverage.TAIR10$TE_type),]
linc_Helitron_only<-as.vector(unique(a$gene [grep(",",a$TE_type,invert = T)] ) )
length(linc_Helitron_only) #586







linc_DNA<-as.vector(unique(lincRNAs_TE_coverage.TAIR10$gene [grep("DNA",lincRNAs_TE_coverage.TAIR10$TE_type)] ) )
length(linc_DNA) #682
linc_DNAother<-as.vector(unique(lincRNAs_TE_coverage.TAIR10$gene [grep("DNA_other",lincRNAs_TE_coverage.TAIR10$TE_type)] ) )
length(linc_DNAother) #333

linc_DNA_mudr<-as.vector(unique(lincRNAs_TE_coverage.TAIR10$gene [grep("DNA_MuDR",lincRNAs_TE_coverage.TAIR10$TE_type)] ) )
length(linc_DNA_mudr) #499

a<-lincRNAs_TE_coverage.TAIR10 [grep("DNA_MuDR",lincRNAs_TE_coverage.TAIR10$TE_type),]
linc_DNA_mudr_only<-as.vector(unique(a$gene [grep(",",a$TE_type,invert = T)] ) )
length(linc_DNA_mudr_only) #323


linc_LTR_copia<-as.vector(unique(lincRNAs_TE_coverage.TAIR10$gene [grep("LTR_Copia",lincRNAs_TE_coverage.TAIR10$TE_type)] ) )
length(linc_LTR_copia)#130
a<-lincRNAs_TE_coverage.TAIR10 [grep("LTR_Copia",lincRNAs_TE_coverage.TAIR10$TE_type),]
linc_LTR_copia_only<-as.vector(unique(a$gene [grep(",",a$TE_type,invert = T)] ) )
length(linc_LTR_copia_only) #79


linc_LTR_gypsy<-as.vector(unique(lincRNAs_TE_coverage.TAIR10$gene [grep("LTR_Gypsy",lincRNAs_TE_coverage.TAIR10$TE_type)] ) )
length(linc_LTR_gypsy) #249
a<-lincRNAs_TE_coverage.TAIR10 [grep("LTR_Gypsy",lincRNAs_TE_coverage.TAIR10$TE_type),]
linc_LTR_gypsy_only<-as.vector(unique(a$gene [grep(",",a$TE_type,invert = T)] ) )
length(linc_LTR_gypsy_only) #203

linc_LTR<-as.vector(unique(lincRNAs_TE_coverage.TAIR10$gene [grep("LTR",lincRNAs_TE_coverage.TAIR10$TE_type)] ) )
length(linc_LTR) #350

linc_SINELINE<-as.vector(unique(lincRNAs_TE_coverage.TAIR10$gene [grep("SINE",lincRNAs_TE_coverage.TAIR10$TE_type)] ) ) 
length(linc_SINELINE) #94


#boxplot CHH in lincs with different TE pieces vs different TAIR10 TE fragments 
##############################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_CHH_6909_1001G_linc_tetypes_and_TEfrags.4types.pdf",height = 4,width =4.5)
##############################################################################
par(mar=c(8,5,4,2))
l1<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% linc_noTE]
l2<- CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% linc_Helitron_only]
l3<-        CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% linc_DNA_mudr_only]
l4<-        CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% linc_LTR_copia_only]
l5<-        CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% linc_LTR_gypsy_only]
t1<-CHH.1001.araport$X6909[CHH.1001.araport$transcript %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="RC_Helitron"]]
t2<-        CHH.1001.araport$X6909[CHH.1001.araport$transcript %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="DNA_MuDR"]]
t3<-        CHH.1001.araport$X6909[CHH.1001.araport$transcript %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="LTR_Copia"]] 
t4<-        CHH.1001.araport$X6909[CHH.1001.araport$transcript %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="LTR_Gypsy"]]
boxplot(l1,l2,l3,l4,l5,t1,t2,t3,t4,
        outline = F,notch = T,
        main="CHH methylation (Col-0 rosette)",cex.main=1.2,ylim=c(0,0.17),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("noTE","Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy","Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy"),las=2)
mtext (text = "lincRNAs",side = 1,at = 3,line = 6)
mtext (text = "TAIR10 TE fragments",side = 1,at =7.5,line = 6)

##############################################################################
#################
#add p values   #
#################
a<-wilcox.test(sample(l1,200),sample(l2,200))
b<-wilcox.test(sample(l1,200),sample(l2,200))
c<-wilcox.test(sample(l1,200),sample(l2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.1)
a<-wilcox.test(sample(l3,200),sample(l2,200))
b<-wilcox.test(sample(l3,200),sample(l2,200))
c<-wilcox.test(sample(l3,200),sample(l2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.1)
a<-wilcox.test(sample(l3,200),l4)
b<-wilcox.test(sample(l3,200),l4)
c<-wilcox.test(sample(l3,200),l4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.1)
a<-wilcox.test(l4,l5)
b<-wilcox.test(l4,l5)
c<-wilcox.test(l4,l5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.1)

a<-wilcox.test(sample(t1,200),sample(t2,200))
b<-wilcox.test(sample(t1,200),sample(t2,200))
c<-wilcox.test(sample(t1,200),sample(t2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=0.1)
a<-wilcox.test(sample(t3,200),sample(t2,200))
b<-wilcox.test(sample(t3,200),sample(t2,200))
c<-wilcox.test(sample(t3,200),sample(t2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.1)
a<-wilcox.test(sample(t3,200),sample(t4,200))
b<-wilcox.test(sample(t3,200),sample(t4,200))
c<-wilcox.test(sample(t3,200),sample(t4,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=0.1)
#################
dev.off()

#boxplot CG in lincs with different TE pieces vs different TAIR10 TE fragments 
##############################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_CG_6909_1001G_linc_tetypes_and_TEfrags.4types.pdf",height = 4,width =4.5)
##############################################################################
par(mar=c(8,5,4,2))
l1<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% linc_noTE]
l2<- CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% linc_Helitron_only]
l3<-        CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% linc_DNA_mudr_only]
l4<-        CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% linc_LTR_copia_only]
l5<-        CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% linc_LTR_gypsy_only]
t1<-CG.1001.araport$X6909[CG.1001.araport$transcript %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="RC_Helitron"]]
t2<-        CG.1001.araport$X6909[CG.1001.araport$transcript %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="DNA_MuDR"]]
t3<-        CG.1001.araport$X6909[CG.1001.araport$transcript %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="LTR_Copia"]] 
t4<-        CG.1001.araport$X6909[CG.1001.araport$transcript %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="LTR_Gypsy"]]
boxplot(l1,l2,l3,l4,l5,t1,t2,t3,t4,
        outline = F,notch = T,
        main="CG methylation (Col-0 rosette)",cex.main=1.2,ylim=c(0,1),               col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("noTE","Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy","Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy"),las=2)
mtext (text = "lincRNAs",side = 1,at = 3,line = 6)
mtext (text = "TAIR10 TE fragments",side = 1,at =7.5,line = 6)

##############################################################################
#################
#add p values   #
#################
a<-wilcox.test(sample(l1,200),sample(l2,200))
b<-wilcox.test(sample(l1,200),sample(l2,200))
c<-wilcox.test(sample(l1,200),sample(l2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
a<-wilcox.test(sample(l3,200),sample(l2,200))
b<-wilcox.test(sample(l3,200),sample(l2,200))
c<-wilcox.test(sample(l3,200),sample(l2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
a<-wilcox.test(sample(l3,200),l4)
b<-wilcox.test(sample(l3,200),l4)
c<-wilcox.test(sample(l3,200),l4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(l4,l5)
b<-wilcox.test(l4,l5)
c<-wilcox.test(l4,l5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1)

a<-wilcox.test(sample(t1,200),sample(t2,200))
b<-wilcox.test(sample(t1,200),sample(t2,200))
c<-wilcox.test(sample(t1,200),sample(t2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=1)
a<-wilcox.test(sample(t3,200),sample(t2,200))
b<-wilcox.test(sample(t3,200),sample(t2,200))
c<-wilcox.test(sample(t3,200),sample(t2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
a<-wilcox.test(sample(t3,200),sample(t4,200))
b<-wilcox.test(sample(t3,200),sample(t4,200))
c<-wilcox.test(sample(t3,200),sample(t4,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=1)
#################
dev.off()



#boxplot K9 in lincs with different TE pieces vs different TAIR10 TE fragments 
##############################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_K9_6909_linc_tetypes_and_TEfrags.4types.pdf",height = 4,width =4)
##############################################################################
par(mar=c(8,5,4,2))
l1<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% linc_noTE]
l2<- chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% linc_Helitron_only]
l3<-        chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% linc_DNA_mudr_only]
l4<-        chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% linc_LTR_copia_only]
l5<-        chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% linc_LTR_gypsy_only]
t1<-chip.araport.log2$K9.6909[chip.araport.log2$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="RC_Helitron"]]
t2<-        chip.araport.log2$K9.6909[chip.araport.log2$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="DNA_MuDR"]]
t3<-        chip.araport.log2$K9.6909[chip.araport.log2$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="LTR_Copia"]] 
t4<-        chip.araport.log2$K9.6909[chip.araport.log2$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="LTR_Gypsy"]]
boxplot(l1,l2,l3,l4,l5,t1,t2,t3,t4,
        outline = F,notch = T,
        main="H3K9me2 level (Col-0 rosette)",cex.main=1.2, ylab="log2(ChIP/INPUT)", col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("noTE","Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy","Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy"),las=2)
mtext (text = "lincRNAs",side = 1,at = 3,line = 6)
mtext (text = "TAIR10 TE fragments",side = 1,at =7.5,line = 6)

##############################################################################
#################
#add p values   #
#################
a<-wilcox.test(sample(l1,200),sample(l2,200))
b<-wilcox.test(sample(l1,200),sample(l2,200))
c<-wilcox.test(sample(l1,200),sample(l2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
a<-wilcox.test(sample(l3,200),sample(l2,200))
b<-wilcox.test(sample(l3,200),sample(l2,200))
c<-wilcox.test(sample(l3,200),sample(l2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
a<-wilcox.test(sample(l3,200),l4)
b<-wilcox.test(sample(l3,200),l4)
c<-wilcox.test(sample(l3,200),l4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(l4,l5)
b<-wilcox.test(l4,l5)
c<-wilcox.test(l4,l5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1)

a<-wilcox.test(sample(t1,200),sample(t2,200))
b<-wilcox.test(sample(t1,200),sample(t2,200))
c<-wilcox.test(sample(t1,200),sample(t2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=1)
a<-wilcox.test(sample(t3,200),sample(t2,200))
b<-wilcox.test(sample(t3,200),sample(t2,200))
c<-wilcox.test(sample(t3,200),sample(t2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
a<-wilcox.test(sample(t3,200),sample(t4,200))
b<-wilcox.test(sample(t3,200),sample(t4,200))
c<-wilcox.test(sample(t3,200),sample(t4,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=1)
#################
dev.off()




#boxplot H1 in lincs with different TE pieces vs different TAIR10 TE fragments 
##############################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_H1_6909_linc_tetypes_and_TEfrags.4types.pdf",height = 4,width =4)
##############################################################################
par(mar=c(8,5,4,2))
l1<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% linc_noTE]
l2<- chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% linc_Helitron_only]
l3<-        chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% linc_DNA_mudr_only]
l4<-        chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% linc_LTR_copia_only]
l5<-        chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% linc_LTR_gypsy_only]
t1<-chip.araport.log2$K9.6909[chip.araport.log2$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="RC_Helitron"]]
t2<-        chip.araport.log2$K9.6909[chip.araport.log2$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="DNA_MuDR"]]
t3<-        chip.araport.log2$K9.6909[chip.araport.log2$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="LTR_Copia"]] 
t4<-        chip.araport.log2$K9.6909[chip.araport.log2$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="LTR_Gypsy"]]
boxplot(l1,l2,l3,l4,l5,t1,t2,t3,t4,
        outline = F,notch = T,
        main="H1 level (Col-0 rosette)",cex.main=1.2, ylab="log2(ChIP/INPUT)", col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("noTE","Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy","Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy"),las=2)
mtext (text = "lincRNAs",side = 1,at = 3,line = 6)
mtext (text = "TAIR10 TE fragments",side = 1,at =7.5,line = 6)

##############################################################################
#################
#add p values   #
#################
a<-wilcox.test(sample(l1,200),sample(l2,200))
b<-wilcox.test(sample(l1,200),sample(l2,200))
c<-wilcox.test(sample(l1,200),sample(l2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
a<-wilcox.test(sample(l3,200),sample(l2,200))
b<-wilcox.test(sample(l3,200),sample(l2,200))
c<-wilcox.test(sample(l3,200),sample(l2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
a<-wilcox.test(sample(l3,200),l4)
b<-wilcox.test(sample(l3,200),l4)
c<-wilcox.test(sample(l3,200),l4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(l4,l5)
b<-wilcox.test(l4,l5)
c<-wilcox.test(l4,l5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1)

a<-wilcox.test(sample(t1,200),sample(t2,200))
b<-wilcox.test(sample(t1,200),sample(t2,200))
c<-wilcox.test(sample(t1,200),sample(t2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=1)
a<-wilcox.test(sample(t3,200),sample(t2,200))
b<-wilcox.test(sample(t3,200),sample(t2,200))
c<-wilcox.test(sample(t3,200),sample(t2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
a<-wilcox.test(sample(t3,200),sample(t4,200))
b<-wilcox.test(sample(t3,200),sample(t4,200))
c<-wilcox.test(sample(t3,200),sample(t4,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=1)
#################
dev.off()

#boxplot _sRNA_24nt in lincs with different TE pieces vs different TAIR10 TE fragments 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_sRNA_24nt_6909_linc_tetypes_and_TEfrags.4types.pdf",height = 4,width =4)
##############################################################################
par(mar=c(8,5,4,2))
l1<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% linc_noTE]
l2<- sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% linc_Helitron_only]
l3<-        sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% linc_DNA_mudr_only]
l4<-        sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% linc_LTR_copia_only]
l5<-        sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% linc_LTR_gypsy_only]
t1<-sRNA.24nt.Ar11.RPM$X6909[sRNA.24nt.Ar11.RPM$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="RC_Helitron"]]
t2<-        sRNA.24nt.Ar11.RPM$X6909[sRNA.24nt.Ar11.RPM$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="DNA_MuDR"]]
t3<-        sRNA.24nt.Ar11.RPM$X6909[sRNA.24nt.Ar11.RPM$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="LTR_Copia"]] 
t4<-        sRNA.24nt.Ar11.RPM$X6909[sRNA.24nt.Ar11.RPM$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="LTR_Gypsy"]]
boxplot(l1,l2,l3,l4,l5,t1,t2,t3,t4,
        outline = F,notch = T,
        main="24nt siRNA coverage (Col-0 flowers)",cex.main=1.2, ylab="log2(ChIP/INPUT)", col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("noTE","Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy","Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy"),las=2)
mtext (text = "lincRNAs",side = 1,at = 3,line = 6)
mtext (text = "TAIR10 TE fragments",side = 1,at =7.5,line = 6)

##############################################################################
#################
#add p values   #
#################
a<-wilcox.test(sample(l1,200),sample(l2,200))
b<-wilcox.test(sample(l1,200),sample(l2,200))
c<-wilcox.test(sample(l1,200),sample(l2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=2)
a<-wilcox.test(sample(l3,200),sample(l2,200))
b<-wilcox.test(sample(l3,200),sample(l2,200))
c<-wilcox.test(sample(l3,200),sample(l2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=2)
a<-wilcox.test(sample(l3,200),l4)
b<-wilcox.test(sample(l3,200),l4)
c<-wilcox.test(sample(l3,200),l4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=2)
a<-wilcox.test(l4,l5)
b<-wilcox.test(l4,l5)
c<-wilcox.test(l4,l5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=2)

a<-wilcox.test(sample(t1,200),sample(t2,200))
b<-wilcox.test(sample(t1,200),sample(t2,200))
c<-wilcox.test(sample(t1,200),sample(t2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=1)
a<-wilcox.test(sample(t3,200),sample(t2,200))
b<-wilcox.test(sample(t3,200),sample(t2,200))
c<-wilcox.test(sample(t3,200),sample(t2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
a<-wilcox.test(sample(t3,200),sample(t4,200))
b<-wilcox.test(sample(t3,200),sample(t4,200))
c<-wilcox.test(sample(t3,200),sample(t4,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=1)
#################
dev.off()




#boxplot _sRNA_24nt in lincs with different TE pieces vs different TAIR10 TE fragments - 6909 early heart Ranj data
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig6/ ",height = 4,width =4)
##############################################################################
par(mar=c(8,5,4,2))
l1<-sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart[sRNA.24nt.denovo2021.RPM.Ranj$gene %in% linc_noTE]
l2<- sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart[sRNA.24nt.denovo2021.RPM.Ranj$gene %in% linc_Helitron_only]
l3<-sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart[sRNA.24nt.denovo2021.RPM.Ranj$gene %in% linc_DNA_mudr_only]
l4<-sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart[sRNA.24nt.denovo2021.RPM.Ranj$gene %in% linc_LTR_copia_only]
l5<-sRNA.24nt.denovo2021.RPM.Ranj$WT.eheart[sRNA.24nt.denovo2021.RPM.Ranj$gene %in% linc_LTR_gypsy_only]
#t1<-sRNA.24nt.Ar11.RPM$X6909[sRNA.24nt.Ar11.RPM$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="RC_Helitron"]]
#t2<-        sRNA.24nt.Ar11.RPM$X6909[sRNA.24nt.Ar11.RPM$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="DNA_MuDR"]]
#t3<-        sRNA.24nt.Ar11.RPM$X6909[sRNA.24nt.Ar11.RPM$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="LTR_Copia"]] 
#t4<-        sRNA.24nt.Ar11.RPM$X6909[sRNA.24nt.Ar11.RPM$gene %in% Ar11_TE_frag_types$TE[Ar11_TE_frag_types$TE_family=="LTR_Gypsy"]]
boxplot(l1,l2,l3,l4,l5,
        outline = F,notch = T,
        main="24nt siRNA coverage (Col-0 early heart)",cex.main=1.2, ylab="RPM", col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("noTE","Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy"),las=2)

boxplot(l1,l2,l3,l4,l5,t1,t2,t3,t4,
        outline = F,notch = T,
        main="H3K9me2 level (Col-0 rosette)",cex.main=1.2, ylab="log2(ChIP/INPUT)", col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("noTE","Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy","Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy"),las=2)
mtext (text = "lincRNAs",side = 1,at = 3,line = 6)
mtext (text = "TAIR10 TE fragments",side = 1,at =7.5,line = 6)

##############################################################################
#################
#add p values   #
#################
a<-wilcox.test(sample(l1,200),sample(l2,200))
b<-wilcox.test(sample(l1,200),sample(l2,200))
c<-wilcox.test(sample(l1,200),sample(l2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=2)
a<-wilcox.test(sample(l3,200),sample(l2,200))
b<-wilcox.test(sample(l3,200),sample(l2,200))
c<-wilcox.test(sample(l3,200),sample(l2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=2)
a<-wilcox.test(sample(l3,200),l4)
b<-wilcox.test(sample(l3,200),l4)
c<-wilcox.test(sample(l3,200),l4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=2)
a<-wilcox.test(l4,l5)
b<-wilcox.test(l4,l5)
c<-wilcox.test(l4,l5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=2)

a<-wilcox.test(sample(t1,200),sample(t2,200))
b<-wilcox.test(sample(t1,200),sample(t2,200))
c<-wilcox.test(sample(t1,200),sample(t2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=1)
a<-wilcox.test(sample(t3,200),sample(t2,200))
b<-wilcox.test(sample(t3,200),sample(t2,200))
c<-wilcox.test(sample(t3,200),sample(t2,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
a<-wilcox.test(sample(t3,200),sample(t4,200))
b<-wilcox.test(sample(t3,200),sample(t4,200))
c<-wilcox.test(sample(t3,200),sample(t4,200))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=8.5,y=1)
#################
dev.off()
















### distribution of TE pieces of different types. 


linc_LTR_gypsy_only
linc_Helitron_only

plot(linc_TE_cov_all_loci_2cols)

a<-merge(linc_TE_cov_all_loci_2cols, lncRNAs.intergenic.loci, by="gene")
a<-merge(a, chip.denovo.log2[,c("gene","H1.6909","K9.6909","K4.6909","K27.6909","K36.6909")], by="gene")
a<-merge(a, chip.denovo.quantstan[,c("gene","variance.hist1", "variance.key4", "variance.key9", "variance.key27", "variance.key36")], by="gene")
a<-merge(a, CG.1001.denovo[c("transcript","X6909","sd")], by.x="gene",by.y="transcript")
a<-merge(a, CG.1001new.denovo[c("transcript","mean.6909","sd_of_means")], by.x="gene",by.y="transcript")
a<-merge(a, sRNA.24nt.denovo2021.RPM[c("gene","X6909")], by.x="gene",by.y="gene")
a<-merge(a, CHH.1001.denovo[c("transcript","X6909","sd")], by.x="gene",by.y="transcript")
a<-merge(a, CHH.1001new.denovo[c("transcript","mean.6909","sd_of_means")], by.x="gene",by.y="transcript")

plot(a$dist_from_centromere,a$coverage)
plot(a$dist_from_centromere,a$K9.6909)
plot(a$dist_from_centromere,a$X6909)



mudr<-a[a$gene %in% linc_DNA_mudr_only,]
ltr_copia<-a[a$gene %in% linc_LTR_copia_only,]

ltr<-a[a$gene %in% linc_LTR_gypsy_only,]
heli<-a[a$gene %in% linc_Helitron_only,]
plot(ltr$dist_from_centromere,ltr$K9.6909,col="blue",pch=18)
points(heli$dist_from_centromere,heli$K9.6909,col="red",pch=1)

#boxplot H3K9me2 in lincs with Helitron or LTR_Gypsy TE pieces vs distance to centromeres
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_H3K9me2_lincs_withHeli_LTR_vs_distance_to_centromeres.pdf",height = 4,width =6)
##############################################################################
par(mar=c(8,5,4,2))
l1<-  ltr$K9.6909[ltr$dist_from_centromere<=1000000]
l2<- heli$K9.6909[heli$dist_from_centromere<=1000000]
l3<- ltr$K9.6909[ltr$dist_from_centromere>1000000&ltr$dist_from_centromere<=2000000]
l4<- heli$K9.6909[heli$dist_from_centromere>1000000&heli$dist_from_centromere<=2000000]
l5<- ltr$K9.6909[ltr$dist_from_centromere>2000000&ltr$dist_from_centromere<=3000000]
l6<- heli$K9.6909[heli$dist_from_centromere>2000000&heli$dist_from_centromere<=3000000]
l7<- ltr$K9.6909[ltr$dist_from_centromere>3000000&ltr$dist_from_centromere<=5000000]
l8<- heli$K9.6909[heli$dist_from_centromere>3000000&heli$dist_from_centromere<=5000000]
l9<- ltr$K9.6909[ltr$dist_from_centromere>5000000&ltr$dist_from_centromere<=8000000]
l10<- heli$K9.6909[heli$dist_from_centromere>5000000&heli$dist_from_centromere<=8000000]
boxplot(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,
        outline = F,notch = T,names=c("LTR_G<1Mb","Heli<1Mb","1Mb<LTR_G<2Mb","1Mb<Heli<2Mb","2Mb<LTR_G<3Mb","2Mb<Heli<3Mb","3Mb<LTR_G<5Mb","3Mb<Heli<5Mb","5Mb<LTR_G<8Mb","5Mb<Heli<8Mb"),las=2, col=c("bisque2","burlywood3"),main="H3K9me2 level\n (Col-0,rosette)",ylab="log2(ChIP/Input)")
mtext(paste("N=",length(l1),sep=""),at=1,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l2),sep=""),at=2,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l3),sep=""),at=3,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l4),sep=""),at=4,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l5),sep=""),at=5,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l6),sep=""),at=6,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l7),sep=""),at=7,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l8),sep=""),at=8,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l9),sep=""),at=9,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l10),sep=""),at=10,line=-1,side=1,cex=0.9)
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
text(b,x=1.5,y=1.5)
a<-wilcox.test(l3,l4)
b<-wilcox.test(l3,l4)
c<-wilcox.test(l3,l4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1.5)
a<-wilcox.test(l5,l6)
b<-wilcox.test(l5,l6)
c<-wilcox.test(l5,l6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1.5)
a<-wilcox.test(l7,l8)
b<-wilcox.test(l7,l8)
c<-wilcox.test(l7,l8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
a<-wilcox.test(l9,l10)
b<-wilcox.test(l9,l10)
c<-wilcox.test(l9,l10)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=9.5,y=1)
#################
dev.off()

#boxplot CG in lincs with Helitron or LTR_Gypsy TE pieces vs distance to centromeres
##############################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_CG_1001_lincs_withHeli_LTR_vs_distance_to_centromeres.pdf",height = 4,width =6)
par(mar=c(8,5,4,2))
l1<-  ltr$X6909.x[ltr$dist_from_centromere<=1000000]
l2<- heli$X6909.x[heli$dist_from_centromere<=1000000]
l3<- ltr$X6909.x[ltr$dist_from_centromere>1000000&ltr$dist_from_centromere<=2000000]
l4<- heli$X6909.x[heli$dist_from_centromere>1000000&heli$dist_from_centromere<=2000000]
l5<- ltr$X6909.x[ltr$dist_from_centromere>2000000&ltr$dist_from_centromere<=3000000]
l6<- heli$X6909.x[heli$dist_from_centromere>2000000&heli$dist_from_centromere<=3000000]
l7<- ltr$X6909.x[ltr$dist_from_centromere>3000000&ltr$dist_from_centromere<=5000000]
l8<- heli$X6909.x[heli$dist_from_centromere>3000000&heli$dist_from_centromere<=5000000]
l9<- ltr$X6909.x[ltr$dist_from_centromere>5000000&ltr$dist_from_centromere<=8000000]
l10<- heli$X6909.x[heli$dist_from_centromere>5000000&heli$dist_from_centromere<=8000000]
boxplot(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,
        outline = F,notch = T,names=c("LTR_G<1Mb","Heli<1Mb","1Mb<LTR_G<2Mb","1Mb<Heli<2Mb","2Mb<LTR_G<3Mb","2Mb<Heli<3Mb","3Mb<LTR_G<5Mb","3Mb<Heli<5Mb","5Mb<LTR_G<8Mb","5Mb<Heli<8Mb"),las=2, col=c("bisque2","burlywood3"),main="CG methylation level\n (Col-0,rosette, 1001G)",ylab="CG")
mtext(paste("N=",length(l1),sep=""),at=1,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l2),sep=""),at=2,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l3),sep=""),at=3,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l4),sep=""),at=4,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l5),sep=""),at=5,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l6),sep=""),at=6,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l7),sep=""),at=7,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l8),sep=""),at=8,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l9),sep=""),at=9,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l10),sep=""),at=10,line=-1,side=1,cex=0.9)
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
text(b,x=1.5,y=1)
a<-wilcox.test(l3,l4)
b<-wilcox.test(l3,l4)
c<-wilcox.test(l3,l4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)
a<-wilcox.test(l5,l6)
b<-wilcox.test(l5,l6)
c<-wilcox.test(l5,l6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1)
a<-wilcox.test(l7,l8)
b<-wilcox.test(l7,l8)
c<-wilcox.test(l7,l8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1)
a<-wilcox.test(l9,l10)
b<-wilcox.test(l9,l10)
c<-wilcox.test(l9,l10)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=9.5,y=1)
#################
dev.off()


#boxplot CHH in lincs with Helitron or LTR_Gypsy TE pieces vs distance to centromeres
##############################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_CHH_1001_lincs_withHeli_LTR_vs_distance_to_centromeres.pdf",height = 4,width =6)
par(mar=c(8,5,4,2))
l1<-  ltr$X6909[ltr$dist_from_centromere<=1000000]
l2<- heli$X6909[heli$dist_from_centromere<=1000000]
l3<- ltr$X6909[ltr$dist_from_centromere>1000000&ltr$dist_from_centromere<=2000000]
l4<- heli$X6909[heli$dist_from_centromere>1000000&heli$dist_from_centromere<=2000000]
l5<- ltr$X6909[ltr$dist_from_centromere>2000000&ltr$dist_from_centromere<=3000000]
l6<- heli$X6909[heli$dist_from_centromere>2000000&heli$dist_from_centromere<=3000000]
l7<- ltr$X6909[ltr$dist_from_centromere>3000000&ltr$dist_from_centromere<=5000000]
l8<- heli$X6909[heli$dist_from_centromere>3000000&heli$dist_from_centromere<=5000000]
l9<- ltr$X6909[ltr$dist_from_centromere>5000000&ltr$dist_from_centromere<=8000000]
l10<- heli$X6909[heli$dist_from_centromere>5000000&heli$dist_from_centromere<=8000000]
boxplot(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,
        outline = F,notch = T,names=c("LTR_G<1Mb","Heli<1Mb","1Mb<LTR_G<2Mb","1Mb<Heli<2Mb","2Mb<LTR_G<3Mb","2Mb<Heli<3Mb","3Mb<LTR_G<5Mb","3Mb<Heli<5Mb","5Mb<LTR_G<8Mb","5Mb<Heli<8Mb"),las=2, col=c("bisque2","burlywood3"),main="CHH methylation level\n (Col-0,rosette, 1001G)",ylab="CHH")
mtext(paste("N=",length(l1),sep=""),at=1,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l2),sep=""),at=2,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l3),sep=""),at=3,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l4),sep=""),at=4,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l5),sep=""),at=5,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l6),sep=""),at=6,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l7),sep=""),at=7,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l8),sep=""),at=8,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l9),sep=""),at=9,line=-1,side=1,cex=0.9)
mtext(paste("N=",length(l10),sep=""),at=10,line=-1,side=1,cex=0.9)
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
text(b,x=1.5,y=0.1)
a<-wilcox.test(l3,l4)
b<-wilcox.test(l3,l4)
c<-wilcox.test(l3,l4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.1)
a<-wilcox.test(l5,l6)
b<-wilcox.test(l5,l6)
c<-wilcox.test(l5,l6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.1)
a<-wilcox.test(l7,l8)
b<-wilcox.test(l7,l8)
c<-wilcox.test(l7,l8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.1)
a<-wilcox.test(l9,l10)
b<-wilcox.test(l9,l10)
c<-wilcox.test(l9,l10)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=9.5,y=0.1)
#################
dev.off()



#boxplot 24nt  coverage lincs with MuDR or LTR_Gypsy TE pieces vs distance to centromeres
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_boxplot_24nt_lincs_withMuDR_LTR_vs_distance_to_centromeres.pdf",height = 5,width =6)
##############################################################################
par(mar=c(8,5,4,2))
l1<-  ltr$X6909.y[ltr$dist_from_centromere<=1000000]
l2<- mudr$X6909.y[mudr$dist_from_centromere<=1000000]
l3<- ltr$X6909.y[ltr$dist_from_centromere>1000000&ltr$dist_from_centromere<=2000000]
l4<- mudr$X6909.y[mudr$dist_from_centromere>1000000&mudr$dist_from_centromere<=2000000]
l5<- ltr$X6909.y[ltr$dist_from_centromere>2000000&ltr$dist_from_centromere<=3000000]
l6<- mudr$X6909.y[mudr$dist_from_centromere>2000000&mudr$dist_from_centromere<=3000000]
l7<- ltr$X6909.y[ltr$dist_from_centromere>3000000&ltr$dist_from_centromere<=5000000]
l8<- mudr$X6909.y[mudr$dist_from_centromere>3000000&mudr$dist_from_centromere<=5000000]
l9<- ltr$X6909.y[ltr$dist_from_centromere>5000000&ltr$dist_from_centromere<=8000000]
l10<- mudr$X6909.y[mudr$dist_from_centromere>5000000&mudr$dist_from_centromere<=8000000]
boxplot(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,
        outline = F,notch = T,names=c("LTR_G<1Mb","MuDR<1Mb","1Mb<LTR_G<2Mb","1Mb<MuDR<2Mb","2Mb<LTR_G<3Mb","2Mb<MuDR<3Mb","3Mb<LTR_G<5Mb","3Mb<MuDR<5Mb","5Mb<LTR_G<8Mb","5Mb<MuDR<8Mb"),las=2, col=c("bisque2","burlywood3"),main="24nt siRNA coverage\n (Col-0,flowers)",ylab="24nt, RPM")
mtext(paste("N=",length(l1),sep=""),at=1,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(l2),sep=""),at=2,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(l3),sep=""),at=3,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(l4),sep=""),at=4,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(l5),sep=""),at=5,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(l6),sep=""),at=6,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(l7),sep=""),at=7,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(l8),sep=""),at=8,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(l9),sep=""),at=9,line=-1,side=1,cex=0.6)
mtext(paste("N=",length(l10),sep=""),at=10,line=-1,side=1,cex=0.6)
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
text(b,x=1.5,y=1.4)
a<-wilcox.test(l3,l4)
b<-wilcox.test(l3,l4)
c<-wilcox.test(l3,l4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1.4)
a<-wilcox.test(l5,l6)
b<-wilcox.test(l5,l6)
c<-wilcox.test(l5,l6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1.4)
a<-wilcox.test(l7,l8)
b<-wilcox.test(l7,l8)
c<-wilcox.test(l7,l8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1.4)
a<-wilcox.test(l9,l10)
b<-wilcox.test(l9,l10)
c<-wilcox.test(l9,l10)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=9.5,y=1.4)
#################
dev.off()



length(ltr$dist_from_centromere)#203
length(heli$dist_from_centromere)#586
length(mudr$dist_from_centromere)#323

# plot the distance from centromeres for lincRNAs with pieces of different types of TEs 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig4/Fig4_SUPPL_lincs_3types_plot_distance_to_centromeres_distib.pdf",height = 4,width =5)
###########################################################
par(mar=c(8,5,4,2))
plot(ltr$dist_from_centromere[order(ltr$dist_from_centromere,decreasing = T)]/1000000,ylim=c(0,16),xlim=c(0,600),pch=".",col="#4ABAB8",ylab="distance from centromere,Mb",xlab="lincRNAs with TE piece",main="")
points(218,15,pch=".",col="#4ABAB8")
text(x =400,y=15,label=paste("lincRNAs with a piece of LTR_Gypsy, N =" ,length(ltr$dist_from_centromere)),cex=0.7)
points(heli$dist_from_centromere[order(heli$dist_from_centromere,decreasing = T)]/1000000,pch=".",col="#C5A49C")
points(230,12.5,pch=".",col="#C5A49C")
text(x =400,y=12.5,label=paste("lincRNAs with a piece of Helitron, N =" ,length(heli$dist_from_centromere)),cex=0.7)
points(mudr$dist_from_centromere[order(mudr$dist_from_centromere,decreasing = T)]/1000000,pch=".",col="#BA77B0")
points(230,10,pch=".",col="#BA77B0")
text(x =400,y=10,label=paste("lincRNAs with a piece of MuDR, N =" ,length(ltr$dist_from_centromere)),cex=0.7)
#points(ltr_copia$dist_from_centromere[order(ltr_copia$dist_from_centromere,decreasing = T)]/1000000,pch=19,col="#4581C2")
abline(h = 1,lty=2,col="orange")
###########################################################
dev.off()













































#########################################################################
### repression signal in pieces vs. the rest of locus 
#########################################################################
#import methylation 
CG.lincs_minus_tepieces.6909 <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/allc_6909.CG.lincs_minus_tepieces.bed", header=FALSE)
names(CG.lincs_minus_tepieces.6909)<-c("Chr","start","end","gene","score","strand","CG_rest")
CG.lincs_minus_tepieces.6909$length<-CG.lincs_minus_tepieces.6909$end-CG.lincs_minus_tepieces.6909$start
CG.tepieces_in_lincs.6909 <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/allc_6909.CG.tepieces_in_lincs.bed", header=FALSE)
names(CG.tepieces_in_lincs.6909)<-c("Chr","start","end","TE_types","gene","strand","CG_tepiece")
CG.tepieces_in_lincs.6909$length_Tepiece<-CG.tepieces_in_lincs.6909$end-CG.tepieces_in_lincs.6909$start

a<-merge(CG.tepieces_in_lincs.6909,CG.lincs_minus_tepieces.6909,by="gene",all.x=T)
a$te_direction[a$strand.x==a$strand.y]<-"forward"
a$te_direction[a$strand.x!=a$strand.y]<-"reverse"
aa_cg<-a[,c("gene","length","TE_types","length_Tepiece","te_direction","CG_rest","CG_tepiece" )]
aa_cg<-aa_cg[aa_cg$length>50,]


CHH.lincs_minus_tepieces.6909 <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/allc_6909.CHH.lincs_minus_tepieces.bed", header=FALSE)
names(CHH.lincs_minus_tepieces.6909)<-c("Chr","start","end","gene","score","strand","CHH_rest")
CHH.lincs_minus_tepieces.6909$length<-CHH.lincs_minus_tepieces.6909$end-CHH.lincs_minus_tepieces.6909$start
CHH.tepieces_in_lincs.6909 <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/allc_6909.CHH.tepieces_in_lincs.bed", header=FALSE)
names(CHH.tepieces_in_lincs.6909)<-c("Chr","start","end","TE_types","gene","strand","CHH_tepiece")
CHH.tepieces_in_lincs.6909$length_Tepiece<-CHH.tepieces_in_lincs.6909$end-CHH.tepieces_in_lincs.6909$start

a<-merge(CHH.tepieces_in_lincs.6909,CHH.lincs_minus_tepieces.6909,by="gene",all.x=T)
a$te_direction[a$strand.x==a$strand.y]<-"forward"
a$te_direction[a$strand.x!=a$strand.y]<-"reverse"
aa_chh<-a[,c("gene","length","TE_types","length_Tepiece","te_direction","CHH_rest","CHH_tepiece" )]
aa_chh<-aa_chh[aa_chh$length>50,]




hist(aa$CG_tepiece[grep("DNA_MuDR",aa$TE_types)]-aa$CG_rest[grep("DNA_MuDR",aa$TE_types)])
hist(aa$CG_tepiece[grep("LTR",aa$TE_types)]-aa$CG_rest[grep("LTR",aa$TE_types)])
hist(aa$CG_tepiece[grep("^(?!.*Heli)^(?!.*DNA).*LTR",aa$TE_types,perl = T)]-aa$CG_rest[grep("^(?!.*Heli)^(?!.*DNA).*LTR",aa$TE_types,perl = T)])


# import miRNAs
sRNA.24nt.lincs_minus_tepieces.6909 <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/6909.24nt.lincs_minus_tepieces.coverage.perbase_calc.bed", header=FALSE)
sRNA.24nt.tepieces_in_lincs.6909 <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/6909.24nt.tepieces_in_lincs.coverage.perbase_calc.bed", header=FALSE)
names(sRNA.24nt.lincs_minus_tepieces.6909)<-c("Chr","start","end","gene","score","strand","sRNA_rest")
sRNA.24nt.lincs_minus_tepieces.6909$length<-sRNA.24nt.lincs_minus_tepieces.6909$end-sRNA.24nt.lincs_minus_tepieces.6909$start
names(sRNA.24nt.tepieces_in_lincs.6909)<-c("Chr","start","end","TE_types","gene","strand","sRNA_tepiece")
sRNA.24nt.tepieces_in_lincs.6909$length_Tepiece<-sRNA.24nt.tepieces_in_lincs.6909$end-sRNA.24nt.tepieces_in_lincs.6909$start


b<-merge(sRNA.24nt.tepieces_in_lincs.6909,sRNA.24nt.lincs_minus_tepieces.6909,by="gene",all.x=T)
b$te_direction[b$strand.x==b$strand.y]<-"forward"
b$te_direction[b$strand.x!=b$strand.y]<-"reverse"
bb<-b[,c("gene","length","TE_types","length_Tepiece","te_direction","sRNA_rest","sRNA_tepiece" )]
bb<-bb[bb$length>50,]

hist(bb$sRNA_tepiece[grep("DNA_MuDR",bb$TE_types)]-bb$sRNA_rest[grep("DNA_MuDR",bb$TE_types)],xlim=c(-3,3),breaks=300,col="grey")
plot(density(bb$sRNA_tepiece[grep("LTR",bb$TE_types)]-bb$sRNA_rest[grep("LTR",bb$TE_types)],n = length(bb$sRNA_tepiece[grep("LTR",bb$TE_types)])),xlim=c(-3,3))
lines(density(bb$sRNA_tepiece[grep("DNA_MuDR",bb$TE_types)]-bb$sRNA_rest[grep("DNA_MuDR",bb$TE_types)],n = length(bb$sRNA_tepiece[grep("DNA_MuDR",bb$TE_types)])),xlim=c(-3,3),col="purple")

boxplot((bb$sRNA_tepiece[grep("DNA_MuDR",bb$TE_types)]+0.01)-(bb$sRNA_rest[grep("DNA_MuDR",bb$TE_types)]+0.01),
        (bb$sRNA_tepiece[grep("LTR",bb$TE_types)]+0.01)-(bb$sRNA_rest[grep("LTR",bb$TE_types)]+0.01),ylim=c(-2,2),outline = F,notch = T)

hist(bb$sRNA_tepiece[grep("LTR",bb$TE_types)]-bb$sRNA_rest[grep("LTR",bb$TE_types)])
boxplot(bb$sRNA_tepiece[grep("^(?!.*Heli)^(?!.*DNA).*LTR",bb$TE_types,perl = T)]-bb$sRNA_rest[grep("^(?!.*Heli)^(?!.*DNA).*LTR",bb$TE_types,perl = T)],
        bb$sRNA_tepiece[grep("^(?!.*LTR)^(?!.*SINE).*DNA_MuDR",bb$TE_types,perl = T)]-bb$sRNA_rest[grep("^(?!.*LTR)^(?!.*SINE).*DNA_MuDR",bb$TE_types,perl = T)],ylim=c(-1.2,1.8),outline = F,notch = T)
boxplot(bb$sRNA_tepiece[grep("^(?!.*Heli)^(?!.*DNA).*LTR",bb$TE_types,perl = T)]-bb$sRNA_rest[grep("^(?!.*Heli)^(?!.*DNA).*LTR",bb$TE_types,perl = T)],
        bb$sRNA_tepiece[grep("^(?!.*LTR)^(?!.*SINE).*DNA_MuDR",bb$TE_types,perl = T)]-bb$sRNA_rest[grep("^(?!.*LTR)^(?!.*SINE).*DNA_MuDR",bb$TE_types,perl = T)],ylim=c(-1.2,1.8),outline = F,notch = T)


#import chipseq

K9_rep2.lincs_minus_tepieces<- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/r14.rep2.H3K9me2.lincs_minus_tepieces.log2.mean_cov.bed", header=FALSE)
names(K9_rep2.lincs_minus_tepieces)<-c("Chr","start","end","gene","score","strand","k9_rest.rep2")
K9_rep2.tepieces_in_lincs<- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/r14.rep2.H3K9me2.tepieces_in_lincs.log2.mean_cov.bed", header=FALSE)
names(K9_rep2.tepieces_in_lincs)<-c("Chr","start","end","TE_types","gene","strand","k9_TE.rep2")
K9_rep2.lincs_minus_tepieces$length<-K9_rep2.lincs_minus_tepieces$end-K9_rep2.lincs_minus_tepieces$start
K9_rep2.tepieces_in_lincs$length_Tepiece<-K9_rep2.tepieces_in_lincs$end-K9_rep2.tepieces_in_lincs$start

K9_rep1.lincs_minus_tepieces<- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/r14.rep1.H3K9me2.lincs_minus_tepieces.log2.mean_cov.bed", header=FALSE)
names(K9_rep1.lincs_minus_tepieces)<-c("Chr","start","end","gene","score","strand","k9_rest.rep1")
K9_rep1.tepieces_in_lincs<- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/05_lncRNAs_vs_TEs/r14.rep1.H3K9me2.tepieces_in_lincs.log2.mean_cov.bed", header=FALSE)
names(K9_rep1.tepieces_in_lincs)<-c("Chr","start","end","TE_types","gene","strand","k9_TE.rep1")
K9_rep1.lincs_minus_tepieces$length<-K9_rep1.lincs_minus_tepieces$end-K9_rep1.lincs_minus_tepieces$start
K9_rep1.tepieces_in_lincs$length_Tepiece<-K9_rep1.tepieces_in_lincs$end-K9_rep1.tepieces_in_lincs$start


b<-merge(K9_rep2.lincs_minus_tepieces,K9_rep2.tepieces_in_lincs,by="gene",all.x=T)
b$te_direction[b$strand.x==b$strand.y]<-"forward"
b$te_direction[b$strand.x!=b$strand.y]<-"reverse"
bb2<-b[,c("gene","length","TE_types","length_Tepiece","te_direction","k9_rest.rep2","k9_TE.rep2" )]
bb2<-bb2[bb2$length>50,]


b<-merge(K9_rep1.lincs_minus_tepieces,K9_rep1.tepieces_in_lincs,by="gene",all.x=T)
b$te_direction[b$strand.x==b$strand.y]<-"forward"
b$te_direction[b$strand.x!=b$strand.y]<-"reverse"
bb1<-b[,c("gene","length","TE_types","length_Tepiece","te_direction","k9_rest.rep1","k9_TE.rep1" )]
bb1<-bb1[bb1$length>50,]
boxplot(bb1$k9_TE.rep1,bb1$k9_rest.rep1)

#make a summary table with CG, CHH, sRNA, H3K9me2

summary_tepiece_rest<-as.data.frame(unique(aa$gene))
summary_tepiece_rest$gene<-summary_tepiece_rest$`unique(aa$gene)`
summary_tepiece_rest$CG_tepiece<-0
summary_tepiece_rest$CG_rest<-0
summary_tepiece_rest$CHH_tepiece<-0
summary_tepiece_rest$CHH_rest<-0
summary_tepiece_rest$sRNA_tepiece<-0
summary_tepiece_rest$sRNA_rest<-0
summary_tepiece_rest$k9_tepiece_rep1<-0
summary_tepiece_rest$k9_rest_rep1<-0
summary_tepiece_rest$k9_tepiece_rep2<-0
summary_tepiece_rest$k9_rest_rep2<-0



for (i in 1:934){
  linc=as.character(summary_tepiece_rest$gene[i])
  df<-aa_cg[aa_cg$gene==linc & !is.na(aa_cg$gene),]
  summary_tepiece_rest$CG_tepiece[i]<-mean(df$CG_tepiece)
  summary_tepiece_rest$CG_rest[i] <-mean(df$CG_rest)
  df<-aa_chh[aa_chh$gene==linc & !is.na(aa_chh$gene),]
  summary_tepiece_rest$CHH_tepiece[i]<-mean(df$CHH_tepiece)
  summary_tepiece_rest$CHH_rest[i] <-mean(df$CHH_rest)
  df<-bb[bb$gene==linc & !is.na(bb$gene),]
  summary_tepiece_rest$sRNA_tepiece[i]<-mean(df$sRNA_tepiece)
  summary_tepiece_rest$sRNA_rest[i] <-mean(df$sRNA_rest)
  df<-bb1[bb1$gene==linc & !is.na(bb1$gene),]
  summary_tepiece_rest$k9_tepiece_rep1[i]<-mean(df$k9_TE.rep1)
  summary_tepiece_rest$k9_rest_rep1[i] <-mean(df$k9_rest.rep1)
  df<-bb2[bb2$gene==linc & !is.na(bb2$gene),]
  summary_tepiece_rest$k9_tepiece_rep2[i]<-mean(df$k9_TE.rep2)
  summary_tepiece_rest$k9_rest_rep2[i] <-mean(df$k9_rest.rep2)
  
}

summary_tepiece_rest$k9_tepiece<-apply(summary_tepiece_rest[,c("k9_tepiece_rep1","k9_tepiece_rep2")],1,mean)
summary_tepiece_rest$k9_rest<-apply(summary_tepiece_rest[,c("k9_rest_rep1","k9_rest_rep2")],1,mean)



hist(summary_tepiece_rest$sRNA_tepiece-summary_tepiece_rest$sRNA_rest,xlim=c(-2,2),breaks=300,col="grey")
hist(summary_tepiece_rest$k9_tepiece_rep1-summary_tepiece_rest$k9_rest_rep1,xlim=c(-2,2),breaks=30,col="grey")
hist(summary_tepiece_rest$k9_tepiece_rep2-summary_tepiece_rest$k9_rest_rep2,xlim=c(-2,2),breaks=30,col="grey")

boxplot(summary_tepiece_rest$k9_tepiece_rep1,summary_tepiece_rest$k9_rest_rep1,outline = F,notch = T)
wilcox.test(summary_tepiece_rest$k9_tepiece_rep1,summary_tepiece_rest$k9_rest_rep1) # p-value = 0.0004473

boxplot(summary_tepiece_rest$k9_tepiece_rep2,summary_tepiece_rest$k9_rest_rep2,outline = F,notch = T)
wilcox.test(summary_tepiece_rest$k9_tepiece_rep2,summary_tepiece_rest$k9_rest_rep2) # p-value = 0.0004473

wilcox.test(summary_tepiece_rest$sRNA_tepiece,summary_tepiece_rest$sRNA_rest) #p-value = 2.604e-05



hist(summary_tepiece_rest$sRNA_tepiece[ summary_tepiece_rest$gene %in% bb$gene[grep("DNA_MuDR",bb$TE_types)]]-summary_tepiece_rest$sRNA_rest[ summary_tepiece_rest$gene %in% bb$gene[grep("DNA_MuDR",bb$TE_types)]] ,xlim=c(-3,3),breaks=300,col="grey")
hist(summary_tepiece_rest$sRNA_tepiece[ summary_tepiece_rest$gene %in% bb$gene[grep("LTR",bb$TE_types)]]-summary_tepiece_rest$sRNA_rest[ summary_tepiece_rest$gene %in% bb$gene[grep("LTR",bb$TE_types)]] ,xlim=c(-3,3),breaks=300,col="grey")
hist(summary_tepiece_rest$sRNA_tepiece[ summary_tepiece_rest$gene %in% bb$gene[grep("Heli",bb$TE_types)]]-summary_tepiece_rest$sRNA_rest[ summary_tepiece_rest$gene %in% bb$gene[grep("Heli",bb$TE_types)]] ,xlim=c(-3,3),breaks=300,col="grey")

boxplot(summary_tepiece_rest$sRNA_tepiece[ summary_tepiece_rest$gene %in% bb$gene[grep("Heli",bb$TE_types)]],summary_tepiece_rest$sRNA_rest[ summary_tepiece_rest$gene %in% bb$gene[grep("Heli",bb$TE_types)]],outline = F,notch = T )
wilcox.test(summary_tepiece_rest$sRNA_tepiece[ summary_tepiece_rest$gene %in% bb$gene[grep("Heli",bb$TE_types)]],summary_tepiece_rest$sRNA_rest[ summary_tepiece_rest$gene %in% bb$gene[grep("Heli",bb$TE_types)]]) #p-value = 0.0001061
boxplot(summary_tepiece_rest$sRNA_tepiece[ summary_tepiece_rest$gene %in% bb$gene[grep("LTR",bb$TE_types)]],summary_tepiece_rest$sRNA_rest[ summary_tepiece_rest$gene %in% bb$gene[grep("LTR",bb$TE_types)]],outline = F,notch = T ) #p-value = 3.43e-05
wilcox.test(summary_tepiece_rest$sRNA_tepiece[ summary_tepiece_rest$gene %in% bb$gene[grep("LTR",bb$TE_types)]],summary_tepiece_rest$sRNA_rest[ summary_tepiece_rest$gene %in% bb$gene[grep("LTR",bb$TE_types)]])
boxplot(summary_tepiece_rest$sRNA_tepiece[ summary_tepiece_rest$gene %in% bb$gene[grep("DNA_MuDR",bb$TE_types)]],summary_tepiece_rest$sRNA_rest[ summary_tepiece_rest$gene %in% bb$gene[grep("DNA_MuDR",bb$TE_types)]],outline = F,notch = T )
wilcox.test(summary_tepiece_rest$sRNA_tepiece[ summary_tepiece_rest$gene %in% bb$gene[grep("DNA_MuDR",bb$TE_types)]],summary_tepiece_rest$sRNA_rest[ summary_tepiece_rest$gene %in% bb$gene[grep("DNA_MuDR",bb$TE_types)]])#p-value = 1.629e-07

#srna - boxplot
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Supple_for_Fig7_24nt_coverage_TEpiece_vs_rest_of_locus_lincs.pdf",height = 4,width = 3)
####################################################
par(mar=c(8,4,3,2)) 
#boxplot(summary_tepiece_rest$sRNA_tepiece,summary_tepiece_rest$sRNA_rest,outline = F,notch = T, col=c("burlywood4","burlywood2"), names=c("TE patches\n inside lincRNA","TE-free part\n of lincRNA") ,las=2, main="24nt coverage\n(flowers, Col-0)",ylab="24nt coverage,RPM")
boxplot(summary_tepiece_rest$sRNA_tepiece,summary_tepiece_rest$sRNA_rest,sRNA.24nt.denovo2021.RPM$X6909[!(sRNA.24nt.denovo2021.RPM$gene %in% linc_TE_cov_loci$gene)],   outline = F,notch = T, col=c("burlywood4","burlywood2","lightgray"), names=c("TE patches\n inside lincRNA","TE-free part\n of lincRNA","lincRNAs wo TE") ,las=2, main="24nt coverage\n(flowers, Col-0)",ylab="24nt coverage,RPM")

a<-wilcox.test(summary_tepiece_rest$sRNA_tepiece,summary_tepiece_rest$sRNA_rest)
a$p.value#2.604364e-05
text(paste("p=",round(a$p.value,6)),x=1.5,y=1.5,cex=0.8)

a<-wilcox.test(summary_tepiece_rest$sRNA_rest,sRNA.24nt.denovo2021.RPM$X6909[!(sRNA.24nt.denovo2021.RPM$gene %in% linc_TE_cov_loci$gene)])
a$p.value#3.40328e-136
text("***",x=2.5,y=1.5,cex=0.8)
####################################################
dev.off()

##srna - histogram
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Suppl_Fig7_histogram_24nt_TEpiece_minus_rest_of_locus.pdf",height = 3,width = 4)
####################################################
vect<-summary_tepiece_rest$sRNA_tepiece-summary_tepiece_rest$sRNA_rest
# Layout to split the screen
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(4,2))
# Draw the boxplot and the histogram
par(mar=c(0, 4, 1.1, 2.1))
hist(vect , col="grey" , border=T , xlab="", xlim=c(-1,1), main="24nt sRNA level difference",ylab="number of lincRNA loci",breaks=300)
par(mar=c(0, 4, 1.1, 2.1))
boxplot(vect , horizontal=TRUE , ylim=c(-1,1), xaxt="n" , col=rgb(0.8,0.8,0,0.5) , frame=F,outline = F,notch=T)
mtext("24nt coverage difference: TE patches - rest of locus",line=-1,side=1)
##################################################
dev.off()      

#CG boxplot  
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Supple_for_Fig7_CG_TEpiece_vs_rest_of_locus_lincs.pdf",height = 4,width = 3)
####################################################
par(mar=c(8,4,3,2)) 
#boxplot(summary_tepiece_rest$CG_tepiece,summary_tepiece_rest$CG_rest,outline = F,notch = T, col=c("burlywood4","burlywood2"), names=c("TE patches\n inside lincRNA","TE-free part\n of lincRNA") ,las=2, main="CG methylation\n(rosette, Col-0)",ylab="methylation level")
boxplot(summary_tepiece_rest$CG_tepiece,summary_tepiece_rest$CG_rest,CG.1001.denovo$X6909[!(CG.1001.denovo$transcript %in% linc_TE_cov_loci$gene)],   outline = F,notch = T, col=c("burlywood4","burlywood2","lightgray"), names=c("TE patches\n inside lincRNA","TE-free part\n of lincRNA","lincRNAs wo TE") ,las=2, main="CG methylation\n(rosette, Col-0)",ylab="methylation level")

wilcox.test(summary_tepiece_rest$CG_tepiece,summary_tepiece_rest$CG_rest) 
#add p values
a<-wilcox.test(summary_tepiece_rest$CG_tepiece,summary_tepiece_rest$CG_rest)
a$p.value#1.965863e-35
text("***",x=1.5,y=1,cex=0.7)
a<-wilcox.test(summary_tepiece_rest$CG_rest,CG.1001.denovo$X6909[!(CG.1001.denovo$transcript %in% linc_TE_cov_loci$gene)])
a$p.value #5.898346e-80
text("***",x=2.5,y=1,cex=0.7)
##################################################
dev.off()

#CG histogram  
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Suppl_Fig7_histogram_CG_TEpiece_minus_rest_of_locus.pdf",height = 3,width = 4)
####################################################
vect<-summary_tepiece_rest$CG_tepiece-summary_tepiece_rest$CG_rest 
#my_variable<-bb$sRNA_tepiece[grep("^(?!.*Heli)^(?!.*DNA).*LTR",bb$TE_types,perl = T)]
# Layout to split the screen
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(4,2))
# Draw the boxplot and the histogram
par(mar=c(0, 4, 1.1, 2.1))
hist(vect , col="grey" , border=T , xlab="", xlim=c(-1,1), main="CG methylation difference",ylab="number of lincRNA loci",breaks=40)
par(mar=c(0, 4, 1.1, 2.1))
boxplot(vect , horizontal=TRUE , ylim=c(-1,1), xaxt="n" , col=rgb(0.8,0.8,0,0.5) , frame=F,outline = F,notch=T)
mtext("meth.level difference: TE patches - rest of locus",line=-1,side=1)
##################################################
dev.off()                                                                                   

#CHH boxplot
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Supple_for_Fig7_CHH_TEpiece_vs_rest_of_locus_lincs.pdf",height = 4,width = 3)
####################################################
par(mar=c(8,4,3,2)) 
#boxplot(summary_tepiece_rest$CHH_tepiece,summary_tepiece_rest$CHH_rest,outline = F,notch = T, col=c("burlywood4","burlywood2"), names=c("TE patches\n inside lincRNA","TE-free part\n of lincRNA") ,las=2, main="CHH methylation\n(rosette, Col-0)",ylab="methylation level")
boxplot(summary_tepiece_rest$CHH_tepiece,summary_tepiece_rest$CHH_rest,CHH.1001.denovo$X6909[!(CHH.1001.denovo$transcript %in% linc_TE_cov_loci$gene)],   outline = F,notch = T, col=c("burlywood4","burlywood2","lightgray"), names=c("TE patches\n inside lincRNA","TE-free part\n of lincRNA","lincRNAs wo TE") ,las=2, main="CHH methylation\n(rosette, Col-0)",ylab="methylation level")

wilcox.test(summary_tepiece_rest$CG_tepiece,summary_tepiece_rest$CG_rest) #p-value = 2.604e-05 
#add p values
a<-wilcox.test(summary_tepiece_rest$CHH_tepiece,summary_tepiece_rest$CHH_rest)
a$p.value #3.189205e-07
text("**",x=1.5,y=0.15,cex=0.8)
a<-wilcox.test(summary_tepiece_rest$CHH_rest,CHH.1001.denovo$X6909[!(CHH.1001.denovo$transcript %in% linc_TE_cov_loci$gene)])
a$p.value #1.118922e-182
text("***",x=2.5,y=0.15,cex=0.8)
##########################################
dev.off()

#CHH histogram difference
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Suppl_Fig7_histogram_CHH_TEpiece_minus_rest_of_locus.pdf",height = 3,width = 4)
####################################################
vect<-summary_tepiece_rest$CHH_tepiece-summary_tepiece_rest$CHH_rest 
# Layout to split the screen
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(4,2))
# Draw the boxplot and the histogram
par(mar=c(0, 4, 1.1, 2.1))
hist(vect , col="grey" , border=T , xlab="", xlim=c(-0.1,0.1), main="CHH methylation difference",ylab="number of lincRNA loci",breaks=100)
par(mar=c(0, 4, 1.1, 2.1))
boxplot(vect , horizontal=TRUE , ylim=c(-0.1,0.1), xaxt="n" , col=rgb(0.8,0.8,0,0.5) , frame=F,outline = F,notch=T)
mtext("meth.level difference: TE patches - rest of locus",line=-1,side=1)
##################################################
dev.off()                                                                                   

#K9 boxplot
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Supple_for_Fig7_K9_reps_averaged_coverage_TEpiece_vs_rest_of_locus_lincs.pdf",height = 4,width = 3)
####################################################
par(mar=c(8,4,3,2)) 
boxplot(summary_tepiece_rest$k9_tepiece,summary_tepiece_rest$k9_rest,chip.denovo.log2$K9.6909[!(chip.denovo.log2$gene %in% linc_TE_cov_loci$gene)],   outline = F,notch = T, col=c("burlywood4","burlywood2","lightgray"), names=c("TE patches\n inside lincRNA","TE-free part\n of lincRNA","lincRNAs wo TE") ,las=2, main="H3K9me2 level\n(rosette, Col-0)",ylab="log2(ChIP/Input)")
a<-wilcox.test(summary_tepiece_rest$k9_tepiece,summary_tepiece_rest$k9_rest)
a$p.value#0.8712916
text("n.s.",x=1.5,y=1.5,cex=0.8)
a<-wilcox.test(summary_tepiece_rest$k9_rest,chip.denovo.log2$K9.6909[!(chip.denovo.log2$gene %in% linc_TE_cov_loci$gene)])
a$p.value#1.375683e-107
text("***",x=2.5,y=1.5,cex=0.8)
####################################################
dev.off()

#K9 histogram difference
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Suppl_Fig7_histogram_K9_TEpiece_minus_rest_of_locus.pdf",height = 3,width = 4)
####################################################
vect<-summary_tepiece_rest$k9_tepiece-summary_tepiece_rest$k9_rest 
# Layout to split the screen
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(4,2))
# Draw the boxplot and the histogram
par(mar=c(0, 4, 1.1, 2.1))
hist(vect , col="grey" , border=T , xlab="", xlim=c(-1,1), main="H3K9me2 level difference",ylab="number of lincRNA loci",breaks=30)
par(mar=c(0, 4, 1.1, 2.1))
boxplot(vect , horizontal=TRUE , ylim=c(-1,1), xaxt="n" , col=rgb(0.8,0.8,0,0.5) , frame=F,outline = F,notch=T)
mtext("H3K9me2 level difference=TE patches-rest of locus",line=-1,side=1)
##################################################
dev.off()  



#K9 rep1 boxplot
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Supple_for_Fig7_K9_rep1_coverage_TEpiece_vs_rest_of_locus_lincs.pdf",height = 4,width = 3)
####################################################
par(mar=c(8,4,3,2)) 
boxplot(summary_tepiece_rest$k9_tepiece_rep1,summary_tepiece_rest$k9_rest_rep1,chip.denovo.log2_all_samples$r14.rep1.H3K9me2[!(chip.denovo.log2_all_samples$gene %in% linc_TE_cov_loci$gene)],   outline = F,notch = T, col=c("burlywood4","burlywood2","lightgray"), names=c("TE patches\n inside lincRNA","TE-free part\n of lincRNA","lincRNAs wo TE") ,las=2, main="H3K9me2 level\n(rosette, Col-0, rep1)",ylab="log2(ChIP/Input)")
a<-wilcox.test(summary_tepiece_rest$k9_tepiece_rep1,summary_tepiece_rest$k9_rest_rep1)
a$p.value#0.0004473028
text("*",x=1.5,y=1.5,cex=0.8)
a<-wilcox.test(summary_tepiece_rest$k9_rest_rep1,chip.denovo.log2_all_samples$r14.rep1.H3K9me2[!(chip.denovo.log2_all_samples$gene %in% linc_TE_cov_loci$gene)])
a$p.value#1.141014e-151
text("***",x=2.5,y=1.5,cex=0.8)
####################################################
dev.off()

#K9 rep2 boxplot
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Supple_for_Fig7_K9_rep2_coverage_TEpiece_vs_rest_of_locus_lincs.pdf",height = 4,width = 3)
####################################################
par(mar=c(8,4,3,2)) 
boxplot(summary_tepiece_rest$k9_tepiece_rep2,summary_tepiece_rest$k9_rest_rep2,chip.denovo.log2_all_samples$r14.rep2.H3K9me2[!(chip.denovo.log2_all_samples$gene %in% linc_TE_cov_loci$gene)],   outline = F,notch = T, col=c("burlywood4","burlywood2","lightgray"), names=c("TE patches\n inside lincRNA","TE-free part\n of lincRNA","lincRNAs wo TE") ,las=2, main="H3K9me2 level\n(rosette, Col-0, rep2)",ylab="log2(ChIP/Input)")
a<-wilcox.test(summary_tepiece_rest$k9_tepiece_rep2,summary_tepiece_rest$k9_rest_rep2)
a$p.value#0.003187825
text("*",x=1.5,y=1.5,cex=0.8)
a<-wilcox.test(summary_tepiece_rest$k9_rest_rep2,chip.denovo.log2_all_samples$r14.rep2.H3K9me2[!(chip.denovo.log2_all_samples$gene %in% linc_TE_cov_loci$gene)])
a$p.value#1.055868e-45
text("***",x=2.5,y=1.5,cex=0.8)
####################################################
dev.off()










#############
# number if TEs expressed vs nr of lincRNAs expressed 
#1001G
a<-as.data.frame(t(denovo2021.TPMs.genes.1001G[,2:462]))

a$N_PC<-apply(a[,colnames(a) %in% denovoPC.loci$gene], 1, function(i) sum(i > 0.5))
a$N_AS<-apply(a[,names(a )%in% lncRNAs.antisense.loci$gene], 1, function(i) sum(i > 0.5))
a$N_TE<-apply(a[,names(a) %in% TE_genes.loci$gene], 1, function(i) sum(i > 0.5))
a$N_TEfrag<-apply(a[,names(a) %in% TE_frags.transcripts$gene], 1, function(i) sum(i > 0.5))

a$N_linc<-apply(a[,names(a) %in% lncRNAs.intergenic.loci$gene], 1, function(i) sum(i > 0.5))
a$N_linc_withTE<-apply(a[,names(a) %in% lincRNAs_TE_coverage.TAIR10$gene], 1, function(i) sum(i > 0.5))
a$N_linc_withoutTE<-apply(a[,names(a) %in% lncRNAs.intergenic.loci$gene & !(names(a) %in% lincRNAs_TE_coverage.TAIR10$gene)], 1, function(i) sum(i > 0.5))

a[a<=0] <- NA
a$TE_expressed_meanTPM<-apply(a[,names(a) %in% lncRNAs.intergenic.loci$gene & !(names(a) %in% lincRNAs_TE_coverage.TAIR10$gene)], 1, function(i) sum(i > 0.5))

b<- a[,names(a) %in% lncRNAs.intergenic.loci$gene & !(names(a) %in% lincRNAs_TE_coverage.TAIR10$gene)]
b$mean_TPM<-lapply(x, function (y) mean(y[y > 0]))

ab<-a[,38155:38161]
# remove outlier 
#X10010 has twice as many lincRNAs and TE genes expressed 
ab<-ab[!(rownames(ab)=="X10010"),]

library(scales)

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Supple_for_Fig7_alllinc_vs_TEgene_number_1001G.pdf",height = 3,width = 3)
####################################################
par(mar=c(4,4,3,2)) 
plot(ab$N_linc,ab$N_TE,ylim=c(40,130),xlim=c(50,140),ylab="number of TE genes expressed",xlab="number of lincRNAs expressed",main="N of TE genes vs lincRNAs expressed\n 460 accessions, rosette",pch=1,cex.main=0.9,col = alpha("black", 0.5),las=2)
text(paste("R=",round(cor(ab$N_linc,ab$N_TE),2),sep=""),x=60,y=120,cex=0.8)
####################################################
dev.off()


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Supple_for_Fig7_TElinc_vs_TEgene_number_1001G.pdf",height = 3,width = 3)
####################################################
par(mar=c(4,4,3,2)) 
plot(ab$N_linc_withTE,ab$N_TE,ylim=c(40,130),xlim=c(20,80),ylab="number of TE genes expressed",xlab="number of lincRNAs expressed",main="lincRNAs with TE patches",pch=1,cex.main=0.9,col = alpha("black", 0.5),las=2)
text(paste("R=",round(cor(ab$N_linc_withTE,ab$N_TE),2),sep=""),x=25,y=120,cex=0.8)
####################################################
dev.off()


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Supple_for_Fig7_noTElinc_vs_TEgene_number_1001G.pdf",height = 3,width = 3)
####################################################
par(mar=c(4,4,3,2)) 
plot(ab$N_linc_withoutTE,ab$N_TE,ylim=c(40,130),xlim=c(20,80),ylab="number of TE genes expressed",xlab="number of lincRNAs expressed",main="lincRNAs w/o TE patches",pch=1,cex.main=0.9,col = alpha("black", 0.5),las=2)
text(paste("R=",round(cor(ab$N_linc_withoutTE,ab$N_TE),2),sep=""),x=25,y=120,cex=0.8)
####################################################
dev.off()


################################ 
### correlation between expression level of expressed lincRNAs and expressed TE genes 

a<-as.data.frame(t(denovo2021.TPMs.genes.1001G[,2:462]))
a$N_PC<-apply(a[,colnames(a) %in% denovoPC.loci$gene], 1, function(i) sum(i > 0.5))
a$N_AS<-apply(a[,names(a )%in% lncRNAs.antisense.loci$gene], 1, function(i) sum(i > 0.5))
a$N_TE<-apply(a[,names(a) %in% TE_genes.loci$gene], 1, function(i) sum(i > 0.5))
a$N_TEfrag<-apply(a[,names(a) %in% TE_frags.transcripts$gene], 1, function(i) sum(i > 0.5))
a$N_linc<-apply(a[,names(a) %in% lncRNAs.intergenic.loci$gene], 1, function(i) sum(i > 0.5))
a$N_linc_withTE<-apply(a[,names(a) %in% lincRNAs_TE_coverage.TAIR10$gene], 1, function(i) sum(i > 0.5))
a$N_linc_withoutTE<-apply(a[,names(a) %in% lncRNAs.intergenic.loci$gene & !(names(a) %in% lincRNAs_TE_coverage.TAIR10$gene)], 1, function(i) sum(i > 0.5))

genenumbers <-a[,38155:38161]
  
  
c<-a[,names(a) %in% TE_genes.loci$gene]
c[c<=0.5] <- NA
c$TE_expressed_meanTPM<-apply(c, 1,mean, na.rm=T) 
b<-a[,names(a) %in% lncRNAs.intergenic.loci$gene]
b[b<=0.5] <- NA
b$linc_expressed_meanTPM<-apply(b, 1,mean, na.rm=T) 
d<-merge(b[,2246:2247],c[,2130:2131], by="row.names")

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Supple_for_Fig7_alllinc_vs_TEgene_mean_expression_of_expressed_genes_1001G.pdf",height = 3,width = 3)
####################################################
par(mar=c(4,4,3,2)) 
plot(d$linc_expressed_meanTPM ,d$TE_expressed_meanTPM,ylab="mean expression of expressed TE genes",xlab="mean expression of expressed lincRNAs",main="N of TE genes vs lincRNAs expressed\n 460 accessions, rosette",pch=1,cex.main=0.9,col = alpha("black", 0.5),las=2)
text(paste("R=",round(cor(d$linc_expressed_meanTPM ,d$TE_expressed_meanTPM),2),sep=""),x=15,y=7,cex=0.8)
####################################################
dev.off()


b<-a[,names(a) %in% lincRNAs_TE_coverage.TAIR10$gene]
b[b<=0.5] <- NA
b$lincwithTE_expressed_meanTPM<-apply(b, 1,mean, na.rm=T) 
d<-merge(b[,1176:1177],c[,2130:2131], by="row.names")


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Supple_for_Fig7_lincwithTE_vs_TEgene_mean_expression_of_expressed_genes_1001G.pdf",height = 3,width = 3)
####################################################
par(mar=c(4,4,3,2)) 
plot(d$lincwithTE_expressed_meanTPM ,d$TE_expressed_meanTPM,ylab="mean expression of expressed TE genes",xlab="mean expression of expressed lincRNAs",main="N of TE genes vs lincRNAs expressed\n 460 accessions, rosette",pch=1,cex.main=0.9,col = alpha("black", 0.5),las=2)
text(paste("R=",round(cor(d$lincwithTE_expressed_meanTPM ,d$TE_expressed_meanTPM),2),sep=""),x=15,y=7,cex=0.8)
####################################################
dev.off()






############################## 
# # heatmaps - accessiosn with many lincRNAs and few lincRNAs 


fewest lincs 
X9756
X7186
X9421
X9702
X6099


most lincs

X10010 - outlier 
X9532
X9944
X9933
X9696
X10016

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Suppl_pheatmap_lincRNAexpression_5mostgenes_5least_genes.pdf",height = 4,width = 4)
####################################################

tmp<-denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene ,c("X9756","X7186","X9421","X9702","X6099","X9944","X9933","X9696","X10016","X9532")]
library(pheatmap)
pheatmap (tmp[apply(tmp,1,max)>0.5,],scale = "row", labels_row = F,cluster_cols = F)
length(tmp$X9532[apply(tmp,1,max)>0.5]) #357
####################################################
dev.off()


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Suppl_pheatmap_TEgeneexpression_5mostgenes_5least_genesMATCHacessions_forlincs.pdf",height = 4,width = 4)
####################################################
tmp<-denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.loci$gene ,c("X9756","X7186","X9421","X9702","X6099","X9944","X9933","X9696","X10016","X9532")]
library(pheatmap)
pheatmap (tmp[apply(tmp,1,max)>0.5,],scale = "row", labels_row = F,cluster_cols = F)
####################################################
dev.off()

least te genes 
X5535
X9780
X9795
X7106
X6138

most te genes 

X10010- outlier 
X9944
X9565
X6917
X9543
X9554

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig7_silenc/Suppl_pheatmap_TEgeneexpression_5mostgenes_5least_genes.pdf",height = 4,width = 4)
####################################################
tmp<-denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.loci$gene ,c("X5535","X9780","X9795","X7106","X6138","X9944","X9565","X6917","X9543","X9554")]
library(pheatmap)
pheatmap (tmp[apply(tmp,1,max)>0.5,],scale = "row", labels_row = F,cluster_cols = F)
length(tmp$X5535[apply(tmp,1,max)>0.5]) #306
####################################################
dev.off()

















































