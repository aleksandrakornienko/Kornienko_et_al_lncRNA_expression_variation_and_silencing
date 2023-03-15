##########################################
###Figure 4 : epigenetic variation
##########################################




# variation between accessions 
#methylation 

# CG variation 1001G
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CGmeth_variation.1001Gdata.allgenes.pdf",height = 3,width = 3.5)
###################################################
par(mar=c(6,6,2,2)) 
pc<-CG.1001.denovo$sd[CG.1001.denovo$transcript %in% denovoPC.loci$gene]
as<- CG.1001.denovo$sd[CG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CG.1001.denovo$sd[CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CG.1001.denovo$sd[CG.1001.denovo$transcript %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE"),main="CG methylation variation among 450 accessions",ylab="standard deviation of CG methylation level",notch = T,outline = F)
mtext("all genes",side = 4)


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
###################################################
dev.off()


# CHH variation 1001G
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CHHmeth_variation.1001Gdata.allgenes.pdf",height = 3,width = 3.5)
###################################################
par(mar=c(6,6,2,2)) 
pc<-CHH.1001.denovo$sd[CG.1001.denovo$transcript %in% denovoPC.loci$gene]
as<- CHH.1001.denovo$sd[CG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CHH.1001.denovo$sd[CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CHH.1001.denovo$sd[CG.1001.denovo$transcript %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE"),main="CHH methylation variation among 450 accessions",ylab="standard deviation of CHH methylation level",notch = T,outline = F)
mtext("all genes",side = 4)
###################################################
##add p values 
###################
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
###################
dev.off()


# CG inter variation 1001G new 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CGmeth_variation.1001GNEWdata.allgenes.pdf",height = 3,width = 3.5)
###################################################
par(mar=c(6,6,2,2)) 
pc<-CG.1001new.denovo$sd_of_means[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]
as<- CG.1001new.denovo$sd_of_means[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CG.1001new.denovo$sd_of_means[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CG.1001new.denovo$sd_of_means[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE"),main="CG methylation variation\n28 accessions",ylab="standard deviation of CG methylation level",notch = T,outline = F)
mtext("all genes",side = 4)

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
###################################################
dev.off()


# CHH inter variation 1001G new
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CHHmeth_variation.1001GNEWdata.allgenes.pdf",height = 3,width = 3.5)
###################################################
par(mar=c(6,6,2,2)) 
pc<-CHH.1001new.denovo$sd[CHH.1001new.denovo$transcript %in% denovoPC.loci$gene]
as<- CHH.1001new.denovo$sd[CHH.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CHH.1001new.denovo$sd[CHH.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CHH.1001new.denovo$sd[CHH.1001new.denovo$transcript %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE"),main="CHH methylation variation\n28 accessions",ylab="standard deviation of CHH methylation level",notch = T,outline = F)
mtext("all genes",side = 4)
###################################################
##add p values 
###################
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
###################
dev.off()


# CG intravariation 1001G new 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CGmeth_intravariation.1001GNEWdata.allgenes.pdf",height = 3,width = 3.5)
###################################################
par(mar=c(6,6,2,2)) 
pc<-CG.1001new.denovo$mean_intra_sd[CG.1001new.denovo$transcript %in% denovoPC.loci$gene]
as<- CG.1001new.denovo$mean_intra_sd[CG.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CG.1001new.denovo$mean_intra_sd[CG.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CG.1001new.denovo$mean_intra_sd[CG.1001new.denovo$transcript %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE"),main="CG methylation variation\n2-4reps",ylab="standard deviation of CG methylation level",notch = T,outline = F,ylim=c(0,0.5))
mtext("all genes",side = 4)

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
###################################################
dev.off()


# CHH intravariation 1001G new
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_CHHmeth_intravariation.1001GNEWdata.allgenes.pdf",height = 3,width = 3.5)
###################################################
par(mar=c(6,6,2,2)) 
pc<-CHH.1001new.denovo$mean_intra_sd[CHH.1001new.denovo$transcript %in% denovoPC.loci$gene]
as<- CHH.1001new.denovo$mean_intra_sd[CHH.1001new.denovo$transcript %in% lncRNAs.antisense.loci$gene]
linc<-CHH.1001new.denovo$mean_intra_sd[CHH.1001new.denovo$transcript %in% lncRNAs.intergenic.loci$gene]
te<-CHH.1001new.denovo$mean_intra_sd[CHH.1001new.denovo$transcript %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE"),main="CHH methylation variation\n2-4reps",ylab="standard deviation of CHH methylation level",notch = T,outline = F,ylim=c(0,0.5))
mtext("all genes",side = 4)
###################################################
##add p values 
###################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.3)

a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.3)

a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.3)
###################
dev.off()
















a<- CHH.1001.denovo[CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene,c("transcript","sd")]
a$gene<-a$transcript
b<-denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene,c("gene","variance")]
a<-merge(a,b,by="gene")
plot(a$sd,a$variance)



# histone marks (normalized by quantiles)
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_H1_variation.allgenes.pdf",height = 3,width = 3.5)
par(mar=c(6,6,2,2)) 
###################################################
pc<-chip.denovo.quantstan$sd.hist1[chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as<-chip.denovo.quantstan$sd.hist1[chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo.quantstan$sd.hist1[chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo.quantstan$sd.hist1[chip.denovo.quantstan$gene %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),main="H1 variation among 14 accessions",ylab="standard deviation of normalized ChIP signal"
        ,notch = T,outline = F)
mtext("all genes",side = 4)
###################################################
##add p values 
###################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.7)

a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.7)

a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.7)
###################
dev.off()

# K9
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K9_variation.allgenes.pdf",height = 3,width = 3.5)
par(mar=c(6,6,2,2)) 
###################################################
pc<-chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as<-chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),main="H3K9me2 variation\namong 13 accessions",ylab="SD of normalized ChIP signal"
        ,notch = T,outline = F)
mtext("all genes",side = 4)
###################################################
##add p values 
###################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.7)

a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.7)

a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.7)
###################
dev.off()


# K27
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K27_variation.allgenes.pdf",height = 3,width = 3.5)
###################################################
par(mar=c(6,6,2,2)) 
pc<-chip.denovo.quantstan$sd.key27[chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as<-chip.denovo.quantstan$sd.key27[chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo.quantstan$sd.key27[chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo.quantstan$sd.key27[chip.denovo.quantstan$gene %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),main="H3K27me3 variation\namong 12 accessions",ylab="SD of normalized ChIP signal"
        ,notch = T,outline = F)
mtext("all genes",side = 4)
###################################################
##add p values 
###################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.7)

a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.7)

a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.7)
###################
dev.off()



# K36
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_K36_variation.allgenes.pdf",height = 3,width = 3.5)
###################################################
par(mar=c(6,6,2,2)) 
pc<-chip.denovo.quantstan$sd.key36[chip.denovo.quantstan$gene %in% denovoPC.loci$gene]
as<-chip.denovo.quantstan$sd.key36[chip.denovo.quantstan$gene %in% lncRNAs.antisense.loci$gene]
linc<-chip.denovo.quantstan$sd.key36[chip.denovo.quantstan$gene %in% lncRNAs.intergenic.loci$gene]
te<-chip.denovo.quantstan$sd.key36[chip.denovo.quantstan$gene %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),main="H3K36me3 variation\namong 14 accessions",ylab="SD of normalized ChIP signal"
        ,notch = T,outline = F)
mtext("all genes",side = 4)
###################################################
##add p values 
###################
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
###################
dev.off()



#small RNAs 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/Fig3_boxplot_24nt_variation.1001Gdata.allgenes.pdf",height = 3,width = 3.5)
###################################################
par(mar=c(6,6,2,2)) 
pc<-sRNA.24nt.denovo2021.RPM$variance[sRNA.24nt.denovo2021.RPM$gene %in% denovoPC.loci$gene]
as<- sRNA.24nt.denovo2021.RPM$variance[sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.antisense.loci$gene]
linc<-sRNA.24nt.denovo2021.RPM$variance[sRNA.24nt.denovo2021.RPM$gene %in% lncRNAs.intergenic.loci$gene]
te<-sRNA.24nt.denovo2021.RPM$variance[sRNA.24nt.denovo2021.RPM$gene %in% TE_genes.loci$gene]
boxplot(pc,as,linc,te,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC","AS","linc","TE"),main="24nt sRNA coverage\nvariation among 14 accessions",ylab="RPM's coefficient of variance",notch = T,outline = F)
mtext("all genes",side = 4)
###################################################
##add p values 
###################
a<-wilcox.test(sample(pc,2000),sample(as,2000))
b<-wilcox.test(sample(pc,2000),sample(as,2000))
c<-wilcox.test(sample(pc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=3.5)

a<-wilcox.test(sample(linc,2000),sample(as,2000))
b<-wilcox.test(sample(linc,2000),sample(as,2000))
c<-wilcox.test(sample(linc,2000),sample(as,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=3.5)

a<-wilcox.test(sample(linc,2000),sample(te,2000))
b<-wilcox.test(sample(linc,2000),sample(te,2000))
c<-wilcox.test(sample(linc,2000),sample(te,2000))
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=3.5)
###################
dev.off()
























############################################
## explaining variation by epigenetics######
############################################

# expression change explained by methylation 

#lincRNAs

LINC_epigenetic_factors<-lncRNAs.intergenic.loci
rownames(LINC_epigenetic_factors)<-LINC_epigenetic_factors$gene
LINC_epigenetic_factors$TE_content_TAIR10<-0
LINC_epigenetic_factors$mean_CN<-0
LINC_epigenetic_factors$sd_CN<-0
LINC_epigenetic_factors$mean_TPM_1001Gnew<-0
LINC_epigenetic_factors$max_TPM_1001Gnew<-0
LINC_epigenetic_factors$variance_1001Gnew<-0

LINC_epigenetic_factors$cor_H1_TPM<-0
LINC_epigenetic_factors$cor_K4_TPM<-0
LINC_epigenetic_factors$cor_K9_TPM<-0
LINC_epigenetic_factors$cor_K27_TPM<-0
LINC_epigenetic_factors$cor_K36_TPM<-0

LINC_epigenetic_factors$mean_TPM_EC_S<-0
LINC_epigenetic_factors$max_TPM_EC_S<-0
LINC_epigenetic_factors$variance_EC_S<-0

LINC_epigenetic_factors$mean_TPM_EC_R<-0
LINC_epigenetic_factors$max_TPM_EC_R<-0
LINC_epigenetic_factors$variance_EC_R<-0

LINC_epigenetic_factors$mean_TPM_EC_F<-0
LINC_epigenetic_factors$max_TPM_EC_F<-0
LINC_epigenetic_factors$variance_EC_F<-0

LINC_epigenetic_factors$mean_TPM_EC_P<-0
LINC_epigenetic_factors$max_TPM_EC_P<-0
LINC_epigenetic_factors$variance_EC_P<-0

LINC_epigenetic_factors$cor_24nt_vs_TPM_pearson<-0
LINC_epigenetic_factors$cor_24nt_vs_TPM_S<-0
LINC_epigenetic_factors$cor_24nt_vs_TPM_R<-0
LINC_epigenetic_factors$cor_24nt_vs_TPM_F<-0
LINC_epigenetic_factors$cor_24nt_vs_TPM_P<-0


LINC_epigenetic_factors$mean_TPM_1001G<-0
LINC_epigenetic_factors$max_TPM_1001G<-0
LINC_epigenetic_factors$variance_1001G<-0

LINC_epigenetic_factors$CG_GB_vs_TPM_pearson<-0
LINC_epigenetic_factors$CG_TSS_vs_TPM_pearson<-0

LINC_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p<-0
LINC_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p<-0

LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p<-0
LINC_epigenetic_factors$CHH_GB_vs_TPM_pearson<-0
LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p<-0
LINC_epigenetic_factors$CHH_TSS_vs_TPM_pearson<-0

LINC_epigenetic_factors$CG_GB_TPMfoldchange<-0
LINC_epigenetic_factors$CG_TSS_TPMfoldchange<-0
LINC_epigenetic_factors$CHH_GB_TPMfoldchange<-0
LINC_epigenetic_factors$CHH_TSS_TPMfoldchange<-0


for ( i in 1:length (LINC_epigenetic_factors$gene)) {
  
  
  linc=as.character( LINC_epigenetic_factors$gene[i])
  print(linc)
  
  LINC_epigenetic_factors$mean_TPM_EC_S[i]<-denovo2021.TPMs.genes.ERACAPS$mean.seedl[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  LINC_epigenetic_factors$max_TPM_EC_S[i]<-denovo2021.TPMs.genes.ERACAPS$max.seedl[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  LINC_epigenetic_factors$variance_EC_S[i]<-denovo2021.TPMs.genes.ERACAPS$var.seedl[denovo2021.TPMs.genes.ERACAPS$gene==linc]
    LINC_epigenetic_factors$mean_TPM_EC_R[i]<-denovo2021.TPMs.genes.ERACAPS$mean.rosette[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  LINC_epigenetic_factors$max_TPM_EC_R[i]<-denovo2021.TPMs.genes.ERACAPS$max.rosette[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  LINC_epigenetic_factors$variance_EC_R[i]<-denovo2021.TPMs.genes.ERACAPS$var.rosette[denovo2021.TPMs.genes.ERACAPS$gene==linc]
    LINC_epigenetic_factors$mean_TPM_EC_F[i]<-denovo2021.TPMs.genes.ERACAPS$mean.flowers[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  LINC_epigenetic_factors$max_TPM_EC_F[i]<-denovo2021.TPMs.genes.ERACAPS$max.flowers[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  LINC_epigenetic_factors$variance_EC_F[i]<-denovo2021.TPMs.genes.ERACAPS$var.flowers[denovo2021.TPMs.genes.ERACAPS$gene==linc]
    LINC_epigenetic_factors$mean_TPM_EC_P[i]<-denovo2021.TPMs.genes.ERACAPS$mean.pollen[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  LINC_epigenetic_factors$max_TPM_EC_P[i]<-denovo2021.TPMs.genes.ERACAPS$max.pollen[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  LINC_epigenetic_factors$variance_EC_P[i]<-denovo2021.TPMs.genes.ERACAPS$var.pollen[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  ###########################################
  ####CHIP-SEQ
  ############ correlation with histone marks 
  # dataset 1001Gnew 13 accessions TPM-CHiP
  
  a<-as.data.frame(t(denovo2021.TPMs.genes.1001Gnew[linc,grep("mean.",names(denovo2021.TPMs.genes.1001Gnew))]))
  a$accession<-unlist(lapply(strsplit(as.character(rownames(a)),".", fixed = T), "[",2))
  a<-a[1:28,]
  b<-as.data.frame(t(chip.denovo.quantstan[linc,2:68]))
  b$accession<-unlist(lapply(strsplit(as.character(rownames(b)),".", fixed = T), "[",2))
  b$chip_sample<-rownames(b)
  b$mark<-unlist(lapply(strsplit(as.character(rownames(b)),".", fixed = T), "[",1))
  aa<-merge(a,b,by="accession")
  
  LINC_epigenetic_factors$TE_content_TAIR10[i]<- linc_TE_cov_all_loci_2cols$coverage[linc_TE_cov_all_loci_2cols$gene==linc]
  LINC_epigenetic_factors$mean_CN[i]<-CN_linc_27genomes$CN_mean[rownames(CN_linc_27genomes)==linc]
  LINC_epigenetic_factors$sd_CN[i]<-CN_linc_27genomes$CN_sd[rownames(CN_linc_27genomes)==linc]
  
  LINC_epigenetic_factors$mean_TPM_1001Gnew[i]<-mean(aa[,2])
  LINC_epigenetic_factors$max_TPM_1001Gnew[i]<-max(aa[,2])
  LINC_epigenetic_factors$variance_1001Gnew[i]<-denovo2021.TPMs.genes.1001Gnew$variance_of_means[denovo2021.TPMs.genes.1001Gnew$gene==linc]
  
  LINC_epigenetic_factors$cor_H1_TPM[i]<-cor(aa[aa$mark=="H1",3],aa[aa$mark=="H1",2])
  LINC_epigenetic_factors$cor_K4_TPM[i]<-cor(aa[aa$mark=="K4",3],aa[aa$mark=="K4",2])
  LINC_epigenetic_factors$cor_K9_TPM[i]<-cor(aa[aa$mark=="K9",3],aa[aa$mark=="K9",2])
  LINC_epigenetic_factors$cor_K27_TPM[i]<-cor(aa[aa$mark=="K27",3],aa[aa$mark=="K27",2])
  LINC_epigenetic_factors$cor_K36_TPM[i]<-cor(aa[aa$mark=="K36",3],aa[aa$mark=="K36",2])
  
  
  #small RNAs 
  
  a<-as.data.frame(t(denovo2021.TPMs.genes.ERACAPS[linc,2:97]))
  a$accession<-unlist(lapply(strsplit(as.character(rownames(a)),".", fixed = T), "[",2))
  a$tissue<-lapply(strsplit(as.character(row.names(a)),".", fixed = T), "[",1)
  b<-as.data.frame(t(sRNA.24nt.denovo2021.RPM[linc,7:20]))
  b$accession<-unlist(lapply(strsplit(as.character(rownames(b)),"X", fixed = T), "[",2))
  aa<-merge(a,b,by="accession")
  
  names(aa)<-c("accession","TPM","tissue","sRNA_RPM")
  
  LINC_epigenetic_factors$cor_24nt_vs_TPM_pearson[i]<-cor(aa[,2],aa[,4])
  LINC_epigenetic_factors$cor_24nt_vs_TPM_S[i]<-cor(aa[aa$tissue=="S",2],aa[aa$tissue=="S",4])
  LINC_epigenetic_factors$cor_24nt_vs_TPM_R[i]<-cor(aa[aa$tissue=="R",2],aa[aa$tissue=="R",4])
  LINC_epigenetic_factors$cor_24nt_vs_TPM_F[i]<-cor(aa[aa$tissue=="F",2],aa[aa$tissue=="F",4])
  LINC_epigenetic_factors$cor_24nt_vs_TPM_P[i]<-cor(aa[aa$tissue=="P",2],aa[aa$tissue=="P",4])
  
  #####################
  ##METHYLATION
  ### CHH and CG
  ## 1001G dataset ~450  accessions
  
  #######
  #CHH  #
  #######
  a<-as.data.frame(t(denovo2021.TPMs.genes.1001G[linc,2:462]))
  a$accession<-rownames(a)
  b<-as.data.frame(t(CHH.linc[linc,7:450]))
  b$accession<-rownames(b)
  a<-merge(a,b,by="accession")
  b<-as.data.frame(t(CHH.linc_TSS[linc,7:450]))
  b$accession<-rownames(b)
  a<-merge(a,b,by="accession")
  b<-as.data.frame(t(CHH.linc_TES[linc,7:450]))
  b$accession<-rownames(b)
  a<-merge(a,b,by="accession")
  names(a)<-c("accession","TPM","CHH","CHH_TSS","CHH_TES")
  
  
  LINC_epigenetic_factors$mean_TPM_1001G[i]<-mean(a$TPM)
  LINC_epigenetic_factors$max_TPM_1001G[i]<-max(a$TPM)
  LINC_epigenetic_factors$variance_1001G[i]<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene==linc]
  
  
  if ( length(a$TPM[a$CHH<0.01 & !is.na(a$CHH)])>10 & length(a$TPM[a$CHH>=0.01&!is.na(a$CHH)])>10 & max(a$TPM>0.5) ) {
    s1<-wilcox.test(a$TPM[a$CHH<0.01],a$TPM[a$CHH>=0.01])
    LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p[i]<-(-log10(s1$p.value))
    LINC_epigenetic_factors$CHH_GB_vs_TPM_pearson[i]<-cor(a$TPM,a$CHH,use="complete.obs")
    LINC_epigenetic_factors$CHH_GB_TPMfoldchange[i]<-median(a$TPM[a$CHH>=0.01]+1)/median(a$TPM[a$CHH<0.01 ]+1)
  } else {if(length(a$TPM[!is.na(a$CHH)])>10) {
    LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    LINC_epigenetic_factors$CHH_GB_vs_TPM_pearson[i]<-cor(a$TPM,a$CHH,use="complete.obs")
  }   else {
    LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    LINC_epigenetic_factors$CHH_GB_vs_TPM_pearson[i]<-NA
  }
    
  }
  
  if ( length(a$TPM[a$CHH_TSS<0.01&!is.na(a$CHH_TSS)])>10 & length(a$TPM[a$CHH_TSS>=0.01&!is.na(a$CHH_TSS)])>10& max(a$TPM>0.5) ) {
    s1<-wilcox.test(a$TPM[a$CHH_TSS<0.01],a$TPM[a$CHH_TSS>=0.01])
    LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p[i]<-(-log10(s1$p.value))
    LINC_epigenetic_factors$CHH_TSS_vs_TPM_pearson[i]<-cor(a$TPM,a$CHH_TSS,use="complete.obs")
    LINC_epigenetic_factors$CHH_TSS_TPMfoldchange[i]<-median(a$TPM[a$CHH_TSS>=0.01]+1)/median(a$TPM[a$CHH_TSS<0.01 ]+1)
  } else { if(length(a$TPM[!is.na(a$CHH_TSS)])>10) {
    LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    LINC_epigenetic_factors$CHH_TSS_vs_TPM_pearson[i]<-cor(a$TPM,a$CHH_TSS,use="complete.obs")
  }   else {
    LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    LINC_epigenetic_factors$CHH_TSS_vs_TPM_pearson[i]<-NA
  }
  }
  
  #######
  #CG   #
  #######
  a<-as.data.frame(t(denovo2021.TPMs.genes.1001G[linc,2:462]))
  a$accession<-rownames(a)
  b<-as.data.frame(t(CG.linc[linc,7:450]))
  b$accession<-rownames(b)
  a<-merge(a,b,by="accession")
  b<-as.data.frame(t(CG.linc_TSS[linc,7:450]))
  b$accession<-rownames(b)
  a<-merge(a,b,by="accession")
  b<-as.data.frame(t(CG.linc_TES[linc,7:450]))
  b$accession<-rownames(b)
  a<-merge(a,b,by="accession")
  names(a)<-c("accession","TPM","CG","CG_TSS","CG_TES")
  if ( length(a$TPM[a$CG<0.5 & !is.na(a$CG)])>10 & length(a$TPM[a$CG>0.5&!is.na(a$CG)])>10 & max(a$TPM>0.5) ) {
    s1<-wilcox.test(a$TPM[a$CG<0.5 ],a$TPM[a$CG>0.5])
    LINC_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p[i]<-(-log10(s1$p.value))
    LINC_epigenetic_factors$CG_GB_vs_TPM_pearson[i]<-cor(a$TPM,a$CG,use="complete.obs")
    LINC_epigenetic_factors$CG_GB_TPMfoldchange[i]<-median(a$TPM[a$CG>0.5]+1)/median(a$TPM[a$CG<0.5 ]+1)
  } else {if(length(a$TPM[!is.na(a$CG)])>10) {
    LINC_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    LINC_epigenetic_factors$CG_GB_vs_TPM_pearson[i]<-cor(a$TPM,a$CG,use="complete.obs")
  }   else {
    LINC_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    LINC_epigenetic_factors$CG_GB_vs_TPM_pearson[i]<-NA
  }
    
  }
  
  if ( length(a$TPM[a$CG_TSS<0.5 &!is.na(a$CG_TSS)])>10 & length(a$TPM[a$CG_TSS>0.5 &!is.na(a$CG_TSS)])>10& max(a$TPM>0.5) ) {
    s1<-wilcox.test(a$TPM[a$CG_TSS<0.5 ],a$TPM[a$CG_TSS>0.5])
    LINC_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p[i]<-(-log10(s1$p.value))
    LINC_epigenetic_factors$CG_TSS_vs_TPM_pearson[i]<-cor(a$TPM,a$CG_TSS,use="complete.obs")
    LINC_epigenetic_factors$CG_TSS_TPMfoldchange[i]<-median(a$TPM[a$CG_TSS>0.5]+1)/median(a$TPM[a$CG_TSS<0.5 ]+1)
  } else { if(length(a$TPM[!is.na(a$CG_TSS)])>10) {
    LINC_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    LINC_epigenetic_factors$CG_TSS_vs_TPM_pearson[i]<-cor(a$TPM,a$CG_TSS,use="complete.obs")
  }   else {
    LINC_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    LINC_epigenetic_factors$CG_TSS_vs_TPM_pearson[i]<-NA
  }
  }
}




######################################
#AS lncRNAs - explaining variation by epigenetics 
#########################################
AS_epigenetic_factors<-lncRNAs.antisense.loci
rownames(AS_epigenetic_factors)<-AS_epigenetic_factors$gene
AS_epigenetic_factors$TE_content_TAIR10<-0
AS_epigenetic_factors$mean_CN<-0
AS_epigenetic_factors$sd_CN<-0
AS_epigenetic_factors$mean_TPM_1001Gnew<-0
AS_epigenetic_factors$max_TPM_1001Gnew<-0
AS_epigenetic_factors$variance_1001Gnew<-0

AS_epigenetic_factors$cor_H1_TPM<-0
AS_epigenetic_factors$cor_K4_TPM<-0
AS_epigenetic_factors$cor_K9_TPM<-0
AS_epigenetic_factors$cor_K27_TPM<-0
AS_epigenetic_factors$cor_K36_TPM<-0

AS_epigenetic_factors$mean_TPM_EC_S<-0
AS_epigenetic_factors$max_TPM_EC_S<-0
AS_epigenetic_factors$variance_EC_S<-0

AS_epigenetic_factors$mean_TPM_EC_R<-0
AS_epigenetic_factors$max_TPM_EC_R<-0
AS_epigenetic_factors$variance_EC_R<-0

AS_epigenetic_factors$mean_TPM_EC_F<-0
AS_epigenetic_factors$max_TPM_EC_F<-0
AS_epigenetic_factors$variance_EC_F<-0

AS_epigenetic_factors$mean_TPM_EC_P<-0
AS_epigenetic_factors$max_TPM_EC_P<-0
AS_epigenetic_factors$variance_EC_P<-0

AS_epigenetic_factors$cor_24nt_vs_TPM_pearson<-0
AS_epigenetic_factors$cor_24nt_vs_TPM_S<-0
AS_epigenetic_factors$cor_24nt_vs_TPM_R<-0
AS_epigenetic_factors$cor_24nt_vs_TPM_F<-0
AS_epigenetic_factors$cor_24nt_vs_TPM_P<-0


AS_epigenetic_factors$mean_TPM_1001G<-0
AS_epigenetic_factors$max_TPM_1001G<-0
AS_epigenetic_factors$variance_1001G<-0

AS_epigenetic_factors$CG_GB_vs_TPM_pearson<-0
AS_epigenetic_factors$CG_TSS_vs_TPM_pearson<-0

AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p<-0
AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p<-0

AS_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p<-0
AS_epigenetic_factors$CHH_GB_vs_TPM_pearson<-0
AS_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p<-0
AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson<-0

AS_epigenetic_factors$CG_GB_TPMfoldchange<-0
AS_epigenetic_factors$CG_TSS_TPMfoldchange<-0
AS_epigenetic_factors$CHH_GB_TPMfoldchange<-0
AS_epigenetic_factors$CHH_TSS_TPMfoldchange<-0

for ( i in 1:length (AS_epigenetic_factors$gene)) {
  
  
  linc=as.character( AS_epigenetic_factors$gene[i])
  print(linc)
  
  AS_epigenetic_factors$mean_TPM_EC_S[i]<-denovo2021.TPMs.genes.ERACAPS$mean.seedl[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  AS_epigenetic_factors$max_TPM_EC_S[i]<-denovo2021.TPMs.genes.ERACAPS$max.seedl[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  AS_epigenetic_factors$variance_EC_S[i]<-denovo2021.TPMs.genes.ERACAPS$var.seedl[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  AS_epigenetic_factors$mean_TPM_EC_R[i]<-denovo2021.TPMs.genes.ERACAPS$mean.rosette[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  AS_epigenetic_factors$max_TPM_EC_R[i]<-denovo2021.TPMs.genes.ERACAPS$max.rosette[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  AS_epigenetic_factors$variance_EC_R[i]<-denovo2021.TPMs.genes.ERACAPS$var.rosette[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  AS_epigenetic_factors$mean_TPM_EC_F[i]<-denovo2021.TPMs.genes.ERACAPS$mean.flowers[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  AS_epigenetic_factors$max_TPM_EC_F[i]<-denovo2021.TPMs.genes.ERACAPS$max.flowers[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  AS_epigenetic_factors$variance_EC_F[i]<-denovo2021.TPMs.genes.ERACAPS$var.flowers[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  AS_epigenetic_factors$mean_TPM_EC_P[i]<-denovo2021.TPMs.genes.ERACAPS$mean.pollen[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  AS_epigenetic_factors$max_TPM_EC_P[i]<-denovo2021.TPMs.genes.ERACAPS$max.pollen[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  AS_epigenetic_factors$variance_EC_P[i]<-denovo2021.TPMs.genes.ERACAPS$var.pollen[denovo2021.TPMs.genes.ERACAPS$gene==linc]
  ###########################################
  ####CHIP-SEQ
  ############ correlation with histone marks 
  # dataset 1001Gnew 13 accessions TPM-CHiP
  
  a<-as.data.frame(t(denovo2021.TPMs.genes.1001Gnew[linc,grep("mean.",names(denovo2021.TPMs.genes.1001Gnew))]))
  a$accession<-unlist(lapply(strsplit(as.character(rownames(a)),".", fixed = T), "[",2))
  a<-a[1:28,]
  b<-as.data.frame(t(chip.denovo.quantstan[linc,2:68]))
  b$accession<-unlist(lapply(strsplit(as.character(rownames(b)),".", fixed = T), "[",2))
  b$chip_sample<-rownames(b)
  b$mark<-unlist(lapply(strsplit(as.character(rownames(b)),".", fixed = T), "[",1))
  aa<-merge(a,b,by="accession")
  
  linc1<-sub("-",".",linc)
  #AS_epigenetic_factors$TE_content_TAIR10[i]<- AS_TE_coverage_anystrand$TAIR10[AS_TE_coverage_anystrand$gene==linc|AS_TE_coverage_anystrand$gene==linc1]
  AS_epigenetic_factors$mean_CN[i]<-CN_as_27genomes$CN_mean[rownames(CN_as_27genomes)==linc|rownames(CN_as_27genomes)==linc1]
  AS_epigenetic_factors$sd_CN[i]<-CN_as_27genomes$CN_sd[rownames(CN_as_27genomes)==linc|rownames(CN_as_27genomes)==linc1]
  
  AS_epigenetic_factors$mean_TPM_1001Gnew[i]<-mean(aa[,2])
  AS_epigenetic_factors$max_TPM_1001Gnew[i]<-max(aa[,2])
  AS_epigenetic_factors$variance_1001Gnew[i]<-denovo2021.TPMs.genes.1001Gnew$variance_of_means[denovo2021.TPMs.genes.1001Gnew$gene==linc]
  
  AS_epigenetic_factors$cor_H1_TPM[i]<-cor(aa[aa$mark=="H1",3],aa[aa$mark=="H1",2])
  AS_epigenetic_factors$cor_K4_TPM[i]<-cor(aa[aa$mark=="K4",3],aa[aa$mark=="K4",2])
  AS_epigenetic_factors$cor_K9_TPM[i]<-cor(aa[aa$mark=="K9",3],aa[aa$mark=="K9",2])
  AS_epigenetic_factors$cor_K27_TPM[i]<-cor(aa[aa$mark=="K27",3],aa[aa$mark=="K27",2])
  AS_epigenetic_factors$cor_K36_TPM[i]<-cor(aa[aa$mark=="K36",3],aa[aa$mark=="K36",2])
  
  
  #small RNAs 
  
  a<-as.data.frame(t(denovo2021.TPMs.genes.ERACAPS[linc,2:97]))
  a$accession<-unlist(lapply(strsplit(as.character(rownames(a)),".", fixed = T), "[",2))
  a$tissue<-lapply(strsplit(as.character(row.names(a)),".", fixed = T), "[",1)
  b<-as.data.frame(t(sRNA.24nt.denovo2021.RPM[linc,7:20]))
  b$accession<-unlist(lapply(strsplit(as.character(rownames(b)),"X", fixed = T), "[",2))
  aa<-merge(a,b,by="accession")
  
  names(aa)<-c("accession","TPM","tissue","sRNA_RPM")
  
  AS_epigenetic_factors$cor_24nt_vs_TPM_pearson[i]<-cor(aa[,2],aa[,4])
  AS_epigenetic_factors$cor_24nt_vs_TPM_S[i]<-cor(aa[aa$tissue=="S",2],aa[aa$tissue=="S",4])
  AS_epigenetic_factors$cor_24nt_vs_TPM_R[i]<-cor(aa[aa$tissue=="R",2],aa[aa$tissue=="R",4])
  AS_epigenetic_factors$cor_24nt_vs_TPM_F[i]<-cor(aa[aa$tissue=="F",2],aa[aa$tissue=="F",4])
  AS_epigenetic_factors$cor_24nt_vs_TPM_P[i]<-cor(aa[aa$tissue=="P",2],aa[aa$tissue=="P",4])
  
  #####################
  ##METHYLATION
  ### CHH and CG
  ## 1001G dataset ~450  accessions
  
  #######
  #CHH  #
  #######
  a<-as.data.frame(t(denovo2021.TPMs.genes.1001G[linc,2:462]))
  a$accession<-rownames(a)
  b<-as.data.frame(t(CHH.as[linc,7:450]))
  b$accession<-rownames(b)
  a<-merge(a,b,by="accession")
  b<-as.data.frame(t(CHH.as_TSS[linc,7:450]))
  b$accession<-rownames(b)
  a<-merge(a,b,by="accession")
  b<-as.data.frame(t(CHH.as_TES[linc,7:450]))
  b$accession<-rownames(b)
  a<-merge(a,b,by="accession")
  names(a)<-c("accession","TPM","CHH","CHH_TSS","CHH_TES")
  
  
  AS_epigenetic_factors$mean_TPM_1001G[i]<-mean(a$TPM)
  AS_epigenetic_factors$max_TPM_1001G[i]<-max(a$TPM)
  AS_epigenetic_factors$variance_1001G[i]<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene==linc]
  
  
  if ( length(a$TPM[a$CHH<0.01 & !is.na(a$CHH)])>10 & length(a$TPM[a$CHH>=0.01&!is.na(a$CHH)])>10 & max(a$TPM>0.5) ) {
    s1<-wilcox.test(a$TPM[a$CHH<0.01],a$TPM[a$CHH>=0.01])
    AS_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p[i]<-(-log10(s1$p.value))
    AS_epigenetic_factors$CHH_GB_vs_TPM_pearson[i]<-cor(a$TPM,a$CHH,use="complete.obs")
    AS_epigenetic_factors$CHH_GB_TPMfoldchange[i]<-median(a$TPM[a$CHH>=0.01]+1)/median(a$TPM[a$CHH<0.01 ]+1)
  } else {if(length(a$TPM[!is.na(a$CHH)])>10) {
    AS_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    AS_epigenetic_factors$CHH_GB_vs_TPM_pearson[i]<-cor(a$TPM,a$CHH,use="complete.obs")
  }   else {
    AS_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    AS_epigenetic_factors$CHH_GB_vs_TPM_pearson[i]<-NA
  }
    
  }
  
  if ( length(a$TPM[a$CHH_TSS<0.01&!is.na(a$CHH_TSS)])>10 & length(a$TPM[a$CHH_TSS>=0.01&!is.na(a$CHH_TSS)])>10& max(a$TPM>0.5) ) {
    s1<-wilcox.test(a$TPM[a$CHH_TSS<0.01],a$TPM[a$CHH_TSS>=0.01])
    AS_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p[i]<-(-log10(s1$p.value))
    AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson[i]<-cor(a$TPM,a$CHH_TSS,use="complete.obs")
    AS_epigenetic_factors$CHH_TSS_TPMfoldchange[i]<-median(a$TPM[a$CHH_TSS>=0.01]+1)/median(a$TPM[a$CHH_TSS<0.01 ]+1)
  } else { if(length(a$TPM[!is.na(a$CHH_TSS)])>10) {
    AS_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson[i]<-cor(a$TPM,a$CHH_TSS,use="complete.obs")
  }   else {
    AS_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson[i]<-NA
  }
  }
  
  #######
  #CG   #
  #######
  a<-as.data.frame(t(denovo2021.TPMs.genes.1001G[linc,2:462]))
  a$accession<-rownames(a)
  b<-as.data.frame(t(CG.as[linc,7:450]))
  b$accession<-rownames(b)
  a<-merge(a,b,by="accession")
  b<-as.data.frame(t(CG.as_TSS[linc,7:450]))
  b$accession<-rownames(b)
  a<-merge(a,b,by="accession")
  b<-as.data.frame(t(CG.as_TES[linc,7:450]))
  b$accession<-rownames(b)
  a<-merge(a,b,by="accession")
  names(a)<-c("accession","TPM","CG","CG_TSS","CG_TES")
  if ( length(a$TPM[a$CG<0.5 & !is.na(a$CG)])>10 & length(a$TPM[a$CG>0.5&!is.na(a$CG)])>10 & max(a$TPM>0.5) ) {
    s1<-wilcox.test(a$TPM[a$CG<0.5 ],a$TPM[a$CG>0.5])
    AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p[i]<-(-log10(s1$p.value))
    AS_epigenetic_factors$CG_GB_vs_TPM_pearson[i]<-cor(a$TPM,a$CG,use="complete.obs")
    AS_epigenetic_factors$CG_GB_TPMfoldchange[i]<-median(a$TPM[a$CG>0.5]+1)/median(a$TPM[a$CG<0.5 ]+1)
    
  } else {if(length(a$TPM[!is.na(a$CG)])>10) {
    AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    AS_epigenetic_factors$CG_GB_vs_TPM_pearson[i]<-cor(a$TPM,a$CG,use="complete.obs")
  }   else {
    AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    AS_epigenetic_factors$CG_GB_vs_TPM_pearson[i]<-NA
  }
    
  }
  
  if ( length(a$TPM[a$CG_TSS<0.5 &!is.na(a$CG_TSS)])>10 & length(a$TPM[a$CG_TSS>0.5 &!is.na(a$CG_TSS)])>10& max(a$TPM>0.5) ) {
    s1<-wilcox.test(a$TPM[a$CG_TSS<0.5 ],a$TPM[a$CG_TSS>0.5])
    AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p[i]<-(-log10(s1$p.value))
    AS_epigenetic_factors$CG_TSS_vs_TPM_pearson[i]<-cor(a$TPM,a$CG_TSS,use="complete.obs")
    AS_epigenetic_factors$CG_TSS_TPMfoldchange[i]<-median(a$TPM[a$CG_TSS>0.5]+1)/median(a$TPM[a$CG_TSS<0.5 ]+1)
  } else { if(length(a$TPM[!is.na(a$CG_TSS)])>10) {
    AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    AS_epigenetic_factors$CG_TSS_vs_TPM_pearson[i]<-cor(a$TPM,a$CG_TSS,use="complete.obs")
  }   else {
    AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p[i]<-NA
    AS_epigenetic_factors$CG_TSS_vs_TPM_pearson[i]<-NA
  }
  }
}



############################################
# analyse the numbers 
###############################################


# lincRNAs that are defined by methylation

# lincRNAs defined by CHH at their GB 
lincs_def_by_CHH_GB<-LINC_epigenetic_factors$gene[(abs(LINC_epigenetic_factors$CHH_GB_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CHH_GB_vs_TPM_pearson)) | (LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p)) ]
length(lincs_def_by_CHH_GB)
#152

# lincRNAs defined by CHH at their promoter 
lincs_def_by_CHH_TSS<-LINC_epigenetic_factors$gene[(abs(LINC_epigenetic_factors$CHH_TSS_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CHH_TSS_vs_TPM_pearson))| (LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p))]
length(lincs_def_by_CHH_TSS)
#168 

# lincRNAs defined by CG at their GB 
lincs_def_by_CG_GB<-LINC_epigenetic_factors$gene[(abs(LINC_epigenetic_factors$CG_GB_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CG_GB_vs_TPM_pearson)) | (LINC_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p))]
length( lincs_def_by_CG_GB)
#326

# lincRNAs defined by CG at their promoter 
lincs_def_by_CG_TSS<-LINC_epigenetic_factors$gene[(abs(LINC_epigenetic_factors$CG_TSS_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CG_TSS_vs_TPM_pearson))| (LINC_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p))]
  length(lincs_def_by_CG_TSS)
#318 

length(intersect(lincs_def_by_CG_GB,lincs_def_by_CG_TSS))
#257 
length(intersect(lincs_def_by_CHH_GB,lincs_def_by_CHH_TSS))
#95 
length(intersect(lincs_def_by_CG_TSS,lincs_def_by_CHH_TSS))
#98
length(intersect(lincs_def_by_CG_GB,lincs_def_by_CHH_GB))
#104

length(union(lincs_def_by_CG_GB,lincs_def_by_CHH_GB))
#374
length(union(lincs_def_by_CHH_GB,lincs_def_by_CHH_TSS))
#225

length(union(lincs_def_by_CG_GB,lincs_def_by_CG_TSS))
#387
length(union(lincs_def_by_CG_TSS,lincs_def_by_CHH_TSS))
#388
length(intersect(intersect(lincs_def_by_CG_GB,lincs_def_by_CHH_GB),intersect(lincs_def_by_CG_TSS,lincs_def_by_CHH_TSS)))
# 54



library(VennDiagram)
venn.diagram(
  x = list(union(lincs_def_by_CG_GB,lincs_def_by_CG_TSS), union(lincs_def_by_CHH_GB,lincs_def_by_CHH_TSS)),
  category.names = c("CG" , "CHH" ),
  filename = 'Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/linc_defined_by_CG_vs_CHH_venn_diagramm.png',
  output=TRUE
)

venn.diagram(
  x = list(union(lincs_def_by_CG_GB,lincs_def_by_CHH_GB), union(lincs_def_by_CG_TSS,lincs_def_by_CHH_TSS)),
  category.names = c("by meth\n at GB" , "by meth\nat promoter" ),
  filename = 'Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/linc_defined_by_GB_vs_TSS_venn_diagramm.png',
  output=TRUE
)



#LINC RNAs defined by methylation CG/CHH GB/promoter 

LINC_defined_by_meth<-LINC_epigenetic_factors$gene[(abs(LINC_epigenetic_factors$CHH_GB_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CHH_GB_vs_TPM_pearson)) |(abs(LINC_epigenetic_factors$CHH_TSS_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CHH_TSS_vs_TPM_pearson))| (LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p)) | (LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p))|
(abs(LINC_epigenetic_factors$CG_GB_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CG_GB_vs_TPM_pearson)) |(abs(LINC_epigenetic_factors$CG_TSS_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CG_TSS_vs_TPM_pearson))| (LINC_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p)) | (LINC_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p))]
length(LINC_defined_by_meth)
#454

length(union(union(lincs_def_by_CG_GB,lincs_def_by_CHH_GB),union(lincs_def_by_CHH_TSS,lincs_def_by_CG_TSS)))
#454



 ##############################
 # how many lincRNAs were informative? 
 ############################## 
  CG.linc$Nacc_highCG<-apply(CG.linc[,7:450], 1, function(i) sum(i >= 0.5))
  CG.linc$Nacc_lowCG<-apply(CG.linc[,7:450], 1, function(i) sum(i < 0.5))
  
  CG.linc_TSS$Nacc_highCG<-apply(CG.linc_TSS[,7:450], 1, function(i) sum(i >= 0.5))
  CG.linc_TSS$Nacc_lowCG<-apply(CG.linc_TSS[,7:450], 1, function(i) sum(i < 0.5))
  
  CHH.linc$Nacc_highCHH<-apply(CHH.linc[,7:450], 1, function(i) sum(i >= 0.01))
  CHH.linc$Nacc_lowCHH<-apply(CHH.linc[,7:450], 1, function(i) sum(i < 0.01))
  
 CHH.linc_TSS$Nacc_highCHH<-apply(CHH.linc_TSS[,7:450], 1, function(i) sum(i >= 0.01))
  CHH.linc_TSS$Nacc_lowCHH<-apply(CHH.linc_TSS[,7:450], 1, function(i) sum(i < 0.01))
  
  
  enough_var_CG_GB<-CG.linc$transcript[CG.linc$Nacc_highCG>10 & CG.linc$Nacc_lowCG>10 ]
  length(enough_var_CG_GB)
  #960 lincRNAs with enough CG diversity  - CG gene body 
  
  informative_CG_GB<-CG.linc$transcript[CG.linc$Nacc_highCG>10 & CG.linc$Nacc_lowCG>10 & CG.linc$transcript %in% linc_1001$gene[apply(linc_1001[,names(linc_1001) %in% names(CG.linc[,7:450])],1,max)>0.5]]
  length(informative_CG_GB)
  #524 lincRNAs with enough CG diversity and expression maximum - CG gene body 
  
  
  enough_var_CG_TSS<-CG.linc_TSS$transcript[CG.linc_TSS$Nacc_highCG>10 & CG.linc_TSS$Nacc_lowCG>10 ]
  length(enough_var_CG_TSS)
  #1027 lincRNAs with enough CG diversity and expression maximum - CG promoter
  
  informative_CG_TSS<-CG.linc_TSS$transcript[CG.linc_TSS$Nacc_highCG>10 & CG.linc_TSS$Nacc_lowCG>10 & CG.linc_TSS$transcript %in% linc_1001$gene[apply(linc_1001[,names(linc_1001) %in% names(CG.linc_TSS[,7:450])],1,max)>0.5]]
  length(informative_CG_TSS)
  #564 lincRNAs with enough CG diversity and expression maximum - CG promoter
 
  enough_var_CHH_GB<-CHH.linc$transcript[CHH.linc$Nacc_highCHH>10 & CHH.linc$Nacc_lowCHH>10 ]
   length(enough_var_CHH_GB)
  #1367 lincRNAs with enough CHH diversity and expression maximum - CHH gene body 
  
    informative_CHH_GB<-CHH.linc$transcript[CHH.linc$Nacc_highCHH>10 & CHH.linc$Nacc_lowCHH>10 & CHH.linc$transcript %in% linc_1001$gene[apply(linc_1001[,names(linc_1001) %in% names(CHH.linc[,7:450])],1,max)>0.5]]
  length(informative_CHH_GB)
  #682 lincRNAs with enough CHH diversity and expression maximum - CHH gene body 
  
  enough_var_CHH_TSS<-CHH.linc_TSS$transcript[CHH.linc_TSS$Nacc_highCHH>10 & CHH.linc_TSS$Nacc_lowCHH>10 ]
  length(enough_var_CHH_TSS)
  #1586 lincRNAs with enough CHH diversity and expression maximum - CHH promoter
  
  informative_CHH_TSS<-CHH.linc_TSS$transcript[CHH.linc_TSS$Nacc_highCHH>10 & CHH.linc_TSS$Nacc_lowCHH>10 & CHH.linc_TSS$transcript %in% linc_1001$gene[apply(linc_1001[,names(linc_1001) %in% names(CHH.linc_TSS[,7:450])],1,max)>0.5]]
  length(informative_CHH_TSS)
  #796 lincRNAs with enough CHH diversity and expression maximum - CHH promoter
  
 ##################################################################################################################### 
  length(union(union(enough_var_CG_GB,enough_var_CG_TSS),union(enough_var_CHH_GB,enough_var_CHH_TSS)))
  #1815
  length(union(union(informative_CG_GB,informative_CG_TSS),union(informative_CHH_GB,informative_CHH_TSS)))
  #895
  
 ################################## 









length(LINC_epigenetic_factors$gene[(abs(LINC_epigenetic_factors$CHH_GB_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CHH_GB_vs_TPM_pearson)) |(abs(LINC_epigenetic_factors$CHH_TSS_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CHH_TSS_vs_TPM_pearson))])
#17

length(LINC_epigenetic_factors$gene[(abs(LINC_epigenetic_factors$CHH_GB_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CHH_GB_vs_TPM_pearson)) |(abs(LINC_epigenetic_factors$CHH_TSS_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CHH_TSS_vs_TPM_pearson))| (LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p)) | (LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p))])
#225
#positive change vs negative change
length(LINC_epigenetic_factors$gene[ (LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p)) &LINC_epigenetic_factors$CHH_GB_TPMfoldchange>1 ])
#15 
length(LINC_epigenetic_factors$gene[ (LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p)) &LINC_epigenetic_factors$CHH_GB_TPMfoldchange<1 ])
#44

length(LINC_epigenetic_factors$gene[ (LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p)) &LINC_epigenetic_factors$CHH_TSS_TPMfoldchange>1 ])
#31 
length(LINC_epigenetic_factors$gene[ (LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p)) &LINC_epigenetic_factors$CHH_TSS_TPMfoldchange<1 ])
#43



length( (LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p)) | (LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p))])






lincs_defined_by_CHH<-LINC_epigenetic_factors$gene[(abs(LINC_epigenetic_factors$CHH_GB_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CHH_GB_vs_TPM_pearson)) |(abs(LINC_epigenetic_factors$CHH_TSS_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CHH_TSS_vs_TPM_pearson))| (LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p)) | (LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p))]

length(LINC_epigenetic_factors$gene[(abs(LINC_epigenetic_factors$CG_GB_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CG_GB_vs_TPM_pearson)) |(abs(LINC_epigenetic_factors$CG_TSS_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CG_TSS_vs_TPM_pearson))| (LINC_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p)) | (LINC_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p))])
#387

linc_defined_by_CG<-LINC_epigenetic_factors$gene[(abs(LINC_epigenetic_factors$CG_GB_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CG_GB_vs_TPM_pearson)) |(abs(LINC_epigenetic_factors$CG_TSS_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CG_TSS_vs_TPM_pearson))| (LINC_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p)) | (LINC_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(LINC_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p))]

length(LINC_epigenetic_factors$gene[(abs(LINC_epigenetic_factors$CG_GB_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CG_GB_vs_TPM_pearson)) |(abs(LINC_epigenetic_factors$CG_TSS_vs_TPM_pearson)>0.5&!is.na(LINC_epigenetic_factors$CG_TSS_vs_TPM_pearson))])
#97 

length(intersect(linc_defined_by_CG,lincs_defined_by_CHH))
#158






##################################
### define by histone marks 
#####################################

# lincRNAs that are defined by histone marks  
#H1
length(LINC_epigenetic_factors$gene[(abs(LINC_epigenetic_factors$cor_H1_TPM)>0.5&!is.na(LINC_epigenetic_factors$cor_H1_TPM) & LINC_epigenetic_factors$max_TPM_1001Gnew>0.5)])
#92
length(LINC_epigenetic_factors$gene[(abs(LINC_epigenetic_factors$cor_K9_TPM)>0.5&!is.na(LINC_epigenetic_factors$cor_K9_TPM) & LINC_epigenetic_factors$max_TPM_1001Gnew>0.5)])
#103
length(LINC_epigenetic_factors$gene[(abs(LINC_epigenetic_factors$cor_K27_TPM)>0.5&!is.na(LINC_epigenetic_factors$cor_K27_TPM) & LINC_epigenetic_factors$max_TPM_1001Gnew>0.5)])
#64



# AS RNAs that are defined by histone marks  
#H1
length(AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$cor_H1_TPM)>0.5&!is.na(AS_epigenetic_factors$cor_H1_TPM) & AS_epigenetic_factors$max_TPM_1001Gnew>0.5)])
#195
AS_by_H1<-AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$cor_H1_TPM)>0.5&!is.na(AS_epigenetic_factors$cor_H1_TPM) & AS_epigenetic_factors$max_TPM_1001Gnew>0.5)]

length(AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$cor_K9_TPM)>0.5&!is.na(AS_epigenetic_factors$cor_K9_TPM) & AS_epigenetic_factors$max_TPM_1001Gnew>0.5)])
#188
AS_by_K9<-AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$cor_K9_TPM)>0.5&!is.na(AS_epigenetic_factors$cor_K9_TPM) & AS_epigenetic_factors$max_TPM_1001Gnew>0.5)]


length(AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$cor_K27_TPM)>0.5&!is.na(AS_epigenetic_factors$cor_K27_TPM) & AS_epigenetic_factors$max_TPM_1001Gnew>0.5)])
#230
AS_by_K27<-AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$cor_K27_TPM)>0.5&!is.na(AS_epigenetic_factors$cor_K27_TPM) & AS_epigenetic_factors$max_TPM_1001Gnew>0.5)]


length(intersect(AS_defined_by_CG,AS_by_H1))
#47



###################################################
# AS lncRNAs that are defined by methylation 
###################################################

##############################
# how many AS lncRNAs were informative? 
############################## 
CG.as$Nacc_highCG<-apply(CG.as[,7:450], 1, function(i) sum(i >= 0.5))
CG.as$Nacc_lowCG<-apply(CG.as[,7:450], 1, function(i) sum(i < 0.5))

CG.as_TSS$Nacc_highCG<-apply(CG.as_TSS[,7:450], 1, function(i) sum(i >= 0.5))
CG.as_TSS$Nacc_lowCG<-apply(CG.as_TSS[,7:450], 1, function(i) sum(i < 0.5))

CHH.as$Nacc_highCHH<-apply(CHH.as[,7:450], 1, function(i) sum(i >= 0.01))
CHH.as$Nacc_lowCHH<-apply(CHH.as[,7:450], 1, function(i) sum(i < 0.01))

CHH.as_TSS$Nacc_highCHH<-apply(CHH.as_TSS[,7:450], 1, function(i) sum(i >= 0.01))
CHH.as_TSS$Nacc_lowCHH<-apply(CHH.as_TSS[,7:450], 1, function(i) sum(i < 0.01))


enough_var_CG_GB<-CG.as$transcript[CG.as$Nacc_highCG>10 & CG.as$Nacc_lowCG>10 ]
length(enough_var_CG_GB)
#1976 asRNAs with enough CG diversity  - CG gene body 

informative_CG_GB<-CG.as$transcript[CG.as$Nacc_highCG>10 & CG.as$Nacc_lowCG>10 & CG.as$transcript %in% as_1001$gene[apply(as_1001[,names(as_1001) %in% names(CG.as[,7:450])],1,max)>0.5]]
length(informative_CG_GB)
#858 asRNAs with enough CG diversity and expression maximum - CG gene body 


enough_var_CG_TSS<-CG.as_TSS$transcript[CG.as_TSS$Nacc_highCG>10 & CG.as_TSS$Nacc_lowCG>10 ]
length(enough_var_CG_TSS)
#2081 asRNAs with enough CG diversity and expression maximum - CG promoter

informative_CG_TSS<-CG.as_TSS$transcript[CG.as_TSS$Nacc_highCG>10 & CG.as_TSS$Nacc_lowCG>10 & CG.as_TSS$transcript %in% as_1001$gene[apply(as_1001[,names(as_1001) %in% names(CG.as_TSS[,7:450])],1,max)>0.5]]
length(informative_CG_TSS)
#887 asRNAs with enough CG diversity and expression maximum - CG promoter

enough_var_CHH_GB<-CHH.as$transcript[CHH.as$Nacc_highCHH>10 & CHH.as$Nacc_lowCHH>10 ]
length(enough_var_CHH_GB)
#1645 asRNAs with enough CHH diversity and expression maximum - CHH gene body 

informative_CHH_GB<-CHH.as$transcript[CHH.as$Nacc_highCHH>10 & CHH.as$Nacc_lowCHH>10 & CHH.as$transcript %in% as_1001$gene[apply(as_1001[,names(as_1001) %in% names(CHH.as[,7:450])],1,max)>0.5]]
length(informative_CHH_GB)
#820 asRNAs with enough CHH diversity and expression maximum - CHH gene body 

enough_var_CHH_TSS<-CHH.as_TSS$transcript[CHH.as_TSS$Nacc_highCHH>10 & CHH.as_TSS$Nacc_lowCHH>10 ]
length(enough_var_CHH_TSS)
#3233 asRNAs with enough CHH diversity and expression maximum - CHH promoter

informative_CHH_TSS<-CHH.as_TSS$transcript[CHH.as_TSS$Nacc_highCHH>10 & CHH.as_TSS$Nacc_lowCHH>10 & CHH.as_TSS$transcript %in% as_1001$gene[apply(as_1001[,names(as_1001) %in% names(CHH.as_TSS[,7:450])],1,max)>0.5]]
length(informative_CHH_TSS)
#1718 asRNAs with enough CHH diversity and expression maximum - CHH promoter


length(union(union(enough_var_CG_GB,enough_var_CG_TSS),union(enough_var_CHH_GB,enough_var_CHH_TSS)))
#4799
length(union(union(informative_CG_GB,informative_CG_TSS),union(informative_CHH_GB,informative_CHH_TSS)))
#2371






# AS RNAs that are defined by methylation

# AS RNAs defined by CHH at their GB 
as_def_by_CHH_GB<-AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$CHH_GB_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CHH_GB_vs_TPM_pearson)) | (AS_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p)) ]
length(as_def_by_CHH_GB)
#232

# ASRNAs defined by CHH at their promoter 
as_def_by_CHH_TSS<-AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson))| (AS_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p))]
length(as_def_by_CHH_TSS)
#292 

# ASRNAs defined by CG at their GB 
as_def_by_CG_GB<-AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$CG_GB_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CG_GB_vs_TPM_pearson)) | (AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p))]
length( as_def_by_CG_GB)
#274

# ASRNAs defined by CG at their promoter 
as_def_by_CG_TSS<-AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$CG_TSS_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CG_TSS_vs_TPM_pearson))| (AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p))]
length(as_def_by_CG_TSS)
#274 

length(intersect(as_def_by_CG_GB,as_def_by_CG_TSS))
#171 
length(intersect(as_def_by_CHH_GB,as_def_by_CHH_TSS))
#160 
length(intersect(as_def_by_CG_TSS,as_def_by_CHH_TSS))
#178
length(intersect(as_def_by_CG_GB,as_def_by_CHH_GB))
#155
length(intersect(intersect(as_def_by_CG_GB,as_def_by_CHH_GB),intersect(as_def_by_CG_TSS,as_def_by_CHH_TSS)))
#112

length(intersect(union(as_def_by_CG_GB,as_def_by_CHH_GB),union(as_def_by_CHH_TSS,as_def_by_CG_TSS)))
# 230
length(union(as_def_by_CG_GB,as_def_by_CHH_GB))
#351
length(union(as_def_by_CHH_TSS,as_def_by_CG_TSS))
#388

length(union(as_def_by_CG_GB,as_def_by_CG_TSS))
#377
length(union(as_def_by_CHH_TSS,as_def_by_CHH_GB))
#364
length(intersect(union(as_def_by_CG_GB,as_def_by_CG_TSS),union(as_def_by_CHH_TSS,as_def_by_CHH_GB)))
# 232


length(union(union(as_def_by_CG_GB,as_def_by_CHH_GB),union(as_def_by_CHH_TSS,as_def_by_CG_TSS)))
#509
length(union(union(as_def_by_CG_GB,as_def_by_CG_TSS),union(as_def_by_CHH_TSS,as_def_by_CHH_GB)))
#509


library(VennDiagram)
venn.diagram(
  x = list(union(as_def_by_CG_GB,as_def_by_CG_TSS), union(as_def_by_CHH_TSS,as_def_by_CHH_GB)),
  category.names = c("by CG" , "by CHH"),
  filename = 'Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/AS_defined_by_CG_vs_CHH_venn_diagramm.png',
  output=TRUE
)


venn.diagram(
  x = list(union(as_def_by_CG_GB,as_def_by_CHH_GB), union(as_def_by_CHH_TSS,as_def_by_CG_TSS)),
  category.names = c("by meth\n at GB" , "by meth\nat promoter"),
  filename = 'Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/AS_defined_by_GB_vs_TSS_venn_diagramm.png',
  output=TRUE
)

length(union(union(as_def_by_CG_GB,as_def_by_CHH_GB),union(as_def_by_CHH_TSS,as_def_by_CG_TSS)))
#509



# how many of the methylation defined RNAs have a TE in TAIR10
as_def_by_CG_GB
as_def_by_CHH_GB
as_def_by_CHH_TSS
as_def_by_CG_TSS

boxplot(as_TE_cov_all_loci_2cols$coverage,
        as_TE_cov_all_loci_2cols$coverage[as_TE_cov_all_loci_2cols$gene %in% as_def_by_CG_GB],
        as_TE_cov_all_loci_2cols$coverage[as_TE_cov_all_loci_2cols$gene %in% as_def_by_CG_TSS],
        as_TE_cov_all_loci_2cols$coverage[as_TE_cov_all_loci_2cols$gene %in% as_def_by_CHH_GB],
        as_TE_cov_all_loci_2cols$coverage[as_TE_cov_all_loci_2cols$gene %in% as_def_by_CHH_TSS],names=c("all","CG GB","CG TSS","CHH GB","CHH TSS"),las=2,ylab="TE content in TAIR10",main="AS lncRNAs explained by methylation",notch = T,outline = F)

boxplot(linc_TE_cov_all_loci_2cols$coverage,
        linc_TE_cov_all_loci_2cols$coverage[linc_TE_cov_all_loci_2cols$gene %in% lincs_def_by_CG_GB],
        linc_TE_cov_all_loci_2cols$coverage[linc_TE_cov_all_loci_2cols$gene %in% lincs_def_by_CG_TSS],
        linc_TE_cov_all_loci_2cols$coverage[linc_TE_cov_all_loci_2cols$gene %in% lincs_def_by_CHH_GB],
        linc_TE_cov_all_loci_2cols$coverage[linc_TE_cov_all_loci_2cols$gene %in% lincs_def_by_CHH_TSS],names=c("all","CG GB","CG TSS","CHH GB","CHH TSS"),las=2,ylab="TE content in TAIR10",main="luncRNAs explained by methylation",notch = T,outline = F)







#AS RNAs defined by methylation CG/CHH GB/promoter 

AS_defined_by_meth<-AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$CHH_GB_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CHH_GB_vs_TPM_pearson)) |(abs(AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson))| (AS_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p)) | (AS_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p))|
(abs(AS_epigenetic_factors$CG_GB_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CG_GB_vs_TPM_pearson)) |(abs(AS_epigenetic_factors$CG_TSS_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CG_TSS_vs_TPM_pearson))| (AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p)) | (AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p))]
length(AS_defined_by_meth)
#509




length(AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$CG_GB_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CG_GB_vs_TPM_pearson)) |(abs(AS_epigenetic_factors$CG_TSS_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CG_TSS_vs_TPM_pearson))| (AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p)) | (AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p))])
#377

AS_defined_by_CG<-AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$CG_GB_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CG_GB_vs_TPM_pearson)) |(abs(AS_epigenetic_factors$CG_TSS_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CG_TSS_vs_TPM_pearson))| (AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p)) | (AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p))]








length(AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$CHH_GB_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CHH_GB_vs_TPM_pearson)) |(abs(AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson))])
#31

length(AS_epigenetic_factors$gene[
  (abs(AS_epigenetic_factors$CHH_GB_vs_TPM_pearson)>0.5
   &!is.na(AS_epigenetic_factors$CHH_GB_vs_TPM_pearson)) |
  (abs(AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson)>0.5&
     !is.na(AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson))|
   (AS_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p>2 
 & !is.na(AS_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p)) | (AS_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p>2 
& !is.na(AS_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p))])
#364

#defined by CHH at TSS 
length(AS_epigenetic_factors$gene[
  (abs(AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson)>0.5&
!is.na(AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson))|
 (AS_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p>2 
& !is.na(AS_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p))])
#292

#defined by CG at TSS 
length(AS_epigenetic_factors$gene[
  (abs(AS_epigenetic_factors$CHH_GB_vs_TPM_pearson)>0.5
   &!is.na(AS_epigenetic_factors$CHH_GB_vs_TPM_pearson)) |
    (AS_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p>2 
     & !is.na(AS_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p))])
#232

AS_defined_by_CHH<-AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$CHH_GB_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CHH_GB_vs_TPM_pearson)) |(abs(AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CHH_TSS_vs_TPM_pearson))| (AS_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CHH_GB_vs_TPM_ManWhittest_minuslog10p)) | (AS_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CHH_TSS_vs_TPM_ManWhittest_minuslog10p))]
length(AS_defined_by_CHH)

length(AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$CG_GB_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CG_GB_vs_TPM_pearson)) |(abs(AS_epigenetic_factors$CG_TSS_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CG_TSS_vs_TPM_pearson))| (AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p)) | (AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p))])
#377

AS_defined_by_CG<-AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$CG_GB_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CG_GB_vs_TPM_pearson)) |(abs(AS_epigenetic_factors$CG_TSS_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CG_TSS_vs_TPM_pearson))| (AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CG_GB_vs_TPM_ManWhittest_minuslog10p)) | (AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p>2 & !is.na(AS_epigenetic_factors$CG_TSS_vs_TPM_ManWhittest_minuslog10p))]

length(AS_epigenetic_factors$gene[(abs(AS_epigenetic_factors$CG_GB_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CG_GB_vs_TPM_pearson)) |(abs(AS_epigenetic_factors$CG_TSS_vs_TPM_pearson)>0.5&!is.na(AS_epigenetic_factors$CG_TSS_vs_TPM_pearson))])
#73 

length(intersect(AS_defined_by_CG,AS_defined_by_CHH))
#157


#scatter plot for CUFF_NC_1515 - anticorrelated with CHH/CG 
linc="CUFF_NC.1515"
a<-as.data.frame(t(denovo2021.TPMs.genes.1001G[linc,2:462]))
a$accession<-rownames(a)
b<-as.data.frame(t(CG.linc[linc,7:450]))
b$accession<-rownames(b)
a<-merge(a,b,by="accession")
b<-as.data.frame(t(CG.linc_TSS[linc,7:450]))
b$accession<-rownames(b)
a<-merge(a,b,by="accession")
b<-as.data.frame(t(CG.linc_TES[linc,7:450]))
b$accession<-rownames(b)
a<-merge(a,b,by="accession")
names(a)<-c("accession","TPM","CG","CG_TSS","CG_TES")


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/CUFF_NC_1515.CG.1001G.scatterplot.pdf",height = 4,width = 4)
par(mar=c(6,4,2,2)) 
plot(a$CG,a$TPM,pch=20,main=linc,ylab="TPM",xlab="CG methylation level",xlim=c(0,1),las=2,col=alpha("darkgreen",alpha=0.4))
wil<-wilcox.test(a$TPM[a$CG>0.5 & !is.na(a$CG)],a$TPM[a$CG<0.5 & !is.na(a$CG)])
text(0.7,18,labels=paste("M-W test:\nTPM(CG<0.5) vs TPM(CG>0.5)"),cex=0.8)

text(0.8,15,labels=paste("p.value=10^",round(log10(wil$p.value),1),sep = ""))
#text(0.8,15,labels=paste("R=",round(cor(a$CG,a$TPM,method = "pearson"),2)))
dev.off()

linc="CUFF_NC.1515"
a<-as.data.frame(t(denovo2021.TPMs.genes.1001G[linc,2:462]))
a$accession<-rownames(a)
b<-as.data.frame(t(CHH.linc[linc,7:450]))
b$accession<-rownames(b)
a<-merge(a,b,by="accession")
b<-as.data.frame(t(CHH.linc_TSS[linc,7:450]))
b$accession<-rownames(b)
a<-merge(a,b,by="accession")
b<-as.data.frame(t(CHH.linc_TES[linc,7:450]))
b$accession<-rownames(b)
a<-merge(a,b,by="accession")
names(a)<-c("accession","TPM","CHH","CHH_TSS","CHH_TES")

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/CUFF_NC_1515.CHH.1001G.scatterplot.pdf",height = 4,width = 4)
par(mar=c(6,4,2,2)) 
plot(a$CHH,a$TPM,pch=20,main=linc,ylab="TPM",xlab="CHH methylation level",xlim=c(0,1),las=2,col=alpha("sienna",alpha=0.4))
#text(aa[aa$mark=="K9",3],aa[aa$mark=="K9",2]+0.3, labels=gsub(pattern = "K9.",replacement = "",aa$chip_sample[aa$mark=="K9"]),cex=0.6,srt=30)
wil<-wilcox.test(a$TPM[a$CHH>0 & !is.na(a$CHH)],a$TPM[a$CHH==0 & !is.na(a$CHH)])
text(0.7,18,labels=paste("M-W test:\nTPM(CHH=0) vs TPM(CHH>0)"),cex=0.8)

text(0.8,15,labels=paste("p.value=10^",round(log10(wil$p.value),1),sep = ""))
dev.off()


#scatter plot for CUFF_NC_1515 - anticorrelated with histone marks
linc="CUFF_NC.1515"
a<-as.data.frame(t(denovo2021.TPMs.genes.1001Gnew[linc,grep("mean.",names(denovo2021.TPMs.genes.1001Gnew))]))
a$accession<-unlist(lapply(strsplit(as.character(rownames(a)),".", fixed = T), "[",2))
a<-a[1:28,]
b<-as.data.frame(t(chip.denovo.quantstan[linc,2:68]))
b$accession<-unlist(lapply(strsplit(as.character(rownames(b)),".", fixed = T), "[",2))
b$chip_sample<-rownames(b)
b$mark<-unlist(lapply(strsplit(as.character(rownames(b)),".", fixed = T), "[",1))
aa<-merge(a,b,by="accession")


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/CUFF_NC_1515.K9.scatterplot.pdf",height = 3,width = 3)
par(mar=c(6,4,2,2)) 
plot(aa[aa$mark=="K9",3],aa[aa$mark=="K9",2],pch=19,main=linc,ylab="TPM",xlab="H3K9me2 normalized signal",ylim=c(0,4.5),xlim=c(-0.7,2.5),las=2,col=alpha("darkmagenta",alpha=0.7))
text(aa[aa$mark=="K9",3],aa[aa$mark=="K9",2]+0.3, labels=gsub(pattern = "K9.",replacement = "",aa$chip_sample[aa$mark=="K9"]),cex=0.6,srt=30)
text(1.5,4,labels=paste("R=",round(cor(aa[aa$mark=="K9",3],aa[aa$mark=="K9",2],method = "pearson"),2)))
dev.off()


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/CUFF_NC_1515.H1.scatterplot.pdf",height = 3,width = 3)
par(mar=c(6,4,2,2)) 
plot(aa[aa$mark=="H1",3],aa[aa$mark=="H1",2],pch=19,main=linc,ylab="TPM",xlab="H1 normalized signal",ylim=c(0,4.8),xlim=c(-1,1.6),las=2,col=alpha("aquamarine4",alpha=0.7))
text(aa[aa$mark=="H1",3],aa[aa$mark=="H1",2]+0.3, labels=gsub(pattern = "H1.",replacement = "",aa$chip_sample[aa$mark=="H1"]),cex=0.6,srt=30)
text(1,4,labels=paste("R=",round(cor(aa[aa$mark=="H1",3],aa[aa$mark=="H1",2],method = "pearson"),2)))
dev.off()





#small RNAs 

a<-as.data.frame(t(denovo2021.TPMs.genes.ERACAPS[linc,2:97]))
a$accession<-unlist(lapply(strsplit(as.character(rownames(a)),".", fixed = T), "[",2))
a$tissue<-lapply(strsplit(as.character(row.names(a)),".", fixed = T), "[",1)
b<-as.data.frame(t(sRNA.24nt.denovo2021.RPM[linc,7:20]))
b$accession<-unlist(lapply(strsplit(as.character(rownames(b)),"X", fixed = T), "[",2))
aa<-merge(a,b,by="accession")

names(aa)<-c("accession","TPM","tissue","sRNA_RPM")


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/CUFF_NC_1515.24nt.sRNA.scatterplot.pdf",height = 3,width = 3)
par(mar=c(6,4,2,2)) ylim=c(0,4.8),xlim=c(-1,1.6),

plot(aa$sRNA_RPM[aa$tissue=="S"],aa$TPM[aa$tissue=="S"],pch=19,main=linc,ylab="TPM",xlab="24nt coverage, RPM",las=2,col=alpha("brown",alpha=0.6),ylim=c(0,5))
#plot(aa$sRNA_RPM[aa$tissue=="R"],aa$TPM[aa$tissue=="R"],pch=19,main=linc,ylab="TPM",xlab="24nt coverage, RPM",las=2,col=alpha("brown",alpha=0.6))
#plot(aa$sRNA_RPM[aa$tissue=="F"],aa$TPM[aa$tissue=="F"],pch=19,main=linc,ylab="TPM",xlab="24nt coverage, RPM",las=2,col=alpha("brown",alpha=0.6))

text(aa$sRNA_RPM[aa$tissue=="S"]+0.05,aa$TPM[aa$tissue=="S"]+0.05, labels=gsub(pattern = "H1.",replacement = "",aa$accession[aa$tissue=="S"]),cex=0.6,srt=30)
#text(1,4,labels=paste("R=",round(cor(aa[aa$mark=="H1",3],aa[aa$mark=="H1",2],method = "pearson"),2)))
dev.off()





linc="CUFF_NC.1960"

a<-as.data.frame(t(denovo2021.TPMs.genes.1001Gnew[linc,grep("mean.",names(denovo2021.TPMs.genes.1001Gnew))]))
a$accession<-unlist(lapply(strsplit(as.character(rownames(a)),".", fixed = T), "[",2))
a<-a[1:28,]
b<-as.data.frame(t(chip.denovo.quantstan[linc,2:68]))
b$accession<-unlist(lapply(strsplit(as.character(rownames(b)),".", fixed = T), "[",2))
b$chip_sample<-rownames(b)
b$mark<-unlist(lapply(strsplit(as.character(rownames(b)),".", fixed = T), "[",1))
aa<-merge(a,b,by="accession")

par(mar=c(6,4,2,2)) 
plot(aa[aa$mark=="K27",3],aa[aa$mark=="K27",2],pch=19,main=linc,ylab="TPM",xlab="H3K27me3 normalized signal",las=2,col=alpha("red",alpha=0.7))
text(aa[aa$mark=="K27",3],aa[aa$mark=="K27",2]+0.03, labels=gsub(pattern = "K27.",replacement = "",aa$chip_sample[aa$mark=="K27"]),cex=0.6,srt=30)
text(3,0.5,labels=paste("R=",round(cor(aa[aa$mark=="K27",3],aa[aa$mark=="K27",2],method = "pearson"),2)))







pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/CUFF_NC_1515.H1.scatterplot.pdf",height = 3,width = 3)
par(mar=c(6,4,2,2)) 
plot(aa[aa$mark=="H1",3],aa[aa$mark=="H1",2],pch=19,main=linc,ylab="TPM",xlab="H1 normalized signal",ylim=c(0,4.8),xlim=c(-1,1.6),las=2,col=alpha("aquamarine4",alpha=0.7))
text(aa[aa$mark=="H1",3],aa[aa$mark=="H1",2]+0.3, labels=gsub(pattern = "H1.",replacement = "",aa$chip_sample[aa$mark=="H1"]),cex=0.6,srt=30)
text(1,4,labels=paste("R=",round(cor(aa[aa$mark=="H1",3],aa[aa$mark=="H1",2],method = "pearson"),2)))
dev.off()

















linc="CUFF_NC.9367"
a<-as.data.frame(t(denovo2021.TPMs.genes.1001Gnew[linc,grep("mean.",names(denovo2021.TPMs.genes.1001Gnew))]))
a$accession<-unlist(lapply(strsplit(as.character(rownames(a)),".", fixed = T), "[",2))
a<-a[1:28,]
b<-as.data.frame(t(chip.denovo.quantstan[linc,2:68]))
b$accession<-unlist(lapply(strsplit(as.character(rownames(b)),".", fixed = T), "[",2))
b$chip_sample<-rownames(b)
b$mark<-unlist(lapply(strsplit(as.character(rownames(b)),".", fixed = T), "[",1))
aa<-merge(a,b,by="accession")


#pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/CUFF_NC.9145.K9.scatterplot.pdf",height = 3,width = 3)
par(mar=c(6,4,2,2)) 
plot(aa[aa$mark=="K9",3],aa[aa$mark=="K9",2],pch=19,main=linc,ylab="TPM",xlab="H3K9me2 normalized signal",las=2)
text(aa[aa$mark=="K9",3],aa[aa$mark=="K9",2]+0.3, labels=gsub(pattern = "K9.",replacement = "",aa$chip_sample[aa$mark=="K9"]),cex=0.6,srt=30)
text(1.5,4,labels=paste("R=",round(cor(aa[aa$mark=="K9",3],aa[aa$mark=="K9",2],method = "pearson"),2)))
#dev.off()

#pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/CUFF_NC.9145.K9.scatterplot.pdf",height = 3,width = 3)
par(mar=c(6,4,2,2)) 
plot(aa[aa$mark=="H1",3],aa[aa$mark=="H1",2],pch=19,main=linc,ylab="TPM",xlab="H1 normalized signal",las=2)
text(aa[aa$mark=="H1",3],aa[aa$mark=="H1",2]+0.03, labels=gsub(pattern = "H1",replacement = "",aa$chip_sample[aa$mark=="H1"]),cex=0.6,srt=30)
text(1.5,4,labels=paste("R=",round(cor(aa[aa$mark=="H1",3],aa[aa$mark=="H1",2],method = "pearson"),2)))
#dev.off()



linc="CUFF_NC.6674"
a<-as.data.frame(t(denovo2021.TPMs.genes.1001G[linc,2:462]))
a$accession<-rownames(a)
b<-as.data.frame(t(CHH.as[linc,7:450]))
b$accession<-rownames(b)
a<-merge(a,b,by="accession")
b<-as.data.frame(t(CHH.as_TSS[linc,7:450]))
b$accession<-rownames(b)
a<-merge(a,b,by="accession")
b<-as.data.frame(t(CHH.as_TES[linc,7:450]))
b$accession<-rownames(b)
a<-merge(a,b,by="accession")
names(a)<-c("accession","TPM","CHH","CHH_TSS","CHH_TES")

#pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/CUFF_NC.9145.CHH.1001G.scatterplot.pdf",height = 4,width = 4)
par(mar=c(6,4,2,2)) 
plot(a$CHH,a$TPM,pch=20,main=linc,ylab="TPM",xlab="CHH methylation level",xlim=c(0,1),las=2,col=alpha("sienna",alpha=0.4))
#text(aa[aa$mark=="K9",3],aa[aa$mark=="K9",2]+0.3, labels=gsub(pattern = "K9.",replacement = "",aa$chip_sample[aa$mark=="K9"]),cex=0.6,srt=30)
wil<-wilcox.test(a$TPM[a$CHH>=0.01 & !is.na(a$CHH)],a$TPM[a$CHH<0.01 & !is.na(a$CHH)])
text(0.7,18,labels=paste("M-W test:\nTPM(CHH=0) vs TPM(CHH>0)"),cex=0.8)

text(0.8,15,labels=paste("p.value=10^",round(log10(wil$p.value),1),sep = ""))
#dev.off()
plot(a$CHH_TSS,a$TPM,pch=20,main=linc,ylab="TPM",xlab="CHH_TSS methylation level",xlim=c(0,1),las=2,col=alpha("sienna",alpha=0.4))
wil<-wilcox.test(a$TPM[a$CHH_TSS>=0.01 & !is.na(a$CHH_TSS)],a$TPM[a$CHH_TSS<0.01 & !is.na(a$CHH_TSS)])
text(0.8,8,labels=paste("p.value=10^",round(log10(wil$p.value),1),sep = ""))




a<-as.data.frame(t(denovo2021.TPMs.genes.1001G[linc,2:462]))
a$accession<-rownames(a)
b<-as.data.frame(t(CG.as[linc,7:450]))
b$accession<-rownames(b)
a<-merge(a,b,by="accession")
b<-as.data.frame(t(CG.as_TSS[linc,7:450]))
b$accession<-rownames(b)
a<-merge(a,b,by="accession")
b<-as.data.frame(t(CG.as_TES[linc,7:450]))
b$accession<-rownames(b)
a<-merge(a,b,by="accession")
names(a)<-c("accession","TPM","CG","CG_TSS","CG_TES")


#pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig3-4/CUFF_NC.9145.CG.1001G.scatterplot.pdf",height = 4,width = 4)
par(mar=c(6,4,2,2)) 
plot(a$CG,a$TPM,pch=20,main=linc,ylab="TPM",xlab="CG methylation level",xlim=c(0,1),las=2,col=alpha("darkgreen",alpha=0.4))
wil<-wilcox.test(a$TPM[a$CG>0.5 & !is.na(a$CG)],a$TPM[a$CG<0.5 & !is.na(a$CG)])
text(0.7,18,labels=paste("M-W test:\nTPM(CG<0.5) vs TPM(CG>0.5)"),cex=0.8)

text(0.8,15,labels=paste("p.value=10^",round(log10(wil$p.value),1),sep = ""))
#text(0.8,15,labels=paste("R=",round(cor(a$CG,a$TPM,method = "pearson"),2)))
#dev.off()


plot(a$CG_TSS,a$TPM,pch=20,main=linc,ylab="TPM",xlab="CG methylation level",xlim=c(0,1),las=2,col=alpha("darkgreen",alpha=0.4))
wil<-wilcox.test(a$TPM[a$CG_TSS>0.5 & !is.na(a$CG)],a$TPM[a$CG_TSS<0.5 & !is.na(a$CG)])
