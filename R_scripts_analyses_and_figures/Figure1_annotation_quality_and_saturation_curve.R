
######################################
# check annotation quality
######################################

#denovo annotation basic stats

#number of exons 
boxplot(denovoPC.transcripts$exon_N, denovo_pseudogene.transcripts$exon_N,
        lncRNAs.antisense.transcripts$exon_N,lncRNAs.intergenic.transcripts$exon_N,
        lncRNAs.AS_to_TE.transcripts$exon_N, TE_genes.transcripts$exon_N,
        TE_frags.transcripts$exon_N,
        main="Exon number for different types of genes\nin the cumulative transcriptome",
        names=c("PC","pseudogenes","AS to PC","lincRNAs","AS to TE genes","TE genes","TE fragments"),las=2,cex.lab=1.3,log='y',
        ylab="number of exons, log scale",notch = T,col=c("#486EB4","#767172","#90C473","#F2AB54","#B294C5","#673A8E","#805FA5"))

#length 
boxplot((denovoPC.transcripts$end -denovoPC.transcripts$start)/1000,(denovo_pseudogene.transcripts$end -denovo_pseudogene.transcripts$start)/1000, 
        (lncRNAs.antisense.transcripts$end -lncRNAs.antisense.transcripts$start)/1000 ,
        (lncRNAs.intergenic.transcripts$end - lncRNAs.intergenic.transcripts$start)/1000,(lncRNAs.AS_to_TE.transcripts$end- lncRNAs.AS_to_TE.transcripts$start)/1000,
        (TE_genes.transcripts$end- TE_genes.transcripts$start)/1000,
        (TE_frags.transcripts$end-TE_frags.transcripts$start)/1000,
        main="locus length for different types of genes\nin the cumulative transcriptome",
        names=c("PC","pseudogenes","AS to PC","lincRNAs","AS to TE genes","TE genes","TE fragments"),
        las=2,cex.lab=1.3,log='y', ylab="locus length,kbp, log scale",notch = T, col=c("#486EB4","#767172","#90C473","#F2AB54","#B294C5","#673A8E","#805FA5"))



#exon coverage / length coverage of the public annotations 

Araport11NC_exon_coverage_by_denovoNC <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/Araport11NC_exon_coverage_by_denovoNC.bed", header=FALSE)
Araport11NC_strict_exon_coverage_by_denovoNC <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/Araport11NC_strict_exon_coverage_by_denovoNC.bed", header=FALSE)
Araport11PC_exon_coverage_by_denovoPC <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/Araport11PC_exon_coverage_by_denovoPC.bed", header=FALSE)
Araport11TEgenes_exon_coverage_by_denovoTEgenes <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/Araport11TEgenes_exon_coverage_by_denovoTEgenes.bed", header=FALSE)

a1<-aggregate(Araport11NC_strict_exon_coverage_by_denovoNC[, 8], list(Araport11NC_strict_exon_coverage_by_denovoNC$V4), sum)
a2<-aggregate(Araport11NC_strict_exon_coverage_by_denovoNC[, 9], list(Araport11NC_strict_exon_coverage_by_denovoNC$V4), sum)
b1<-aggregate(Araport11PC_exon_coverage_by_denovoPC[, 8], list(Araport11PC_exon_coverage_by_denovoPC$V4), mean)
b2<-aggregate(Araport11PC_exon_coverage_by_denovoPC[, 9], list(Araport11PC_exon_coverage_by_denovoPC$V4), mean)
c1<-aggregate(Araport11TEgenes_exon_coverage_by_denovoTEgenes[,8], list(Araport11TEgenes_exon_coverage_by_denovoTEgenes$V4), mean)
c2<-aggregate(Araport11TEgenes_exon_coverage_by_denovoTEgenes[, 9], list(Araport11TEgenes_exon_coverage_by_denovoTEgenes$V4), mean)

c1$coverage<-c1$x/c2$x


col=c("#486EB4","#90C473","#F2AB54","#673A8E")

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig1/Fig1_SUPPL_boxplot_exon_coverage_annotated_by_denovo.pdf",height = 3,width = 3.5)
###########################################################################
par(mar=c(6,5,4,2))
boxplot(b1$x/b2$x,a1$x/a2$x,c1$x/c2$x,
        names=c("PC genes","lncRNAs","TE genes"),outline = F,notch = T,ylab="exon coverage",main="exon coverage of annotated genes by our cumulative annotation",col=c("#486EB4","#90C473","#673A8E"),las=2)
dev.off()



##################################################################
#saturation curve analysis                         ###############
##################################################################

#Z:\01_POSTDOC\03_Projects\2018_lncRNA_variation_paper\01_lncRNA_identification\02_Identification_Saturation_Curve\2021\20211013_cuffmerge_filtering_and_classification_automated_pipeline_for_satur_curve_rep1.bash



satur_curves_gene_numbers_1001G <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/02_Identification_Saturation_Curve/2021/gene_numbers_1001G_sat_curve.txt")

plot(satur_curves_gene_numbers_1001G$accN,satur_curves_gene_numbers_1001G$AS_loci,ylim=c(0,4000))

points(satur_curves_gene_numbers_1001G$accN,satur_curves_gene_numbers_1001G$linc_loci,ylim=c(0,4000))

plot(satur_curves_gene_numbers_1001G$accN,satur_curves_gene_numbers_1001G$lncRNA_6types_loci,ylim=c(0,6000))
pcmax<-max(satur_curves_gene_numbers_1001G$PC_loci)
pseudomax<-max(satur_curves_gene_numbers_1001G$pseudo_loci)
asmax<-max(satur_curves_gene_numbers_1001G$AS_loci)
lincmax<-max(satur_curves_gene_numbers_1001G$linc_loci)
tegenemax<-max(satur_curves_gene_numbers_1001G$TEgene_loci)
tefragmax<-max(satur_curves_gene_numbers_1001G$TEs_loci)
lncmax<-max(satur_curves_gene_numbers_1001G$lncRNA_6types_loci)
  
seq_n<-seq(10, 460, by=10)
satur_curve1001<-as.data.frame(seq_n)
satur_curve1001$accN<-satur_curve1001[,1]
satur_curve1001$PC_percent<-0
satur_curve1001$PC_pecent_sd<-0
satur_curve1001$PC<-0
satur_curve1001$pseudo_percent<-0
satur_curve1001$pseudo_pecent_sd<-0
satur_curve1001$pseudo<-0
satur_curve1001$AS_percent<-0
satur_curve1001$AS_pecent_sd<-0
satur_curve1001$AS<-0
satur_curve1001$linc_percent<-0
satur_curve1001$linc_pecent_sd<-0
satur_curve1001$linc<-0
satur_curve1001$TEgene_percent<-0
satur_curve1001$TEgene_pecent_sd<-0
satur_curve1001$TEgene<-0
satur_curve1001$TEfrag_percent<-0
satur_curve1001$TEfrag_pecent_sd<-0
satur_curve1001$TEfrag<-0
satur_curve1001$lnc_percent<-0
satur_curve1001$lnc_pecent_sd<-0
satur_curve1001$lnc<-0

for ( i in 1:46){
  accN=i*10
  a<-satur_curves_gene_numbers_1001G[satur_curves_gene_numbers_1001G$accN==accN,]
  satur_curve1001$PC_percent[i]<-mean(a$PC_loci)*100/pcmax
  satur_curve1001$PC_pecent_sd[i]<-sd(a$PC_loci*100/pcmax)
  satur_curve1001$PC[i]<-mean(a$PC_loci)
  satur_curve1001$pseudo_percent[i]<-mean(a$pseudo_loci)*100/pseudomax
  satur_curve1001$pseudo_pecent_sd[i]<-sd(a$pseudo_loci*100/pseudomax)
  satur_curve1001$pseudo[i]<-mean(a$pseudo_loci)
  satur_curve1001$AS_percent[i]<-mean(a$AS_loci)*100/asmax
  satur_curve1001$AS_pecent_sd[i]<-sd(a$AS_loci*100/asmax)
  satur_curve1001$AS[i]<-mean(a$AS_loci)
  satur_curve1001$linc_percent[i]<-mean(a$linc_loci)*100/lincmax
  satur_curve1001$linc_pecent_sd[i]<-sd(a$linc_loci*100/lincmax)
  satur_curve1001$linc[i]<-mean(a$linc_loci)
  satur_curve1001$lnc_percent[i]<-mean(a$lncRNA_6types_loci)*100/lncmax
  satur_curve1001$lnc_pecent_sd[i]<-sd(a$lncRNA_6types_loci*100/lncmax)
  satur_curve1001$lnc[i]<-mean(a$lncRNA_6types_loci)
  satur_curve1001$TEgene_percent[i]<-mean(a$TEgene_loci)*100/tegenemax
  satur_curve1001$TEgene_pecent_sd[i]<-sd(a$TEgene_loci*100/tegenemax)
  satur_curve1001$TEgene[i]<-mean(a$TEgene_loci)
  satur_curve1001$TEfrag_percent[i]<-mean(a$TEs_loci)*100/tefragmax
  satur_curve1001$TEfrag_pecent_sd[i]<-sd(a$TEs_loci*100/tefragmax)
  satur_curve1001$TEfrag[i]<-mean(a$TEs_loci)
  
  
}


#fit polynomial curve 


x <- satur_curve1001$accN
y <- satur_curve1001$lnc_percent 
y.sd<-satur_curve1001$lnc_pecent_sd
z <- satur_curve1001$PC_percent 
z.sd<-satur_curve1001$PC_pecent_sd       


x1=as.data.frame(seq(10, 1000, by=10))
x1seq<-seq(10, 1000, by=10)
names(x1)<-"x"
# plot of x and y :

model_NC <- lm(y ~ x + I(log2(x)) )

model_PC <- lm(z ~ x + I(log2(x)) )

# I can get the features of this model :
summary(model_NC)
summary(model_PC)
model_NC$coefficients
(Intercept)           x  I(log2(x)) 
8.40967409  0.02370699  9.11054079 

model_PC$coefficients
summary(model_NC)$adj.r.squared
[1] 0.9989852
summary(model_PC)$adj.r.squared
[1] 0.9765819
# For each value of x, I can get the value of y estimated by the model, and add it to the current plot !
myPredict_NC <- predict( model_NC, data=x ) 
myPredict_PC <- predict( model_PC, data=x ) 

myPredict_NC_extra <- predict( model_NC, newdata = x1 ) 
myPredict_PC_extra <- predict( model_PC, newdata = x1 ) 


ix <- sort(x,index.return=T)$ix
ix1 <- sort(x1,index.return=T)$ix1





pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig1/Supplem_sat_curve_PC_NC_model_prediction.pdf",height = 4,width =5)
par(mar=c(4,4,2,1)) 
par(mfrow=c(1,1))
plot(x,y,col="green",pch=16 , cex=0.1, ylim=c(30,120),xlim=c(0,1000), main ="Gene identification saturation curve:\n fitting non-linear model ", ylab="percent of genes (loci) identified", xlab="number of accessions used for gene annotation",xaxt="n",yaxt="n") 
axis(1, at = seq(0, 1000, by = 100), las=2)
axis(2, at = seq(30, 120, by = 10), las=1)
abline(v=(seq(0,1000,100)), col="lightgray", lty="dotted")
abline(h=(seq(30,100,10)), col="lightgray", lty="dotted")

lines(x1seq, myPredict_NC_extra, col="darkgreen", lwd=1,ylim=c(30,200),xlim=c(0,1000) ) 
lines(x1seq, myPredict_PC_extra,  lwd=1,ylim=c(30,200),xlim=c(0,1000) , col="darkblue")  
#arrows(x0=x, y0=z-z.sd, x1=x, y1=z+z.sd, code=3, angle=90, length=0.02,col="darkgrey")
#arrows(x0=x, y0=y-y.sd, x1=x, y1=y+y.sd, code=3, angle=90, length=0.02,col="darkgrey")

text(600, 50 , paste("Model lncRNAs: ",round(model_NC$coefficients[1],2)," + ",round(model_NC$coefficients[2],2),"*x + ",round(model_NC$coefficients[3],2),"*log2(x)",sep=""), col="darkgreen",cex=0.8)
text(600, 45 , paste("R squared:",round(summary(model_NC)$adj.r.squared,3) ), col="darkgreen",cex=0.8)

text(600, 65, paste("Model PC genes: ",round(model_PC$coefficients[1],2)," + ",round(model_PC$coefficients[2],2),"*x + ",round(model_PC$coefficients[3],2),"*log2(x)",sep=""), col="darkblue",cex=0.8)
text(600, 60 , paste("R squared:",round(summary(model_PC)$adj.r.squared,3) ), col="darkblue",cex=0.8)
arrows(x0=x, y0=z-z.sd, x1=x, y1=z+z.sd, code=3, angle=90, length=0.02,col="#486EB4")
arrows(x0=x, y0=y-y.sd, x1=x, y1=y+y.sd, code=3, angle=90, length=0.02,col="darkgreen")
points(x,z,col="#486EB4",pch=16 , cex=0.5, ylim=c(30,120),xlim=c(0,1000)) 
points(x,y,col="green",pch=16 , cex=0.5, ylim=c(30,120),xlim=c(0,1000)) 
dev.off()


col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5")

col=c("#486EB4","#90C473","#F2AB54","#673A8E")





x <- satur_curve1001$accN
y_lnc <- satur_curve1001$lnc_percent 
y_lnc.sd<-satur_curve1001$lnc_pecent_sd
y_pc <- satur_curve1001$PC_percent 
y_pc.sd<-satur_curve1001$PC_pecent_sd       
y_linc <- satur_curve1001$linc_percent 
y_linc.sd<-satur_curve1001$linc_pecent_sd
y_as <- satur_curve1001$AS_percent
y_as.sd<-satur_curve1001$AS_pecent_sd
y_tegen <- satur_curve1001$TEgene_percent 
y_tegen.sd<-satur_curve1001$TEgene_pecent_sd
y_tefrag <- satur_curve1001$TEfrag_percent
y_tefrag.sd<-satur_curve1001$TEfrag_pecent_sd

model_NC <- lm(y_lnc ~ x + I(log2(x)) )
model_PC <- lm(y_pc ~ x + I(log2(x)) )
model_linc <- lm(y_linc ~ x + I(log2(x)) )
model_as <- lm(y_as ~ x + I(log2(x)) )
model_tegen <- lm(y_tegen ~ x + I(log2(x)) )
model_tefrag<- lm(y_tefrag ~ x + I(log2(x)) )


# For each value of x, I can get the value of y estimated by the model, and add it to the current plot !
myPredict_NC <- predict( model_NC, data=x ) 
myPredict_PC <- predict( model_PC, data=x ) 
myPredict_linc<- predict( model_linc, data=x ) 
myPredict_as <- predict( model_as, data=x ) 
myPredict_tegen <- predict( model_tegen, data=x ) 
myPredict_tefrag<- predict( model_tefrag, data=x ) 

#myPredict_NC_extra <- predict( model_NC, newdata = x1 ) 
#myPredict_PC_extra <- predict( model_PC, newdata = x1 ) 





col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5")




pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig1/Supplem_sat_curve_5genetypes.pdf",height = 4.5,width =5)
par(mar=c(4,4,2,1)) 
par(mfrow=c(1,1))
plot(x,y,col="white",pch=16 , cex=1, ylim=c(10,100),xlim=c(0,460), main ="Gene identification saturation curve", ylab="percent of genes (loci) identified", xlab="number of accessions used for gene annotation", xaxt="n",yaxt="n") 
axis(1, at = seq(0, 460, by = 20), las=2)
axis(2, at = seq(10, 100, by = 10), las=1)
abline(v=(seq(0,460,20)), col="lightgray", lty="dotted")
abline(h=(seq(10,100,10)), col="lightgray", lty="dotted")

lines(x, myPredict_PC,  lwd=1.5,ylim=c(30,200),xlim=c(0,1000) , col="#486EB4")  
lines(x, myPredict_as,  lwd=1.5,ylim=c(30,200),xlim=c(0,1000) , col="#90C473")  
lines(x, myPredict_linc,  lwd=1.5,ylim=c(30,200),xlim=c(0,1000) , col="#F2AB54")  
lines(x, myPredict_tefrag,  lwd=1.5,ylim=c(30,200),xlim=c(0,1000) , col="#805FA5")  
lines(x, myPredict_tegen,  lwd=1.5,ylim=c(30,200),xlim=c(0,1000) , col="#673A8E")  

arrows(x0=x, y0=y_pc-y_pc.sd, x1=x, y1=y_pc+y_pc.sd, code=3, angle=90, length=0.02,col="#486EB4")
arrows(x0=x, y0=y_linc-y_linc.sd, x1=x, y1=y_linc+y_linc.sd, code=3, angle=90, length=0.02,col="#F2AB54")
arrows(x0=x, y0=y_as-y_as.sd, x1=x, y1=y_as+y_as.sd, code=3, angle=90, length=0.02,col="#90C473")
arrows(x0=x, y0=y_tegen-y_tegen.sd, x1=x, y1=y_tegen+y_tegen.sd, code=3, angle=90, length=0.02,col="#673A8E")
arrows(x0=x, y0=y_tefrag-y_tefrag.sd, x1=x, y1=y_tefrag+y_tefrag.sd, code=3, angle=90, length=0.02,col="#805FA5")

points(satur_curve1001$accN,satur_curve1001$PC_percent,col="#486EB4",pch=16 , cex=0.5, ylim=c(30,100),xlim=c(0,460)) 
points(satur_curve1001$accN,satur_curve1001$linc_percent,col="#F2AB54",pch=16 , cex=0.5, ylim=c(30,100),xlim=c(0,460)) 
points(satur_curve1001$accN,satur_curve1001$AS_percent,col="#90C473",pch=16 , cex=0.5, ylim=c(30,100),xlim=c(0,460)) 
points(satur_curve1001$accN,satur_curve1001$TEgene_percent,col="#673A8E",pch=16 , cex=0.5, ylim=c(30,100),xlim=c(0,460)) 
points(satur_curve1001$accN,satur_curve1001$TEfrag_percent,col="#805FA5",pch=16 , cex=0.5, ylim=c(30,100),xlim=c(0,460)) 

text(310, 60, paste("AS lncRNA: ",round(model_as$coefficients[1],2)," + ",round(model_as$coefficients[2],2),"*x + ",round(model_as$coefficients[3],2),"*log2(x)",sep=""), col="#90C473",cex=0.8)
text(310, 55 , paste("R squared:",round(summary(model_as)$adj.r.squared,3) ), col="#90C473",cex=0.8)

text(310, 45 , paste("lincRNA: ",round(model_linc$coefficients[1],2)," + ",round(model_linc$coefficients[2],2),"*x + ",round(model_linc$coefficients[3],2),"*log2(x)",sep=""), col="#F2AB54",cex=0.8)
text(310, 40 , paste("R squared:",round(summary(model_linc)$adj.r.squared,3) ), col="#F2AB54",cex=0.8)

text(310, 30 , paste("TE genes: ",round(model_tegen$coefficients[1],2)," + ",round(model_tegen$coefficients[2],2),"*x + ",round(model_tegen$coefficients[3],2),"*log2(x)",sep=""), col="#673A8E",cex=0.8)
text(310, 25 , paste("R squared:",round(summary(model_tegen)$adj.r.squared,3) ), col="#673A8E",cex=0.8)

text(310, 15 , paste("TE frag.: ",round(model_tefrag$coefficients[1],2)," + ",round(model_tefrag$coefficients[2],2),"*x + ",round(model_tefrag$coefficients[3],2),"*log2(x)",sep=""), col="#805FA5",cex=0.8)
text(310, 10 , paste("R squared:",round(summary(model_tefrag)$adj.r.squared,3) ), col="#805FA5",cex=0.8)

dev.off()


# cortijo read N control 

#Z:\01_POSTDOC\03_Projects\2018_lncRNA_variation_paper\01_lncRNA_identification\02_Identification_Saturation_Curve\2021\20221109_merge_sat_curve_gene_numbers.bash
#Z:\01_POSTDOC\03_Projects\2018_lncRNA_variation_paper\01_lncRNA_identification\02_Identification_Saturation_Curve\2021\20221109_calculate_readN_for_1001G_satcurve_and_Cortijo_control.bash
#Z:\01_POSTDOC\03_Projects\2018_lncRNA_variation_paper\01_lncRNA_identification\02_Identification_Saturation_Curve\2021\20211013_cuffmerge_filtering_and_classification_automated_pipeline_for_satur_curve_rep1.bash
#Z:\01_POSTDOC\03_Projects\2018_lncRNA_variation_paper\01_lncRNA_identification\02_Identification_Saturation_Curve\2021\20221109_CortijoContorl_cuffmerge_filtering_and_classification_automated_pipeline_for_satur_curve_rep1.bash


PC_NC_genenumber_readnumber_1001G <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/02_Identification_Saturation_Curve/2021/1001G_sat_curve_PC_NC_genenumber_readnumber.txt")

PC_NC_genenumber_readnumber_Cortijo <- read.delim("Z:/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/02_Identification_Saturation_Curve/2021/Cortijo_sat_curve_PC_NC_genenumber_readnumber.txt")





x_1001 <- PC_NC_genenumber_readnumber_1001G$readN
z_1001<- PC_NC_genenumber_readnumber_1001G$accN
y_PC_1001  <- PC_NC_genenumber_readnumber_1001G$PC_loci 
y_NC_1001 <-PC_NC_genenumber_readnumber_1001G$lncRNA_6types_loci

x_Cort <- PC_NC_genenumber_readnumber_Cortijo$readN
y_PC_Cort  <- PC_NC_genenumber_readnumber_Cortijo$PC_loci 
y_NC_Cort <-PC_NC_genenumber_readnumber_Cortijo$lncRNA_6types_loci


model_NC_1001 <- lm(y_NC_1001 ~ x_1001 + I(log2(x_1001)))
model_PC_1001 <- lm(y_PC_1001 ~ x_1001 + I(log2(x_1001)) )

model_NC_Cort<- lm(y_NC_Cort ~ x_Cort + I(log2(x_Cort)) )
model_PC_Cort <- lm(y_PC_Cort ~ x_Cort + I(log2(x_Cort)) )
summary(model_NC_1001)
summary(model_NC_Cort)
model_NC_1001$coefficients
model_NC_Cort$coefficients
summary(model_NC_1001)$adj.r.squared
[1] 0.8452668
summary(model_NC_Cort)$adj.r.squared
[1] 0.8273996




plot(PC_NC_genenumber_readnumber_Cortijo$readN,PC_NC_genenumber_readnumber_Cortijo$lncRNA_6types_loci)
plot(PC_NC_genenumber_readnumber_Cortijo$readN,PC_NC_genenumber_readnumber_Cortijo$PC_loci)


# For each value of x, I can get the value of y estimated by the model, and add it to the current plot !
myPredict_NC_1001 <- predict( model_NC_1001, data=x_1001 ) 
myPredict_PC_1001 <- predict( model_PC_1001, data=x_1001 ) 
myPredict_NC_Cort <- predict( model_NC_Cort, data=x_Cort ) 
myPredict_PC_Cort <- predict( model_PC_Cort, data=x_Cort ) 


col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5")


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig1/Supplem_sat_curve_Cortijo_control.pdf",height = 4.5,width =5)
par(mar=c(4,4,2,1)) 
par(mfrow=c(1,1))
plot(x_1001,y_NC_1001,col="white",pch=16 , cex=1, ylim=c(10,10000),xlim=c(0,1000000000), main ="Gene identification saturation curve", ylab="percent of genes (loci) identified", xlab="number of accessions used for gene annotation", xaxt="n",yaxt="n") 
#axis(1, at = seq(0, 460, by = 20), las=2)
#axis(2, at = seq(10, 100, by = 10), las=1)
abline(v=(seq(0,460,20)), col="lightgray", lty="dotted")
abline(h=(seq(10,100,10)), col="lightgray", lty="dotted")

lines(x_1001, myPredict_NC_1001,  lwd=1.5, col="#486EB4")  
lines(x_Cort, myPredict_NC_Cort,  lwd=1.5, col="#486EB4")  

lines(x, myPredict_as,  lwd=1.5,ylim=c(30,200),xlim=c(0,1000) , col="#90C473")  
lines(x, myPredict_linc,  lwd=1.5,ylim=c(30,200),xlim=c(0,1000) , col="#F2AB54")  
lines(x, myPredict_tefrag,  lwd=1.5,ylim=c(30,200),xlim=c(0,1000) , col="#805FA5")  
lines(x, myPredict_tegen,  lwd=1.5,ylim=c(30,200),xlim=c(0,1000) , col="#673A8E")  

arrows(x0=x, y0=y_pc-y_pc.sd, x1=x, y1=y_pc+y_pc.sd, code=3, angle=90, length=0.02,col="#486EB4")
arrows(x0=x, y0=y_linc-y_linc.sd, x1=x, y1=y_linc+y_linc.sd, code=3, angle=90, length=0.02,col="#F2AB54")
arrows(x0=x, y0=y_as-y_as.sd, x1=x, y1=y_as+y_as.sd, code=3, angle=90, length=0.02,col="#90C473")
arrows(x0=x, y0=y_tegen-y_tegen.sd, x1=x, y1=y_tegen+y_tegen.sd, code=3, angle=90, length=0.02,col="#673A8E")
arrows(x0=x, y0=y_tefrag-y_tefrag.sd, x1=x, y1=y_tefrag+y_tefrag.sd, code=3, angle=90, length=0.02,col="#805FA5")

points(satur_curve1001$accN,satur_curve1001$PC_percent,col="#486EB4",pch=16 , cex=0.5, ylim=c(30,100),xlim=c(0,460)) 
points(satur_curve1001$accN,satur_curve1001$linc_percent,col="#F2AB54",pch=16 , cex=0.5, ylim=c(30,100),xlim=c(0,460)) 
points(satur_curve1001$accN,satur_curve1001$AS_percent,col="#90C473",pch=16 , cex=0.5, ylim=c(30,100),xlim=c(0,460)) 
points(satur_curve1001$accN,satur_curve1001$TEgene_percent,col="#673A8E",pch=16 , cex=0.5, ylim=c(30,100),xlim=c(0,460)) 
points(satur_curve1001$accN,satur_curve1001$TEfrag_percent,col="#805FA5",pch=16 , cex=0.5, ylim=c(30,100),xlim=c(0,460)) 

text(310, 60, paste("AS lncRNA: ",round(model_as$coefficients[1],2)," + ",round(model_as$coefficients[2],2),"*x + ",round(model_as$coefficients[3],2),"*log2(x)",sep=""), col="#90C473",cex=0.8)
text(310, 55 , paste("R squared:",round(summary(model_as)$adj.r.squared,3) ), col="#90C473",cex=0.8)

text(310, 45 , paste("lincRNA: ",round(model_linc$coefficients[1],2)," + ",round(model_linc$coefficients[2],2),"*x + ",round(model_linc$coefficients[3],2),"*log2(x)",sep=""), col="#F2AB54",cex=0.8)
text(310, 40 , paste("R squared:",round(summary(model_linc)$adj.r.squared,3) ), col="#F2AB54",cex=0.8)

text(310, 30 , paste("TE genes: ",round(model_tegen$coefficients[1],2)," + ",round(model_tegen$coefficients[2],2),"*x + ",round(model_tegen$coefficients[3],2),"*log2(x)",sep=""), col="#673A8E",cex=0.8)
text(310, 25 , paste("R squared:",round(summary(model_tegen)$adj.r.squared,3) ), col="#673A8E",cex=0.8)

text(310, 15 , paste("TE frag.: ",round(model_tefrag$coefficients[1],2)," + ",round(model_tefrag$coefficients[2],2),"*x + ",round(model_tefrag$coefficients[3],2),"*log2(x)",sep=""), col="#805FA5",cex=0.8)
text(310, 10 , paste("R squared:",round(summary(model_tefrag)$adj.r.squared,3) ), col="#805FA5",cex=0.8)

dev.off()


