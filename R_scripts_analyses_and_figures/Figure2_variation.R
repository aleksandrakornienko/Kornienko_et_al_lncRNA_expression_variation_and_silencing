install.packages("viridis")
install.packages("ggpubr")
install.packages("tidyverse")
install.packages("hrbrthemes")
library(MASS)
library(ggplot2)

library(viridis)
library(ggfortify)

library(ggpubr)
library(pheatmap)
library(tidyverse)
library(hrbrthemes)
library(viridis)

#> Loading required package: viridisLite
theme_set(theme_bw(base_size = 16))



########### intraindividual variation (in my replicated dataset 1001Gnew)

  
########### analyses and figures for FIGURE 2 
# Fig.2: lncRNAs display high inter-individual variability between A. thaliana natural accessions
###########
  

# in how many accessions are lncRNAs expressed? 

# figure 2A 
#1001G %

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_boxplot_how_many_accessions_1001G_percent.pdf",height = 5,width = 3)
par(mar=c(10,4,4,2)) 
boxplot(pc$Nacc_where_expressed05[pc$max>0.5]/4.61,as$Nacc_where_expressed05[as$max>0.5]/4.61,linc$Nacc_where_expressed05[linc$max>0.5]/4.61,te$Nacc_where_expressed05[te$max>0.5]/4.61,
        as_te$Nacc_where_expressed05[as_te$max>0.5]/4.61,te_frag$Nacc_where_expressed05[te_frag$max>0.5]/4.61,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "lncRNAs\nAS to TE genes","TE fragments"),las=2, notch = T, outline = F, ylab="% accessions where expressed, TPM>0.5")
dev.off()
#1001G % zoom in 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_boxplot_how_many_accessions_1001G_ylim_5percent.pdf",height = 5,width = 3)
par(mar=c(10,4,4,2)) 
boxplot(pc$Nacc_where_expressed05[pc$max>0.5]/4.61,as$Nacc_where_expressed05[as$max>0.5]/4.61,linc$Nacc_where_expressed05[linc$max>0.5]/4.61,te$Nacc_where_expressed05[te$max>0.5]/4.61,
        as_te$Nacc_where_expressed05[as_te$max>0.5]/4.61,te_frag$Nacc_where_expressed05[te_frag$max>0.5]/4.61,ylim=c(0,5),
        col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "lncRNAs\nAS to TE genes","TE fragments"),las=2, notch = T, outline = F, ylab="% accessions where expressed, TPM>0.5")
dev.off()


#############################
#expression frequency histogram plot FIG2 
###############################
  # pc genes
  "#486EB4"
  #lincs 
  "#F2AB54"
  #AS
  "#90C473"
  #AS to TEs
  "#B294C5"
  #AS to pseudogenes
  "#7F7F7F"
  #TE genes
  "#673A8E"
  #TE fragments 
  "#805FA5"
  
  
  
  #expression frequency 461 accessions - 1001G rosete 
  
  pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/expression_freq_461acc_1001G_FIG2.pdf",height = 5.6,width = 2.2)
 par(mfrow=c(4,1))
 par(mar=c(3,5,2,2))
 hist(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]/461,breaks=30, main="PC genes",col="#486EB4",xlim=c(0,1),xlab=" ",ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
 text(x=0.5,y=10000,adj=c(0.5,0.5),
      paste(
        "singletons",
        round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==1 &denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]),1)
        ,"%","\n",
        "all lines",
        round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==461 &denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]),1),"%"))

 
 hist(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]/461, main="AS lncRNAs",xlab=" ",col="#90C473",breaks=30,ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
 text(x=0.5,y=5000,adj=c(0.5,0.5),
      paste(
        "singletons",
        round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==1 &denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]),1)
        ,"%","\n",
        "all lines",
        round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==461 &denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]),1),"%"))

  hist(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]/461, main="lincRNAs",xlab=" ",col="#F2AB54",breaks=30,ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
  text(x=0.5,y=1500,
       paste(
         "singletons",
         round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==1 &denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]),1)
         ,"%","\n",
         "all lines",
         round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==461 &denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]),1),"%"))
  
  hist(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.transcripts$gene]/461, main="TE genes",xlab=" ",col="#673A8E",breaks=30,ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
  text(x=0.5,y=1500,
       paste(
         "singletons",
         round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==1 &denovo2021.TPMs.genes.1001G$gene %in% TE_genes.transcripts$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.transcripts$gene]),1)
         ,"%","\n",
         "all lines",
         round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==461 &denovo2021.TPMs.genes.1001G$gene %in% TE_genes.transcripts$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.transcripts$gene]),1),"%"))
  
  
  dev.off()
  
  
  
  
  
  pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/expression_freq_27acc_eracaps_rosette.pdf",height = 5.6,width = 2.2)
  par(mfrow=c(4,1))
  par(mar=c(3,5,2,2))
  hist(pc_EC$rosette.expr.frequency,breaks=30, main="PC genes",col="#486EB4",xlim=c(0,1),xlab=" ",ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
  text(x=0.5,y=10000,adj=c(0.5,0.5),
       paste(
         "singletons",
         round(100*length(pc_EC$rosette.expr.frequency[pc_EC$rosette.expr.frequency<0.05 & pc_EC$rosette.expr.frequency>0])/length(pc_EC$rosette.expr.frequency),1)
         ,"%","\n",
         "all lines",
         round(100*length(pc_EC$rosette.expr.frequency[pc_EC$rosette.expr.frequency==1])/length(pc_EC$rosette.expr.frequency),1),"%"))
  
  
  hist(as_EC$rosette.expr.frequency, main="AS lncRNAs",xlab=" ",col="#90C473",breaks=30,ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
  text(x=0.5,y=5000,adj=c(0.5,0.5),
       paste(
         "singletons",
         round(100*length(as_EC$rosette.expr.frequency[as_EC$rosette.expr.frequency<0.05 & as_EC$rosette.expr.frequency>0])/length(as_EC$rosette.expr.frequency),1)
         ,"%","\n",
         "all lines",
         round(100*length(as_EC$rosette.expr.frequency[as_EC$rosette.expr.frequency==1])/length(as_EC$rosette.expr.frequency),1),"%"))
  
  
  hist(linc_EC$rosette.expr.frequency, main="lincRNAs",xlab=" ",col="#F2AB54",breaks=30,ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
  text(x=0.5,y=1500,
       paste(
         "singletons",
         round(100*length(linc_EC$rosette.expr.frequency[linc_EC$rosette.expr.frequency<0.05 & linc_EC$rosette.expr.frequency>0])/length(linc_EC$rosette.expr.frequency),1)
         ,"%","\n",
         "all lines",
         round(100*length(linc_EC$rosette.expr.frequency[linc_EC$rosette.expr.frequency==1])/length(linc_EC$rosette.expr.frequency),1),"%"))
  
  
  hist(te_EC$rosette.expr.frequency, main="TE genes",xlab=" ",col="#673A8E",breaks=30,ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
  text(x=0.5,y=1500,
       paste(
         "singletons",
         round(100*length(te_EC$rosette.expr.frequency[te_EC$rosette.expr.frequency<0.05 & te_EC$rosette.expr.frequency>0])/length(te_EC$rosette.expr.frequency),1)
         ,"%","\n",
         "all lines",
         round(100*length(te_EC$rosette.expr.frequency[te_EC$rosette.expr.frequency==1])/length(te_EC$rosette.expr.frequency),1),"%"))
  
  dev.off()
  
  
  
  
  
  
  
  
  pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/expression_freq_27acc_eracaps_pollen.pdf",height = 5.6,width = 2.2)
  par(mfrow=c(4,1))
  par(mar=c(3,5,2,2))
  hist(pc_EC$pollen.expr.frequency,breaks=30, main="PC genes",col="#486EB4",xlim=c(0,1),xlab=" ",ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
  text(x=0.5,y=10000,adj=c(0.5,0.5),
       paste(
         "singletons",
         round(100*length(pc_EC$pollen.expr.frequency[pc_EC$pollen.expr.frequency<0.05 & pc_EC$pollen.expr.frequency>0])/length(pc_EC$pollen.expr.frequency),1)
         ,"%","\n",
         "all lines",
         round(100*length(pc_EC$pollen.expr.frequency[pc_EC$pollen.expr.frequency==1])/length(pc_EC$pollen.expr.frequency),1),"%"))
  
  
  hist(as_EC$pollen.expr.frequency, main="AS lncRNAs",xlab=" ",col="#90C473",breaks=30,ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
  text(x=0.5,y=5000,adj=c(0.5,0.5),
       paste(
         "singletons",
         round(100*length(as_EC$pollen.expr.frequency[as_EC$pollen.expr.frequency<0.05 & as_EC$pollen.expr.frequency>0])/length(as_EC$pollen.expr.frequency),1)
         ,"%","\n",
         "all lines",
         round(100*length(as_EC$pollen.expr.frequency[as_EC$pollen.expr.frequency==1])/length(as_EC$pollen.expr.frequency),1),"%"))
  
  
  hist(linc_EC$pollen.expr.frequency, main="lincRNAs",xlab=" ",col="#F2AB54",breaks=30,ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
  text(x=0.5,y=1500,
       paste(
         "singletons",
         round(100*length(linc_EC$pollen.expr.frequency[linc_EC$pollen.expr.frequency<0.05 & linc_EC$pollen.expr.frequency>0])/length(linc_EC$pollen.expr.frequency),1)
         ,"%","\n",
         "all lines",
         round(100*length(linc_EC$pollen.expr.frequency[linc_EC$pollen.expr.frequency==1])/length(linc_EC$pollen.expr.frequency),1),"%"))
  
  
  hist(te_EC$pollen.expr.frequency, main="TE genes",xlab=" ",col="#673A8E",breaks=30,ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
  text(x=0.5,y=1500,
       paste(
         "singletons",
         round(100*length(te_EC$pollen.expr.frequency[te_EC$pollen.expr.frequency<0.05 & te_EC$pollen.expr.frequency>0])/length(te_EC$pollen.expr.frequency),1)
         ,"%","\n",
         "all lines",
         round(100*length(te_EC$pollen.expr.frequency[te_EC$pollen.expr.frequency==1])/length(te_EC$pollen.expr.frequency),1),"%"))
  
  dev.off()
  
  
  
  
  
  pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/expression_freq_1001G_FIG2_2row2col.pdf",height = 4,width = 5.5)
  par(mfrow=c(2,2))
  par(mar=c(3,5,2,2))
  hist(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]/461,breaks=30, main="PC genes",col="#486EB4",xlim=c(0,1),xlab=" ",ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
  text(x=0.5,y=10000,adj=c(0.5,0.5),
       paste(
         "singletons",
         round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==1 &denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]),1)
         ,"%","\n",
         "all lines",
         round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==461 &denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]),1),"%"))
  
  
  hist(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]/461, main="AS lncRNAs",xlab=" ",col="#90C473",breaks=30,ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
  text(x=0.5,y=5000,adj=c(0.5,0.5),
       paste(
         "singletons",
         round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==1 &denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]),1)
         ,"%","\n",
         "all lines",
         round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==461 &denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]),1),"%"))
  
  hist(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]/461, main="lincRNAs",xlab=" ",col="#F2AB54",breaks=30,ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
  text(x=0.5,y=1500,
       paste(
         "singletons",
         round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==1 &denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]),1)
         ,"%","\n",
         "all lines",
         round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==461 &denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]),1),"%"))
  
  hist(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.transcripts$gene]/461, main="TE genes",xlab=" ",col="#673A8E",breaks=30,ylab="number of loci",cex.lab=1.3,cex.axis=1,las=2,mgp=c(3,1,0))
  text(x=0.5,y=1500,
       paste(
         "singletons",
         round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==1 &denovo2021.TPMs.genes.1001G$gene %in% TE_genes.transcripts$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.transcripts$gene]),1)
         ,"%","\n",
         "all lines",
         round(100*length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$Nacc_where_expressed05==461 &denovo2021.TPMs.genes.1001G$gene %in% TE_genes.transcripts$gene])/length(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.transcripts$gene]),1),"%"))
  
  
  dev.off()
  
 
 write.table(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene],"PC_nACC_expressed_05.txt",col.names = F,row.names = F)
 write.table(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene],"AS_nACC_expressed_05.txt",col.names = F,row.names = F)
 
 write.table(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene],"LINC_nACC_expressed_05.txt",col.names = F,row.names = F)
 
 write.table(denovo2021.TPMs.genes.1001G$Nacc_where_expressed05[denovo2021.TPMs.genes.1001G$gene %in% transcripts.expressed_TEs$gene],"TE_nACC_expressed_05.txt",col.names = F,row.names = F)
 
 
 #############################################################
 # in how many accessions on average are genes expressed? 
 ###############################################################
 
 
 #1001G - Figure 2A
 pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_boxplot_how_many_accessions_1001G.pdf",height = 5,width = 3)
 ######################################################
 par(mar=c(10,4,4,2)) 
 boxplot(pc$Nacc_where_expressed05[pc$max>0.5],as$Nacc_where_expressed05[as$max>0.5],linc$Nacc_where_expressed05[linc$max>0.5],te$Nacc_where_expressed05[te$max>0.5],
         as_te$Nacc_where_expressed05[as_te$max>0.5],te_frag$Nacc_where_expressed05[te_frag$max>0.5],
         col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "lncRNAs\nAS to TE genes","TE fragments"),las=2, notch = T, outline = F, ylab="# of accessions where expressed")
 ######################################################
 dev.off()
 
 

  
  #ERACAPS %
  pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_boxplot_how_many_accessions_EC_percent.pdf",height = 5,width = 3)
  ######################################################
  par(mar=c(10,4,4,2)) 
  boxplot(pc_EC$rosette.expr.frequency[pc_EC$max.rosette>0.5]*100,
          as_EC$rosette.expr.frequency[as_EC$max.rosette>0.5]*100,
          linc_EC$rosette.expr.frequency[linc_EC$max.rosette>0.5]*100,
          te_EC$rosette.expr.frequency[te_EC$max.rosette>0.5]*100,
          as_te_EC$rosette.expr.frequency[as_te_EC$max.rosette>0.5]*100,
          te_frag_EC$rosette.expr.frequency[te_frag_EC$max.rosette>0.5]*100,
          col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "lncRNAs\nAS to TE genes","TE fragments"),las=2, notch = T, outline = F, ylab="% accessions where expressed, TPM>0.5")
  ######################################################
  dev.off()
  
  #ERACAPS % zoom in 
  pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_boxplot_how_many_accessions_expresed_EC_ylim_5percent_tpmmax05.pdf",height = 5,width = 3)
  ######################################################
  par(mar=c(10,4,4,2)) 
  boxplot(pc_EC$rosette.expr.frequency[pc_EC$max.rosette>0.5]*100,
          as_EC$rosette.expr.frequency[as_EC$max.rosette>0.5]*100,
          linc_EC$rosette.expr.frequency[linc_EC$max.rosette>0.5]*100,
          te_EC$rosette.expr.frequency[te_EC$max.rosette>0.5]*100,
          as_te_EC$rosette.expr.frequency[as_te_EC$max.rosette>0.5]*100,
          te_frag_EC$rosette.expr.frequency[te_frag_EC$max.rosette>0.5]*100,ylim=c(0,100),
          col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "lncRNAs\nAS to TE genes","TE fragments"),las=2, notch = T, outline = F, ylab="% accessions where expressed, TPM>0.5")
  ######################################################
  dev.off()
  
  pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_boxplot_how_many_accessions_EC_4genetypes.pdf",height = 4,width = 3)
  ######################################################
  par(mar=c(10,4,4,2)) 
  boxplot(pc_EC$rosette.expr.frequency*100[pc_EC$max.rosette>0.5],
          as_EC$rosette.expr.frequency*100[as_EC$max.rosette>0.5],
          linc_EC$rosette.expr.frequency*100[linc_EC$max.rosette>0.5],
          te_EC$rosette.expr.frequency*100[te_EC$max.rosette>0.5],
   ylim=c(0,100),
          col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="% accessions where expressed, TPM>0.5")
  ######################################################
  dev.off()
  
  
  
  
  # number of accessions where the gene is present 
  
    pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_boxplot_how_many_accessions_locus_present_EC_percent.pdf",height = 4,width = 3)
    ######################################################
  par(mar=c(10,4,4,2)) 
  boxplot(pc_Oct2021_genetic_var_stats$inhowmanyaccessions[pc_Oct2021_genetic_var_stats$gene %in% pc_EC$gene[pc_EC$max.rosette>0.5]]*100/27,
          as_Oct2021_genetic_var_stats$inhowmanyaccessions[as_Oct2021_genetic_var_stats$gene %in% as_EC$gene[as_EC$max.rosette>0.5]]*100/27,
          linc_Oct2021_genetic_var_stats$inhowmanyaccessions[linc_Oct2021_genetic_var_stats$gene %in% linc_EC$gene[linc_EC$max.rosette>0.5]]*100/27,
          te_Oct2021_genetic_var_stats$inhowmanyaccessions[te_Oct2021_genetic_var_stats$gene %in% te_EC$gene[te_EC$max.rosette>0.5]]*100/27,
         
          col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="% accessions where expressed, TPM>0.5")
  ######################################################
  dev.off()
  
  
  
 
 
 
 
 #########################
 # Fig2 
 #expression variation boxplot - 1001G data 461 accession 
 par(mar=c(10,4,4,2)) 
 
 pc<-denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene,]
 as<-denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene,]
 linc<-denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene,]
 as_te<-denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.AS_to_TE.loci$gene,]
 te<-denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.transcripts$gene,]
 te_frag<-denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene %in% TE_frags.transcripts$gene,]
 
 pc_EC<-denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene %in% denovoPC.loci$gene,]
 as_EC<-denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.antisense.loci$gene,]
 linc_EC<-denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.intergenic.loci$gene,]
 as_te_EC<-denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.AS_to_TE.loci$gene,]
 te_EC<-denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene %in% TE_genes.transcripts$gene,]
 te_frag_EC<-denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene %in% TE_frags.transcripts$gene,]
 
 
 
 # pc genes  "#486EB4"
  #lincs  "#F2AB54"
 #AS  "#90C473"
  #AS to TEs  "#B294C5"
  #AS to pseudogenes  "#7F7F7F"
 #TE genes  "#673A8E"
 #TE fragments   "#805FA5"
 dev.off()
 
 
 pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2_boxplot_intervariance_1001G.pdf",height = 5,width = 5)
 par(mar=c(10,4,4,2)) 
 boxplot(pc$variance[pc$max>1],as$variance[as$max>1],linc$variance[linc$max>1],te$variance[te$max>1],
         as_te$variance[as_te$max>1],te_frag$variance[te_frag$max>1],
         col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "lncRNAs\nAS to TE genes","TE fragments"),las=2, notch = T, outline = F, ylab="coefficient of variance")
 dev.off()
 
 
 # check N in each box 
 
 length( pc$variance[pc$max>1]) #19253
 length(as$variance[as$max>1]) #2344
 length(linc$variance[linc$max>1]) #762 
 length(te$variance[te$max>1]) #900 
 length(as_te$variance[as_te$max>1]) # 195 
 length(te_frag$variance[te_frag$max>1]) #203 

 
 # check median value 
 median ( pc$variance[pc$max>1]) #0.4304806
 median(as$variance[as$max>1]) #1.256707
 median(linc$variance[linc$max>1]) #3.584102 
 median(te$variance[te$max>1]) #6.858215 
 median(as_te$variance[as_te$max>1]) # 6.197651 
 median(te_frag$variance[te_frag$max>1]) #6.976937 
 
 
 # check p values 
 
 wilcox.test(sample(pc$variance[pc$max>1],2344),as$variance[as$max>1]) # p-value < 2.2e-16
 wilcox.test(sample(as$variance[as$max>1],762),linc$variance[linc$max>1])  # p-value < 2.2e-16
 wilcox.test(linc$variance[linc$max>1],te$variance[te$max>1]) # p-value < 2.2e-16
 wilcox.test(te$variance[te$max>1],as_te$variance[as_te$max>1]) #p-value = 0.9761
 wilcox.test( as_te$variance[as_te$max>1],te_frag$variance[te_frag$max>1]) #p-value = 0.6598
 
 wilcox.test(sample(linc$variance[linc$max>1],195),as_te$variance[as_te$max>1]) # p-value = 3.313e-07 p-value = 1.873e-09 p-value = 1.433e-07
 
 
 
 pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2_boxplot_intervariance_1001G_max1-2.pdf",height = 5,width = 5)
 par(mar=c(10,4,4,2)) 
 boxplot(pc$variance[pc$max>1&pc$max<2],as$variance[as$max>1&as$max<2],linc$variance[linc$max>1&linc$max<2],te$variance[te$max>1&te$max<2],
         as_te$variance[as_te$max>1&as_te$max<2],te_frag$variance[te_frag$max>1&te_frag$max<2],
         col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "lncRNAs\nAS to TE genes","TE fragments"),las=2, notch = T, outline = F, ylab="coefficient of variance")
 dev.off()
 
 
 
 
 
 
 
 
 
 
 
 
 ##################################################
 ## Expression variability controls 
 ##################################################
 
 
 
 
 
 ###############################
 ### Controls 
 ################################### 
 install.packages("ggplot2")
 install.packages("Z:/programs/Rpackages/rlang_1.0.6.tar.gz", repos = NULL,lib ="Z:/programs/Rpackages" )
 install.packages("Z:/programs/Rpackages/vctrs_0.3.8.tar.gz", repos = NULL,lib ="Z:/programs/Rpackages" )
 install.packages ("magrittr",lib ="Z:/programs/Rpackages")
 install.packages ("ggthemes",lib ="Z:/programs/Rpackages")
 install.packages("pillar",lib ="Z:/programs/Rpackages")
 install.packages("ggplot2",lib ="Z:/programs/Rpackages")
 install.packages("ggfortify",lib ="Z:/programs/Rpackages")
 install.packages("ggbiplot",lib ="Z:/programs/Rpackages")
 install.packages("plotly",lib ="Z:/programs/Rpackages")
 install.packages("ggpubr",lib ="Z:/programs/Rpackages")
 
 install.packages("Z:/programs/Rpackages/vctrs_0.3.8.tar.gz" , repos=NULL, lib ="Z:/programs/Rpackages")
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
 #absolute expression level 
 library (scales)
 library (ggpubr)
 library(hexbin)
 
 
 
 # variance vs absolute expression level (max )
 
 pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/max_expression_1001G_4genetypes_supplem_tpmmax05.pdf",height = 4.5,width = 3)
 ######################################################
 par(mar=c(10,4,4,2)) 
 boxplot(pc$max[pc$max>0.5],as$max[as$max>0.5],linc$max[linc$max>0.5],te$max[te$max>0.5],log="y",
          col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="TPMmax among 461acc, log scale")
 ######################################################
 dev.off()
 
 
 pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/length_4genetypes_supplem_tpmmax05.pdf",height = 4.5,width = 3)
 ####################################################
 par(mar=c(10,4,4,2)) 
 boxplot(denovoPC.loci$length[denovoPC.loci$gene %in% pc$gene[pc$max>0.5]],
         lncRNAs.antisense.loci$length[lncRNAs.antisense.loci$gene %in% as$gene[as$max>0.5]],
         lncRNAs.intergenic.loci$length[lncRNAs.intergenic.loci$gene %in% linc$gene[linc$max>0.5]],
TE_genes.loci$length[TE_genes.loci$gene %in% te$gene[te$max>0.5]],log="y",
         col=c("#486EB4","#90C473","#F2AB54","#673A8E"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes"),las=2, notch = T, outline = F, ylab="gene length in bp, log scale")
 ######################################################
 dev.off()
 
 
 pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Control_variance_vs_MAX.1001G.PC.pdf",height = 2,width = 3)
 ######################################################
 par(mar=c(4,4,4,2)) 
  a<-merge( denovo2021.TPMs.genes.1001G[,c("max","mean","variance","gene","gene_type")],denovoPC.loci[,c("gene","length")])
 a$log_var<-log2(a$variance)
 a$log_len<-log2(a$length)
 a$log_max<-log2(a$max)
 a<-a[!is.na(a$variance),]
 a<-a[sample(1:length(a$gene),2000),]
 
ggplot(a[a$gene_type=="pc",], aes(x=log_max, y=log_var) ) +
   geom_hex(bins = 70) +
   scale_fill_continuous(type = "viridis") +
   geom_smooth(method='loess', formula= y~x,col="red") +
   ylim(-3,5) +
  xlim(-8, 15) +
   labs(title="PC genes",
        x ="TPMmax (461acc), log2", y = "coef. variance,log2")+
   theme_bw()
######################################################
dev.off() 

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Control_variance_vs_MAX.1001G.AS.pdf",height = 2,width = 3)
######################################################
par(mar=c(4,4,4,2)) 

 a<-merge( denovo2021.TPMs.genes.1001G[,c("max","mean","variance","gene","gene_type")],lncRNAs.antisense.loci[,c("gene","length")])
 a$log_var<-log2(a$variance)
 a$log_len<-log2(a$length)
 a$log_max<-log2(a$max)
 a<-a[!is.na(a$variance),]
 a<-a[sample(1:length(a$gene),2000),]
 
 ggplot(a[a$gene_type=="as",], aes(x=log_max, y=log_var) ) +
   geom_hex(bins = 70) +
   scale_fill_continuous(type = "viridis") +
   geom_smooth(method='loess', formula= y~x,col="red") +
   ylim(-3,5) +
   xlim(-8, 15) +
   labs(title="AS lncRNAs",
        x ="TPMmax (461acc), log2", y = "coef. variance,log2")+
   theme_bw()
 ######################################################
 dev.off() 
 
 pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Control_variance_vs_MAX.1001G.LINC.pdf",height = 2,width = 3)
 ######################################################
 par(mar=c(4,4,4,2)) 
 
 a<-merge( denovo2021.TPMs.genes.1001G[,c("max","mean","variance","gene","gene_type")],lncRNAs.intergenic.loci[,c("gene","length")])
 a$log_var<-log2(a$variance)
 a$log_len<-log2(a$length)
 a$log_max<-log2(a$max)
 a<-a[!is.na(a$variance),]
 #a<-a[sample(1:length(a$gene),2000),]
 
 ggplot(a[a$gene_type=="linc",], aes(x=log_max, y=log_var) ) +
   geom_hex(bins = 70) +
   scale_fill_continuous(type = "viridis") +
   geom_smooth(method='loess', formula= y~x,col="red") +
   ylim(-3,5) +
   xlim(-8, 15) +
   labs(title="lincRNAs",
        x ="TPMmax (461acc), log2", y = "coef. variance,log2")+
   theme_bw()
 ######################################################
 dev.off() 
 
 pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Control_variance_vs_MAX.1001G.TE.pdf",height = 2,width = 3)
 ######################################################
 par(mar=c(4,4,4,2)) 
 
 a<-merge( denovo2021.TPMs.genes.1001G[,c("max","mean","variance","gene","gene_type")],TE_genes.loci[,c("gene","length")])
 a$log_var<-log2(a$variance)
 a$log_len<-log2(a$length)
 a$log_max<-log2(a$max)
 a<-a[!is.na(a$variance),]
 ggplot(a, aes(x=log_max, y=log_var) ) +
   geom_hex(bins = 70) +
   scale_fill_continuous(type = "viridis") +
   geom_smooth(method='loess', formula= y~x,col="red") +
   ylim(-3,5) +
   xlim(-8, 15) +
   labs(title="TE genes",
        x ="TPMmax (461acc), log2", y = "coef. variance,log2")+
   theme_bw()
#####################################################

 dev.off() 
 
 
 
 
 
 
 
 
# variance vs length


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Control_variance_vs_LEN.1001G.PC.pdf",height = 2,width = 3)
###################################################
par(mar=c(4,4,4,2)) 
a<-merge( denovo2021.TPMs.genes.1001G[,c("max","mean","variance","gene","gene_type")],denovoPC.loci[,c("gene","length")])
a$log_var<-log2(a$variance)
a$log_len<-log2(a$length)
a$log_max<-log2(a$max)
a<-a[!is.na(a$variance),]
a<-a[sample(1:length(a$gene),2000),]

ggplot(a[a$gene_type=="pc",], aes(x=log_len, y=log_var) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method='loess', formula= y~x,col="red") +
  ylim(-3,5) +
  xlim(7, 15) +
  labs(title="PC genes",
       x ="gene length in bp, log2", y = "coef. variance,log2")+
  theme_bw()
###################################################
dev.off() 

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Control_variance_vs_LEN.1001G.AS.pdf",height = 2,width = 3)
par(mar=c(4,4,4,2)) 
###################################################
a<-merge( denovo2021.TPMs.genes.1001G[,c("max","mean","variance","gene","gene_type")],lncRNAs.antisense.loci[,c("gene","length")])
a$log_var<-log2(a$variance)
a$log_len<-log2(a$length)
a$log_max<-log2(a$max)
a<-a[!is.na(a$variance),]
a<-a[sample(1:length(a$gene),2000),]

ggplot(a[a$gene_type=="as",], aes(x=log_len, y=log_var) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method='loess', formula= y~x,col="red") +
  ylim(-3,5) +
  xlim(7, 15) +
  labs(title="AS lncRNAs",
       x ="gene length in bp, log2", y = "coef. variance,log2")+
  theme_bw()
###################################################
dev.off() 

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Control_variance_vs_LEN.1001G.LINC.pdf",height = 2,width = 3)
###################################################
par(mar=c(4,4,4,2)) 

a<-merge( denovo2021.TPMs.genes.1001G[,c("max","mean","variance","gene","gene_type")],lncRNAs.intergenic.loci[,c("gene","length")])
a$log_var<-log2(a$variance)
a$log_len<-log2(a$length)
a$log_max<-log2(a$max)
a<-a[!is.na(a$variance),]
#a<-a[sample(1:length(a$gene),2000),]

ggplot(a[a$gene_type=="linc",], aes(x=log_len, y=log_var) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method='loess', formula= y~x,col="red") +
  ylim(-3,5) +
  xlim(7, 15) +
  labs(title="lincRNAs",
       x ="gene length in bp, log2", y = "coef. variance,log2")+
  theme_bw()
###################################################
dev.off() 

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Control_variance_vs_LEN.1001G.TE.pdf",height = 2,width = 3)
###################################################
par(mar=c(4,4,4,2)) 

a<-merge( denovo2021.TPMs.genes.1001G[,c("max","mean","variance","gene","gene_type")],TE_genes.loci[,c("gene","length")])
a$log_var<-log2(a$variance)
a$log_len<-log2(a$length)
a$log_max<-log2(a$max)
a<-a[!is.na(a$variance),]
ggplot(a, aes(x=log_len, y=log_var) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method='loess', formula= y~x,col="red") +
  ylim(-3,5) +
  xlim(7, 15) +
  labs(title="TE genes",
       x ="gene length in bp, log2", y = "coef. variance,log2")+
  theme_bw()
###################################################
dev.off() 



 
 

##################################################
#control boxplots - certain lengh/expression range
#################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/SUPPLEM_Fig2_boxplot_intervariance_1001G_max1-2_length800-1500.pdf",height = 4,width = 3)
###################################################
par(mar=c(10,4,4,2)) 
a1<-pc$variance[pc$max>1 & pc$max<2 & pc$gene %in% denovoPC.loci$gene[denovoPC.loci$length<1500 &denovoPC.loci$length>500]]
a2<-      as$variance[as$max>1&as$max<2 & as$gene %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$length<1500 &lncRNAs.antisense.loci$length>500]]
a3<-     linc$variance[linc$max>1&linc$max<2  &   linc$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$length<1500 &lncRNAs.intergenic.loci$length>500]]
a4<-        te$variance[te$max>1&te$max<2& te$gene %in%TE_genes.loci$gene[TE_genes.loci$length<1500 &TE_genes.loci$length>500]]
a5<-      as_te$variance[as_te$max>1&as_te$max<2& as_te$gene %in%lncRNAs.AS_to_TE.loci$gene[lncRNAs.AS_to_TE.loci$length<1500 &lncRNAs.AS_to_TE.loci$length>500]]
a6<-      te_frag$variance[te_frag$max>1&te_frag$max<2& te_frag$gene %in%TE_frags.transcripts$gene[TE_frags.transcripts$length<1500 &TE_frags.transcripts$length>500]]
boxplot(a1,a2,a3,a4,a5,a6,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"),main="Expression variability between 461 acc\nonly genes with 1<TPMmax<2 and 800bp<Length<1500bp" ,names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE frag-s"),las=2, notch = T, outline = F, ylab="coefficient of variance")

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
text(b,x=1.5,y=15)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=15)
a<-wilcox.test(a4,a3)
b<-wilcox.test(a4,a3)
c<-wilcox.test(a4,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=15)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=15)
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=15)
#################
###################################################
dev.off()



pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/SUPPLEM_Fig2_boxplot_intervariance_1001G_max1-2.pdf",height = 4,width = 3)
###################################################
par(mar=c(10,4,4,2)) 
a1<-pc$variance[pc$max>1 & pc$max<2 ]
a2<-      as$variance[as$max>1&as$max<2 ]
a3<-     linc$variance[linc$max>1&linc$max<2 ]
a4<-        te$variance[te$max>1&te$max<2]
a5<-      as_te$variance[as_te$max>1&as_te$max<2]
a6<-      te_frag$variance[te_frag$max>1&te_frag$max<2]
boxplot(a1,a2,a3,a4,a5,a6,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"),main="Expression variability between 461 acc\nonly genes with 1<TPMmax<2 and 800bp<Length<1500bp" ,names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE frag-s"),las=2, notch = T, outline = F, ylab="coefficient of variance")

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
text(b,x=1.5,y=15)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=15)
a<-wilcox.test(a4,a3)
b<-wilcox.test(a4,a3)
c<-wilcox.test(a4,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=15)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=15)
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=15)
#################

dev.off()

length(pc$variance[pc$max>1 & pc$max<2 & pc$gene %in% denovoPC.loci$gene[denovoPC.loci$length<1500 &denovoPC.loci$length>500]])

length(as$variance[as$max>1&as$max<2 & as$gene %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$length<1500 &lncRNAs.antisense.loci$length>500]])

length(linc$variance[linc$max>1&linc$max<2  &   linc$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$length<1500 &lncRNAs.intergenic.loci$length>500]])
     
length(te$variance[te$max>1&te$max<2& te$gene %in%TE_genes.loci$gene[TE_genes.loci$length<1500 &TE_genes.loci$length>500]])

length(as_te$variance[as_te$max>1&as_te$max<2& as_te$gene %in%lncRNAs.AS_to_TE.loci$gene[lncRNAs.AS_to_TE.loci$length<1500 &lncRNAs.AS_to_TE.loci$length>500]])

length(te_frag$variance[te_frag$max>1&te_frag$max<2& te_frag$gene %in%TE_frags.transcripts$gene[TE_frags.transcripts$length<1500 &TE_frags.transcripts$length>500]])




pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/SUPPLEM_Fig2_boxplot_intervariance_1001G_max2-5_length500-1000.pdf",height = 4,width = 3)
###################################################
par(mar=c(10,4,4,2)) 
a1<-pc$variance[pc$max>2 & pc$max<5 & pc$gene %in% denovoPC.loci$gene[denovoPC.loci$length<1000 &denovoPC.loci$length>500]]
a2<-as$variance[as$max>2&as$max<5 & as$gene %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$length<1000 &lncRNAs.antisense.loci$length>500]]
a3<-linc$variance[linc$max>2&linc$max<5   &   linc$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$length<1000 &lncRNAs.intergenic.loci$length>500]]
a4<-te$variance[te$max>2&te$max<5& te$gene %in%TE_genes.loci$gene[TE_genes.loci$length<1000 &TE_genes.loci$length>500]]
a5<-as_te$variance[as_te$max>2&as_te$max<5& as_te$gene %in%lncRNAs.AS_to_TE.loci$gene[lncRNAs.AS_to_TE.loci$length<1000 &lncRNAs.AS_to_TE.loci$length>500]]
a6<-te_frag$variance[te_frag$max>2&te_frag$max<5& te_frag$gene %in%TE_frags.transcripts$gene[TE_frags.transcripts$length<1000 &TE_frags.transcripts$length>500]]
boxplot(a1,a2,a3,a4,a5,a6,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"),main="Expression variability between 461 acc\nonly genes with 2<TPMmax<5 and 500bp<Length<1000bp" , names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE frag-s"),las=2, notch = T, outline = F, ylab="coefficient of variance")

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
text(b,x=1.5,y=15)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=15)
a<-wilcox.test(a4,a3)
b<-wilcox.test(a4,a3)
c<-wilcox.test(a4,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=15)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=15)
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=15)
#################
dev.off()

 
length(pc$variance[pc$max>2 & pc$max<5 & pc$gene %in% denovoPC.loci$gene[denovoPC.loci$length<1000 &denovoPC.loci$length>500]])
length(as$variance[as$max>2&as$max<5 & as$gene %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$length<1000 &lncRNAs.antisense.loci$length>500]])
length(linc$variance[linc$max>2&linc$max<5   &   linc$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$length<1000 &lncRNAs.intergenic.loci$length>500]])
length(te$variance[te$max>2&te$max<5& te$gene %in%TE_genes.loci$gene[TE_genes.loci$length<1000 &TE_genes.loci$length>500]]) 
       length(as_te$variance[as_te$max>2&as_te$max<5& as_te$gene %in%lncRNAs.AS_to_TE.loci$gene[lncRNAs.AS_to_TE.loci$length<1000 &lncRNAs.AS_to_TE.loci$length>500]])
       length(te_frag$variance[te_frag$max>2&te_frag$max<5& te_frag$gene %in%TE_frags.transcripts$gene[TE_frags.transcripts$length<1000 &TE_frags.transcripts$length>500]])
 
 
##############################################################################
##############################################################################
# Intravariance
##############################################################################
       
       
 ########### intraindividual variation (in my replicated dataset 1001Gnew)
 pc_sasa<-denovo2021.TPMs.genes.1001Gnew[denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene,]
 as_sasa<-denovo2021.TPMs.genes.1001Gnew[denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene,]
 linc_sasa<-denovo2021.TPMs.genes.1001Gnew[denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene,]
 as_te_sasa<-denovo2021.TPMs.genes.1001Gnew[denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.AS_to_TE.loci$gene,]
 te_sasa<-denovo2021.TPMs.genes.1001Gnew[denovo2021.TPMs.genes.1001Gnew$gene %in% TE_genes.transcripts$gene,]
 te_frag_sasa<-denovo2021.TPMs.genes.1001Gnew[denovo2021.TPMs.genes.1001Gnew$gene %in% TE_frags.transcripts$gene,]
 
 
 
 # pc genes  "#486EB4"
 #lincs  "#F2AB54"
 #AS  "#90C473"
 #AS to TEs  "#B294C5"
 #AS to pseudogenes  "#7F7F7F"
 #TE genes  "#673A8E"
 #TE fragments   "#805FA5"
 
 
###################################################
 #boxplot inter vs intra variability 1001G new data TPMmax>1
 ###################################################
 
 pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/SUPPLEM_Fig2_boxplot_inter_and_intravariance_1001Gnew_max1.pdf",height = 4,width = 6)
  ###################################################
 
 a1<-pc_sasa$variance_of_means[pc_sasa$ma_x>1]
 a2<-as_sasa$variance_of_means[as_sasa$ma_x>1]
 a3<-linc_sasa$variance_of_means[linc_sasa$ma_x>1]
 a4<-te_sasa$variance_of_means[te_sasa$ma_x>1]
 a5<- as_te_sasa$variance_of_means[as_te_sasa$ma_x>1]
 a6<-te_frag_sasa$variance_of_means[te_frag_sasa$ma_x>1]
 
 b1<-pc_sasa$mean_intravariance[pc_sasa$ma_x>1]
 b2<- as_sasa$mean_intravariance[as_sasa$ma_x>1]
 b3<-linc_sasa$mean_intravariance[linc_sasa$ma_x>1]
 b4<- te_sasa$mean_intravariance[te_sasa$ma_x>1]
 b5<- as_te_sasa$mean_intravariance[as_te_sasa$ma_x>1]
 b6<- te_frag_sasa$mean_intravariance[te_frag_sasa$ma_x>1]
 
 
 boxplot(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,
          col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE frags","PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE frags"),las=2, notch = T, outline = F, ylab="coefficient of variance",main="inter- and intra-individual expression variability\n 28 accessions, 2-4 replicates")
 #################
 #add p values   #
 #################
 a<-wilcox.test(a1,b1)
 b<-wilcox.test(a1,b1)
 c<-wilcox.test(a1,b1)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=7,y=0)
 a<-wilcox.test(a2,b2)
 b<-wilcox.test(a2,b2)
 c<-wilcox.test(a2,b2)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=8,y=0)
 a<-wilcox.test(b3,a3)
 b<-wilcox.test(b3,a3)
 c<-wilcox.test(b3,a3)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=9,y=0)
 a<-wilcox.test(a4,b4)
 b<-wilcox.test(a4,b4)
 c<-wilcox.test(a4,b4)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=10,y=0)
 a<-wilcox.test(b5,a5)
 b<-wilcox.test(b5,a5)
 c<-wilcox.test(b5,a5)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=11,y=0)
 a<-wilcox.test(a6,b6)
 b<-wilcox.test(a6,b6)
 c<-wilcox.test(a6,b6)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=12,y=0)
 a<-wilcox.test(b1,b2)
 b<-wilcox.test(b1,b2)
 c<-wilcox.test(b1,b2)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=7.5,y=1.8)
 a<-wilcox.test(b2,b3)
 b<-wilcox.test(b2,b3)
 c<-wilcox.test(b2,b3)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=8.5,y=1.8)
 a<-wilcox.test(b4,b3)
 b<-wilcox.test(b4,b3)
 c<-wilcox.test(b4,b3)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=9.5,y=1.8)
 a<-wilcox.test(b4,b5)
 b<-wilcox.test(b4,b5)
 c<-wilcox.test(b4,b5)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=10.5,y=1.8)
 a<-wilcox.test(b6,b5)
 b<-wilcox.test(b6,b5)
 c<-wilcox.test(b6,b5)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=11.5,y=1.8)
 #################
 dev.off()
 
 
 ###################################################
 #boxplot intra variability 1001G new data TPMmax>1
 ###################################################
 pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/SUPPLEM_Fig2_boxplot_intravariance_1001Gnew_max1-2.pdf",height = 4,width = 4)
 ######################################################
 a1<-pc_sasa$mean_intravariance[pc_sasa$ma_x>1&pc_sasa$ma_x<2]
 a2<- as_sasa$mean_intravariance[as_sasa$ma_x>1&as_sasa$ma_x<2]
 a3<-linc_sasa$mean_intravariance[linc_sasa$ma_x>1&linc_sasa$ma_x<2]
 a4<- te_sasa$mean_intravariance[te_sasa$ma_x>1&te_sasa$ma_x<2]
 a5<- as_te_sasa$mean_intravariance[as_te_sasa$ma_x>1&as_te_sasa$ma_x<2]
 a6<- te_frag_sasa$mean_intravariance[te_frag_sasa$ma_x>1&te_frag_sasa$ma_x<2]
 boxplot(a1,a2,a3,a4,a5,a6,
  col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE frags"),las=2, notch = T, outline = F, ylab="coefficient of variance",main="intraindividual expression variability (2-4 reps)")
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
 a<-wilcox.test(a4,a3)
 b<-wilcox.test(a4,a3)
 c<-wilcox.test(a4,a3)
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
 a<-wilcox.test(a6,a5)
 b<-wilcox.test(a6,a5)
 c<-wilcox.test(a6,a5)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=5.5,y=2)
 #################
 dev.off()
 
 
 
 ###################################################
 #boxplot intra variability 1001G new data TPMmax>1 0.8kb<length<1.5kb
 ###################################################
  pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/SUPPLEM_Fig2_boxplot_intravariance_1001Gnew_max1-2_length800-1500.pdf",height = 4,width = 4)
 ######################################################
  a1<-pc_sasa$mean_intravariance[pc_sasa$ma_x>1&pc_sasa$ma_x<2 & pc_sasa$gene %in% denovoPC.loci$gene[denovoPC.loci$length<1500 &denovoPC.loci$length>800]]
  a2<- as_sasa$mean_intravariance[as_sasa$ma_x>1&as_sasa$ma_x<2 & as_sasa$gene %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$length<1500 &lncRNAs.antisense.loci$length>800]]
 a3<-linc_sasa$mean_intravariance[linc_sasa$ma_x>1&linc_sasa$ma_x<2 & linc_sasa$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$length<1500 &lncRNAs.intergenic.loci$length>800]]
 a4<- te_sasa$mean_intravariance[te_sasa$ma_x>1&te_sasa$ma_x<2 & te_sasa$gene %in% TE_genes.loci$gene[TE_genes.loci$length<1500 &TE_genes.loci$length>800]]
 a5<- as_te_sasa$mean_intravariance[as_te_sasa$ma_x>1&as_te_sasa$ma_x<2 & as_te_sasa$gene %in% lncRNAs.AS_to_TE.loci$gene[lncRNAs.AS_to_TE.loci$length<1500 &lncRNAs.AS_to_TE.loci$length>800]]
 a6<- te_frag_sasa$mean_intravariance[te_frag_sasa$ma_x>1&te_frag_sasa$ma_x<2 & te_frag_sasa$gene %in% TE_frags.transcripts$gene[TE_frags.transcripts$length<1500 &TE_frags.transcripts$length>800]]
  boxplot(a1,a2,a3,a4,a5,a6,
         col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE frags"),las=2, notch = T, outline = F, ylab="coefficient of variance",main="intraindividual expression variability (2-4 reps)\n 1<TPMmax<2 and 0.8kb<length<1.5kb")
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
 text(b,x=1.5,y=1.7)
 a<-wilcox.test(a2,a3)
 b<-wilcox.test(a2,a3)
 c<-wilcox.test(a2,a3)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=2.5,y=1.7)
 a<-wilcox.test(a4,a3)
 b<-wilcox.test(a4,a3)
 c<-wilcox.test(a4,a3)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=3.5,y=1.7)
 a<-wilcox.test(a4,a5)
 b<-wilcox.test(a4,a5)
 c<-wilcox.test(a4,a5)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=4.5,y=1.7)
 a<-wilcox.test(a6,a5)
 b<-wilcox.test(a6,a5)
 c<-wilcox.test(a6,a5)
 d<-mean(c(a$p.value,b$p.value,c$p.value))
 if (d<0.0000000001){b="***"} 
 if (d>=0.0000000001 &d<0.00001 ){b="**"} 
 if (d>=0.00001  &d<0.01){b="*"} 
 if (d>=0.01){b="n.s."} 
 text(b,x=5.5,y=1.7)
 #################
 dev.off()
 
 
 
 
 
 
 boxplot(pc_sasa$var_3random_8[pc_sasa$ma_x>1],as_sasa$var_3random_mean[as_sasa$ma_x>1],linc_sasa$var_3random_mean[linc_sasa$ma_x>1],te_sasa$var_3random_mean[te_sasa$ma_x>1],
         as_te_sasa$var_3random_mean[as_te_sasa$ma_x>1],te_frag_sasa$var_3random_mean[te_frag_sasa$ma_x>1],
         
         pc_sasa$mean_intravariance[pc_sasa$ma_x>1],as_sasa$mean_intravariance[as_sasa$ma_x>1],linc_sasa$mean_intravariance[linc_sasa$ma_x>1],te_sasa$mean_intravariance[te_sasa$ma_x>1],
         as_te_sasa$mean_intravariance[as_te_sasa$ma_x>1],te_frag_sasa$mean_intravariance[te_frag_sasa$ma_x>1],
         
         col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE frags","PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE frags"),las=2, notch = T, outline = F, ylab="coefficient of variance",main="inter- and intra-individual expression variability")
 
 
 
 
 dev.off()

  boxplot(pc_sasa$mean_intravariance[pc_sasa$mean_of_means>0.5],as_sasa$mean_intravariance[as_sasa$mean_of_means>0.5],linc_sasa$mean_intravariance[linc_sasa$mean_of_means>0.5],te_sasa$mean_intravariance[te_sasa$mean_of_means>0.5],
         as_te_sasa$mean_intravariance[as_te_sasa$mean_of_means>0.5],te_frag_sasa$mean_intravariance[te_frag_sasa$mean_of_means>0.5],
         col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "lncRNAs\nAS to TE genes","TE fragments"),las=2, notch = T, outline = F, ylab="coefficient of variance",main="intra-individual expression variability")
 
 
 
 
 
 
 
 
 #anova test with replicates
  
  
 
 
 
  
  
  
#####################
# how many lincRNAs are shared between two accessions? 
##############################
  
a<-as.vector(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene & denovo2021.TPMs.genes.1001G$X6906>0.5]) #93
b<-as.vector(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene & denovo2021.TPMs.genes.1001G$X6390>0.5])#77
  
length(intersect(a,b))    #41

a<-as.vector(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene & denovo2021.TPMs.genes.1001G$X6906>0.5]) #545
b<-as.vector(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene & denovo2021.TPMs.genes.1001G$X6390>0.5])#530

length(intersect(a,b))    #323


a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.intergenic.loci$gene & denovo2021.TPMs.genes.ERACAPS$R.10015>0.5]) #100
b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.intergenic.loci$gene & denovo2021.TPMs.genes.ERACAPS$R.6024>0.5]) #63
length(intersect(a,b)) #40

a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.antisense.loci$gene & denovo2021.TPMs.genes.ERACAPS$R.10015>0.5]) #479
b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.antisense.loci$gene & denovo2021.TPMs.genes.ERACAPS$R.6024>0.5]) #471
length(intersect(a,b)) #276
#####  

# make a table of shared genes 
tab_linc_shared<-data.frame(1:461)
tab_as_shared<-data.frame(1:461)
tab_pc_shared<-data.frame(1:461)
tab_te_shared<-data.frame(1:461)
tab_linc_shared$acc1<-"X"
tab_linc_shared$acc2<-"X"
tab_linc_shared$N_acc1<-0
tab_linc_shared$N_acc2<-0
tab_linc_shared$N_shared<-0

tab_as_shared$acc1<-"X"
tab_as_shared$acc2<-"X"
tab_as_shared$N_acc1<-0
tab_as_shared$N_acc2<-0
tab_as_shared$N_shared<-0

tab_pc_shared$acc1<-"X"
tab_pc_shared$acc2<-"X"
tab_pc_shared$N_acc1<-0
tab_pc_shared$N_acc2<-0
tab_pc_shared$N_shared<-0

tab_te_shared$acc1<-"X"
tab_te_shared$acc2<-"X"
tab_te_shared$N_acc1<-0
tab_te_shared$N_acc2<-0
tab_te_shared$N_shared<-0



  for (i in 1:461) {
    acc1<-as.character(names(denovo2021.TPMs.genes.1001G)[i+1])
    acc2<-as.character(names(denovo2021.TPMs.genes.1001G)[sample(2:462,1)])
    tab_linc_shared$acc1[i]<-acc1
    tab_linc_shared$acc2[i]<-acc2

    tab_as_shared$acc1[i]<-acc1
    tab_as_shared$acc2[i]<-acc2
    
    tab_pc_shared$acc1[i]<-acc1
    tab_pc_shared$acc2[i]<-acc2
    
    tab_te_shared$acc1[i]<-acc1
    tab_te_shared$acc2[i]<-acc2
    
        a<-as.vector(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene & denovo2021.TPMs.genes.1001G[,as.character(acc1)]>0.5]) 
    b<-as.vector(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene & denovo2021.TPMs.genes.1001G[,as.character(acc2)]>0.5])
    tab_linc_shared$N_acc1[i]<-length(a)
    tab_linc_shared$N_acc2[i]<-length(b)
    tab_linc_shared$N_shared[i]<- length(intersect(a,b))
   
    a<-as.vector(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene & denovo2021.TPMs.genes.1001G[,as.character(acc1)]>0.5]) 
    b<-as.vector(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene & denovo2021.TPMs.genes.1001G[,as.character(acc2)]>0.5])
    tab_pc_shared$N_acc1[i]<-length(a)
    tab_pc_shared$N_acc2[i]<-length(b)
    tab_pc_shared$N_shared[i]<- length(intersect(a,b))
    
    a<-as.vector(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.loci$gene & denovo2021.TPMs.genes.1001G[,as.character(acc1)]>0.5]) 
    b<-as.vector(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.loci$gene & denovo2021.TPMs.genes.1001G[,as.character(acc2)]>0.5])
    tab_te_shared$N_acc1[i]<-length(a)
    tab_te_shared$N_acc2[i]<-length(b)
    tab_te_shared$N_shared[i]<- length(intersect(a,b))
    
    a<-as.vector(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene & denovo2021.TPMs.genes.1001G[,as.character(acc1)]>0.5]) 
    b<-as.vector(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene & denovo2021.TPMs.genes.1001G[,as.character(acc2)]>0.5])
    tab_as_shared$N_acc1[i]<-length(a)
    tab_as_shared$N_acc2[i]<-length(b)
    tab_as_shared$N_shared[i]<- length(intersect(a,b))
    
  }


tab_linc_shared<-tab_linc_shared[tab_linc_shared$acc1!=tab_linc_shared$acc2,2:6]
tab_as_shared<-tab_as_shared[tab_as_shared$acc1!=tab_as_shared$acc2,2:6]
tab_pc_shared<-tab_pc_shared[tab_pc_shared$acc1!=tab_pc_shared$acc2,2:6]
tab_te_shared<-tab_te_shared[tab_te_shared$acc1!=tab_te_shared$acc2,2:6]


mean(tab_linc_shared$N_acc1) #96.76356
mean(tab_linc_shared$N_acc2) #97.16486
mean(tab_linc_shared$N_shared) #48.54664

mean(tab_as_shared$N_acc1) #543.9566
mean(tab_as_shared$N_acc2) #540.872
mean(tab_as_shared$N_shared)#313.7202

mean(tab_pc_shared$N_acc1) #15918.68
mean(tab_pc_shared$N_acc2) #15932.11
mean(tab_pc_shared$N_shared) #15228.41

mean(tab_te_shared$N_acc1) #62.73753
mean(tab_te_shared$N_acc2) #63.10846
mean(tab_te_shared$N_shared) #29.10846

par(mfrow=c(2,2))
par(mar=c(10,4,4,2)) 


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/venndiag_PC_shared_2_acc.pdf",height = 2,width = 2)
par(mar=c(10,4,4,2)) 

a<-round(mean(tab_pc_shared$N_acc1),digits = 1)
b<-round(mean(tab_pc_shared$N_acc2),digits = 1)
c<-round(mean(tab_pc_shared$N_shared),digits = 1)
draw.pairwise.venn(area1=a, area2=b,cross.area=c,
                   category=c("acc1","acc2"),fill=c("#486EB4","#486EB4"),)
round(100*mean(tab_pc_shared$N_shared)/mean(tab_pc_shared$N_acc1),digits = 1)
dev.off()

dev.off()
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/venndiag_LINC_shared_2_acc.pdf",height = 2,width = 2)
par(mar=c(10,4,4,2)) 

a<-round(mean(tab_linc_shared$N_acc1),digits = 1)
b<-round(mean(tab_linc_shared$N_acc2),digits = 1)
c<-round(mean(tab_linc_shared$N_shared),digits = 1)
draw.pairwise.venn(area1=a, area2=b,cross.area=c,
                   category=c("acc1","acc2"),fill=c("#F2AB54","#F2AB54"))
round(100*mean(tab_linc_shared$N_shared)/mean(tab_linc_shared$N_acc1),digits = 1)

dev.off()

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/venndiag_AS_shared_2_acc.pdf",height = 2,width = 2)
par(mar=c(10,4,4,2)) 

a<-round(mean(tab_as_shared$N_acc1),digits = 1)
b<-round(mean(tab_as_shared$N_acc2),digits = 1)
c<-round(mean(tab_as_shared$N_shared),digits = 1)
draw.pairwise.venn(area1=a, area2=b,cross.area=c,
                   category=c("acc1","acc2"),fill=c("#90C473","#90C473"))
round(100*mean(tab_as_shared$N_shared)/mean(tab_as_shared$N_acc1),digits = 1)

dev.off()

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/venndiag_TE_shared_2_acc.pdf",height = 2,width = 2)
par(mar=c(10,4,4,2)) 

a<-round(mean(tab_te_shared$N_acc1),digits = 1)
b<-round(mean(tab_te_shared$N_acc2),digits = 1)
c<-round(mean(tab_te_shared$N_shared),digits = 1)
draw.pairwise.venn(area1=a, area2=b,cross.area=c,
                   category=c("acc1","acc2"),fill=c("#673A8E","#673A8E"))
round(100*mean(tab_te_shared$N_shared)/mean(tab_te_shared$N_acc1),digits = 1)

dev.off()






# make a table of shared genes 
tab_linc_shared_EC<-data.frame(1:23)
tab_as_shared_EC<-data.frame(1:23)
tab_pc_shared_EC<-data.frame(1:23)
tab_te_shared_EC<-data.frame(1:23)
tab_linc_shared_EC$acc1<-"X"
tab_linc_shared_EC$acc2<-"X"
tab_linc_shared_EC$N_acc1_S<-0
tab_linc_shared_EC$N_acc2_S<-0
tab_linc_shared_EC$N_shared_S<-0
tab_linc_shared_EC$N_acc1_R<-0
tab_linc_shared_EC$N_acc2_R<-0
tab_linc_shared_EC$N_shared_R<-0
tab_linc_shared_EC$N_acc1_F<-0
tab_linc_shared_EC$N_acc2_F<-0
tab_linc_shared_EC$N_shared_F<-0
tab_linc_shared_EC$N_acc1_P<-0
tab_linc_shared_EC$N_acc2_P<-0
tab_linc_shared_EC$N_shared_P<-0

tab_as_shared_EC$acc1<-"X"
tab_as_shared_EC$acc2<-"X"
tab_as_shared_EC$N_acc1_S<-0
tab_as_shared_EC$N_acc2_S<-0
tab_as_shared_EC$N_shared_S<-0
tab_as_shared_EC$N_acc1_R<-0
tab_as_shared_EC$N_acc2_R<-0
tab_as_shared_EC$N_shared_R<-0
tab_as_shared_EC$N_acc1_F<-0
tab_as_shared_EC$N_acc2_F<-0
tab_as_shared_EC$N_shared_F<-0
tab_as_shared_EC$N_acc1_P<-0
tab_as_shared_EC$N_acc2_P<-0
tab_as_shared_EC$N_shared_P<-0

tab_pc_shared_EC$acc1<-"X"
tab_pc_shared_EC$acc2<-"X"
tab_pc_shared_EC$N_acc1_S<-0
tab_pc_shared_EC$N_acc2_S<-0
tab_pc_shared_EC$N_shared_S<-0
tab_pc_shared_EC$N_acc1_R<-0
tab_pc_shared_EC$N_acc2_R<-0
tab_pc_shared_EC$N_shared_R<-0
tab_pc_shared_EC$N_acc1_F<-0
tab_pc_shared_EC$N_acc2_F<-0
tab_pc_shared_EC$N_shared_F<-0
tab_pc_shared_EC$N_acc1_P<-0
tab_pc_shared_EC$N_acc2_P<-0
tab_pc_shared_EC$N_shared_P<-0

tab_te_shared_EC$acc1<-"X"
tab_te_shared_EC$acc2<-"X"
tab_te_shared_EC$N_acc1_S<-0
tab_te_shared_EC$N_acc2_S<-0
tab_te_shared_EC$N_shared_S<-0
tab_te_shared_EC$N_acc1_R<-0
tab_te_shared_EC$N_acc2_R<-0
tab_te_shared_EC$N_shared_R<-0
tab_te_shared_EC$N_acc1_F<-0
tab_te_shared_EC$N_acc2_F<-0
tab_te_shared_EC$N_shared_F<-0
tab_te_shared_EC$N_acc1_P<-0
tab_te_shared_EC$N_acc2_P<-0
tab_te_shared_EC$N_shared_P<-0



for (i in 1:23) {
  acc1<-as.character(ERACAPS_accession_list [i,1])
  acc2<-as.character(ERACAPS_accession_list [sample(1:27,1),1])
  tab_linc_shared_EC$acc1[i]<-acc1
  tab_linc_shared_EC$acc2[i]<-acc2
  
  tab_as_shared_EC$acc1[i]<-acc1
  tab_as_shared_EC$acc2[i]<-acc2
  
  tab_pc_shared_EC$acc1[i]<-acc1
  tab_pc_shared_EC$acc2[i]<-acc2
  
  tab_te_shared_EC$acc1[i]<-acc1
  tab_te_shared_EC$acc2[i]<-acc2
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.intergenic.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("S.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.intergenic.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("S.",as.character(acc2),sep = '')]>0.5])
  tab_linc_shared_EC$N_acc1_S[i]<-length(a)
  tab_linc_shared_EC$N_acc2_S[i]<-length(b)
  tab_linc_shared_EC$N_shared_S[i]<- length(intersect(a,b))
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.intergenic.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("R.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.intergenic.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("R.",as.character(acc2),sep = '')]>0.5])
  tab_linc_shared_EC$N_acc1_R[i]<-length(a)
  tab_linc_shared_EC$N_acc2_R[i]<-length(b)
  tab_linc_shared_EC$N_shared_R[i]<- length(intersect(a,b))
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.intergenic.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("F.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.intergenic.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("F.",as.character(acc2),sep = '')]>0.5])
  tab_linc_shared_EC$N_acc1_F[i]<-length(a)
  tab_linc_shared_EC$N_acc2_F[i]<-length(b)
  tab_linc_shared_EC$N_shared_F[i]<- length(intersect(a,b))
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.intergenic.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("P.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.intergenic.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("P.",as.character(acc2),sep = '')]>0.5])
  tab_linc_shared_EC$N_acc1_P[i]<-length(a)
  tab_linc_shared_EC$N_acc2_P[i]<-length(b)
  tab_linc_shared_EC$N_shared_P[i]<- length(intersect(a,b))
  
  
  
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.antisense.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("S.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.antisense.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("S.",as.character(acc2),sep = '')]>0.5])
  tab_as_shared_EC$N_acc1_S[i]<-length(a)
  tab_as_shared_EC$N_acc2_S[i]<-length(b)
  tab_as_shared_EC$N_shared_S[i]<- length(intersect(a,b))
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.antisense.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("R.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.antisense.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("R.",as.character(acc2),sep = '')]>0.5])
  tab_as_shared_EC$N_acc1_R[i]<-length(a)
  tab_as_shared_EC$N_acc2_R[i]<-length(b)
  tab_as_shared_EC$N_shared_R[i]<- length(intersect(a,b))
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.antisense.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("F.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.antisense.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("F.",as.character(acc2),sep = '')]>0.5])
  tab_as_shared_EC$N_acc1_F[i]<-length(a)
  tab_as_shared_EC$N_acc2_F[i]<-length(b)
  tab_as_shared_EC$N_shared_F[i]<- length(intersect(a,b))
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.antisense.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("P.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.antisense.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("P.",as.character(acc2),sep = '')]>0.5])
  tab_as_shared_EC$N_acc1_P[i]<-length(a)
  tab_as_shared_EC$N_acc2_P[i]<-length(b)
  tab_as_shared_EC$N_shared_P[i]<- length(intersect(a,b))
  
  
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% denovoPC.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("S.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% denovoPC.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("S.",as.character(acc2),sep = '')]>0.5])
  tab_pc_shared_EC$N_acc1_S[i]<-length(a)
  tab_pc_shared_EC$N_acc2_S[i]<-length(b)
  tab_pc_shared_EC$N_shared_S[i]<- length(intersect(a,b))
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% denovoPC.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("R.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% denovoPC.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("R.",as.character(acc2),sep = '')]>0.5])
  tab_pc_shared_EC$N_acc1_R[i]<-length(a)
  tab_pc_shared_EC$N_acc2_R[i]<-length(b)
  tab_pc_shared_EC$N_shared_R[i]<- length(intersect(a,b))
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% denovoPC.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("F.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% denovoPC.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("F.",as.character(acc2),sep = '')]>0.5])
  tab_pc_shared_EC$N_acc1_F[i]<-length(a)
  tab_pc_shared_EC$N_acc2_F[i]<-length(b)
  tab_pc_shared_EC$N_shared_F[i]<- length(intersect(a,b))
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% denovoPC.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("P.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% denovoPC.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("P.",as.character(acc2),sep = '')]>0.5])
  tab_pc_shared_EC$N_acc1_P[i]<-length(a)
  tab_pc_shared_EC$N_acc2_P[i]<-length(b)
  tab_pc_shared_EC$N_shared_P[i]<- length(intersect(a,b))
  
  
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% TE_genes.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("S.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% TE_genes.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("S.",as.character(acc2),sep = '')]>0.5])
  tab_te_shared_EC$N_acc1_S[i]<-length(a)
  tab_te_shared_EC$N_acc2_S[i]<-length(b)
  tab_te_shared_EC$N_shared_S[i]<- length(intersect(a,b))
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% TE_genes.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("R.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% TE_genes.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("R.",as.character(acc2),sep = '')]>0.5])
  tab_te_shared_EC$N_acc1_R[i]<-length(a)
  tab_te_shared_EC$N_acc2_R[i]<-length(b)
  tab_te_shared_EC$N_shared_R[i]<- length(intersect(a,b))
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% TE_genes.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("F.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% TE_genes.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("F.",as.character(acc2),sep = '')]>0.5])
  tab_te_shared_EC$N_acc1_F[i]<-length(a)
  tab_te_shared_EC$N_acc2_F[i]<-length(b)
  tab_te_shared_EC$N_shared_F[i]<- length(intersect(a,b))
  
  a<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% TE_genes.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("P.",as.character(acc1),sep = '')]>0.5]) 
  b<-as.vector(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$gene %in% TE_genes.loci$gene & denovo2021.TPMs.genes.ERACAPS[,paste("P.",as.character(acc2),sep = '')]>0.5])
  tab_te_shared_EC$N_acc1_P[i]<-length(a)
  tab_te_shared_EC$N_acc2_P[i]<-length(b)
  tab_te_shared_EC$N_shared_P[i]<- length(intersect(a,b))
  
  
  
}



tab_linc_shared_EC<-tab_linc_shared_EC[tab_linc_shared_EC$acc1!=tab_linc_shared_EC$acc2,]
tab_as_shared_EC<-tab_as_shared_EC[tab_as_shared_EC$acc1!=tab_as_shared_EC$acc2,]
tab_pc_shared_EC<-tab_pc_shared_EC[tab_pc_shared_EC$acc1!=tab_pc_shared_EC$acc2,]
tab_te_shared_EC<-tab_te_shared_EC[tab_te_shared_EC$acc1!=tab_te_shared_EC$acc2,]


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/venndiag_PC_shared_2_acc_EC_S.pdf",height = 2,width = 2)
par(mar=c(10,4,4,2)) 

a<-round(mean(tab_pc_shared_EC$N_acc1_S),digits = 1)
b<-round(mean(tab_pc_shared_EC$N_acc2_S),digits = 1)
c<-round(mean(tab_pc_shared_EC$N_shared_S),digits = 1)
draw.pairwise.venn(area1=a, area2=b,cross.area=c,
                   category=c("acc1","acc2"),fill=c("#486EB4","#486EB4"),)
round(100*mean(tab_pc_shared_EC$N_shared_S)/mean(tab_pc_shared_EC$N_acc1_S),digits = 1)
dev.off()

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/venndiag_PC_shared_2_acc_EC_F.pdf",height = 2,width = 2)
par(mar=c(10,4,4,2)) 

a<-round(mean(tab_pc_shared_EC$N_acc1_F),digits = 1)
b<-round(mean(tab_pc_shared_EC$N_acc2_F),digits = 1)
c<-round(mean(tab_pc_shared_EC$N_shared_F),digits = 1)
draw.pairwise.venn(area1=a, area2=b,cross.area=c,
                   category=c("acc1","acc2"),fill=c("#486EB4","#486EB4"),)
round(100*mean(tab_pc_shared_EC$N_shared_F)/mean(tab_pc_shared_EC$N_acc1_F),digits = 1)
dev.off()

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/venndiag_PC_shared_2_acc_EC_P.pdf",height = 2,width = 2)
par(mar=c(10,4,4,2)) 

a<-round(mean(tab_pc_shared_EC$N_acc1_P),digits = 1)
b<-round(mean(tab_pc_shared_EC$N_acc2_P),digits = 1)
c<-round(mean(tab_pc_shared_EC$N_shared_P),digits = 1)
draw.pairwise.venn(area1=a, area2=b,cross.area=c,
                   category=c("acc1","acc2"),fill=c("#486EB4","#486EB4"),)
round(100*mean(tab_pc_shared_EC$N_shared_P)/mean(tab_pc_shared_EC$N_acc1_P),digits = 1)
dev.off()






pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/venndiag_LINC_shared_2_acc_EC.pdf",height = 2,width = 2)
par(mar=c(10,4,4,2)) 

a<-round(mean(tab_linc_shared$N_acc1),digits = 1)
b<-round(mean(tab_linc_shared$N_acc2),digits = 1)
c<-round(mean(tab_linc_shared$N_shared),digits = 1)
draw.pairwise.venn(area1=a, area2=b,cross.area=c,
                   category=c("acc1","acc2"),fill=c("#F2AB54","#F2AB54"))
round(100*mean(tab_linc_shared$N_shared)/mean(tab_linc_shared$N_acc1),digits = 1)

dev.off()

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/venndiag_AS_shared_2_acc_EC.pdf",height = 2,width = 2)
par(mar=c(10,4,4,2)) 

a<-round(mean(tab_as_shared$N_acc1),digits = 1)
b<-round(mean(tab_as_shared$N_acc2),digits = 1)
c<-round(mean(tab_as_shared$N_shared),digits = 1)
draw.pairwise.venn(area1=a, area2=b,cross.area=c,
                   category=c("acc1","acc2"),fill=c("#90C473","#90C473"))
round(100*mean(tab_as_shared$N_shared)/mean(tab_as_shared$N_acc1),digits = 1)

dev.off()

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/venndiag_TE_shared_2_acc_EC.pdf",height = 2,width = 2)
par(mar=c(10,4,4,2)) 

a<-round(mean(tab_te_shared$N_acc1),digits = 1)
b<-round(mean(tab_te_shared$N_acc2),digits = 1)
c<-round(mean(tab_te_shared$N_shared),digits = 1)
draw.pairwise.venn(area1=a, area2=b,cross.area=c,
                   category=c("acc1","acc2"),fill=c("#673A8E","#673A8E"))
round(100*mean(tab_te_shared$N_shared)/mean(tab_te_shared$N_acc1),digits = 1)

dev.off()



#############################################################################
#############################################################################

# expression variation in different tissues 


#1001Gnew rosette 
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_SUPPL_boxplot_intervariance_rosette1001Gnew_all.pdf",height = 3.5,width = 3.5)
######################################################
par(mar=c(7,3,3,2)) 
a1<-pc_sasa$variance_of_means[pc_sasa$ma_x>1]
a2<-as_sasa$variance_of_means[as_sasa$ma_x>1]
a3<-linc_sasa$variance_of_means[linc_sasa$ma_x>1]
a4<-te_sasa$variance_of_means[te_sasa$ma_x>1]
a5<-as_te_sasa$variance_of_means[as_te_sasa$ma_x>1]
a6<-te_frag_sasa$variance_of_means[te_frag_sasa$ma_x>1]
boxplot(a1, a2,  a3    ,  a4      ,   a5     ,  a6      ,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Expression variability\n28 accessions, rosette")
title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
mtext('TPMmax>1', side=1, line=6, at=0,cex=0.7)

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
a<-wilcox.test(a4,a3)
b<-wilcox.test(a4,a3)
c<-wilcox.test(a4,a3)
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
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=5)
#################

dev.off()


pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_SUPPL_boxplot_intervariance_rosette1001Gnew_max1-2.pdf",height = 3.5,width = 3.5)
######################################################
par(mar=c(7,3,3,2)) 
a1<-pc_sasa$variance_of_means[pc_sasa$ma_x>1 &pc_sasa$ma_x<2]
a2<-as_sasa$variance_of_means[as_sasa$ma_x>1 &as_sasa$ma_x<2]
a3<-linc_sasa$variance_of_means[linc_sasa$ma_x>1 &linc_sasa$ma_x<2]
a4<-te_sasa$variance_of_means[te_sasa$ma_x>1 &te_sasa$ma_x<2]
a5<-as_te_sasa$variance_of_means[as_te_sasa$ma_x>1 &as_te_sasa$ma_x<2]
a6<-te_frag_sasa$variance_of_means[te_frag_sasa$ma_x>1 &te_frag_sasa$ma_x<2]
boxplot(a1, a2,  a3    ,  a4      ,   a5     ,  a6      ,
        col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Expression variability\n28 accessions, rosette")
title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
mtext('1<TPMmax<2', side=1, line=6, at=0,cex=0.7)

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
text(b,x=1.5,y=4)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=4)
a<-wilcox.test(a4,a3)
b<-wilcox.test(a4,a3)
c<-wilcox.test(a4,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=4)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=4)
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=4)
#################

dev.off()




 ### expression variation in different tissues 
  
  
#seedling eracaps  
  pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_SUPPL_boxplot_intervariance_EC_seedling_all.pdf",height = 3.5,width = 3.5)
  ######################################################
  par(mar=c(7,3,3,2)) 
  a1<-pc_EC$var.seedl[pc_EC$max.seedl>1]
  a2<-as_EC$var.seedl[as_EC$max.seedl>1]
  a3<-  linc_EC$var.seedl[linc_EC$max.seedl>1]
  a4<-te_EC$var.seedl[te_EC$max.seedl>1]
  a5<-as_te_EC$var.seedl[as_te_EC$max.seedl>1]
  a6<- te_frag_EC$var.seedl[te_frag_EC$max.seedl>1]
  boxplot(a1,a2,a3,a4,a5,a6,
          col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Expression variability\n27 accessions, seedlings")
  title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
  mtext('TPMmax>1', side=1, line=6, at=0,cex=0.7)
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
  text(b,x=1.5,y=4)
  a<-wilcox.test(a2,a3)
  b<-wilcox.test(a2,a3)
  c<-wilcox.test(a2,a3)
  d<-mean(c(a$p.value,b$p.value,c$p.value))
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=4)
  a<-wilcox.test(a4,a3)
  b<-wilcox.test(a4,a3)
  c<-wilcox.test(a4,a3)
  d<-mean(c(a$p.value,b$p.value,c$p.value))
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=4)
  a<-wilcox.test(a4,a5)
  b<-wilcox.test(a4,a5)
  c<-wilcox.test(a4,a5)
  d<-mean(c(a$p.value,b$p.value,c$p.value))
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=4.5,y=4)
  a<-wilcox.test(a6,a5)
  b<-wilcox.test(a6,a5)
  c<-wilcox.test(a6,a5)
  d<-mean(c(a$p.value,b$p.value,c$p.value))
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=5.5,y=4)
  #################
  dev.off()
  
  pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_SUPPL_boxplot_intervariance_EC_seedling_max1-2.pdf",height = 3.5,width = 3.5)
  ######################################################
  par(mar=c(7,3,3,2))
  a1<-pc_EC$var.seedl[pc_EC$max.seedl>1&pc_EC$max.seedl<2]
    a2<-as_EC$var.seedl[as_EC$max.seedl>1&as_EC$max.seedl<2]
    a3<-  linc_EC$var.seedl[linc_EC$max.seedl>1&linc_EC$max.seedl<2]
    a4<-te_EC$var.seedl[te_EC$max.seedl>1&te_EC$max.seedl<2]
    a5<-as_te_EC$var.seedl[as_te_EC$max.seedl>1&as_te_EC$max.seedl<2]
    a6<- te_frag_EC$var.seedl[te_frag_EC$max.seedl>1&te_frag_EC$max.seedl<2]
      boxplot(a1,a2,a3,a4,a5,a6,
                    col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Expression variability\n27 accessions, seedlings")
    title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
    mtext('1<TPMmax<2', side=1, line=6, at=0,cex=0.7)
     #text(x= 1, y= max(te_EC$var.seedl[te_EC$max.seedl>1&te_EC$max.seedl<2]), labels= "1<TPMmax<2",cex = 0.8)
    
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
    text(b,x=1.5,y=4)
    a<-wilcox.test(a2,a3)
    b<-wilcox.test(a2,a3)
    c<-wilcox.test(a2,a3)
    d<-mean(c(a$p.value,b$p.value,c$p.value))
    if (d<0.0000000001){b="***"} 
    if (d>=0.0000000001 &d<0.00001 ){b="**"} 
    if (d>=0.00001  &d<0.01){b="*"} 
    if (d>=0.01){b="n.s."} 
    text(b,x=2.5,y=4)
    a<-wilcox.test(a4,a3)
    b<-wilcox.test(a4,a3)
    c<-wilcox.test(a4,a3)
    d<-mean(c(a$p.value,b$p.value,c$p.value))
    if (d<0.0000000001){b="***"} 
    if (d>=0.0000000001 &d<0.00001 ){b="**"} 
    if (d>=0.00001  &d<0.01){b="*"} 
    if (d>=0.01){b="n.s."} 
    text(b,x=3.5,y=4)
    a<-wilcox.test(a4,a5)
    b<-wilcox.test(a4,a5)
    c<-wilcox.test(a4,a5)
    d<-mean(c(a$p.value,b$p.value,c$p.value))
    if (d<0.0000000001){b="***"} 
    if (d>=0.0000000001 &d<0.00001 ){b="**"} 
    if (d>=0.00001  &d<0.01){b="*"} 
    if (d>=0.01){b="n.s."} 
    text(b,x=4.5,y=4)
    a<-wilcox.test(a6,a5)
    b<-wilcox.test(a6,a5)
    c<-wilcox.test(a6,a5)
    d<-mean(c(a$p.value,b$p.value,c$p.value))
    if (d<0.0000000001){b="***"} 
    if (d>=0.0000000001 &d<0.00001 ){b="**"} 
    if (d>=0.00001  &d<0.01){b="*"} 
    if (d>=0.01){b="n.s."} 
    text(b,x=5.5,y=4)
    #################
    
dev.off()
  
#rosette eracaps
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_SUPPL_boxplot_intervariance_EC_rosette_all.pdf",height = 3.5,width = 3.5)
######################################################
      par(mar=c(7,3,3,2)) 
      a1<-     pc_EC$var.rosette[pc_EC$max.rosette>1]
      a2<-      as_EC$var.rosette[as_EC$max.rosette>1]
      a3<-    linc_EC$var.rosette[linc_EC$max.rosette>1]
      a4<-               te_EC$var.rosette[te_EC$max.rosette>1]
      a5<-               as_te_EC$var.rosette[as_te_EC$max.rosette>1]
      a6<-              te_frag_EC$var.rosette[te_frag_EC$max.rosette>1]
                
    boxplot(a1,a2,a3,a4,a5,a6,col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Expression variability\n27 accessions, rosette")
      title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
      mtext('TPMmax>1', side=1, line=6, at=0,cex=0.7)
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
      text(b,x=1.5,y=4)
      a<-wilcox.test(a2,a3)
      b<-wilcox.test(a2,a3)
      c<-wilcox.test(a2,a3)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=2.5,y=4)
      a<-wilcox.test(a4,a3)
      b<-wilcox.test(a4,a3)
      c<-wilcox.test(a4,a3)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=3.5,y=4)
      a<-wilcox.test(a4,a5)
      b<-wilcox.test(a4,a5)
      c<-wilcox.test(a4,a5)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=4.5,y=4)
      a<-wilcox.test(a6,a5)
      b<-wilcox.test(a6,a5)
      c<-wilcox.test(a6,a5)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=5.5,y=4)
#################
      dev.off()
      
      
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_SUPPL_boxplot_intervariance_EC_rosette_max1-2.pdf",height = 3.5,width = 3.5)
######################################################
      par(mar=c(7,3,3,2)) 
      a1<-         pc_EC$var.rosette[pc_EC$max.rosette>1&pc_EC$max.rosette<2]
      a2<-               as_EC$var.rosette[as_EC$max.rosette>1&as_EC$max.rosette<2]
      a3<-               linc_EC$var.rosette[linc_EC$max.rosette>1&linc_EC$max.rosette<2]
      a4<-              te_EC$var.rosette[te_EC$max.rosette>1&te_EC$max.rosette<2]
      a5<-              as_te_EC$var.rosette[as_te_EC$max.rosette>1&as_te_EC$max.rosette<2]
      a6<-              te_frag_EC$var.rosette[te_frag_EC$max.rosette>1]
               
     boxplot(a1,a2,a3,a4,a5,a6,col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Expression variability\n27 accessions, rosette")
      title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
      mtext('1<TPMmax<2', side=1, line=6, at=0,cex=0.7)
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
      text(b,x=1.5,y=4)
      a<-wilcox.test(a2,a3)
      b<-wilcox.test(a2,a3)
      c<-wilcox.test(a2,a3)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=2.5,y=4)
      a<-wilcox.test(a4,a3)
      b<-wilcox.test(a4,a3)
      c<-wilcox.test(a4,a3)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=3.5,y=4)
      a<-wilcox.test(a4,a5)
      b<-wilcox.test(a4,a5)
      c<-wilcox.test(a4,a5)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=4.5,y=4)
      a<-wilcox.test(a6,a5)
      b<-wilcox.test(a6,a5)
      c<-wilcox.test(a6,a5)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=5.5,y=4)
      #################
      dev.off()
      

#flowers eracaps  
      pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_SUPPL_boxplot_intervariance_EC_flowers_all.pdf",height = 3.5,width = 3.5)
      ######################################################
      par(mar=c(7,3,3,2)) 
      a1<-       pc_EC$var.flowers[pc_EC$max.flowers>1]
      a2<-       as_EC$var.flowers[as_EC$max.flowers>1]
      a3<-      linc_EC$var.flowers[linc_EC$max.flowers>1]
      a4<-    te_EC$var.flowers[te_EC$max.flowers>1]
      a5<-    as_te_EC$var.flowers[as_te_EC$max.flowers>1]
      a6<-    te_frag_EC$var.flowers[te_frag_EC$max.flowers>1] 
    boxplot(a1,a2,a3,a4,a5,a6,
              col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Expression variability\n27 accessions, flowers")
      title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
      mtext('TPMmax>1', side=1, line=6, at=0,cex=0.7)
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
      text(b,x=1.5,y=4)
      a<-wilcox.test(a2,a3)
      b<-wilcox.test(a2,a3)
      c<-wilcox.test(a2,a3)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=2.5,y=4)
      a<-wilcox.test(a4,a3)
      b<-wilcox.test(a4,a3)
      c<-wilcox.test(a4,a3)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=3.5,y=4)
      a<-wilcox.test(a4,a5)
      b<-wilcox.test(a4,a5)
      c<-wilcox.test(a4,a5)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=4.5,y=4)
      a<-wilcox.test(a6,a5)
      b<-wilcox.test(a6,a5)
      c<-wilcox.test(a6,a5)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=5.5,y=4)
      #################
      dev.off()
      
      pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_SUPPL_boxplot_intervariance_EC_flowers_max1-2.pdf",height = 3.5,width = 3.5)
      ######################################################
      par(mar=c(7,3,3,2)) 
      a1<-       pc_EC$var.flowers[pc_EC$max.flowers>1&pc_EC$max.flowers<2]
      a2<-               as_EC$var.flowers[as_EC$max.flowers>1&as_EC$max.flowers<2]
      a3<-                linc_EC$var.flowers[linc_EC$max.flowers>1&linc_EC$max.flowers<2]
      a4<-               te_EC$var.flowers[te_EC$max.flowers>1&te_EC$max.flowers<2]
      a5<-              as_te_EC$var.flowers[as_te_EC$max.flowers>1&as_te_EC$max.flowers<2]
      a6<-    te_frag_EC$var.flowers[te_frag_EC$max.flowers>1&te_frag_EC$max.flowers<2]
               
     boxplot(a1,a2,a3,a4,a5,a6,col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Expression variability\n27 accessions, flowers")
      title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
      mtext('1<TPMmax<2', side=1, line=6, at=0,cex=0.7)
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
      text(b,x=1.5,y=4)
      a<-wilcox.test(a2,a3)
      b<-wilcox.test(a2,a3)
      c<-wilcox.test(a2,a3)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=2.5,y=4)
      a<-wilcox.test(a4,a3)
      b<-wilcox.test(a4,a3)
      c<-wilcox.test(a4,a3)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=3.5,y=4)
      a<-wilcox.test(a4,a5)
      b<-wilcox.test(a4,a5)
      c<-wilcox.test(a4,a5)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=4.5,y=4)
      a<-wilcox.test(a6,a5)
      b<-wilcox.test(a6,a5)
      c<-wilcox.test(a6,a5)
      d<-mean(c(a$p.value,b$p.value,c$p.value))
      if (d<0.0000000001){b="***"} 
      if (d>=0.0000000001 &d<0.00001 ){b="**"} 
      if (d>=0.00001  &d<0.01){b="*"} 
      if (d>=0.01){b="n.s."} 
      text(b,x=5.5,y=4)
      #################
       dev.off()
      
#pollen eracaps  
       pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_SUPPL_boxplot_intervariance_EC_pollen_all.pdf",height = 3.5,width = 3.5)
       ######################################################
       par(mar=c(7,3,3,2)) 
       a1<-           pc_EC$var.pollen[pc_EC$max.pollen>1]
       a2<-           as_EC$var.pollen[as_EC$max.pollen>1]
       a3<-          linc_EC$var.pollen[linc_EC$max.pollen>1]
       a4<-    te_EC$var.pollen[te_EC$max.pollen>1]
       a5<-    as_te_EC$var.pollen[as_te_EC$max.pollen>1]
       a6<-    te_frag_EC$var.pollen[te_frag_EC$max.pollen>1]
               
       boxplot(a1,a2,a3,a4,a5,a6,
               col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Expression variability\n27 accessions, pollen")
       title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
       mtext('TPMmax>1', side=1, line=6, at=0,cex=0.7)
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
       text(b,x=1.5,y=4)
       a<-wilcox.test(a2,a3)
       b<-wilcox.test(a2,a3)
       c<-wilcox.test(a2,a3)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=2.5,y=4)
       a<-wilcox.test(a4,a3)
       b<-wilcox.test(a4,a3)
       c<-wilcox.test(a4,a3)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=3.5,y=4)
       a<-wilcox.test(a4,a5)
       b<-wilcox.test(a4,a5)
       c<-wilcox.test(a4,a5)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=4.5,y=4)
       a<-wilcox.test(a6,a5)
       b<-wilcox.test(a6,a5)
       c<-wilcox.test(a6,a5)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=5.5,y=4)
       #################
       
       dev.off()
       
       pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_SUPPL_boxplot_intervariance_EC_pollen_max1-2.pdf",height = 3.5,width = 3.5)
       ######################################################
       par(mar=c(7,3,3,2)) 
       a1<-      pc_EC$var.pollen[pc_EC$max.pollen>1&pc_EC$max.pollen<2]
       a2<-               as_EC$var.pollen[as_EC$max.pollen>1&as_EC$max.pollen<2]
       a3<-               linc_EC$var.pollen[linc_EC$max.pollen>1&linc_EC$max.pollen<2]
       a4<-               te_EC$var.pollen[te_EC$max.pollen>1&te_EC$max.pollen<2]
       a5<-               as_te_EC$var.pollen[as_te_EC$max.pollen>1&as_te_EC$max.pollen<2]
       a6<-              te_frag_EC$var.pollen[te_frag_EC$max.pollen>1&te_frag_EC$max.pollen<2]
               
      boxplot( a1,a2,a3,a4,a5,a6,col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Expression variability\n27 accessions, pollen")
       title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
       mtext('1<TPMmax<2', side=1, line=6, at=0,cex=0.7)
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
       text(b,x=1.5,y=4)
       a<-wilcox.test(a2,a3)
       b<-wilcox.test(a2,a3)
       c<-wilcox.test(a2,a3)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=2.5,y=4)
       a<-wilcox.test(a4,a3)
       b<-wilcox.test(a4,a3)
       c<-wilcox.test(a4,a3)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=3.5,y=4)
       a<-wilcox.test(a4,a5)
       b<-wilcox.test(a4,a5)
       c<-wilcox.test(a4,a5)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=4.5,y=4)
       a<-wilcox.test(a6,a5)
       b<-wilcox.test(a6,a5)
       c<-wilcox.test(a6,a5)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=5.5,y=4)
       #################
       dev.off()     
      
      
  
  # expression noise (Cortijo)
  
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_SUPPL_boxplot_exp_noise_Cortijo_all.pdf",height = 3.5,width = 3.5)
######################################################
       par(mar=c(7,3,3,2)) 
     a1<-  pc_cortijo$noise_average[pc_cortijo$max>1]
     a2<-         as_cortijo$noise_average[as_cortijo$max>1]
     a3<-              linc_cortijo$noise_average[linc_cortijo$max>1]
     a4<-              te_gene_cortijo$noise_average[te_gene_cortijo$max>1]
     a5<-              as_te_cortijo$noise_average[as_te_cortijo$max>1]
     a6<-              te_frag_cortijo$noise_average[te_frag_cortijo$max>1]
               
       boxplot(a1,a2,a3,a4,a5,a6,col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Expression noise\nCortijo et al, Col-0 seedlings")
       title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
       mtext('TPMmax>1', side=1, line=6, at=0,cex=0.7)
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
       text(b,x=1.5,y=4)
       a<-wilcox.test(a2,a3)
       b<-wilcox.test(a2,a3)
       c<-wilcox.test(a2,a3)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=2.5,y=4)
       a<-wilcox.test(a4,a3)
       b<-wilcox.test(a4,a3)
       c<-wilcox.test(a4,a3)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=3.5,y=4)
       a<-wilcox.test(a4,a5)
       b<-wilcox.test(a4,a5)
       c<-wilcox.test(a4,a5)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=4.5,y=4)
       a<-wilcox.test(a6,a5)
       b<-wilcox.test(a6,a5)
       c<-wilcox.test(a6,a5)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=5.5,y=4)
       #################
       dev.off()
       
       pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_SUPPL_boxplot_exp_noise_Cortijo_max1-2.pdf",height = 3.5,width = 3.5)
       ######################################################
       par(mar=c(7,3,3,2)) 
       a1<-     pc_cortijo$noise_average[pc_cortijo$max>1 &pc_cortijo$max<2]
       a2<-              as_cortijo$noise_average[as_cortijo$max>1 &as_cortijo$max<2]
       a3<-              linc_cortijo$noise_average[linc_cortijo$max>1 &linc_cortijo$max<2]
       a4<-              te_gene_cortijo$noise_average[te_gene_cortijo$max>1 &te_gene_cortijo$max<2]
       a5<-            as_te_cortijo$noise_average[as_te_cortijo$max>1 &as_te_cortijo$max<2]
       a6<-        te_frag_cortijo$noise_average[te_frag_cortijo$max>1 &te_frag_cortijo$max<3]
                
      boxplot(a1,a2,a3,a4,a5,a6,col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Expression noise\nCortijo et al, Col-0 seedlings")
       title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
       mtext('1<TPMmax<2', side=1, line=6, at=0,cex=0.7)
       mtext(paste("n=",length(a1),sep=""), side=1, line=-1, at=1,cex=0.5)
       mtext(paste("n=",length(a2),sep=""), side=1, line=-1, at=2,cex=0.5)
       mtext(paste("n=",length(a3),sep=""), side=1, line=-1, at=3,cex=0.5)
       mtext(paste("n=",length(a4),sep=""), side=1, line=-1, at=4,cex=0.5)
       mtext(paste("n=",length(a5),sep=""), side=1, line=-1, at=5,cex=0.5)
       mtext(paste("n=",length(a6),sep=""), side=1, line=-1, at=6,cex=0.5)
       
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
       text(b,x=1.5,y=4)
       a<-wilcox.test(a2,a3)
       b<-wilcox.test(a2,a3)
       c<-wilcox.test(a2,a3)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=2.5,y=4)
       a<-wilcox.test(a4,a3)
       b<-wilcox.test(a4,a3)
       c<-wilcox.test(a4,a3)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=3.5,y=4)
       a<-wilcox.test(a4,a5)
       b<-wilcox.test(a4,a5)
       c<-wilcox.test(a4,a5)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=4.5,y=4)
       a<-wilcox.test(a6,a5)
       b<-wilcox.test(a6,a5)
       c<-wilcox.test(a6,a5)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=5.5,y=4)
       #################
       dev.off()    
  
       
       
       pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_SUPPL_boxplot_exp_noise_Cortijo_max1-2_length800-1500.pdf",height = 3.5,width = 3.5)
       ######################################################
       par(mar=c(7,3,3,2)) 
       a1<-   pc_cortijo$noise_average[pc_cortijo$max>1 &pc_cortijo$max<2 &pc_cortijo$gene %in% denovoPC.loci$gene[denovoPC.loci$length<1500 &denovoPC.loci$length>800 ]]
a2<-              as_cortijo$noise_average[as_cortijo$max>1 &as_cortijo$max<2&as_cortijo$gene %in% lncRNAs.antisense.loci$gene[lncRNAs.antisense.loci$length<1500 &lncRNAs.antisense.loci$length>800]]
a3<-     linc_cortijo$noise_average[linc_cortijo$max>1 &linc_cortijo$max<2&linc_cortijo$gene %in% lncRNAs.intergenic.loci$gene[lncRNAs.intergenic.loci$length<1500 &lncRNAs.intergenic.loci$length>800]]
a4<-        te_gene_cortijo$noise_average[te_gene_cortijo$max>1 &te_gene_cortijo$max<2&te_gene_cortijo$gene %in% TE_genes.loci$gene[TE_genes.loci$length<1500 &TE_genes.loci$length>800]]
a5<-      as_te_cortijo$noise_average[as_te_cortijo$max>1 &as_te_cortijo$max<2&as_te_cortijo$gene %in% lncRNAs.AS_to_TE.loci$gene[lncRNAs.AS_to_TE.loci$length<1500 &lncRNAs.AS_to_TE.loci$length>800]]
a6<-     te_frag_cortijo$noise_average[te_frag_cortijo$max>1 &te_frag_cortijo$max<2&te_frag_cortijo$gene %in% TE_frags.transcripts$gene[TE_frags.transcripts$length<1500 &TE_frags.transcripts$length>800]]

               boxplot( a1,a2,a3,a4,a5,a6, col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE frags"),las=2, notch = T, outline = F, main="Expression noise\nCortijo et al, Col-0 seedlings")
       title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
       mtext('1<TPMmax<2 + length 0.8-1.5kb', side=1, line=6, at=0,cex=0.7)
       mtext(paste("n=",length(a1),sep=""), side=1, line=-2, at=1,cex=0.7)
       mtext(paste("n=",length(a2),sep=""), side=1, line=-2, at=2,cex=0.7)
       mtext(paste("n=",length(a3),sep=""), side=1, line=-2, at=3,cex=0.7)
       mtext(paste("n=",length(a4),sep=""), side=1, line=-2, at=4,cex=0.7)
       mtext(paste("n=",length(a5),sep=""), side=1, line=-2, at=5,cex=0.7)
       mtext(paste("n=",length(a6),sep=""), side=1, line=-2, at=6,cex=0.7)
       
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
       a<-wilcox.test(a4,a3)
       b<-wilcox.test(a4,a3)
       c<-wilcox.test(a4,a3)
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
       a<-wilcox.test(a6,a5)
       b<-wilcox.test(a6,a5)
       c<-wilcox.test(a6,a5)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=5.5,y=2)
       #################
       dev.off()    
       
       
      
       
       
       
       
# noise 
       pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2_SUPPL_boxplot_exp_noise_Cortijo_max1-2.pdf",height = 3.5,width = 3.5)
######################################################
       par(mar=c(7,3,3,2)) 
       
       boxplot(pc_cortijo$noise_average[pc_cortijo$max>1 &pc_cortijo$max<2],
               as_cortijo$noise_average[as_cortijo$max>1 &as_cortijo$max<2],
               linc_cortijo$noise_average[linc_cortijo$max>1 &linc_cortijo$max<2&  linc_cortijo$gene %in% linc_TE_cov_all_loci_2cols$gene[linc_TE_cov_all_loci_2cols$coverage==0] ],
               linc_cortijo$noise_average[linc_cortijo$max>1&linc_cortijo$max<2  & linc_cortijo$gene %in% linc_TE_cov_all_loci_2cols$gene[linc_TE_cov_all_loci_2cols$coverage>0
                                                                                                                                          &linc_cortijo$max<2&linc_TE_cov_all_loci_2cols$coverage<=0.5] ],
               linc_cortijo$noise_average[linc_cortijo$max>1  & linc_cortijo$gene %in% linc_TE_cov_all_loci_2cols$gene[linc_TE_cov_all_loci_2cols$coverage>0.5
                                                                                                                       &linc_TE_cov_all_loci_2cols$coverage<=0.8] ],
               linc_cortijo$noise_average[linc_cortijo$max>1 &linc_cortijo$max<2 & linc_cortijo$gene %in% linc_TE_cov_all_loci_2cols$gene[linc_TE_cov_all_loci_2cols$coverage>0.8]],
               te_gene_cortijo$noise_average[te_gene_cortijo$max>1 &te_gene_cortijo$max<2],
               as_te_cortijo$noise_average[as_te_cortijo$max>1 &as_te_cortijo$max<2],
               te_frag_cortijo$noise_average[te_frag_cortijo$max>1 &te_frag_cortijo$max<2],
               col=c("#486EB4","#90C473","#F2AB54","#F2AB54","#F2AB54","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","lincRNAs","lincRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Expression noise\nCortijo et al, Col-0 seedlings")
       
       title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
       mtext('1<TPMmax<2', side=1, line=6, at=0,cex=0.7)
       ######################################################
       dev.off()    
       
       
       
       
       
       
 
 # circadian variance Cortijo 
       
       
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2/Fig2_SUPPL_boxplot_circadian_variance_Cortijo_tpmmax1.pdf",height = 3.5,width = 3.5)
######################################################
       par(mar=c(7,3,3,2)) 
     a1<- pc_cortijo$variance_circadian[pc_cortijo$max>1]
     a2<-         as_cortijo$variance_circadian[as_cortijo$max>1]
     a3<-           linc_cortijo$variance_circadian[linc_cortijo$max>1]
     a4<-           te_gene_cortijo$variance_circadian[te_gene_cortijo$max>1]
     a5<-             as_te_cortijo$variance_circadian[as_te_cortijo$max>1]
     a6<-              te_frag_cortijo$variance_circadian[te_frag_cortijo$max>1]
      
      boxplot(a1,a2,a3,a4,a5,a6,
               col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Circadian expression variance\nCortijo et al, Col-0 seedlings")
       title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
       mtext('TPMmax>1', side=1, line=6, at=0,cex=0.7)
       mtext(paste("n=",length(a1),sep=""), side=1, line=-1, at=1,cex=0.5)
       mtext(paste("n=",length(a2),sep=""), side=1, line=-1, at=2,cex=0.5)
       mtext(paste("n=",length(a3),sep=""), side=1, line=-1, at=3,cex=0.5)
       mtext(paste("n=",length(a4),sep=""), side=1, line=-1, at=4,cex=0.5)
       mtext(paste("n=",length(a5),sep=""), side=1, line=-1, at=5,cex=0.5)
       mtext(paste("n=",length(a6),sep=""), side=1, line=-1, at=6,cex=0.5)
       
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
       text(b,x=1.5,y=1.5)
       a<-wilcox.test(a2,a3)
       b<-wilcox.test(a2,a3)
       c<-wilcox.test(a2,a3)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=2.5,y=1.5)
       a<-wilcox.test(a4,a3)
       b<-wilcox.test(a4,a3)
       c<-wilcox.test(a4,a3)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=3.5,y=1.5)
       a<-wilcox.test(a4,a5)
       b<-wilcox.test(a4,a5)
       c<-wilcox.test(a4,a5)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=4.5,y=1.5)
       a<-wilcox.test(a6,a5)
       b<-wilcox.test(a6,a5)
       c<-wilcox.test(a6,a5)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=5.5,y=1.5)
       #################      
        dev.off()
       
       
       
       pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2_SUPPL_boxplot_circadian_variance_Cortijo_max1-2.pdf",height = 3.5,width = 3.5)
       par(mar=c(7,3,3,2)) 
 a1<-     pc_cortijo$variance_circadian[pc_cortijo$max>1 &pc_cortijo$max<2]
 a2<-          as_cortijo$variance_circadian[as_cortijo$max>1 &as_cortijo$max<2]
 a3<-            linc_cortijo$variance_circadian[linc_cortijo$max>1 &linc_cortijo$max<2]
 a4<-             te_gene_cortijo$variance_circadian[te_gene_cortijo$max>1 &te_gene_cortijo$max<2]
 a5<-             as_te_cortijo$variance_circadian[as_te_cortijo$max>1 &as_te_cortijo$max<2]
 a6<-             te_frag_cortijo$variance_circadian[te_frag_cortijo$max>1 &te_frag_cortijo$max<2]
 
               boxplot( a1,a2,a3,a4,a5,a6,col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Circadian expression variance\nCortijo et al, Col-0 seedlings")
       title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
       mtext('1<TPMmax<2', side=1, line=6, at=0,cex=0.7)
       mtext(paste("n=",length(a1),sep=""), side=1, line=-1, at=1,cex=0.5)
       mtext(paste("n=",length(a2),sep=""), side=1, line=-1, at=2,cex=0.5)
       mtext(paste("n=",length(a3),sep=""), side=1, line=-1, at=3,cex=0.5)
       mtext(paste("n=",length(a4),sep=""), side=1, line=-1, at=4,cex=0.5)
       mtext(paste("n=",length(a5),sep=""), side=1, line=-1, at=5,cex=0.5)
       mtext(paste("n=",length(a6),sep=""), side=1, line=-1, at=6,cex=0.5)
       
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
       text(b,x=1.5,y=1.5)
       a<-wilcox.test(a2,a3)
       b<-wilcox.test(a2,a3)
       c<-wilcox.test(a2,a3)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=2.5,y=1.5)
       a<-wilcox.test(a4,a3)
       b<-wilcox.test(a4,a3)
       c<-wilcox.test(a4,a3)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=3.5,y=1.5)
       a<-wilcox.test(a4,a5)
       b<-wilcox.test(a4,a5)
       c<-wilcox.test(a4,a5)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=4.5,y=1.5)
       a<-wilcox.test(a6,a5)
       b<-wilcox.test(a6,a5)
       c<-wilcox.test(a6,a5)
       d<-mean(c(a$p.value,b$p.value,c$p.value))
       if (d<0.0000000001){b="***"} 
       if (d>=0.0000000001 &d<0.00001 ){b="**"} 
       if (d>=0.00001  &d<0.01){b="*"} 
       if (d>=0.01){b="n.s."} 
       text(b,x=5.5,y=1.5)
       #################      
       
       dev.off()      
 
 
 
       
# intraindividual variance sasa 
 
       par(mar=c(7,3,3,2)) 
       boxplot(pc_sasa$mean_intravariance[pc_sasa$ma_x>1 &pc_sasa$ma_x<2],
               as_sasa$mean_intravariance[as_sasa$ma_x>1 &as_sasa$ma_x<2],
               linc_sasa$mean_intravariance[linc_sasa$ma_x>1 &linc_sasa$ma_x<2],
               te_sasa$mean_intravariance[te_sasa$ma_x>1 &te_sasa$ma_x<2],
               as_te_sasa$mean_intravariance[as_te_sasa$ma_x>1 &as_te_sasa$ma_x<2],
               te_frag_sasa$mean_intravariance[te_frag_sasa$ma_x>1 &te_frag_sasa$ma_x<2],
               col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Inter-replicate (same-accession) expression variance\n rosette")
       title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
       mtext('1<TPMmax<2', side=1, line=6, at=0,cex=0.7)
       
       
       
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2_SUPPL_boxplot_1001Gnew_inter_vs_intra_variability.pdf",height = 3.5,width = 8)
       par(mar=c(7,3,3,2)) 
       boxplot(pc_sasa$variance_of_means[pc_sasa$ma_x>1],
           pc_sasa$var_3random_mean[pc_sasa$ma_x>1],     pc_sasa$mean_intravariance[pc_sasa$ma_x>1],
              
         as_sasa$variance_of_means[as_sasa$ma_x>1],
          as_sasa$var_3random_mean[as_sasa$ma_x>1],   as_sasa$mean_intravariance[as_sasa$ma_x>1],

                  linc_sasa$variance_of_means[linc_sasa$ma_x>1],
        linc_sasa$var_3random_mean[linc_sasa$ma_x>1], linc_sasa$mean_intravariance[linc_sasa$ma_x>1],

                 te_sasa$variance_of_means[te_sasa$ma_x>1],
         te_sasa$var_3random_mean[te_sasa$ma_x>1], te_sasa$mean_intravariance[te_sasa$ma_x>1],

as_te_sasa$variance_of_means[as_te_sasa$ma_x>1],
as_te_sasa$var_3random_mean[as_te_sasa$ma_x>1],as_te_sasa$mean_intravariance[as_te_sasa$ma_x>1],

         te_frag_sasa$variance_of_means[te_frag_sasa$ma_x>1],
te_frag_sasa$var_3random_mean[te_frag_sasa$ma_x>1],te_frag_sasa$mean_intravariance[te_frag_sasa$ma_x>1],


col=c("#486EB4","#486EB4","#486EB4","#90C473","#90C473","#90C473","#F2AB54","#F2AB54","#F2AB54","#673A8E","#673A8E","#673A8E","#B294C5","#B294C5","#B294C5","#805FA5","#805FA5","#805FA5"),las=2, notch = T, outline = F, main="Inter- vs intra-individual expression variability\n rosette")
       title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
       mtext('TPMmax>1', side=1, line=6, at=0,cex=0.7)
       , names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"))
       
       dev.off()
       
       
       pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2_SUPPL_boxplot_circadian_variance_Cortijo_max1-2.pdf",height = 3.5,width = 3.5)
       par(mar=c(7,3,3,2)) 
       boxplot(pc_cortijo$variance_circadian[pc_cortijo$max>1 &pc_cortijo$max<2],
               as_cortijo$variance_circadian[as_cortijo$max>1 &as_cortijo$max<2],
               linc_cortijo$variance_circadian[linc_cortijo$max>1 &linc_cortijo$max<2],
               te_gene_cortijo$variance_circadian[te_gene_cortijo$max>1 &te_gene_cortijo$max<2],
               as_te_cortijo$variance_circadian[as_te_cortijo$max>1 &as_te_cortijo$max<2],
               te_frag_cortijo$variance_circadian[te_frag_cortijo$max>1 &te_frag_cortijo$max<2],
               col=c("#486EB4","#90C473","#F2AB54","#673A8E","#B294C5","#805FA5"), names=c("PC genes","AS lncRNAs","lincRNAs","TE genes", "AS_to_TE","TE fragments"),las=2, notch = T, outline = F, main="Circadian expression variance\nCortijo et al, Col-0 seedlings")
       title(ylab="coefficient of variance", mgp=c(2,1,0), cex.lab=1)
       mtext('1<TPMmax<2', side=1, line=6, at=0,cex=0.7)
       
       dev.off()      
       
       
 
 
       
       
       
       
       
       
       
       
       
       
       
       
       
 
 
 
##########################################################
    ### FIGURE 2 : pheatmap to illustrate striking variability 
##########################################################
       
library(pheatmap) 
       
   
       
     length(denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene_type=="linc" & apply(denovo2021.TPMs.genes.ERACAPS[,grep("6909|1741|22005",names(denovo2021.TPMs.genes.ERACAPS))],1,max),1])
       #1621
       length(denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene_type=="linc" & apply(denovo2021.TPMs.genes.ERACAPS[,grep("6909|9888|9905",names(denovo2021.TPMs.genes.ERACAPS))],1,max),1])
       #1483
       length(denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene_type=="linc" & apply(denovo2021.TPMs.genes.ERACAPS[,grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))],1,max),1])
       #1515
       a<-denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene_type=="pc" & apply(denovo2021.TPMs.genes.ERACAPS[,grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))],1,max),grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))]
       pheatmap (a[sample(1:length(a$R.6244),1515),],scale = "row",labels_row = F)
       
       a<-denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene_type=="te_gene" & apply(denovo2021.TPMs.genes.ERACAPS[,grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))],1,max),grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))]
       pheatmap (a,scale = "row",labels_row = F)
       
       a<-denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene_type=="as" & apply(denovo2021.TPMs.genes.ERACAPS[,grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))],1,max),grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))]
       pheatmap (a[sample(1:length(a$R.6244),1515),],scale = "row",labels_row = F)
       
       
        pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2_heatmap_3_accessions_4tissues_linc.pdf",height = 3,width =2.5)
       par(mar=c(3,3,3,2)) 
       pheatmap (denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene_type=="linc" & apply(denovo2021.TPMs.genes.ERACAPS[,grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))],1,max),grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))],scale = "row",labels_row = F,treeheight_row = 1)
              dev.off()      
       
              pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2_heatmap_3_accessions_4tissues_PC.pdf",height = 3,width =2.5)
              par(mar=c(3,3,3,2)) 
              a<-denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene_type=="pc" & apply(denovo2021.TPMs.genes.ERACAPS[,grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))],1,max),grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))]
              pheatmap (a[sample(1:length(a$R.6244),1515),],scale = "row",labels_row = F,treeheight_row = 1)
              dev.off()      
              
              
              pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2_heatmap_3_accessions_4tissues_AS.pdf",height = 3,width =2.5)
              par(mar=c(3,3,3,2)) 
              a<-denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene_type=="as" & apply(denovo2021.TPMs.genes.ERACAPS[,grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))],1,max),grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))]
              pheatmap (a[sample(1:length(a$R.6244),1515),],scale = "row",labels_row = F,treeheight_row = 1)
              dev.off()      
              
              
              pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig2_heatmap_3_accessions_4tissues_TE.pdf",height = 3,width =2.5)
              par(mar=c(3,3,3,2)) 
              a<-denovo2021.TPMs.genes.ERACAPS[denovo2021.TPMs.genes.ERACAPS$gene_type=="te_gene" & apply(denovo2021.TPMs.genes.ERACAPS[,grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))],1,max),grep("6909|6244|9888",names(denovo2021.TPMs.genes.ERACAPS))]
              pheatmap (a,scale = "row",labels_row = F,treeheight_row = 1)
              dev.off()      
              
              
              
##########################################################
##########################################################
  
              
              
              
                        
       
              
              
       
       

