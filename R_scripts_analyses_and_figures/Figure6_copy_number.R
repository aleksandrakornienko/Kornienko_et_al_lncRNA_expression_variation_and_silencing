####################################################
#Kornienko et. al 
#analyses for Figure 5 and supplements

###################################################
#copy number in TAIR10

#copy number distribution 

length(CN_pc_27genomes$gene[CN_pc_27genomes$TAIR10==1]) #22080
length(CN_pc_27genomes$gene[CN_pc_27genomes$TAIR10==2]) #1061
length(CN_pc_27genomes$gene[CN_pc_27genomes$TAIR10>2&CN_pc_27genomes$TAIR10<=10])#482
length(CN_pc_27genomes$gene[CN_pc_27genomes$TAIR10>10])#53

length(CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10==1])#1625
length(CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10==2])#212
length(CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10>2&CN_linc_27genomes$TAIR10<=10])#229
length(CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10>10])#180

length(CN_as_27genomes$gene[CN_as_27genomes$TAIR10==1])#7515
length(CN_as_27genomes$gene[CN_as_27genomes$TAIR10==2])#468
length(CN_as_27genomes$gene[CN_as_27genomes$TAIR10>2&CN_as_27genomes$TAIR10<=10])#198
length(CN_as_27genomes$gene[CN_as_27genomes$TAIR10>10])#14

length(CN_te_27genomes$gene[CN_te_27genomes$TAIR10==1])#668
length(CN_te_27genomes$gene[CN_te_27genomes$TAIR10==2])#288
length(CN_te_27genomes$gene[CN_te_27genomes$TAIR10>2&CN_te_27genomes$TAIR10<=10])#575
length(CN_te_27genomes$gene[CN_te_27genomes$TAIR10>10])#599


#CN and TEcontent 

length(CN_pc_27genomes$gene[CN_pc_27genomes$TAIR10==1& !(CN_pc_27genomes$gene %in% pc_TE_cov_loci$gene) ]) #17460 
length(CN_pc_27genomes$gene[CN_pc_27genomes$TAIR10==1& CN_pc_27genomes$gene %in% pc_TE_cov_loci$gene ]) #4620
length(CN_pc_27genomes$gene[CN_pc_27genomes$TAIR10>1 & !(CN_pc_27genomes$gene  %in% pc_TE_cov_loci$gene) ])#1004
length(CN_pc_27genomes$gene[CN_pc_27genomes$TAIR10>1 & CN_pc_27genomes$gene  %in% pc_TE_cov_loci$gene ])#592

length(CN_as_27genomes$gene[CN_as_27genomes$TAIR10==1& !(CN_as_27genomes$gene %in% as_TE_cov_loci$gene) ]) #6344 
length(CN_as_27genomes$gene[CN_as_27genomes$TAIR10==1& CN_as_27genomes$gene %in% as_TE_cov_loci$gene ]) #1171
length(CN_as_27genomes$gene[CN_as_27genomes$TAIR10>1 & !(CN_as_27genomes$gene  %in% as_TE_cov_loci$gene) ])#494
length(CN_as_27genomes$gene[CN_as_27genomes$TAIR10>1 & CN_as_27genomes$gene  %in% as_TE_cov_loci$gene ])#186

length(CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10==1& !(CN_linc_27genomes$gene %in% linc_TE_cov_loci$gene) ]) #919 
length(CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10==1& CN_linc_27genomes$gene %in% linc_TE_cov_loci$gene ]) #706
length(CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10>1 & !(CN_linc_27genomes$gene  %in% linc_TE_cov_loci$gene) ])#151
length(CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10>1 & CN_linc_27genomes$gene  %in% linc_TE_cov_loci$gene ])#470


length(CN_te_27genomes$gene[CN_te_27genomes$TAIR10==1& !(CN_te_27genomes$gene %in% te_TE_cov_loci$gene) ]) #16 
length(CN_te_27genomes$gene[CN_te_27genomes$TAIR10==1& CN_te_27genomes$gene %in% te_TE_cov_loci$gene ]) #652
length(CN_te_27genomes$gene[CN_te_27genomes$TAIR10>1 & !(CN_te_27genomes$gene  %in% te_TE_cov_loci$gene) ])#15
length(CN_te_27genomes$gene[CN_te_27genomes$TAIR10>1 & CN_te_27genomes$gene  %in% te_TE_cov_loci$gene ])#1447



#expressed silent

lincs_1copy<-CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10==1]#1625
lincs_2copies<-CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10==2]#212
lincs3_10_copies<-CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10>=3 &CN_linc_27genomes$TAIR10<=10]#229
lincs_more10copies<-CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10>10]#180

lincs_singlecopy_noTE<-CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10==1& !(CN_linc_27genomes$gene %in% linc_TE_cov_loci$gene)  ] 
lincs_singlecopy_withTE<-CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10==1& CN_linc_27genomes$gene %in% linc_TE_cov_loci$gene ] 
lincs_multicopy_noTE<-CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10>1 & !(CN_linc_27genomes$gene  %in% linc_TE_cov_loci$gene) ]
lincs_multicopy_withTE<-CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10>1 & CN_linc_27genomes$gene  %in% linc_TE_cov_loci$gene ]


tegene_singlecopy<-CN_te_27genomes$gene[CN_te_27genomes$TAIR10==1] #668
tegene_multicopy<-CN_te_27genomes$gene[CN_te_27genomes$TAIR10>1 ]# 1462


tegene_singlecopy<-CN_te_27genomes$gene[CN_te_27genomes$TAIR10==1] #668
tegene_2copies<-CN_te_27genomes$gene[CN_te_27genomes$TAIR10==2 ]# 288
tegene3_10_copies<-CN_te_27genomes$gene[CN_te_27genomes$TAIR10>=3 &CN_te_27genomes$TAIR10<=10]# 575
tegene_more10copies<-CN_te_27genomes$gene[CN_te_27genomes$TAIR10>10 ]# 599


#1001G - Col0

length(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lincs_singlecopy_noTE])#60
length(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<=0.5 & denovo2021.TPMs.genes.1001G$gene %in% lincs_singlecopy_noTE])#859
length(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lincs_singlecopy_withTE])#40
length(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<=0.5 & denovo2021.TPMs.genes.1001G$gene %in% lincs_singlecopy_withTE])#666
length(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lincs_multicopy_noTE])#2
length(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<=0.5 & denovo2021.TPMs.genes.1001G$gene %in% lincs_multicopy_noTE])#149
length(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% lincs_multicopy_withTE])#10
length(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<=0.5 & denovo2021.TPMs.genes.1001G$gene %in% lincs_multicopy_withTE])#460
length(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% tegene_singlecopy])#63
length(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<=0.5 & denovo2021.TPMs.genes.1001G$gene %in% tegene_singlecopy])#605
length(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909>0.5 & denovo2021.TPMs.genes.1001G$gene %in% tegene_multicopy])#11
length(denovo2021.TPMs.genes.1001G$gene[denovo2021.TPMs.genes.1001G$X6909<=0.5 & denovo2021.TPMs.genes.1001G$gene %in% tegene_multicopy])#1451


#1001Gnew rosette col-0
length(denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lincs_singlecopy_noTE])#42
length(denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<=0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lincs_singlecopy_noTE])#877
length(denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lincs_singlecopy_withTE])#35
length(denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<=0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lincs_singlecopy_withTE])#671

length(denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lincs_multicopy_noTE])#1
length(denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<=0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lincs_multicopy_noTE])#150
length(denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lincs_multicopy_withTE])#11
length(denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<=0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% lincs_multicopy_withTE])#459

length(denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% tegene_singlecopy])#53
length(denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<=0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% tegene_singlecopy])#615
length(denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909>0.5 & denovo2021.TPMs.genes.1001Gnew$gene %in% tegene_multicopy])#6
length(denovo2021.TPMs.genes.1001Gnew$gene[denovo2021.TPMs.genes.1001Gnew$mean.6909<=0.5 & denovo2021.TPMs.genes.1001G$gene %in% tegene_multicopy])#1456



#eracaps pollen col-0
length(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% lincs_singlecopy_noTE])#88
length(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$P.6909<=0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% lincs_singlecopy_noTE])#831
length(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% lincs_singlecopy_withTE])#103
length(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$P.6909<=0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% lincs_singlecopy_withTE])#603

length(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% lincs_multicopy_noTE])#3
length(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$P.6909<=0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% lincs_multicopy_noTE])#148
length(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% lincs_multicopy_withTE])#44
length(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$P.6909<=0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% lincs_multicopy_withTE])#426

length(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% tegene_singlecopy])#71
length(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$P.6909<=0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% tegene_singlecopy])#597
length(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$P.6909>0.5 & denovo2021.TPMs.genes.ERACAPS$gene %in% tegene_multicopy])#  
length(denovo2021.TPMs.genes.ERACAPS$gene[denovo2021.TPMs.genes.ERACAPS$P.6909<=0.5 & denovo2021.TPMs.genes.1001G$gene %in% tegene_multicopy])#1416




# TE content of singlecopy and multicopy lincRNAs 
#1001G rosette expression
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Supplem_Fig5CN_boxplot_linc_TE_content_anystrand_vs_CN.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-linc_TE_cov_all_loci_2cols$coverage[ linc_TE_cov_all_loci_2cols$gene %in% lincs_1copy]
a2<-linc_TE_cov_all_loci_2cols$coverage[ linc_TE_cov_all_loci_2cols$gene %in% lincs_2copies]
a3<-linc_TE_cov_all_loci_2cols$coverage[linc_TE_cov_all_loci_2cols$gene %in% lincs3_10_copies]
a4<-linc_TE_cov_all_loci_2cols$coverage[ linc_TE_cov_all_loci_2cols$gene %in% lincs_more10copies]

boxplot(        a1,a2,a3,a4, main="TE content",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E"), names=c("CN=1","CN=2","CN:3-10","CN>10"),las=2, notch = T, outline = F, ylab="% of loci covered by TE patches")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)
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
text(b,x=1.5,y=1.05)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1.05)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1.05)
#################
dev.off()

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Supplem_Fig5CN_boxplot_linc_TE_content_Same_strand_vs_CN.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-lincRNAs_TE_coverage_forward$TAIR10[ lincRNAs_TE_coverage_forward$gene %in% lincs_1copy]
a2<-lincRNAs_TE_coverage_forward$TAIR10[ lincRNAs_TE_coverage_forward$gene %in% lincs_2copies]
a3<-lincRNAs_TE_coverage_forward$TAIR10[lincRNAs_TE_coverage_forward$gene %in% lincs3_10_copies]
a4<-lincRNAs_TE_coverage_forward$TAIR10[ lincRNAs_TE_coverage_forward$gene %in% lincs_more10copies]

boxplot(        a1,a2,a3,a4, main="TE content - TE-linc parallel ",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E"), names=c("CN=1","CN=2","CN:3-10","CN>10"),las=2, notch = T, outline = F, ylab="% of loci covered by TE patches")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)
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
#################
dev.off()

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Supplem_Fig5CN_boxplot_linc_TE_content_opposite_strand_vs_CN.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-lincRNAs_TE_coverage_reverse$TAIR10[ lincRNAs_TE_coverage_reverse$gene %in% lincs_1copy]
a2<-lincRNAs_TE_coverage_reverse$TAIR10[ lincRNAs_TE_coverage_reverse$gene %in% lincs_2copies]
a3<-lincRNAs_TE_coverage_reverse$TAIR10[lincRNAs_TE_coverage_reverse$gene %in% lincs3_10_copies]
a4<-lincRNAs_TE_coverage_reverse$TAIR10[ lincRNAs_TE_coverage_reverse$gene %in% lincs_more10copies]

boxplot(        a1,a2,a3,a4, main="TE content - TE-linc antiparallel ",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E"), names=c("CN=1","CN=2","CN:3-10","CN>10"),las=2, notch = T, outline = F, ylab="% of loci covered by TE patches")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)
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
#################
dev.off()


#how many multicopy lincRNAs without TE pieces? 

length(CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10>1 & CN_linc_27genomes$gene %in% linc_TE_cov_all_loci_2cols$gene [linc_TE_cov_all_loci_2cols$coverage==0]])
#151
length(CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10>2& CN_linc_27genomes$gene %in% linc_TE_cov_all_loci_2cols$gene [linc_TE_cov_all_loci_2cols$coverage==0]])
#68

# what is the family of TE pieces for multicopy lincRNAs? 

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Supplem_Fig5CN_boxplot_CN_TEpiece_family_TEcont_0.2.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))

a1<-CN_linc_27genomes$TAIR10[CN_linc_27genomes$gene %in% linc_Helitron_only & linc_TE_cov_all_loci_2cols$coverage<0.2]
a2<-CN_linc_27genomes$TAIR10[CN_linc_27genomes$gene %in% linc_DNA_mudr_only & linc_TE_cov_all_loci_2cols$coverage<0.2]
a3<-CN_linc_27genomes$TAIR10[CN_linc_27genomes$gene %in% linc_LTR_copia_only & linc_TE_cov_all_loci_2cols$coverage<0.2]
a4<-CN_linc_27genomes$TAIR10[CN_linc_27genomes$gene %in% linc_LTR_gypsy_only & linc_TE_cov_all_loci_2cols$coverage<0.2]

boxplot(        a1,a2,a3,a4, main="Copy number",cex.main=1.2,     col=c("#F2AB54"), names=c("Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy"),las=2, notch = T, outline = F, ylab="lincRNA copy number in TAIR10")
mtext (text = "TE-piece-family within lincRNA",side = 1,at = 2.5,line = 3)
mtext (text = "only this family, TE content<20%",side = 1,at = 2.5,line = 4)
mtext (text = paste("N=",length(a1)),side = 1,at = 1,line = -1,cex=0.8)
mtext (text = paste("N=",length(a2)),side = 1,at = 2,line = -1,cex=0.8)
mtext (text = paste("N=",length(a3)),side = 1,at = 3,line = -1,cex=0.8)
mtext (text = paste("N=",length(a4)),side = 1,at = 4,line = -1,cex=0.8)

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
text(b,x=1.5,y=10)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=10)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=25)
#################
dev.off()

pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Supplem_Fig5CN_boxplot_CN_TEpiece_family_TEcont_any.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-CN_linc_27genomes$TAIR10[CN_linc_27genomes$gene %in% linc_Helitron_only ]
a2<-CN_linc_27genomes$TAIR10[CN_linc_27genomes$gene %in% linc_DNA_mudr_only ]
a3<-CN_linc_27genomes$TAIR10[CN_linc_27genomes$gene %in% linc_LTR_copia_only ]
a4<-CN_linc_27genomes$TAIR10[CN_linc_27genomes$gene %in% linc_LTR_gypsy_only]


boxplot(        a1,a2,a3,a4, main="Copy number",cex.main=1.2,     col=c("#F2AB54"), names=c("Helitron","DNA_MuDR","LTR_Copia","LTR_Gypsy"),las=2, notch = T, outline = F, ylab="lincRNA copy number in TAIR10")
mtext (text = "TE-piece-family within lincRNA",side = 1,at = 2.5,line = 3)
mtext (text = "only this family, TE content any",side = 1,at = 2.5,line = 4)
mtext (text = paste("N=",length(a1)),side = 1,at = 1,line = -1,cex=0.8)
mtext (text = paste("N=",length(a2)),side = 1,at = 2,line = -1,cex=0.8)
mtext (text = paste("N=",length(a3)),side = 1,at = 3,line = -1,cex=0.8)
mtext (text = paste("N=",length(a4)),side = 1,at = 4,line = -1,cex=0.8)

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
text(b,x=1.5,y=10)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=10)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=25)
#################
dev.off()



# do multicopy lincRNAs have LTRs on their borders? 
length(CN_linc_27genomes$gene[CN_linc_27genomes$TAIR10>1 & CN_linc_27genomes$gene %in% linc_TSS_Oct2021_genetic_var_stats])










# do lincRNAs with with different copy numbers have different expression? 

#1001G rosette expression
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Supplem_Fig5CN_boxplot_linc_CN_1001G.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-denovo2021.TPMs.genes.1001G$X6909[ denovo2021.TPMs.genes.1001G$gene %in% lincs_1copy]
a2<-denovo2021.TPMs.genes.1001G$X6909[ denovo2021.TPMs.genes.1001G$gene %in% lincs_2copies]
a3<-denovo2021.TPMs.genes.1001G$X6909[denovo2021.TPMs.genes.1001G$gene %in% lincs3_10_copies]
a4<-denovo2021.TPMs.genes.1001G$X6909[ denovo2021.TPMs.genes.1001G$gene %in% lincs_more10copies]

a5<-denovo2021.TPMs.genes.1001G$X6909[chip.denovo.log2$gene %in% tegene_singlecopy]
a6<-denovo2021.TPMs.genes.1001G$X6909[chip.denovo.log2$gene %in% tegene_multicopy]
boxplot(        a1,a2,a3,a4,a5,a6, main="Expression level\n(Col-0 rosette)",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E"), names=c("CN=1","CN=2","CN:3-10","CN>10","TE genes:CN=1","TE genes:CN>1"),las=2, notch = T, outline = F, ylab="expression, TPM")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 2)
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

#H3K9me2 CN 1,>1, with TE without TE
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Fig5CN_boxplot_linc_TEgene_CN-TE_H3K9me2_Col0.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% lincs_singlecopy_noTE]
a2<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% lincs_singlecopy_withTE]
a3<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% lincs_multicopy_noTE]
a4<-chip.denovo.log2$K9.6909[ chip.denovo.log2$gene %in% lincs_multicopy_withTE]
a5<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% tegene_singlecopy]
a6<-chip.denovo.log2$K9.6909[chip.denovo.log2$gene %in% tegene_multicopy]

boxplot(        a1,a2,a3,a4,a5,a6,main="H3K9me2 level\n(Col-0, rosette)",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("1 copy,no TE","1 copy,TE",">1 copy,no TE",">1 copy,TE","1 copy",">1 copy"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 5.5,line = 3)

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





#H3K9me2  CN 1,2,3-10,>10
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Supplem_Fig5CN_boxplot_linc_TEgene_CN-1_2_3_10_H3K9me2.pdf",height = 4,width = 4)
###########################################################################
par(mar=c(6,4,4,2))
a1<-chip.denovo.log2$K9.6909[ chip.denovo.quantstan$gene %in% lincs_1copy]
a2<-chip.denovo.log2$K9.6909[ chip.denovo.quantstan$gene %in% lincs_2copies]
a3<-chip.denovo.log2$K9.6909[chip.denovo.quantstan$gene %in% lincs3_10_copies]
a4<-chip.denovo.log2$K9.6909[ chip.denovo.quantstan$gene %in% lincs_more10copies]

a5<-chip.denovo.log2$K9.6909[chip.denovo.quantstan$gene %in% tegene_singlecopy ]
a6<-chip.denovo.log2$K9.6909[chip.denovo.quantstan$gene %in% tegene_2copies]
a7<-chip.denovo.log2$K9.6909[chip.denovo.quantstan$gene %in% tegene3_10_copies ]
a8<-chip.denovo.log2$K9.6909[chip.denovo.quantstan$gene %in% tegene_more10copies ]

boxplot(        a1,a2,a3,a4,a5,a6,a7,a8, main="H3K9me2 level\n(Col-0, rosette)",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("CN=1","CN=2","CN:3-10","CN>10","CN=1","CN>1","CN:3-10","CN>10"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 6.5,line = 3)
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
text(b,x=1.5,y=1.2)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1.2)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1.2)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1.2)
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1.2)
a<-wilcox.test(a6,a7)
b<-wilcox.test(a6,a7)
c<-wilcox.test(a6,a7)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=1.2)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1.2)
a<-wilcox.test(a1,a5)
b<-wilcox.test(a1,a5)
c<-wilcox.test(a1,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=2.27)
a<-wilcox.test(a2,a6)
b<-wilcox.test(a2,a6)
c<-wilcox.test(a2,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=2.27)
a<-wilcox.test(a3,a7)
b<-wilcox.test(a3,a7)
c<-wilcox.test(a3,a7)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=2.27)
a<-wilcox.test(a4,a8)
b<-wilcox.test(a4,a8)
c<-wilcox.test(a4,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=2.27)
#################
dev.off()




#H1 CN 1,>1, with TE without TE
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Supplem_Fig5CN_boxplot_linc_TEgene_CN-TE_H1_Col0.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% lincs_singlecopy_noTE]
a2<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% lincs_singlecopy_withTE]
a3<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% lincs_multicopy_noTE]
a4<-chip.denovo.log2$H1.6909[ chip.denovo.log2$gene %in% lincs_multicopy_withTE]
a5<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% tegene_singlecopy]
a6<-chip.denovo.log2$H1.6909[chip.denovo.log2$gene %in% tegene_multicopy]

boxplot(        a1,a2,a3,a4,a5,a6,main="H1 level\n(Col-0, rosette)",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("1 copy,no TE","1 copy,TE",">1 copy,no TE",">1 copy,TE","1 copy",">1 copy"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 5.5,line = 3)

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
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1)
#################
dev.off()


#H1  CN 1,2,3-10,>10
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Supplem_Fig5CN_boxplot_linc_TEgene_CN-1_2_3_10_H1.pdf",height = 4,width = 4)
###########################################################################
par(mar=c(6,4,4,2))
a1<-chip.denovo.log2$H1.6909[ chip.denovo.quantstan$gene %in% lincs_1copy]
a2<-chip.denovo.log2$H1.6909[ chip.denovo.quantstan$gene %in% lincs_2copies]
a3<-chip.denovo.log2$H1.6909[chip.denovo.quantstan$gene %in% lincs3_10_copies]
a4<-chip.denovo.log2$H1.6909[ chip.denovo.quantstan$gene %in% lincs_more10copies]

a5<-chip.denovo.log2$H1.6909[chip.denovo.quantstan$gene %in% tegene_singlecopy ]
a6<-chip.denovo.log2$H1.6909[chip.denovo.quantstan$gene %in% tegene_2copies]
a7<-chip.denovo.log2$H1.6909[chip.denovo.quantstan$gene %in% tegene3_10_copies ]
a8<-chip.denovo.log2$H1.6909[chip.denovo.quantstan$gene %in% tegene_more10copies ]

boxplot(        a1,a2,a3,a4,a5,a6,a7,a8, main="H1 level\n(Col-0, rosette)",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("CN=1","CN=2","CN:3-10","CN>10","CN=1","CN>1","CN:3-10","CN>10"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 6.5,line = 3)
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
text(b,x=1.5,y=1.2)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1.2)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1.2)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1.2)
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1.2)
a<-wilcox.test(a6,a7)
b<-wilcox.test(a6,a7)
c<-wilcox.test(a6,a7)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=1.2)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1.2)
a<-wilcox.test(a1,a5)
b<-wilcox.test(a1,a5)
c<-wilcox.test(a1,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=1.3)
a<-wilcox.test(a2,a6)
b<-wilcox.test(a2,a6)
c<-wilcox.test(a2,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=1.3)
a<-wilcox.test(a3,a7)
b<-wilcox.test(a3,a7)
c<-wilcox.test(a3,a7)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=1.3)
a<-wilcox.test(a4,a8)
b<-wilcox.test(a4,a8)
c<-wilcox.test(a4,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=1.3)
#################
dev.off()



#H3K27me3 CN 1,>1, with TE without TE
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Supplem_Fig5CN_boxplot_linc_TEgene_CN-TE_H3K27me3_Col0.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-chip.denovo.log2$K27.6909[ chip.denovo.log2$gene %in% lincs_singlecopy_noTE]
a2<-chip.denovo.log2$K27.6909[ chip.denovo.log2$gene %in% lincs_singlecopy_withTE]
a3<-chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% lincs_multicopy_noTE]
a4<-chip.denovo.log2$K27.6909[ chip.denovo.log2$gene %in% lincs_multicopy_withTE]
a5<-chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% tegene_singlecopy]
a6<-chip.denovo.log2$K27.6909[chip.denovo.log2$gene %in% tegene_multicopy]

boxplot(        a1,a2,a3,a4,a5,a6,main="H3K27me3 level\n(Col-0, rosette)",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("1 copy,no TE","1 copy,TE",">1 copy,no TE",">1 copy,TE","1 copy",">1 copy"),las=2, notch = T, outline = F, ylab="log2(ChIP/Input)")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 5.5,line = 3)

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
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1)
#################
dev.off()




#CG 1001G    CN 1,>1, with TE without TE
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Fig5CN_boxplot_linc_TEgene_CN-TE_CG_1001G.pdf",height = 3,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% lincs_singlecopy_noTE]
a2<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% lincs_singlecopy_withTE]
a3<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lincs_multicopy_noTE]
a4<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% lincs_multicopy_withTE]
a5<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% tegene_singlecopy]
a6<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% tegene_multicopy]

boxplot(        a1,a2,a3,a4,a5,a6,main="CG methylation level\n(Col-0, rosette)",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("1 copy,no TE","1 copy,TE",">1 copy,no TE",">1 copy,TE","1 copy",">1 copy"),las=2, notch = T, outline = F, ylab="methylation level")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 5.5,line = 3)

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
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1)
#################
dev.off()

#CG 1001G   linc with CN 1,2,3-10,>10
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Fig5CN_boxplot_linc_TEgene_CG_vs_CN_1_2_3_10.pdf",height = 4,width = 4)
###########################################################################
par(mar=c(6,4,4,2))
a1<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% lincs_1copy]
a2<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% lincs_2copies]
a3<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% lincs3_10_copies]
a4<-CG.1001.denovo$X6909[ CG.1001.denovo$transcript %in% lincs_more10copies]
a5<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% tegene_singlecopy]
a6<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% tegene_2copies]
a7<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% tegene3_10_copies]
a8<-CG.1001.denovo$X6909[CG.1001.denovo$transcript %in% tegene_more10copies]

boxplot(        a1,a2,a3,a4,a5,a6,a7,a8,main="CG methylation",cex.main=1.2,      col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("CN=1","CN=2","CN:3-10","CN>10","CN=1","CN>1","CN:3-10","CN>10"),las=2, notch = T, outline = F, ylab="methylation level",ylim=c(0,1.1))
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 5.5,line = 3)

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
text(b,x=1.5,y=0.98)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.98)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.98)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.98)
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.98)
a<-wilcox.test(a6,a7)
b<-wilcox.test(a6,a7)
c<-wilcox.test(a6,a7)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=0.98)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.98)
a<-wilcox.test(a1,a5)
b<-wilcox.test(a1,a5)
c<-wilcox.test(a1,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=1.1)
a<-wilcox.test(a2,a6)
b<-wilcox.test(a2,a6)
c<-wilcox.test(a2,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=1.1)
a<-wilcox.test(a3,a7)
b<-wilcox.test(a3,a7)
c<-wilcox.test(a3,a7)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=1.1)
a<-wilcox.test(a4,a8)
b<-wilcox.test(a4,a8)
c<-wilcox.test(a4,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=1.1)
#################
dev.off()



#CHH 1001G    CN 1,>1, with TE without TE
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Fig5CN_boxplot_linc_TEgene_CN-TE_CHH_1001G.pdf",height = 4,width = 4)
###########################################################################
par(mar=c(6,4,4,2))
a1<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% lincs_singlecopy_noTE]
a2<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% lincs_singlecopy_withTE]
a3<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% lincs_multicopy_noTE]
a4<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% lincs_multicopy_withTE]
a5<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% tegene_singlecopy]
a6<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% tegene_multicopy]

boxplot(        a1,a2,a3,a4,a5,a6,main="CHH methylation level\n(Col-0, rosette)",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("1 copy,no TE","1 copy,TE",">1 copy,no TE",">1 copy,TE","1 copy",">1 copy"),las=2, notch = T, outline = F, ylab="methylation level")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 5.5,line = 3)

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
text(b,x=1.5,y=0.11)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.11)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.11)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.11)
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.11)
#################
dev.off()



#CHH 1001G   linc with CN 1,2,3-10,>10
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Fig5CN_boxplot_linc_TEgene_CHH_vs_CN_1_2_3_10.pdf",height = 4,width = 4)
###########################################################################
par(mar=c(6,4,4,2))
a1<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% lincs_1copy]
a2<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% lincs_2copies]
a3<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% lincs3_10_copies]
a4<-CHH.1001.denovo$X6909[ CHH.1001.denovo$transcript %in% lincs_more10copies]
a5<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% tegene_singlecopy]
a6<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% tegene_2copies]
a7<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% tegene3_10_copies]
a8<-CHH.1001.denovo$X6909[CHH.1001.denovo$transcript %in% tegene_more10copies]

boxplot(        a1,a2,a3,a4,a5,a6,a7,a8,main="CHH methylation",cex.main=1.2,      col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("CN=1","CN=2","CN:3-10","CN>10","CN=1","CN>1","CN:3-10","CN>10"),las=2, notch = T, outline = F, ylab="methylation level")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 5.5,line = 3)

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
text(b,x=1.5,y=0.10)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.10)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.10)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.10)
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.10)
a<-wilcox.test(a6,a7)
b<-wilcox.test(a6,a7)
c<-wilcox.test(a6,a7)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=0.10)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=0.10)
a<-wilcox.test(a1,a5)
b<-wilcox.test(a1,a5)
c<-wilcox.test(a1,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=0.14)
a<-wilcox.test(a2,a6)
b<-wilcox.test(a2,a6)
c<-wilcox.test(a2,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=0.14)
a<-wilcox.test(a3,a7)
b<-wilcox.test(a3,a7)
c<-wilcox.test(a3,a7)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=0.14)
a<-wilcox.test(a4,a8)
b<-wilcox.test(a4,a8)
c<-wilcox.test(a4,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=0.14)
#################
dev.off()



# small RNA coverage

#sRNA 24nt    CN 1,>1, with TE without TE
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Fig5CN_boxplot_linc_TEgene_CN-TE_24nt.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% lincs_singlecopy_noTE]
a2<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% lincs_singlecopy_withTE]
a3<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% lincs_multicopy_noTE]
a4<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% lincs_multicopy_withTE]
a5<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% tegene_singlecopy]
a6<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% tegene_multicopy]

boxplot(        a1,a2,a3,a4,a5,a6,main="24nt sRNA coverage\n(Col-0, flowers)",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("1 copy,no TE","1 copy,TE",">1 copy,no TE",">1 copy,TE","1 copy",">1 copy"),las=2, notch = T, outline = F, ylab="methylation level")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 5.5,line = 3)

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
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1)
#################

a<-wilcox.test(a4,a6)
b<-wilcox.test(a4,a6)
c<-wilcox.test(a4,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5,y=2)

dev.off()

#sRNA 24nt  linc with CN 1,2,3-10,>10
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Fig5CN_boxplot_linc_TEgene_24ntcoverage_vs_CN.pdf",height = 4,width = 4)
###########################################################################
par(mar=c(6,4,4,2))
a1<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% lincs_1copy]
a2<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% lincs_2copies]
a3<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% lincs3_10_copies]
a4<-sRNA.24nt.denovo2021.RPM$X6909[ sRNA.24nt.denovo2021.RPM$gene %in% lincs_more10copies]
a5<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% tegene_singlecopy]
a6<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% tegene_2copies]
a7<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% tegene3_10_copies]
a8<-sRNA.24nt.denovo2021.RPM$X6909[sRNA.24nt.denovo2021.RPM$gene %in% tegene_more10copies]

boxplot(        a1,a2,a3,a4,a5,a6,a7,a8,main="24nt sRNA coverage\n(Col-0, flowers)",cex.main=1.2,      col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("CN=1","CN=2","CN:3-10","CN>10","CN=1","CN>1","CN:3-10","CN>10"),las=2, notch = T, outline = F, ylab="methylation level")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 5.5,line = 3)

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
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
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
a<-wilcox.test(a6,a7)
b<-wilcox.test(a6,a7)
c<-wilcox.test(a6,a7)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=1.5)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1.5)
a<-wilcox.test(a1,a5)
b<-wilcox.test(a1,a5)
c<-wilcox.test(a1,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=2.2)
a<-wilcox.test(a2,a6)
b<-wilcox.test(a2,a6)
c<-wilcox.test(a2,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=2.2)
a<-wilcox.test(a3,a7)
b<-wilcox.test(a3,a7)
c<-wilcox.test(a3,a7)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=2.2)
a<-wilcox.test(a4,a8)
b<-wilcox.test(a4,a8)
c<-wilcox.test(a4,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=2.2)
#################

dev.off()


#expression variation


#1001G rosette expression variation CN 1,2,3-10,>10
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Supplem_Fig5CN_boxplot_linc_TEgene_expression_var_CN_1_2_3_10_1001G.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% lincs_1copy&denovo2021.TPMs.genes.1001G$max>1 ]
a2<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% lincs_2copies&denovo2021.TPMs.genes.1001G$max>1]
a3<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lincs3_10_copies&denovo2021.TPMs.genes.1001G$max>1]
a4<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% lincs_more10copies&denovo2021.TPMs.genes.1001G$max>1]

a5<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% tegene_singlecopy&denovo2021.TPMs.genes.1001G$max>1  ]
a6<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% tegene_2copies&denovo2021.TPMs.genes.1001G$max>1 ]
a7<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% tegene3_10_copies&denovo2021.TPMs.genes.1001G$max>1  ]
a8<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% tegene_more10copies&denovo2021.TPMs.genes.1001G$max>1 ]

boxplot(        a1,a2,a3,a4,a5,a6,a7,a8, main="Expression variation\n(461 accessions, rosette)",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("CN=1","CN=2","CN:3-10","CN>10","CN=1","CN>1","CN:3-10","CN>10"),las=2, notch = T, outline = F, ylab="coefficient of variance")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 6.5,line = 3)
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
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
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
a<-wilcox.test(a6,a7)
b<-wilcox.test(a6,a7)
c<-wilcox.test(a6,a7)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=15)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=15)
a<-wilcox.test(a1,a5)
b<-wilcox.test(a1,a5)
c<-wilcox.test(a1,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=22)
a<-wilcox.test(a2,a6)
b<-wilcox.test(a2,a6)
c<-wilcox.test(a2,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=22)
a<-wilcox.test(a3,a7)
b<-wilcox.test(a3,a7)
c<-wilcox.test(a3,a7)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=22)
a<-wilcox.test(a4,a8)
b<-wilcox.test(a4,a8)
c<-wilcox.test(a4,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=22)
#################
dev.off()



#1001G rosette expression variation CN 1,>1, with TE without TE
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Fig5CN_boxplot_linc_TEgene_CN_expression_variation_1001G.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% lincs_singlecopy_noTE&denovo2021.TPMs.genes.1001G$max>1 ]
a2<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% lincs_singlecopy_withTE&denovo2021.TPMs.genes.1001G$max>1]
a3<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lincs_multicopy_noTE&denovo2021.TPMs.genes.1001G$max>1]
a4<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% lincs_multicopy_withTE&denovo2021.TPMs.genes.1001G$max>1]
a5<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% tegene_singlecopy&denovo2021.TPMs.genes.1001G$max>1  ]
a6<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% tegene_multicopy&denovo2021.TPMs.genes.1001G$max>1 ]

boxplot(        a1,a2,a3,a4,a5,a6,main="Expression variation\n(461 accessions, rosette)",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("1 copy,no TE","1 copy,TE",">1 copy,no TE",">1 copy,TE","1 copy",">1 copy"),las=2, notch = T, outline = F, ylab="coefficient of variance")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 5.5,line = 3)
mtext (text = "TPMmax>1",side = 3,at = 1.5,line = -0.9)

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
a<-wilcox.test(a1,a3)
b<-wilcox.test(a1,a3)
c<-wilcox.test(a1,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=20)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=15)
a<-wilcox.test(a2,a4)
b<-wilcox.test(a2,a4)
c<-wilcox.test(a2,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=22)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
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



#1001G rosette expression variation CN 1,>1, with TE without TE 1<TPM<2
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Supplem_Fig5CN_boxplot_linc_TEgene_CN_expression_1001G_max_1_2.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% lincs_singlecopy_noTE&denovo2021.TPMs.genes.1001G$max>1  &denovo2021.TPMs.genes.1001G$max<2]
a2<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% lincs_singlecopy_withTE&denovo2021.TPMs.genes.1001G$max>1&denovo2021.TPMs.genes.1001G$max<2]
a3<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% lincs_multicopy_noTE&denovo2021.TPMs.genes.1001G$max>1&denovo2021.TPMs.genes.1001G$max<2]
a4<-denovo2021.TPMs.genes.1001G$variance[ denovo2021.TPMs.genes.1001G$gene %in% lincs_multicopy_withTE&denovo2021.TPMs.genes.1001G$max>1&denovo2021.TPMs.genes.1001G$max<2]
a5<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% tegene_singlecopy&denovo2021.TPMs.genes.1001G$max>1 &denovo2021.TPMs.genes.1001G$max<2 ]
a6<-denovo2021.TPMs.genes.1001G$variance[denovo2021.TPMs.genes.1001G$gene %in% tegene_multicopy&denovo2021.TPMs.genes.1001G$max>1 &denovo2021.TPMs.genes.1001G$max<2]

boxplot(        a1,a2,a3,a4,a5,a6,main="Expression variation\n(461 accessions, rosette)",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("1 copy,no TE","1 copy,TE",">1 copy,no TE",">1 copy,TE","1 copy",">1 copy"),las=2, notch = T, outline = F, ylab="coefficient of variance")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 5.5,line = 3)
mtext (text = "1<TPMmax<2",side = 3,at = 1.5,line = -0.9)

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
a<-wilcox.test(a1,a3)
b<-wilcox.test(a1,a3)
c<-wilcox.test(a1,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=20)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=15)
a<-wilcox.test(a2,a4)
b<-wilcox.test(a2,a4)
c<-wilcox.test(a2,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=22)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
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


# epigenetic variation 




#H3K9me2 variation CN 1,2,3-10,>10
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Supplem_Fig5CN_boxplot_linc_TEgene_CN_K9variation.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-chip.denovo.quantstan$sd.key9[ chip.denovo.quantstan$gene %in% lincs_1copy]
a2<-chip.denovo.quantstan$sd.key9[ chip.denovo.quantstan$gene %in% lincs_2copies]
a3<-chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% lincs3_10_copies]
a4<-chip.denovo.quantstan$sd.key9[ chip.denovo.quantstan$gene %in% lincs_more10copies]

a5<-chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% tegene_singlecopy ]
a6<-chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% tegene_2copies]
a7<-chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% tegene3_10_copies ]
a8<-chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% tegene_more10copies ]

boxplot(        a1,a2,a3,a4,a5,a6,a7,a8, main="H3K9me2 variation\n(13 accessions, rosette)",cex.main=1.2,     col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("CN=1","CN=2","CN:3-10","CN>10","CN=1","CN>1","CN:3-10","CN>10"),las=2, notch = T, outline = F, ylab="standard deviation of normalized ChIP signal")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 6.5,line = 3)
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
text(b,x=1.5,y=1.2)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1.2)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1.2)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1.2)
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1.2)
a<-wilcox.test(a6,a7)
b<-wilcox.test(a6,a7)
c<-wilcox.test(a6,a7)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=6.5,y=1.2)
a<-wilcox.test(a7,a8)
b<-wilcox.test(a7,a8)
c<-wilcox.test(a7,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=7.5,y=1.2)
a<-wilcox.test(a1,a5)
b<-wilcox.test(a1,a5)
c<-wilcox.test(a1,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1,y=1.27)
a<-wilcox.test(a2,a6)
b<-wilcox.test(a2,a6)
c<-wilcox.test(a2,a6)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2,y=1.27)
a<-wilcox.test(a3,a7)
b<-wilcox.test(a3,a7)
c<-wilcox.test(a3,a7)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3,y=1.27)
a<-wilcox.test(a4,a8)
b<-wilcox.test(a4,a8)
c<-wilcox.test(a4,a8)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4,y=1.27)
#################
dev.off()




#H3K9me2 variation CN=1 , >1 with TE without TE 
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Fig5CN_boxplot_linc_TEgene_K9variation_CN_TE_noTE.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-chip.denovo.quantstan$sd.key9[ chip.denovo.quantstan$gene %in% lincs_singlecopy_noTE]
a2<-chip.denovo.quantstan$sd.key9[ chip.denovo.quantstan$gene %in% lincs_singlecopy_withTE]
a3<-chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% lincs_multicopy_noTE]
a4<-chip.denovo.quantstan$sd.key9[ chip.denovo.quantstan$gene %in% lincs_multicopy_withTE]

a5<-chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% tegene_singlecopy ]
a6<-chip.denovo.quantstan$sd.key9[chip.denovo.quantstan$gene %in% tegene_multicopy]

boxplot(        a1,a2,a3,a4,a5,a6,main="H3K9me2 variation\n(13 accessions, rosette)",cex.main=1.2,      col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("1 copy,no TE","1 copy,TE",">1 copy,no TE",">1 copy,TE","1 copy",">1 copy"),las=2, notch = T, outline = F, ylab="standard deviation of normalized ChIP signal")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 6.5,line = 3)
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
text(b,x=1.5,y=1.2)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1.2)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1.2)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1.2)
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=1.2)

#################
dev.off()


#H1 variation CN=1 , >1 with TE without TE 
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Fig5CN_boxplot_linc_TEgene_H1variation_CN_TE_noTE.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-chip.denovo.quantstan$sd.hist1[ chip.denovo.quantstan$gene %in% lincs_singlecopy_noTE]
a2<-chip.denovo.quantstan$sd.hist1[ chip.denovo.quantstan$gene %in% lincs_singlecopy_withTE]
a3<-chip.denovo.quantstan$sd.hist1[chip.denovo.quantstan$gene %in% lincs_multicopy_noTE]
a4<-chip.denovo.quantstan$sd.hist1[ chip.denovo.quantstan$gene %in% lincs_multicopy_withTE]

a5<-chip.denovo.quantstan$sd.hist1[chip.denovo.quantstan$gene %in% tegene_singlecopy ]
a6<-chip.denovo.quantstan$sd.hist1[chip.denovo.quantstan$gene %in% tegene_multicopy]

boxplot(        a1,a2,a3,a4,a5,a6,main="H1 variation\n(13 accessions, rosette)",cex.main=1.2,      col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("1 copy,no TE","1 copy,TE",">1 copy,no TE",">1 copy,TE","1 copy",">1 copy"),las=2, notch = T, outline = F, ylab="standard deviation of normalized ChIP signal")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 6.5,line = 3)
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
text(b,x=1.5,y=0.8)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.8)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.8)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=1)
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.8)

#################
dev.off()



#CG variation CN=1 , >1 with TE without TE 
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Fig5CN_boxplot_linc_TEgene_CG_1001G_variation_CN_TE_noTE.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-CG.1001.denovo$sd[ CG.1001.denovo$transcript %in% lincs_singlecopy_noTE]
a2<-CG.1001.denovo$sd[ CG.1001.denovo$transcript %in% lincs_singlecopy_withTE]
a3<-CG.1001.denovo$sd[CG.1001.denovo$transcript %in% lincs_multicopy_noTE]
a4<-CG.1001.denovo$sd[ CG.1001.denovo$transcript %in% lincs_multicopy_withTE]

a5<-CG.1001.denovo$sd[CG.1001.denovo$transcript %in% tegene_singlecopy ]
a6<-CG.1001.denovo$sd[CG.1001.denovo$transcript %in% tegene_multicopy]

boxplot(        a1,a2,a3,a4,a5,a6,main="CG methylation variation\n(450 accessions, rosette)",cex.main=1.2,      col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("1 copy,no TE","1 copy,TE",">1 copy,no TE",">1 copy,TE","1 copy",">1 copy"),las=2, notch = T, outline = F, ylab="sd of methylation level")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 6.5,line = 3)
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
text(b,x=1.5,y=0.5)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.5)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.5)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.5)
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.5)

#################
dev.off()


#CHH variation CN=1 , >1 with TE without TE 
###########################################################################
pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Fig5CN_boxplot_linc_TEgene_CHH_1001G_variation_CN_TE_noTE.pdf",height = 4,width = 3)
###########################################################################
par(mar=c(6,4,4,2))
a1<-CHH.1001.denovo$sd[ CHH.1001.denovo$transcript %in% lincs_singlecopy_noTE]
a2<-CHH.1001.denovo$sd[ CHH.1001.denovo$transcript %in% lincs_singlecopy_withTE]
a3<-CHH.1001.denovo$sd[CHH.1001.denovo$transcript %in% lincs_multicopy_noTE]
a4<-CHH.1001.denovo$sd[ CHH.1001.denovo$transcript %in% lincs_multicopy_withTE]

a5<-CHH.1001.denovo$sd[CHH.1001.denovo$transcript %in% tegene_singlecopy ]
a6<-CHH.1001.denovo$sd[CHH.1001.denovo$transcript %in% tegene_multicopy]

boxplot(        a1,a2,a3,a4,a5,a6,main="CHH methylation variation\n(450 accessions, rosette)",cex.main=1.2,      col=c("#F2AB54","#C55B16","#C55B16","#C55B16","#673A8E","#673A8E","#673A8E","#673A8E"), names=c("1 copy,no TE","1 copy,TE",">1 copy,no TE",">1 copy,TE","1 copy",">1 copy"),las=2, notch = T, outline = F, ylab="sd of methylation level")
mtext (text = "lincRNAs",side = 1,at = 2.5,line = 3)
mtext (text = "TE genes",side = 1,at = 6.5,line = 3)
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
text(b,x=1.5,y=0.5)
a<-wilcox.test(a2,a3)
b<-wilcox.test(a2,a3)
c<-wilcox.test(a2,a3)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.5)
a<-wilcox.test(a3,a4)
b<-wilcox.test(a3,a4)
c<-wilcox.test(a3,a4)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.5)
a<-wilcox.test(a4,a5)
b<-wilcox.test(a4,a5)
c<-wilcox.test(a4,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=4.5,y=0.5)
a<-wilcox.test(a6,a5)
b<-wilcox.test(a6,a5)
c<-wilcox.test(a6,a5)
d<-mean(c(a$p.value,b$p.value,c$p.value))
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=5.5,y=0.5)

#################
dev.off()





# do single/multicopy lincRNAs with and without TEs have different expression? 







pdf(file="Z:/01_POSTDOC/04_writing/lncRNA variation paper/figures/Fig5_CN/Fig5_lincs_CN_inTAIR10_hist.pdf",height = 3,width = 3)
par(mar=c(4,4,4,2)) 
hist (CN_linc_27genomes$TAIR10[CN_linc_27genomes$TAIR10>0], breaks=10, main="Copy number\n of lincRNA loci", ylab="number of loci", xlab='CN in TAIR10',xaxt='n', col="#F2AB54",las=2,panel.first=grid(lty = 1, col = "darkgray"),ylim=c(0,950))
axis(side=1, at=seq(0,100, 10), labels=seq(0,100,10),las=2)
dev.off()

hist (ifelse(CN_linc_27genomes$TAIR10>=8,8,CN_linc_27genomes$TAIR10), breaks=12,main="length of TE patches\n within lincRNA loci", ylab="number of TE patches", xlab='TE patch length,bp',xaxt='n', col="#F2AB54",las=2,,panel.first=grid(lty = 1, col = "darkgray",lwd=0.6), xlim=c(0,12))
axis(side=1, at=seq(0,500, 50), labels=seq(0,500,50),las=2,lwd=0.6)
axis(side=1, at=575, labels=">500",las=2)
text(paste("median =", median(lincRNAs_TE_coverage.TAIR10$TE_length),"bp"),x=400,y=1500)
text(paste("min =", min(lincRNAs_TE_coverage.TAIR10$TE_length),"bp"),x=400,y=1300)


































