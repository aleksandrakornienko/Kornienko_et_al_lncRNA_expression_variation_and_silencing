setwd("/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/")
setwd("Z:/01_POSTDOC/")

#upload public annotations 

#methylation 1001G

CG.pc<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.pc.1001G_data.bed")
CG.pc_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.pc_TES.1001G_data.bed")
CG.pc_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.pc_TSS.1001G_data.bed")
rownames(CG.pc)<-CG.pc$transcript
rownames(CG.pc_TES)<-CG.pc_TES$transcript
rownames(CG.pc_TSS)<-CG.pc_TSS$transcript


CG.pc[CG.pc == 22] <- NA
CG.pc_TES[CG.pc_TES == 22] <- NA
CG.pc_TSS[CG.pc_TSS == 22] <- NA

CG.linc<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.linc.1001G_data.bed")
CG.linc_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.linc_TES.1001G_data.bed")
CG.linc_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.linc_TSS.1001G_data.bed")
rownames(CG.linc)<-CG.linc$transcript
rownames(CG.linc_TES)<-CG.linc_TES$transcript
rownames(CG.linc_TSS)<-CG.linc_TSS$transcript


CG.linc[CG.linc == 22] <- NA
CG.linc_TES[CG.linc_TES == 22] <- NA
CG.linc_TSS[CG.linc_TSS == 22] <- NA

CG.te<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.te.1001G_data.bed")
CG.te_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.te_TES.1001G_data.bed")
CG.te_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.te_TSS.1001G_data.bed")
rownames(CG.te)<-CG.te$transcript
rownames(CG.te_TES)<-CG.te_TES$transcript
rownames(CG.te_TSS)<-CG.te_TSS$transcript

CG.te[CG.te == 22] <- NA
CG.te_TES[CG.te_TES == 22] <- NA
CG.te_TSS[CG.te_TSS == 22] <- NA

CG.as<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.as.1001G_data.bed")
CG.as_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.as_TES.1001G_data.bed")
CG.as_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.as_TSS.1001G_data.bed")
rownames(CG.as)<-CG.as$transcript
rownames(CG.as_TES)<-CG.as_TES$transcript
rownames(CG.as_TSS)<-CG.as_TSS$transcript

CG.as[CG.as == 22] <- NA
CG.as_TES[CG.as_TES == 22] <- NA
CG.as_TSS[CG.as_TSS == 22] <- NA



CG.ar11_pc<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.ar11_pc.1001G_data.bed")
CG.ar11_pc_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.ar11_pc_TES.1001G_data.bed")
CG.ar11_pc_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CG.coverage.ar11_pc_TSS.1001G_data.bed")
rownames(CG.ar11_pc)<-CG.ar11_pc$transcript
rownames(CG.ar11_pc_TES)<-CG.ar11_pc_TES$transcript
rownames(CG.ar11_pc_TSS)<-CG.ar11_pc_TSS$transcript

CG.ar11_pc[CG.ar11_pc == 22] <- NA
CG.ar11_pc[CG.ar11_pc == 22] <- NA
CG.ar11_pc[CG.ar11_pc == 22] <- NA

CHH.pc<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHH.coverage.pc.1001G_data.bed")
CHH.pc_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHH.coverage.pc_TES.1001G_data.bed")
CHH.pc_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHH.coverage.pc_TSS.1001G_data.bed")
rownames(CHH.pc)<-CHH.pc$transcript
rownames(CHH.pc_TES)<-CHH.pc_TES$transcript
rownames(CHH.pc_TSS)<-CHH.pc_TSS$transcript

CHH.pc[CHH.pc == 22] <- NA
CHH.pc_TES[CHH.pc_TES == 22] <- NA
CHH.pc_TSS[CHH.pc_TSS == 22] <- NA

CHH.linc<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHH.coverage.linc.1001G_data.bed")
CHH.linc_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHH.coverage.linc_TES.1001G_data.bed")
CHH.linc_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHH.coverage.linc_TSS.1001G_data.bed")
rownames(CHH.linc)<-CHH.linc$transcript
rownames(CHH.linc_TES)<-CHH.linc_TES$transcript
rownames(CHH.linc_TSS)<-CHH.linc_TSS$transcript

CHH.linc[CHH.linc == 22] <- NA
CHH.linc_TES[CHH.linc_TES == 22] <- NA
CHH.linc_TSS[CHH.linc_TSS == 22] <- NA

CHH.te<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHH.coverage.te.1001G_data.bed")
CHH.te_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHH.coverage.te_TES.1001G_data.bed")
CHH.te_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHH.coverage.te_TSS.1001G_data.bed")
rownames(CHH.te)<-CHH.te$transcript
rownames(CHH.te_TES)<-CHH.te_TES$transcript
rownames(CHH.te_TSS)<-CHH.te_TSS$transcript

CHH.te[CHH.te == 22] <- NA
CHH.te_TES[CHH.te_TES == 22] <- NA
CHH.te_TSS[CHH.te_TSS == 22] <- NA

CHH.as<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHH.coverage.as.1001G_data.bed")
CHH.as_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHH.coverage.as_TES.1001G_data.bed")
CHH.as_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHH.coverage.as_TSS.1001G_data.bed")
rownames(CHH.as)<-CHH.as$transcript
rownames(CHH.as_TES)<-CHH.as_TES$transcript
rownames(CHH.as_TSS)<-CHH.as_TSS$transcript

CHH.as[CHH.as == 22] <- NA
CHH.as_TES[CHH.as_TES == 22] <- NA
CHH.as_TSS[CHH.as_TSS == 22] <- NA




CHG.pc<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHG.coverage.pc.1001G_data.bed")
CHG.pc_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHG.coverage.pc_TES.1001G_data.bed")
CHG.pc_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHG.coverage.pc_TSS.1001G_data.bed")
rownames(CHG.pc)<-CHG.pc$transcript
rownames(CHG.pc_TES)<-CHG.pc_TES$transcript
rownames(CHG.pc_TSS)<-CHG.pc_TSS$transcript

CHG.pc[CHG.pc == 22] <- NA
CHG.pc_TES[CHG.pc_TES == 22] <- NA
CHG.pc_TSS[CHG.pc_TSS == 22] <- NA

CHG.linc<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHG.coverage.linc.1001G_data.bed")
CHG.linc_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHG.coverage.linc_TES.1001G_data.bed")
CHG.linc_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHG.coverage.linc_TSS.1001G_data.bed")
rownames(CHG.linc)<-CHG.linc$transcript
rownames(CHG.linc_TES)<-CHG.linc_TES$transcript
rownames(CHG.linc_TSS)<-CHG.linc_TSS$transcript

CHG.linc[CHG.linc == 22] <- NA
CHG.linc_TES[CHG.linc_TES == 22] <- NA
CHG.linc_TSS[CHG.linc_TSS == 22] <- NA

CHG.te<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHG.coverage.te.1001G_data.bed")
CHG.te_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHG.coverage.te_TES.1001G_data.bed")
CHG.te_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHG.coverage.te_TSS.1001G_data.bed")
rownames(CHG.te)<-CHG.te$transcript
rownames(CHG.te_TES)<-CHG.te_TES$transcript
rownames(CHG.te_TSS)<-CHG.te_TSS$transcript

CHG.te[CHG.te == 22] <- NA
CHG.te_TES[CHG.te_TES == 22] <- NA
CHG.te_TSS[CHG.te_TSS == 22] <- NA

CHG.as<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHG.coverage.as.1001G_data.bed")
CHG.as_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHG.coverage.as_TES.1001G_data.bed")
CHG.as_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/CHG.coverage.as_TSS.1001G_data.bed")
rownames(CHG.as)<-CHG.as$transcript
rownames(CHG.as_TES)<-CHG.as_TES$transcript
rownames(CHG.as_TSS)<-CHG.as_TSS$transcript

CHG.as[CHG.as == 22] <- NA
CHG.as_TES[CHG.as_TES == 22] <- NA
CHG.as_TSS[CHG.as_TSS == 22] <- NA


CG.pc$mean<-apply(CG.pc[,7:450],1,mean)
CG.pc$max<-apply(CG.pc[,7:450],1,max)
CG.pc$range<-apply(CG.pc[,7:450],1,max)-apply(CG.pc[,7:450],1,min)
CG.pc$range_percent_of_mean<-100*(apply(CG.pc[,7:450],1,max)-apply(CG.pc[,7:450],1,min))/apply(CG.pc[,7:450],1,mean)
CG.pc$sd<-apply(CG.pc[,7:450],1,sd)
CG.pc$variance<-apply(CG.pc[,7:450],1,sd)/apply(CG.pc[,7:450],1,mean)

CG.linc$mean<-apply(CG.linc[,7:450],1,mean)
CG.linc$max<-apply(CG.linc[,7:450],1,max)
CG.linc$range<-apply(CG.linc[,7:450],1,max)-apply(CG.linc[,7:450],1,min)
CG.linc$range_percent_of_mean<-100*(apply(CG.linc[,7:450],1,max)-apply(CG.linc[,7:450],1,min))/apply(CG.linc[,7:450],1,mean)
CG.linc$sd<-apply(CG.linc[,7:450],1,sd)
CG.linc$variance<-apply(CG.linc[,7:450],1,sd)/apply(CG.linc[,7:450],1,mean)

CG.as$mean<-apply(CG.as[,7:450],1,mean)
CG.as$max<-apply(CG.as[,7:450],1,max)
CG.as$range<-apply(CG.as[,7:450],1,max)-apply(CG.as[,7:450],1,min)
CG.as$range_percent_of_mean<-100*(apply(CG.as[,7:450],1,max)-apply(CG.as[,7:450],1,min))/apply(CG.as[,7:450],1,mean)
CG.as$sd<-apply(CG.as[,7:450],1,sd)
CG.as$variance<-apply(CG.as[,7:450],1,sd)/apply(CG.as[,7:450],1,mean)

CG.te$mean<-apply(CG.te[,7:450],1,mean)
CG.te$max<-apply(CG.te[,7:450],1,max)
CG.te$range<-apply(CG.te[,7:450],1,max)-apply(CG.te[,7:450],1,min)
CG.te$range_percent_of_mean<-100*(apply(CG.te[,7:450],1,max)-apply(CG.te[,7:450],1,min))/apply(CG.te[,7:450],1,mean)
CG.te$sd<-apply(CG.te[,7:450],1,sd)
CG.te$variance<-apply(CG.te[,7:450],1,sd)/apply(CG.te[,7:450],1,mean)

CG.ar11_pc$mean<-apply(CG.ar11_pc[,7:450],1,mean)
CG.ar11_pc$max<-apply(CG.ar11_pc[,7:450],1,max)
CG.ar11_pc$range<-apply(CG.ar11_pc[,7:450],1,max)-apply(CG.ar11_pc[,7:450],1,min)
CG.ar11_pc$range_percent_of_mean<-100*(apply(CG.ar11_pc[,7:450],1,max)-apply(CG.ar11_pc[,7:450],1,min))/apply(CG.ar11_pc[,7:450],1,mean)
CG.ar11_pc$sd<-apply(CG.ar11_pc[,7:450],1,sd)
CG.ar11_pc$variance<-apply(CG.ar11_pc[,7:450],1,sd)/apply(CG.ar11_pc[,7:450],1,mean)



CHH.pc$mean<-apply(CHH.pc[,7:450],1,mean)
CHH.pc$max<-apply(CHH.pc[,7:450],1,max)
CHH.pc$range<-apply(CHH.pc[,7:450],1,max)-apply(CHH.pc[,7:450],1,min)
CHH.pc$range_percent_of_mean<-100*(apply(CHH.pc[,7:450],1,max)-apply(CHH.pc[,7:450],1,min))/apply(CHH.pc[,7:450],1,mean)
CHH.pc$sd<-apply(CHH.pc[,7:450],1,sd)
CHH.pc$variance<-apply(CHH.pc[,7:450],1,sd)/apply(CHH.pc[,7:450],1,mean)

CHH.linc$mean<-apply(CHH.linc[,7:450],1,mean)
CHH.linc$max<-apply(CHH.linc[,7:450],1,max)
CHH.linc$range<-apply(CHH.linc[,7:450],1,max)-apply(CHH.linc[,7:450],1,min)
CHH.linc$range_percent_of_mean<-100*(apply(CHH.linc[,7:450],1,max)-apply(CHH.linc[,7:450],1,min))/apply(CHH.linc[,7:450],1,mean)
CHH.linc$sd<-apply(CHH.linc[,7:450],1,sd)
CHH.linc$variance<-apply(CHH.linc[,7:450],1,sd)/apply(CHH.linc[,7:450],1,mean)

CHH.as$mean<-apply(CHH.as[,7:450],1,mean)
CHH.as$max<-apply(CHH.as[,7:450],1,max)
CHH.as$range<-apply(CHH.as[,7:450],1,max)-apply(CHH.as[,7:450],1,min)
CHH.as$range_percent_of_mean<-100*(apply(CHH.as[,7:450],1,max)-apply(CHH.as[,7:450],1,min))/apply(CHH.as[,7:450],1,mean)
CHH.as$sd<-apply(CHH.as[,7:450],1,sd)
CHH.as$variance<-apply(CHH.as[,7:450],1,sd)/apply(CHH.as[,7:450],1,mean)

CHH.te$mean<-apply(CHH.te[,7:450],1,mean)
CHH.te$max<-apply(CHH.te[,7:450],1,max)
CHH.te$range<-apply(CHH.te[,7:450],1,max)-apply(CHH.te[,7:450],1,min)
CHH.te$range_percent_of_mean<-100*(apply(CHH.te[,7:450],1,max)-apply(CHH.te[,7:450],1,min))/apply(CHH.te[,7:450],1,mean)
CHH.te$sd<-apply(CHH.te[,7:450],1,sd)
CHH.te$variance<-apply(CHH.te[,7:450],1,sd)/apply(CHH.te[,7:450],1,mean)

CHH.ar11_pc$mean<-apply(CHH.ar11_pc[,7:450],1,mean)
CHH.ar11_pc$max<-apply(CHH.ar11_pc[,7:450],1,max)
CHH.ar11_pc$range<-apply(CHH.ar11_pc[,7:450],1,max)-apply(CHH.ar11_pc[,7:450],1,min)
CHH.ar11_pc$range_percent_of_mean<-100*(apply(CHH.ar11_pc[,7:450],1,max)-apply(CHH.ar11_pc[,7:450],1,min))/apply(CHH.ar11_pc[,7:450],1,mean)
CHH.ar11_pc$sd<-apply(CHH.ar11_pc[,7:450],1,sd)
CHH.ar11_pc$variance<-apply(CHH.ar11_pc[,7:450],1,sd)/apply(CHH.ar11_pc[,7:450],1,mean)




CHG.pc$mean<-apply(CHG.pc[,7:450],1,mean)
CHG.pc$max<-apply(CHG.pc[,7:450],1,max)
CHG.pc$range<-apply(CHG.pc[,7:450],1,max)-apply(CHG.pc[,7:450],1,min)
CHG.pc$range_percent_of_mean<-100*(apply(CHG.pc[,7:450],1,max)-apply(CHG.pc[,7:450],1,min))/apply(CHG.pc[,7:450],1,mean)
CHG.pc$sd<-apply(CHG.pc[,7:450],1,sd)
CHG.pc$variance<-apply(CHG.pc[,7:450],1,sd)/apply(CHG.pc[,7:450],1,mean)

CHG.linc$mean<-apply(CHG.linc[,7:450],1,mean)
CHG.linc$max<-apply(CHG.linc[,7:450],1,max)
CHG.linc$range<-apply(CHG.linc[,7:450],1,max)-apply(CHG.linc[,7:450],1,min)
CHG.linc$range_percent_of_mean<-100*(apply(CHG.linc[,7:450],1,max)-apply(CHG.linc[,7:450],1,min))/apply(CHG.linc[,7:450],1,mean)
CHG.linc$sd<-apply(CHG.linc[,7:450],1,sd)
CHG.linc$variance<-apply(CHG.linc[,7:450],1,sd)/apply(CHG.linc[,7:450],1,mean)

CHG.as$mean<-apply(CHG.as[,7:450],1,mean)
CHG.as$max<-apply(CHG.as[,7:450],1,max)
CHG.as$range<-apply(CHG.as[,7:450],1,max)-apply(CHG.as[,7:450],1,min)
CHG.as$range_percent_of_mean<-100*(apply(CHG.as[,7:450],1,max)-apply(CHG.as[,7:450],1,min))/apply(CHG.as[,7:450],1,mean)
CHG.as$sd<-apply(CHG.as[,7:450],1,sd)
CHG.as$variance<-apply(CHG.as[,7:450],1,sd)/apply(CHG.as[,7:450],1,mean)

CHG.te$mean<-apply(CHG.te[,7:450],1,mean)
CHG.te$max<-apply(CHG.te[,7:450],1,max)
CHG.te$range<-apply(CHG.te[,7:450],1,max)-apply(CHG.te[,7:450],1,min)
CHG.te$range_percent_of_mean<-100*(apply(CHG.te[,7:450],1,max)-apply(CHG.te[,7:450],1,min))/apply(CHG.te[,7:450],1,mean)
CHG.te$sd<-apply(CHG.te[,7:450],1,sd)
CHG.te$variance<-apply(CHG.te[,7:450],1,sd)/apply(CHG.te[,7:450],1,mean)

CHG.ar11_pc$mean<-apply(CHG.ar11_pc[,7:450],1,mean)
CHG.ar11_pc$max<-apply(CHG.ar11_pc[,7:450],1,max)
CHG.ar11_pc$range<-apply(CHG.ar11_pc[,7:450],1,max)-apply(CHG.ar11_pc[,7:450],1,min)
CHG.ar11_pc$range_percent_of_mean<-100*(apply(CHG.ar11_pc[,7:450],1,max)-apply(CHG.ar11_pc[,7:450],1,min))/apply(CHG.ar11_pc[,7:450],1,mean)
CHG.ar11_pc$sd<-apply(CHG.ar11_pc[,7:450],1,sd)
CHG.ar11_pc$variance<-apply(CHG.ar11_pc[,7:450],1,sd)/apply(CHG.ar11_pc[,7:450],1,mean)






