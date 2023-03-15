setwd("/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/")
setwd("Z:/01_POSTDOC/")

#methylation 1001G

CG.1001.denovo<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CG.coverage.denovo.1001G_data.bed")
CG.1001.denovo_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CG.coverage.denovo_TES.1001G_data.bed")
CG.1001.denovo_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CG.coverage.denovo_TSS.1001G_data.bed")
rownames(CG.1001.denovo)<-CG.1001.denovo$transcript
rownames(CG.1001.denovo_TES)<-CG.1001.denovo_TES$transcript
rownames(CG.1001.denovo_TSS)<-CG.1001.denovo_TSS$transcript
CG.1001.denovo[CG.1001.denovo == 22] <- NA
CG.1001.denovo_TES[CG.1001.denovo_TES == 22] <- NA
CG.1001.denovo_TSS[CG.1001.denovo_TSS == 22] <- NA



CHG.1001.denovo<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CHG.coverage.denovo.1001G_data.bed")
CHG.1001.denovo_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CHG.coverage.denovo_TES.1001G_data.bed")
CHG.1001.denovo_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CHG.coverage.denovo_TSS.1001G_data.bed")
rownames(CHG.1001.denovo)<-CHG.1001.denovo$transcript
rownames(CHG.1001.denovo_TES)<-CHG.1001.denovo_TES$transcript
rownames(CHG.1001.denovo_TSS)<-CHG.1001.denovo_TSS$transcript
CHG.1001.denovo[CHG.1001.denovo == 22] <- NA
CHG.1001.denovo_TES[CHG.1001.denovo_TES == 22] <- NA
CHG.1001.denovo_TSS[CHG.1001.denovo_TSS == 22] <- NA



CHH.1001.denovo<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CHH.coverage.denovo.1001G_data.bed")
CHH.1001.denovo_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CHH.coverage.denovo_TES.1001G_data.bed")
CHH.1001.denovo_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CHH.coverage.denovo_TSS.1001G_data.bed")
rownames(CHH.1001.denovo)<-CHH.1001.denovo$transcript
rownames(CHH.1001.denovo_TES)<-CHH.1001.denovo_TES$transcript
rownames(CHH.1001.denovo_TSS)<-CHH.1001.denovo_TSS$transcript
CHH.1001.denovo[CHH.1001.denovo == 22] <- NA
CHH.1001.denovo_TES[CHH.1001.denovo_TES == 22] <- NA
CHH.1001.denovo_TSS[CHH.1001.denovo_TSS == 22] <- NA





CG.1001.araport<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CG.coverage.araport.1001G_data.bed")
CG.1001.araport_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CG.coverage.araport_TES.1001G_data.bed")
CG.1001.araport_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CG.coverage.araport_TSS.1001G_data.bed")
CG.1001.araport<-CG.1001.araport[!duplicated(CG.1001.araport$transcript),]
#CG.1001.araport_TES<-CG.1001.araport_TES[!duplicated(CG.1001.araport_TES$transcript),]
CG.1001.araport_TSS<-CG.1001.araport_TSS[!duplicated(CG.1001.araport_TSS$transcript),]



rownames(CG.1001.araport)<-CG.1001.araport$transcript
rownames(CG.1001.araport_TES)<-CG.1001.araport_TES$transcript
rownames(CG.1001.araport_TSS)<-CG.1001.araport_TSS$transcript
CG.1001.araport[CG.1001.araport == 22] <- NA
CG.1001.araport_TES[CG.1001.araport_TES == 22] <- NA
CG.1001.araport_TSS[CG.1001.araport_TSS == 22] <- NA

CHG.1001.araport<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CHG.coverage.araport.1001G_data.bed")
CHG.1001.araport_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CHG.coverage.araport_TES.1001G_data.bed")
CHG.1001.araport_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CHG.coverage.araport_TSS.1001G_data.bed")
CHG.1001.araport<-CHG.1001.araport[!duplicated(CHG.1001.araport$transcript),]
#CHG.1001.araport_TES<-CHG.1001.araport_TES[!duplicated(CHG.1001.araport_TES$transcript),]
CHG.1001.araport_TSS<-CHG.1001.araport_TSS[!duplicated(CHG.1001.araport_TSS$transcript),]


rownames(CHG.1001.araport)<-CHG.1001.araport$transcript
rownames(CHG.1001.araport_TES)<-CHG.1001.araport_TES$transcript
rownames(CHG.1001.araport_TSS)<-CHG.1001.araport_TSS$transcript
CHG.1001.araport[CHG.1001.araport == 22] <- NA
CHG.1001.araport_TES[CHG.1001.araport_TES == 22] <- NA
CHG.1001.araport_TSS[CHG.1001.araport_TSS == 22] <- NA



CHH.1001.araport<- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CHH.coverage.araport.1001G_data.bed")
CHH.1001.araport_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CHH.coverage.araport_TES.1001G_data.bed")
CHH.1001.araport_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2021_methylation/coverage/CHH.coverage.araport_TSS.1001G_data.bed")
CHH.1001.araport<-CHH.1001.araport[!duplicated(CHH.1001.araport$transcript),]
#CHH.1001.araport_TES<-CHH.1001.araport_TES[!duplicated(CHH.1001.araport_TES$transcript),]
CHH.1001.araport_TSS<-CHH.1001.araport_TSS[!duplicated(CHH.1001.araport_TSS$transcript),]

rownames(CHH.1001.araport)<-CHH.1001.araport$transcript
rownames(CHH.1001.araport_TES)<-CHH.1001.araport_TES$transcript
rownames(CHH.1001.araport_TSS)<-CHH.1001.araport_TSS$transcript
CHH.1001.araport[CHH.1001.araport == 22] <- NA
CHH.1001.araport_TES[CHH.1001.araport_TES == 22] <- NA
CHH.1001.araport_TSS[CHH.1001.araport_TSS == 22] <- NA




a<-CG.1001.denovo
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CG.1001.denovo<-a

a<-CG.1001.denovo_TES
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CG.1001.denovo_TES<-a

a<-CG.1001.denovo_TSS
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CG.1001.denovo_TSS<-a


a<-CHG.1001.denovo
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CHG.1001.denovo<-a

a<-CHG.1001.denovo_TES
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CHG.1001.denovo_TES<-a

a<-CHG.1001.denovo_TSS
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CHG.1001.denovo_TSS<-a


a<-CHH.1001.denovo
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CHH.1001.denovo<-a

a<-CHH.1001.denovo_TES
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CHH.1001.denovo_TES<-a

a<-CHH.1001.denovo_TSS
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CHH.1001.denovo_TSS<-a



a<-CG.1001.araport
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CG.1001.araport<-a

a<-CG.1001.araport_TES
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CG.1001.araport_TES<-a

a<-CG.1001.araport_TSS
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CG.1001.araport_TSS<-a


a<-CHG.1001.araport
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CHG.1001.araport<-a

a<-CHG.1001.araport_TES
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CHG.1001.araport_TES<-a

a<-CHG.1001.araport_TSS
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CHG.1001.araport_TSS<-a


a<-CHH.1001.araport
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CHH.1001.araport<-a

a<-CHH.1001.araport_TES
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CHH.1001.araport_TES<-a

a<-CHH.1001.araport_TSS
a$mean<-apply(a[,7:450],1,mean)
a$max<-apply(a[,7:450],1,max)
a$range<-apply(a[,7:450],1,max)-apply(a[,7:450],1,min)
a$range_percent_of_mean<-100*(apply(a[,7:450],1,max)-apply(a[,7:450],1,min))/apply(a[,7:450],1,mean)
a$sd<-apply(a[,7:450],1,sd)
a$variance<-apply(a[,7:450],1,sd)/apply(a[,7:450],1,mean)
CHH.1001.araport_TSS<-a





#######################

#1001G additional samples

CG.1001new.denovo<- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CG.coverage.denovo.additional1001G_data.bed")
CG.1001new.denovo_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CG.coverage.denovo_TES.additional1001G_data.bed")
CG.1001new.denovo_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CG.coverage.denovo_TSS.additional1001G_data.bed")
CG.1001new.denovo<-CG.1001new.denovo[!duplicated(CG.1001new.denovo$transcript),]
CG.1001new.denovo_TSS<-CG.1001new.denovo_TSS[!duplicated(CG.1001new.denovo_TSS$transcript),]

rownames(CG.1001new.denovo)<-CG.1001new.denovo$transcript
rownames(CG.1001new.denovo_TES)<-CG.1001new.denovo_TES$transcript
rownames(CG.1001new.denovo_TSS)<-CG.1001new.denovo_TSS$transcript
CG.1001new.denovo[CG.1001new.denovo == 22] <- NA
CG.1001new.denovo_TES[CG.1001new.denovo_TES == 22] <- NA
CG.1001new.denovo_TSS[CG.1001new.denovo_TSS == 22] <- NA



CHG.1001new.denovo<- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CHG.coverage.denovo.additional1001G_data.bed")
CHG.1001new.denovo_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CHG.coverage.denovo_TES.additional1001G_data.bed")
CHG.1001new.denovo_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CHG.coverage.denovo_TSS.additional1001G_data.bed")
CHG.1001new.denovo<-CHG.1001new.denovo[!duplicated(CHG.1001new.denovo$transcript),]
CHG.1001new.denovo_TSS<-CHG.1001new.denovo_TSS[!duplicated(CHG.1001new.denovo_TSS$transcript),]

rownames(CHG.1001new.denovo)<-CHG.1001new.denovo$transcript
rownames(CHG.1001new.denovo_TES)<-CHG.1001new.denovo_TES$transcript
rownames(CHG.1001new.denovo_TSS)<-CHG.1001new.denovo_TSS$transcript
CHG.1001new.denovo[CHG.1001new.denovo == 22] <- NA
CHG.1001new.denovo_TES[CHG.1001new.denovo_TES == 22] <- NA
CHG.1001new.denovo_TSS[CHG.1001new.denovo_TSS == 22] <- NA



CHH.1001new.denovo<- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CHH.coverage.denovo.additional1001G_data.bed")
CHH.1001new.denovo_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CHH.coverage.denovo_TES.additional1001G_data.bed")
CHH.1001new.denovo_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CHH.coverage.denovo_TSS.additional1001G_data.bed")
CHH.1001new.denovo<-CHH.1001new.denovo[!duplicated(CHH.1001new.denovo$transcript),]
CHH.1001new.denovo_TSS<-CHH.1001new.denovo_TSS[!duplicated(CHH.1001new.denovo_TSS$transcript),]

rownames(CHH.1001new.denovo)<-CHH.1001new.denovo$transcript
rownames(CHH.1001new.denovo_TES)<-CHH.1001new.denovo_TES$transcript
rownames(CHH.1001new.denovo_TSS)<-CHH.1001new.denovo_TSS$transcript
CHH.1001new.denovo[CHH.1001new.denovo == 22] <- NA
CHH.1001new.denovo_TES[CHH.1001new.denovo_TES == 22] <- NA
CHH.1001new.denovo_TSS[CHH.1001new.denovo_TSS == 22] <- NA





CG.1001new.araport<- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CG.coverage.araport.additional1001G_data.bed")
CG.1001new.araport_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CG.coverage.araport_TES.additional1001G_data.bed")
CG.1001new.araport_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CG.coverage.araport_TSS.additional1001G_data.bed")
CG.1001new.araport<-CG.1001new.araport[!duplicated(CG.1001new.araport$transcript),]
CG.1001new.araport_TSS<-CG.1001new.araport_TSS[!duplicated(CG.1001new.araport_TSS$transcript),]

rownames(CG.1001new.araport)<-CG.1001new.araport$transcript
rownames(CG.1001new.araport_TES)<-CG.1001new.araport_TES$transcript
rownames(CG.1001new.araport_TSS)<-CG.1001new.araport_TSS$transcript
CG.1001new.araport[CG.1001new.araport == 22] <- NA
CG.1001new.araport_TES[CG.1001new.araport_TES == 22] <- NA
CG.1001new.araport_TSS[CG.1001new.araport_TSS == 22] <- NA

CHG.1001new.araport<- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CHG.coverage.araport.additional1001G_data.bed")
CHG.1001new.araport_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CHG.coverage.araport_TES.additional1001G_data.bed")
CHG.1001new.araport_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CHG.coverage.araport_TSS.additional1001G_data.bed")
CHG.1001new.araport<-CHG.1001new.araport[!duplicated(CHG.1001new.araport$transcript),]
CHG.1001new.araport_TSS<-CHG.1001new.araport_TSS[!duplicated(CHG.1001new.araport_TSS$transcript),]

rownames(CHG.1001new.araport)<-CHG.1001new.araport$transcript
rownames(CHG.1001new.araport_TES)<-CHG.1001new.araport_TES$transcript
rownames(CHG.1001new.araport_TSS)<-CHG.1001new.araport_TSS$transcript
CHG.1001new.araport[CHG.1001new.araport == 22] <- NA
CHG.1001new.araport_TES[CHG.1001new.araport_TES == 22] <- NA
CHG.1001new.araport_TSS[CHG.1001new.araport_TSS == 22] <- NA



CHH.1001new.araport<- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CHH.coverage.araport.additional1001G_data.bed")
CHH.1001new.araport_TES <- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CHH.coverage.araport_TES.additional1001G_data.bed")
CHH.1001new.araport_TSS <- read.delim("03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/coverage/CHH.coverage.araport_TSS.additional1001G_data.bed")
CHH.1001new.araport<-CHH.1001new.araport[!duplicated(CHH.1001new.araport$transcript),]
CHH.1001new.araport_TSS<-CHH.1001new.araport_TSS[!duplicated(CHH.1001new.araport_TSS$transcript),]

rownames(CHH.1001new.araport)<-CHH.1001new.araport$transcript
rownames(CHH.1001new.araport_TES)<-CHH.1001new.araport_TES$transcript
rownames(CHH.1001new.araport_TSS)<-CHH.1001new.araport_TSS$transcript
CHH.1001new.araport[CHH.1001new.araport == 22] <- NA
CHH.1001new.araport_TES[CHH.1001new.araport_TES == 22] <- NA
CHH.1001new.araport_TSS[CHH.1001new.araport_TSS == 22] <- NA




a<-as.data.frame(CG.1001new.denovo)

#a[is.nan(a)]<-NA

a$mean.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)/apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)

a$mean.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)/apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)

a$mean.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)/apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)

a$mean.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)/apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)

a$mean.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)/apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)

a$mean.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)/apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)

a$mean.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)/apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)

a$mean.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)/apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)

a$mean.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)/apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)

a$mean.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)/apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)

a$mean.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)/apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)

a$mean.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)/apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)

a$mean.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)/apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)

a$mean.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)/apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)

a$mean.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)/apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)

a$mean.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)/apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)

a$mean.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)/apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)

a$mean.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)/apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)

a$mean.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)/apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)

a$mean.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)/apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)

a$mean.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)/apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)

a$mean.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)/apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)

a$mean.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)/apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)

a$mean.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)/apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)

a$mean.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)/apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)

a$mean.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)/apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)

a$mean.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)

a$mean.9470<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)

a$mean_intravariance<-apply(a[,grep("var.",names(a))],1,mean,na.rm = TRUE)
a$mean_intra_sd<-apply(a[,grep("sd.",names(a))],1,mean,na.rm = TRUE)

b<-a[,100:182]
a$mean_of_means<-apply(b[,grep("mean.",names(b))],1,mean,na.rm = TRUE)
a$variance_of_means<-apply(b[,grep("mean.",names(b))],1,sd)/a$mean_of_means
a$sd_of_means<-apply(b[,grep("mean.",names(b))],1,sd)

a$range_all<-apply(a[,7:99],1,max)-apply(a[,7:99],1,min)
a$range_all_percent_of_mean<-100*(apply(a[,7:99],1,max)-apply(a[,7:99],1,min))/apply(a[,7:99],1,mean)
a$range_means<-apply(b[,grep("mean.",names(b))],1,max)-apply(b[,grep("mean.",names(b))],1,min)
a$sd<-apply(a[,7:99],1,sd)
a$variance<-apply(a[,7:99],1,sd)/apply(a[,7:99],1,mean)
CG.1001new.denovo<-a


a<-as.data.frame(CG.1001new.denovo_TES)

#a[is.nan(a)]<-NA

a$mean.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)/apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)

a$mean.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)/apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)

a$mean.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)/apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)

a$mean.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)/apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)

a$mean.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)/apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)

a$mean.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)/apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)

a$mean.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)/apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)

a$mean.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)/apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)

a$mean.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)/apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)

a$mean.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)/apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)

a$mean.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)/apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)

a$mean.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)/apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)

a$mean.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)/apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)

a$mean.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)/apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)

a$mean.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)/apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)

a$mean.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)/apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)

a$mean.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)/apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)

a$mean.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)/apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)

a$mean.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)/apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)

a$mean.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)/apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)

a$mean.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)/apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)

a$mean.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)/apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)

a$mean.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)/apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)

a$mean.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)/apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)

a$mean.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)/apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)

a$mean.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)/apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)

a$mean.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)

a$mean.9470<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)

a$mean_intravariance<-apply(a[,grep("var.",names(a))],1,mean,na.rm = TRUE)
a$mean_intra_sd<-apply(a[,grep("sd.",names(a))],1,mean,na.rm = TRUE)

b<-a[,100:182]
a$mean_of_means<-apply(b[,grep("mean.",names(b))],1,mean,na.rm = TRUE)
a$variance_of_means<-apply(b[,grep("mean.",names(b))],1,sd)/a$mean_of_means
a$sd_of_means<-apply(b[,grep("mean.",names(b))],1,sd)

a$range_all<-apply(a[,7:99],1,max)-apply(a[,7:99],1,min)
a$range_all_percent_of_mean<-100*(apply(a[,7:99],1,max)-apply(a[,7:99],1,min))/apply(a[,7:99],1,mean)
a$range_means<-apply(b[,grep("mean.",names(b))],1,max)-apply(b[,grep("mean.",names(b))],1,min)
a$sd<-apply(a[,7:99],1,sd)
a$variance<-apply(a[,7:99],1,sd)/apply(a[,7:99],1,mean)
CG.1001new.denovo_TES<-a


a<-as.data.frame(CG.1001new.denovo_TSS)

#a[is.nan(a)]<-NA

a$mean.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)/apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)

a$mean.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)/apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)

a$mean.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)/apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)

a$mean.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)/apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)

a$mean.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)/apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)

a$mean.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)/apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)

a$mean.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)/apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)

a$mean.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)/apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)

a$mean.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)/apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)

a$mean.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)/apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)

a$mean.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)/apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)

a$mean.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)/apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)

a$mean.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)/apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)

a$mean.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)/apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)

a$mean.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)/apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)

a$mean.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)/apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)

a$mean.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)/apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)

a$mean.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)/apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)

a$mean.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)/apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)

a$mean.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)/apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)

a$mean.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)/apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)

a$mean.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)/apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)

a$mean.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)/apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)

a$mean.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)/apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)

a$mean.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)/apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)

a$mean.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)/apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)

a$mean.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)

a$mean.9470<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)

a$mean_intravariance<-apply(a[,grep("var.",names(a))],1,mean,na.rm = TRUE)
a$mean_intra_sd<-apply(a[,grep("sd.",names(a))],1,mean,na.rm = TRUE)

b<-a[,100:182]
a$mean_of_means<-apply(b[,grep("mean.",names(b))],1,mean,na.rm = TRUE)
a$variance_of_means<-apply(b[,grep("mean.",names(b))],1,sd)/a$mean_of_means
a$sd_of_means<-apply(b[,grep("mean.",names(b))],1,sd)

a$range_all<-apply(a[,7:99],1,max)-apply(a[,7:99],1,min)
a$range_all_percent_of_mean<-100*(apply(a[,7:99],1,max)-apply(a[,7:99],1,min))/apply(a[,7:99],1,mean)
a$range_means<-apply(b[,grep("mean.",names(b))],1,max)-apply(b[,grep("mean.",names(b))],1,min)
a$sd<-apply(a[,7:99],1,sd)
a$variance<-apply(a[,7:99],1,sd)/apply(a[,7:99],1,mean)
CG.1001new.denovo_TSS<-a





a<-as.data.frame(CG.1001new.araport)

#a[is.nan(a)]<-NA

a$mean.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)/apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)

a$mean.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)/apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)

a$mean.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)/apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)

a$mean.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)/apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)

a$mean.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)/apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)

a$mean.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)/apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)

a$mean.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)/apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)

a$mean.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)/apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)

a$mean.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)/apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)

a$mean.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)/apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)

a$mean.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)/apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)

a$mean.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)/apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)

a$mean.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)/apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)

a$mean.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)/apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)

a$mean.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)/apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)

a$mean.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)/apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)

a$mean.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)/apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)

a$mean.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)/apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)

a$mean.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)/apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)

a$mean.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)/apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)

a$mean.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)/apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)

a$mean.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)/apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)

a$mean.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)/apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)

a$mean.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)/apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)

a$mean.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)/apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)

a$mean.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)/apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)

a$mean.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)

a$mean.9470<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)

a$mean_intravariance<-apply(a[,grep("var.",names(a))],1,mean,na.rm = TRUE)
a$mean_intra_sd<-apply(a[,grep("sd.",names(a))],1,mean,na.rm = TRUE)

b<-a[,100:182]
a$mean_of_means<-apply(b[,grep("mean.",names(b))],1,mean,na.rm = TRUE)
a$variance_of_means<-apply(b[,grep("mean.",names(b))],1,sd)/a$mean_of_means
a$sd_of_means<-apply(b[,grep("mean.",names(b))],1,sd)

a$range_all<-apply(a[,7:99],1,max)-apply(a[,7:99],1,min)
a$range_all_percent_of_mean<-100*(apply(a[,7:99],1,max)-apply(a[,7:99],1,min))/apply(a[,7:99],1,mean)
a$range_means<-apply(b[,grep("mean.",names(b))],1,max)-apply(b[,grep("mean.",names(b))],1,min)
a$sd<-apply(a[,7:99],1,sd)
a$variance<-apply(a[,7:99],1,sd)/apply(a[,7:99],1,mean)
CG.1001new.araport<-a




a<-as.data.frame(CHG.1001new.denovo)

#a[is.nan(a)]<-NA

a$mean.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)/apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)

a$mean.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)/apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)

a$mean.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)/apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)

a$mean.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)/apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)

a$mean.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)/apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)

a$mean.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)/apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)

a$mean.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)/apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)

a$mean.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)/apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)

a$mean.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)/apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)

a$mean.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)/apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)

a$mean.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)/apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)

a$mean.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)/apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)

a$mean.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)/apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)

a$mean.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)/apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)

a$mean.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)/apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)

a$mean.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)/apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)

a$mean.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)/apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)

a$mean.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)/apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)

a$mean.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)/apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)

a$mean.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)/apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)

a$mean.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)/apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)

a$mean.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)/apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)

a$mean.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)/apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)

a$mean.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)/apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)

a$mean.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)/apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)

a$mean.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)/apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)

a$mean.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)

a$mean.9470<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)

a$mean_intravariance<-apply(a[,grep("var.",names(a))],1,mean,na.rm = TRUE)
a$mean_intra_sd<-apply(a[,grep("sd.",names(a))],1,mean,na.rm = TRUE)

b<-a[,100:182]
a$mean_of_means<-apply(b[,grep("mean.",names(b))],1,mean,na.rm = TRUE)
a$variance_of_means<-apply(b[,grep("mean.",names(b))],1,sd)/a$mean_of_means
a$sd_of_means<-apply(b[,grep("mean.",names(b))],1,sd)

a$range_all<-apply(a[,7:99],1,max)-apply(a[,7:99],1,min)
a$range_all_percent_of_mean<-100*(apply(a[,7:99],1,max)-apply(a[,7:99],1,min))/apply(a[,7:99],1,mean)
a$range_means<-apply(b[,grep("mean.",names(b))],1,max)-apply(b[,grep("mean.",names(b))],1,min)
a$sd<-apply(a[,7:99],1,sd)
a$variance<-apply(a[,7:99],1,sd)/apply(a[,7:99],1,mean)
CHG.1001new.denovo<-a



a<-as.data.frame(CHG.1001new.denovo_TES)

#a[is.nan(a)]<-NA

a$mean.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)/apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)

a$mean.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)/apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)

a$mean.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)/apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)

a$mean.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)/apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)

a$mean.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)/apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)

a$mean.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)/apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)

a$mean.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)/apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)

a$mean.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)/apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)

a$mean.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)/apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)

a$mean.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)/apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)

a$mean.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)/apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)

a$mean.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)/apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)

a$mean.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)/apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)

a$mean.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)/apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)

a$mean.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)/apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)

a$mean.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)/apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)

a$mean.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)/apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)

a$mean.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)/apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)

a$mean.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)/apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)

a$mean.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)/apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)

a$mean.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)/apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)

a$mean.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)/apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)

a$mean.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)/apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)

a$mean.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)/apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)

a$mean.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)/apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)

a$mean.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)/apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)

a$mean.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)

a$mean.9470<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)

a$mean_intravariance<-apply(a[,grep("var.",names(a))],1,mean,na.rm = TRUE)
a$mean_intra_sd<-apply(a[,grep("sd.",names(a))],1,mean,na.rm = TRUE)

b<-a[,100:182]
a$mean_of_means<-apply(b[,grep("mean.",names(b))],1,mean,na.rm = TRUE)
a$variance_of_means<-apply(b[,grep("mean.",names(b))],1,sd)/a$mean_of_means
a$sd_of_means<-apply(b[,grep("mean.",names(b))],1,sd)

a$range_all<-apply(a[,7:99],1,max)-apply(a[,7:99],1,min)
a$range_all_percent_of_mean<-100*(apply(a[,7:99],1,max)-apply(a[,7:99],1,min))/apply(a[,7:99],1,mean)
a$range_means<-apply(b[,grep("mean.",names(b))],1,max)-apply(b[,grep("mean.",names(b))],1,min)
a$sd<-apply(a[,7:99],1,sd)
a$variance<-apply(a[,7:99],1,sd)/apply(a[,7:99],1,mean)
CHG.1001new.denovo_TES<-a


a<-as.data.frame(CHG.1001new.denovo_TSS)

#a[is.nan(a)]<-NA

a$mean.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)/apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)

a$mean.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)/apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)

a$mean.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)/apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)

a$mean.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)/apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)

a$mean.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)/apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)

a$mean.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)/apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)

a$mean.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)/apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)

a$mean.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)/apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)

a$mean.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)/apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)

a$mean.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)/apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)

a$mean.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)/apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)

a$mean.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)/apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)

a$mean.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)/apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)

a$mean.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)/apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)

a$mean.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)/apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)

a$mean.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)/apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)

a$mean.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)/apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)

a$mean.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)/apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)

a$mean.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)/apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)

a$mean.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)/apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)

a$mean.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)/apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)

a$mean.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)/apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)

a$mean.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)/apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)

a$mean.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)/apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)

a$mean.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)/apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)

a$mean.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)/apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)

a$mean.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)

a$mean.9470<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)

a$mean_intravariance<-apply(a[,grep("var.",names(a))],1,mean,na.rm = TRUE)
a$mean_intra_sd<-apply(a[,grep("sd.",names(a))],1,mean,na.rm = TRUE)

b<-a[,100:182]
a$mean_of_means<-apply(b[,grep("mean.",names(b))],1,mean,na.rm = TRUE)
a$variance_of_means<-apply(b[,grep("mean.",names(b))],1,sd)/a$mean_of_means
a$sd_of_means<-apply(b[,grep("mean.",names(b))],1,sd)

a$range_all<-apply(a[,7:99],1,max)-apply(a[,7:99],1,min)
a$range_all_percent_of_mean<-100*(apply(a[,7:99],1,max)-apply(a[,7:99],1,min))/apply(a[,7:99],1,mean)
a$range_means<-apply(b[,grep("mean.",names(b))],1,max)-apply(b[,grep("mean.",names(b))],1,min)
a$sd<-apply(a[,7:99],1,sd)
a$variance<-apply(a[,7:99],1,sd)/apply(a[,7:99],1,mean)
CHG.1001new.denovo_TSS<-a



a<-as.data.frame(CHH.1001new.denovo)

#a[is.nan(a)]<-NA

a$mean.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)/apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)

a$mean.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)/apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)

a$mean.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)/apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)

a$mean.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)/apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)

a$mean.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)/apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)

a$mean.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)/apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)

a$mean.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)/apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)

a$mean.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)/apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)

a$mean.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)/apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)

a$mean.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)/apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)

a$mean.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)/apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)

a$mean.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)/apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)

a$mean.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)/apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)

a$mean.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)/apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)

a$mean.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)/apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)

a$mean.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)/apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)

a$mean.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)/apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)

a$mean.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)/apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)

a$mean.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)/apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)

a$mean.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)/apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)

a$mean.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)/apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)

a$mean.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)/apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)

a$mean.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)/apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)

a$mean.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)/apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)

a$mean.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)/apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)

a$mean.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)/apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)

a$mean.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)

a$mean.9470<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)

a$mean_intravariance<-apply(a[,grep("var.",names(a))],1,mean,na.rm = TRUE)
a$mean_intra_sd<-apply(a[,grep("sd.",names(a))],1,mean,na.rm = TRUE)

b<-a[,100:182]
a$mean_of_means<-apply(b[,grep("mean.",names(b))],1,mean,na.rm = TRUE)
a$variance_of_means<-apply(b[,grep("mean.",names(b))],1,sd)/a$mean_of_means
a$sd_of_means<-apply(b[,grep("mean.",names(b))],1,sd)

a$range_all<-apply(a[,7:99],1,max)-apply(a[,7:99],1,min)
a$range_all_percent_of_mean<-100*(apply(a[,7:99],1,max)-apply(a[,7:99],1,min))/apply(a[,7:99],1,mean)
a$range_means<-apply(b[,grep("mean.",names(b))],1,max)-apply(b[,grep("mean.",names(b))],1,min)
a$sd<-apply(a[,7:99],1,sd)
a$variance<-apply(a[,7:99],1,sd)/apply(a[,7:99],1,mean)
CHH.1001new.denovo<-a


a<-as.data.frame(CHH.1001new.denovo_TES)

#a[is.nan(a)]<-NA

a$mean.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)/apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)

a$mean.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)/apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)

a$mean.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)/apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)

a$mean.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)/apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)

a$mean.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)/apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)

a$mean.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)/apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)

a$mean.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)/apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)

a$mean.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)/apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)

a$mean.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)/apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)

a$mean.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)/apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)

a$mean.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)/apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)

a$mean.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)/apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)

a$mean.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)/apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)

a$mean.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)/apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)

a$mean.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)/apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)

a$mean.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)/apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)

a$mean.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)/apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)

a$mean.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)/apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)

a$mean.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)/apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)

a$mean.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)/apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)

a$mean.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)/apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)

a$mean.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)/apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)

a$mean.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)/apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)

a$mean.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)/apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)

a$mean.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)/apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)

a$mean.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)/apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)

a$mean.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)

a$mean.9470<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)

a$mean_intravariance<-apply(a[,grep("var.",names(a))],1,mean,na.rm = TRUE)
a$mean_intra_sd<-apply(a[,grep("sd.",names(a))],1,mean,na.rm = TRUE)

b<-a[,100:182]
a$mean_of_means<-apply(b[,grep("mean.",names(b))],1,mean,na.rm = TRUE)
a$variance_of_means<-apply(b[,grep("mean.",names(b))],1,sd)/a$mean_of_means
a$sd_of_means<-apply(b[,grep("mean.",names(b))],1,sd)

a$range_all<-apply(a[,7:99],1,max)-apply(a[,7:99],1,min)
a$range_all_percent_of_mean<-100*(apply(a[,7:99],1,max)-apply(a[,7:99],1,min))/apply(a[,7:99],1,mean)
a$range_means<-apply(b[,grep("mean.",names(b))],1,max)-apply(b[,grep("mean.",names(b))],1,min)
a$sd<-apply(a[,7:99],1,sd)
a$variance<-apply(a[,7:99],1,sd)/apply(a[,7:99],1,mean)
CHH.1001new.denovo_TES<-a



a<-as.data.frame(CHH.1001new.denovo_TSS)

#a[is.nan(a)]<-NA

a$mean.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)/apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)

a$mean.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)/apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)

a$mean.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)/apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)

a$mean.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)/apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)

a$mean.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)/apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)

a$mean.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)/apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)

a$mean.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)/apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)

a$mean.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)/apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)

a$mean.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)/apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)

a$mean.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)/apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)

a$mean.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)/apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)

a$mean.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)/apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)

a$mean.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)/apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)

a$mean.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)/apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)

a$mean.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)/apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)

a$mean.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)/apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)

a$mean.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)/apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)

a$mean.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)/apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)

a$mean.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)/apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)

a$mean.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)/apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)

a$mean.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)/apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)

a$mean.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)/apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)

a$mean.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)/apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)

a$mean.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)/apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)

a$mean.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)/apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)

a$mean.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)/apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)

a$mean.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)

a$mean.9470<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)

a$mean_intravariance<-apply(a[,grep("var.",names(a))],1,mean,na.rm = TRUE)
a$mean_intra_sd<-apply(a[,grep("sd.",names(a))],1,mean,na.rm = TRUE)

b<-a[,100:182]
a$mean_of_means<-apply(b[,grep("mean.",names(b))],1,mean,na.rm = TRUE)
a$variance_of_means<-apply(b[,grep("mean.",names(b))],1,sd)/a$mean_of_means
a$sd_of_means<-apply(b[,grep("mean.",names(b))],1,sd)

a$range_all<-apply(a[,7:99],1,max)-apply(a[,7:99],1,min)
a$range_all_percent_of_mean<-100*(apply(a[,7:99],1,max)-apply(a[,7:99],1,min))/apply(a[,7:99],1,mean)
a$range_means<-apply(b[,grep("mean.",names(b))],1,max)-apply(b[,grep("mean.",names(b))],1,min)
a$sd<-apply(a[,7:99],1,sd)
a$variance<-apply(a[,7:99],1,sd)/apply(a[,7:99],1,mean)
CHH.1001new.denovo_TSS<-a



a<-as.data.frame(CHH.1001new.araport)

#a[is.nan(a)]<-NA

a$mean.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)/apply(a[,grep("1741",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1741<-apply(a[,grep("1741",names(a[,1:99]))],1,sd)

a$mean.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)/apply(a[,grep("4807",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.4807<-apply(a[,grep("4807",names(a[,1:99]))],1,sd)

a$mean.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)/apply(a[,grep("5210",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5210<-apply(a[,grep("5210",names(a[,1:99]))],1,sd)

a$mean.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)/apply(a[,grep("5772",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5772<-apply(a[,grep("5772",names(a[,1:99]))],1,sd)

a$mean.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)/apply(a[,grep("5784",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5784<-apply(a[,grep("5784",names(a[,1:99]))],1,sd)

a$mean.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)/apply(a[,grep("5856",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.5856<-apply(a[,grep("5856",names(a[,1:99]))],1,sd)

a$mean.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)/apply(a[,grep("6021",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6021<-apply(a[,grep("6021",names(a[,1:99]))],1,sd)

a$mean.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)/apply(a[,grep("6220",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6220<-apply(a[,grep("6220",names(a[,1:99]))],1,sd)

a$mean.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)/apply(a[,grep("6909",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6909<-apply(a[,grep("6909",names(a[,1:99]))],1,sd)

a$mean.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)/apply(a[,grep("6911",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6911<-apply(a[,grep("6911",names(a[,1:99]))],1,sd)

a$mean.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)/apply(a[,grep("6966",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6966<-apply(a[,grep("6966",names(a[,1:99]))],1,sd)

a$mean.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)/apply(a[,grep("8244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8244<-apply(a[,grep("8244",names(a[,1:99]))],1,sd)

a$mean.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)/apply(a[,grep("8366",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.8366<-apply(a[,grep("8366",names(a[,1:99]))],1,sd)

a$mean.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)/apply(a[,grep("9518",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9518<-apply(a[,grep("9518",names(a[,1:99]))],1,sd)

a$mean.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)/apply(a[,grep("9588",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9588<-apply(a[,grep("9588",names(a[,1:99]))],1,sd)

a$mean.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)/apply(a[,grep("9888",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9888<-apply(a[,grep("9888",names(a[,1:99]))],1,sd)

a$mean.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)/apply(a[,grep("9905",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9905<-apply(a[,grep("9905",names(a[,1:99]))],1,sd)

a$mean.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)/apply(a[,grep("10012",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.10012<-apply(a[,grep("10012",names(a[,1:99]))],1,sd)

a$mean.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)/apply(a[,grep("1254",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.1254<-apply(a[,grep("1254",names(a[,1:99]))],1,sd)

a$mean.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)/apply(a[,grep("6024",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6024<-apply(a[,grep("6024",names(a[,1:99]))],1,sd)

a$mean.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)/apply(a[,grep("6069",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6069<-apply(a[,grep("6069",names(a[,1:99]))],1,sd)

a$mean.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)/apply(a[,grep("6076",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6076<-apply(a[,grep("6076",names(a[,1:99]))],1,sd)

a$mean.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)/apply(a[,grep("6184",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6184<-apply(a[,grep("6184",names(a[,1:99]))],1,sd)

a$mean.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)/apply(a[,grep("6189",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6189<-apply(a[,grep("6189",names(a[,1:99]))],1,sd)

a$mean.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)/apply(a[,grep("6244",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.6244<-apply(a[,grep("6244",names(a[,1:99]))],1,sd)

a$mean.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)/apply(a[,grep("9057",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9057<-apply(a[,grep("9057",names(a[,1:99]))],1,sd)

a$mean.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9412<-apply(a[,grep("9412",names(a[,1:99]))],1,sd)

a$mean.9470<-apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$var.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)/apply(a[,grep("9412",names(a[,1:99]))],1,mean,na.rm = TRUE)
a$sd.9470<-apply(a[,grep("9470",names(a[,1:99]))],1,sd,na.rm = TRUE)

a$mean_intravariance<-apply(a[,grep("var.",names(a))],1,mean,na.rm = TRUE)
a$mean_intra_sd<-apply(a[,grep("sd.",names(a))],1,mean,na.rm = TRUE)

b<-a[,100:182]
a$mean_of_means<-apply(b[,grep("mean.",names(b))],1,mean,na.rm = TRUE)
a$variance_of_means<-apply(b[,grep("mean.",names(b))],1,sd)/a$mean_of_means
a$sd_of_means<-apply(b[,grep("mean.",names(b))],1,sd)

a$range_all<-apply(a[,7:99],1,max)-apply(a[,7:99],1,min)
a$range_all_percent_of_mean<-100*(apply(a[,7:99],1,max)-apply(a[,7:99],1,min))/apply(a[,7:99],1,mean)
a$range_means<-apply(b[,grep("mean.",names(b))],1,max)-apply(b[,grep("mean.",names(b))],1,min)
a$sd<-apply(a[,7:99],1,sd)
a$variance<-apply(a[,7:99],1,sd)/apply(a[,7:99],1,mean)
CHH.1001new.araport<-a



a<-CG.1001.denovo
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CG.1001.denovo<-a



a<-CG.1001.denovo_TSS
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CG.1001.denovo_TSS<-a

a<-CG.1001.denovo_TES
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CG.1001.denovo_TES<-a


a<-CHG.1001.denovo
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CHG.1001.denovo<-a



a<-CHG.1001.denovo_TSS
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CHG.1001.denovo_TSS<-a

a<-CHG.1001.denovo_TES
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CHG.1001.denovo_TES<-a


a<-CHH.1001.denovo
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CHH.1001.denovo<-a



a<-CHH.1001.denovo_TSS
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CHH.1001.denovo_TSS<-a

a<-CHH.1001.denovo_TES
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CHH.1001.denovo_TES<-a



a<-CG.1001new.denovo
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CG.1001new.denovo<-a




a<-CG.1001new.denovo_TSS
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CG.1001new.denovo_TSS<-a





a<-CG.1001new.denovo_TES
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CG.1001new.denovo_TES<-a





a<-CHG.1001new.denovo
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CHG.1001new.denovo<-a


a<-CHG.1001new.denovo_TSS
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CHG.1001new.denovo_TSS<-a

a<-CHG.1001new.denovo_TES
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CHG.1001new.denovo_TES<-a





a<-CHH.1001new.denovo
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CHH.1001new.denovo<-a


a<-CHH.1001new.denovo_TSS
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CHH.1001new.denovo_TSS<-a

a<-CHH.1001new.denovo_TES
a$genetype<-"other"
a$genetype[a$transcript %in% denovoPC.loci$gene]<-"PC"
a$genetype[a$transcript %in% lncRNAs.antisense.loci$gene]<-"AS"
a$genetype[a$transcript %in% lncRNAs.intergenic.loci$gene]<-"linc"
a$genetype[a$transcript %in% TE_frags.transcripts$gene]<-"TEfrag"
a$genetype[a$transcript %in% TE_genes.loci$gene]<-"TEgen"
a$genetype[a$transcript %in% lncRNAs.AS_to_pseudo.loci$gene]<-"ASpseudo"
a$genetype[a$transcript %in% lncRNAs.AS_to_TE.loci$gene]<-"AS_to_TE"
CHH.1001new.denovo_TES<-a









#find correlates of lncRNA expression
CHH.linc<-CHH.1001.denovo[CHH.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene,]
CHG.linc<-CHG.1001.denovo[CHG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene,]
CG.linc<-CG.1001.denovo[CG.1001.denovo$transcript %in% lncRNAs.intergenic.loci$gene,]

CHH.linc_TSS<-CHH.1001.denovo_TSS[CHH.1001.denovo_TSS$transcript %in% lncRNAs.intergenic.loci$gene,]
CHG.linc_TSS<-CHG.1001.denovo_TSS[CHG.1001.denovo_TSS$transcript %in% lncRNAs.intergenic.loci$gene,]
CG.linc_TSS<-CG.1001.denovo_TSS[CG.1001.denovo_TSS$transcript %in% lncRNAs.intergenic.loci$gene,]

CHH.linc_TES<-CHH.1001.denovo_TES[CHH.1001.denovo_TES$gene %in% lncRNAs.intergenic.loci$gene,]
CHG.linc_TES<-CHG.1001.denovo_TES[CHG.1001.denovo_TES$gene %in% lncRNAs.intergenic.loci$gene,]
CG.linc_TES<-CG.1001.denovo_TES[CG.1001.denovo_TES$gene %in% lncRNAs.intergenic.loci$gene,]

CHH.as<-CHH.1001.denovo[CHH.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene,]
#CHG.as<-CHG.1001.denovo[CHG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene,]
CG.as<-CG.1001.denovo[CG.1001.denovo$transcript %in% lncRNAs.antisense.loci$gene,]

CHH.as_TSS<-CHH.1001.denovo_TSS[CHH.1001.denovo_TSS$transcript %in% lncRNAs.antisense.loci$gene,]
#CHG.as_TSS<-CHG.1001.denovo_TSS[CHG.1001.denovo_TSS$transcript %in% lncRNAs.antisense.loci$gene,]
CG.as_TSS<-CG.1001.denovo_TSS[CG.1001.denovo_TSS$transcript %in% lncRNAs.antisense.loci$gene,]

CHH.as_TES<-CHH.1001.denovo_TES[CHH.1001.denovo_TES$gene %in% lncRNAs.antisense.loci$gene,]
#CHG.as_TES<-CHG.1001.denovo_TES[CHG.1001.denovo_TES$gene %in% lncRNAs.antisense.loci$gene,]
CG.as_TES<-CG.1001.denovo_TES[CG.1001.denovo_TES$gene %in% lncRNAs.antisense.loci$gene,]





