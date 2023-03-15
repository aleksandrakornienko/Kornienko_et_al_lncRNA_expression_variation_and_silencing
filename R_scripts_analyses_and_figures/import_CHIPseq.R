setwd("/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/")
setwd("Z:/01_POSTDOC/")

names_chipseqsamples<-read.delim("03_Projects/ChIP-seq_data/names_usable_chipsamples.txt",header = F)
names_chip_r<-read.delim("03_Projects/ChIP-seq_data/names_usable_chipsamples_rCodes.txt",header = F)

All_ChIP_data_samples_readN_duplicN <- read.delim("03_Projects/ChIP-seq_data/All_ChIP_data_samples_readN_duplicN.txt")
All_ChIP_data_samples_readN_duplicN$mln_usable_reads<-All_ChIP_data_samples_readN_duplicN$number_of_usable_reads_unique_nondup/1000000

chip.denovo.log2_all_samples<- read.delim("03_Projects/ChIP-seq_data/2021/coverage/20220408_Chipseq_coverage.denovo.log2.bed")

chip.araport.log2_all_samples<- read.delim("03_Projects/ChIP-seq_data/2021/coverage/20220408_Chipseq_coverage.araport.log2.bed")

chip.denovo_TSS.log2_prelim<- read.delim("03_Projects/ChIP-seq_data/2021/coverage/20220408_Chipseq_coverage.denovo_TSS.log2.bed")
chip.denovo_TES.log2_prelim<- read.delim("03_Projects/ChIP-seq_data/2021/coverage/20220408_Chipseq_coverage.denovo_TES.log2.bed")
chip.denovo.log2_prelim<- read.delim("03_Projects/ChIP-seq_data/2021/coverage/20220408_Chipseq_coverage.denovo.log2.bed")
chip.araport_TSS.log2_prelim<- read.delim("03_Projects/ChIP-seq_data/2021/coverage/20220408_Chipseq_coverage.araport_TSS.log2.bed")
chip.araport_TES.log2_prelim<- read.delim("03_Projects/ChIP-seq_data/2021/coverage/20220408_Chipseq_coverage.araport_TES.log2.bed")
chip.araport.log2_prelim<- read.delim("03_Projects/ChIP-seq_data/2021/coverage/20220408_Chipseq_coverage.araport.log2.bed")




#chip.denovo_TSS.subtr<- read.delim("03_Projects/ChIP-seq_data/2021/coverage/20220408_Chipseq_coverage.denovo_TSS.subtr.bed")
#chip.denovo_TES.subtr<- read.delim("03_Projects/ChIP-seq_data/2021/coverage/20220408_Chipseq_coverage.denovo_TES.subtr.bed")
#chip.denovo.subtr<- read.delim("03_Projects/ChIP-seq_data/2021/coverage/20220408_Chipseq_coverage.denovo.subtr.bed")

#chip.araport_TSS.subtr<- read.delim("03_Projects/ChIP-seq_data/2021/coverage/20220408_Chipseq_coverage.araport_TSS.subtr.bed")
#chip.araport_TES.subtr<- read.delim("03_Projects/ChIP-seq_data/2021/coverage/20220408_Chipseq_coverage.araport_TES.subtr.bed")
#chip.araport.subtr<- read.delim("03_Projects/ChIP-seq_data/2021/coverage/20220408_Chipseq_coverage.araport.subtr.bed")


#check quality and correlations 

chipcor<-cor(chip.denovo.log2_all_samples[,7:131])




plot(k27$r10.rep2.H3K27me3,k27$r11.rep2.H3K27me3)


k27<-chip.denovo.log2_all_samples[,grep("K27",names(chip.denovo.log2_all_samples))]
k4<-chip.denovo.log2_all_samples[,grep("K4",names(chip.denovo.log2_all_samples))]
h1<-chip.denovo.log2_all_samples[,grep("H1",names(chip.denovo.log2_all_samples))]
k9<-chip.denovo.log2_all_samples[,grep("K9",names(chip.denovo.log2_all_samples))]
k36<-chip.denovo.log2_all_samples[,grep("K36",names(chip.denovo.log2_all_samples))]

cor27<-cor(k27)
corh1<-cor(h1)
cor9<-cor(k9)
cork4<-cor(k4)
cork36<-cor(k36)


cor27[apply(cor27,1,mean)<0.7,] 
#r24.rep1.H3K27me3
#r24.rep2.H3K27me3
cor27[apply(cor27,1,mean)<0.5,]
#r24.rep2.H3K27me3

corh1[apply(corh1,1,mean)<0.7,]
#r10.rep2.H1
#r21.rep1.H1
corh1[apply(corh1,1,mean)<0.5,]
#none

cor9[apply(cor9,1,mean)<0.7,]
#r10.rep2.H3K9me2
#r15.rep1.H3K9me2
cor9[apply(cor9,1,mean)<0.5,]
#r15.rep1.H3K9me2
cork4[apply(cork4,1,mean)<0.7,]
#many
cork4[apply(cork4,1,mean)<0.5,]
#r24.rep2.H3K4me3
#r28.rep2.H3K4me3

cork36[apply(cork36,1,mean)<0.7,]
cork36[apply(cork36,1,mean)<0.5,]
#none
#none


# exclude: 
#r24.rep2.H3K4me3
#r24.rep2.H3K27me3
#r15.rep1.H3K9me2
#r28.rep2.H3K4me3
#r15.rep1.H3K27me3

### 
#r16.rep1.H2K27me3 - seems to be contaminated with K9


# custom function to implement min max scaling
#################################################################
minMax <- function(x) {  (x - min(x)) / (max(x) - min(x))}
quantile90minmax <- function(x) {  (x - quantile(x,.20)) / (quantile(x,.80) - quantile(x,.20))}
##################################################################

#####################################################################
processchipseq_table <- function(a){
  a<-a[!duplicated(a$gene),]
  rownames(a)<-a$gene
  a$H1.1741<-apply(a[,c( "r2.rep1.H1" ,"r2.rep2.H1"  )],1,mean)
  a$K4.1741<-apply(a[,c( "r2.rep1.H3K4me3" ,"r2.rep2.H3K4me3"  )],1,mean)
  a$K9.1741  <-apply(a[,c( "r2.rep1.H3K9me2" ,"r2.rep2.H3K9me2"  )],1,mean)
  a$K27.1741<-apply(a[,c( "r2.rep1.H3K27me3" ,"r2.rep2.H3K27me3"  )],1,mean)
  a$K36.1741<-apply(a[,c( "r2.rep1.H3K36me3" ,"r2.rep2.H3K36me3"  )],1,mean)
  #r8
  a$H1.5784<-apply(a[,c( "r8.rep1.H1" ,"r8.rep2.H1"  )],1,mean)
  a$K4.5784<-apply(a[,c( "r8.rep1.H3K4me3" ,"r8.rep2.H3K4me3"  )],1,mean)
  a$K9.5784  <-apply(a[,c( "r8.rep1.H3K9me2" ,"r8.rep2.H3K9me2"  )],1,mean)
  a$K27.5784<-apply(a[,c( "r8.rep1.H3K27me3" ,"r8.rep2.H3K27me3"  )],1,mean)
  a$K36.5784<-apply(a[,c( "r8.rep1.H3K36me3" ,"r8.rep2.H3K36me3"  )],1,mean)
  #r10
  a$H1.5856<-apply(a[,c( "r10.rep1.H1" ,"r10.rep2.H1"  )],1,mean)
  a$K4.5856<-apply(a[,c( "r10.rep1.H3K4me3" ,"r10.rep2.H3K4me3"  )],1,mean)
  a$K9.5856  <-a[,c( "r10.rep2.H3K9me2"  )]
  a$K27.5856<-a[,c( "r10.rep2.H3K27me3"  )]
  a$K36.5856<-apply(a[,c( "r10.rep1.H3K36me3" ,"r10.rep2.H3K36me3"  )],1,mean)
  #r11
  a$H1.6021<-a[,c( "r11.rep2.H1"  )]
  a$K4.6021<-a[,c( "r11.rep2.H3K4me3"  )]
  a$K9.6021  <-a[,c( "r11.rep2.H3K9me2"  )]
  a$K27.6021<-a[,c( "r11.rep2.H3K27me3"  )]
  a$K36.6021<-a[,c( "r11.rep2.H3K36me3"  )]
  #r14
  a$H1.6909<-apply(a[,c( "r14.rep1.H1" ,"r14.rep2.H1"  )],1,mean)
  a$K4.6909<-apply(a[,c( "r14.rep1.H3K4me3" ,"r14.rep2.H3K4me3"  )],1,mean)
  a$K9.6909  <-apply(a[,c( "r14.rep1.H3K9me2" ,"r14.rep2.H3K9me2"  )],1,mean)
  a$K27.6909<-apply(a[,c( "r14.rep1.H3K27me3" ,"r14.rep2.H3K27me3"  )],1,mean)
  a$K36.6909<-apply(a[,c( "r14.rep1.H3K36me3" ,"r14.rep2.H3K36me3"  )],1,mean)
  #r15
  a$H1.6911<-a[,c( "r15.rep1.H1"  )]
  a$K4.6911<-a[,c( "r15.rep1.H3K4me3"  )]
  ##a$K9.6911  <-a[,c( "r15.rep1.H3K9me2"  )]
  ##a$K27.6911<-a[,c( "r15.rep1.H3K27me3"  )]
  a$K36.6911<-a[,c( "r15.rep1.H3K36me3"  )]
  #r16
  a$H1.6966<-apply(a[,c( "r16.rep1.H1" ,"r16.rep2.H1"  )],1,mean)
  a$K4.6966<-apply(a[,c( "r16.rep1.H3K4me3" ,"r16.rep2.H3K4me3"  )],1,mean)
  a$K9.6966  <-apply(a[,c( "r16.rep1.H3K9me2" ,"r16.rep2.H3K9me2"  )],1,mean)
  a$K27.6966<-apply(a[,c( "r16.rep1.H3K27me3" ,"r16.rep2.H3K27me3"  )],1,mean)
  a$K36.6966<-apply(a[,c( "r16.rep1.H3K36me3" ,"r16.rep2.H3K36me3"  )],1,mean)
  #r21
  a$H1.9518<-a[,c( "r21.rep1.H1"  )]
  a$K4.9518<-a[,c( "r21.rep1.H3K4me3"  )]
  a$K9.9518  <-a[,c( "r21.rep1.H3K9me2"  )]
  a$K36.9518<-a[,c( "r21.rep1.H3K36me3"  )]
  #r24
  a$H1.9888<-apply(a[,c( "r24.rep1.H1" ,"r24.rep2.H1"  )],1,mean)
  a$K4.9888<-a[,c( "r24.rep1.H3K4me3" )]
  a$K9.9888  <-apply(a[,c( "r24.rep1.H3K9me2" ,"r24.rep2.H3K9me2"  )],1,mean)
  a$K27.9888<-a[,c( "r24.rep1.H3K27me3"  )]
  a$K36.9888<-apply(a[,c( "r24.rep1.H3K36me3" ,"r24.rep2.H3K36me3"  )],1,mean)
  #r25
  a$H1.9905<-apply(a[,c( "r25.rep1.H1" ,"r25.rep2.H1"  )],1,mean)
  a$K4.9905<-apply(a[,c( "r25.rep1.H3K4me3" ,"r25.rep2.H3K4me3"  )],1,mean)
  a$K9.9905  <-apply(a[,c( "r25.rep1.H3K9me2" ,"r25.rep2.H3K9me2"  )],1,mean)
  a$K27.9905<-apply(a[,c( "r25.rep1.H3K27me3" ,"r25.rep2.H3K27me3"  )],1,mean)
  a$K36.9905<-apply(a[,c( "r25.rep1.H3K36me3" ,"r25.rep2.H3K36me3"  )],1,mean)
  #r26
  a$H1.10012<-apply(a[,c( "r26.rep1.H1" ,"r26.rep2.H1"  )],1,mean)
  a$K4.10012<-apply(a[,c( "r26.rep1.H3K4me3" ,"r26.rep2.H3K4me3"  )],1,mean)
  a$K9.10012  <-apply(a[,c( "r26.rep1.H3K9me2" ,"r26.rep2.H3K9me2"  )],1,mean)
  a$K27.10012<-a[,c( "r26.rep1.H3K27me3" )]
  a$K36.10012<-apply(a[,c( "r26.rep1.H3K36me3" ,"r26.rep2.H3K36me3"  )],1,mean)
  #r27
  a$H1.1254<-apply(a[,c( "r27.rep1.H1" ,"r27.rep2.H1"  )],1,mean)
  a$K4.1254<-apply(a[,c( "r27.rep1.H3K4me3" ,"r27.rep2.H3K4me3"  )],1,mean)
  a$K9.1254  <-a[,c( "r27.rep2.H3K9me2"  )]
  a$K27.1254<-apply(a[,c( "r27.rep1.H3K27me3" ,"r27.rep2.H3K27me3"  )],1,mean)
  a$K36.1254<-apply(a[,c( "r27.rep1.H3K36me3" ,"r27.rep2.H3K36me3"  )],1,mean)
  #r28
  a$H1.6024<-apply(a[,c( "r28.rep1.H1" ,"r28.rep2.H1"  )],1,mean)
  a$K4.6024<-apply(a[,c( "r28.rep1.H3K4me3" ,"r28.rep2.H3K4me3"  )],1,mean)
  a$K9.6024  <-apply(a[,c( "r28.rep1.H3K9me2" ,"r28.rep2.H3K9me2"  )],1,mean)
  a$K27.6024<-apply(a[,c( "r28.rep1.H3K27me3" ,"r28.rep2.H3K27me3"  )],1,mean)
  a$K36.6024<-apply(a[,c( "r28.rep1.H3K36me3" ,"r28.rep2.H3K36me3"  )],1,mean)
  #r35
  a$H1.9057 <-apply(a[,c( "r35.rep1.H1" ,"r35.rep2.H1"  )],1,mean)
  a$K4.9057 <-apply(a[,c( "r35.rep1.H3K4me3" ,"r35.rep2.H3K4me3"  )],1,mean)
  a$K9.9057 <-apply(a[,c( "r35.rep1.H3K9me2" ,"r35.rep2.H3K9me2"  )],1,mean)
  a$K27.9057 <-apply(a[,c( "r35.rep1.H3K27me3" ,"r35.rep2.H3K27me3"  )],1,mean)
  a$K36.9057 <-a[,c( "r35.rep2.H3K36me3"  )]
  a<-a[,c(1, 127:193)]
  return(a)
}
###################################################

#chip.denovo.minmax <- cbind(chip.denovo.log2$gene,as.data.frame(lapply(chip.denovo.log2[,2:68], minMax)))
#chip.denovo.minmax$gene<-chip.denovo.minmax[,1]
#rownames(chip.denovo.minmax)<-chip.denovo.minmax$gene


####################
####process data####
####################


chip.denovo.allsamples.quantstan <- cbind(chip.denovo.log2_all_samples$gene,as.data.frame(lapply(chip.denovo.log2_all_samples[,7:length(chip.denovo.log2_all_samples)], quantile90minmax)))
chip.denovo.allsamples.quantstan$gene<-chip.denovo.allsamples.quantstan[,1]
rownames(chip.denovo.allsamples.quantstan)<-chip.denovo.allsamples.quantstan$gene
chip.denovo.allsamples.quantstan<-chip.denovo.allsamples.quantstan[,c(127,2:126)]



chip.araport.allsamples.quantstan <- cbind(chip.araport.log2_all_samples$gene,as.data.frame(lapply(chip.araport.log2_all_samples[,7:length(chip.araport.log2_all_samples)], quantile90minmax)))
chip.araport.allsamples.quantstan<-chip.araport.allsamples.quantstan[chip.araport.allsamples.quantstan[,1] %in% Araport11_protein_coding.201606.genes$gene,]
chip.araport.allsamples.quantstan$gene<-chip.araport.allsamples.quantstan[,1]
chip.araport.allsamples.quantstan<-chip.araport.allsamples.quantstan[!duplicated(chip.araport.allsamples.quantstan$gene),]
rownames(chip.araport.allsamples.quantstan)<-chip.araport.allsamples.quantstan$gene
chip.araport.allsamples.quantstan<-chip.araport.allsamples.quantstan[,c(127,2:126)]
rm(chip.araport.log2_all_samples)



chip.denovo.quantstan <- cbind(chip.denovo.log2_prelim$gene,as.data.frame(lapply(chip.denovo.log2_prelim[,7:131], quantile90minmax)))
chip.denovo.quantstan$gene<-chip.denovo.quantstan[,1]
rownames(chip.denovo.quantstan)<-chip.denovo.quantstan$gene
chip.denovo.quantstan<-processchipseq_table(chip.denovo.quantstan[,c(127,2:126)])


chip.denovo_TES.quantstan <- cbind(chip.denovo_TES.log2_prelim$gene,as.data.frame(lapply(chip.denovo_TES.log2_prelim[,7:131], quantile90minmax)))
chip.denovo_TES.quantstan$gene<-chip.denovo_TES.quantstan[,1]
rownames(chip.denovo_TES.quantstan)<-chip.denovo_TES.quantstan$gene
chip.denovo_TES.quantstan<-processchipseq_table(chip.denovo_TES.quantstan[,c(127,2:126)])


chip.denovo_TSS.quantstan <- cbind(chip.denovo_TSS.log2_prelim$gene,as.data.frame(lapply(chip.denovo_TSS.log2_prelim[,7:131], quantile90minmax)))
chip.denovo_TSS.quantstan$gene<-chip.denovo_TSS.quantstan[,1]
rownames(chip.denovo_TSS.quantstan)<-chip.denovo_TSS.quantstan$gene
chip.denovo_TSS.quantstan<-processchipseq_table(chip.denovo_TSS.quantstan[,c(127,2:126)])


chip.denovo.log2<-processchipseq_table(chip.denovo.log2_prelim[,c(4,7:131)])

chip.denovo_TSS.log2<-processchipseq_table(chip.denovo_TSS.log2_prelim[,c(4,7:131)])

chip.denovo_TES.log2<-processchipseq_table(chip.denovo_TES.log2_prelim[,c(4,7:131)])



chip.araport.quantstan <- cbind(chip.araport.log2_prelim$gene,as.data.frame(lapply(chip.araport.log2_prelim[,7:131], quantile90minmax)))
chip.araport.quantstan$gene<-chip.araport.quantstan[,1]
chip.araport.quantstan<-processchipseq_table(chip.araport.quantstan[,c(127,2:126)])

chip.araport_TES.quantstan <- cbind(chip.araport_TES.log2_prelim$gene,as.data.frame(lapply(chip.araport_TES.log2_prelim[,7:131], quantile90minmax)))
chip.araport_TES.quantstan$gene<-chip.araport_TES.quantstan[,1]
chip.araport_TES.quantstan<-processchipseq_table(chip.araport_TES.quantstan[,c(127,2:126)])


chip.araport_TSS.quantstan <- cbind(chip.araport_TSS.log2_prelim$gene,as.data.frame(lapply(chip.araport_TSS.log2_prelim[,7:131], quantile90minmax)))
chip.araport_TSS.quantstan$gene<-chip.araport_TSS.quantstan[,1]
chip.araport_TSS.quantstan<-processchipseq_table(chip.araport_TSS.quantstan[,c(127,2:126)])

chip.araport.log2<-processchipseq_table(chip.araport.log2_prelim[,c(4,7:131)])

chip.araport_TSS.log2<-processchipseq_table(chip.araport_TSS.log2_prelim[,c(4,7:131)])

chip.araport_TES.log2<-processchipseq_table(chip.araport_TES.log2_prelim[,c(4,7:131)])














# standardize chipseq data 

# min-max transformation 

a<-chip.denovo.log2


# custom function to implement min max scaling
minMax <- function(x) {  (x - min(x)) / (max(x) - min(x))}
quantile90minmax <- function(x) {  (x - quantile(x,.10)) / (quantile(x,.80) - quantile(x,.20))}

quantile(x,.90)
quantile(chip.denovo.log2$K27.1741,probs = seq(0, 1, 0.25))
#normalise data using custom function

chip.denovo.minmax <- cbind(chip.denovo.log2$gene,as.data.frame(lapply(chip.denovo.log2[,2:68], minMax)))
chip.denovo.minmax$gene<-chip.denovo.minmax[,1]
rownames(chip.denovo.minmax)<-chip.denovo.minmax$gene

chip.denovo.quantstan <- cbind(chip.denovo.log2$gene,as.data.frame(lapply(chip.denovo.log2[,2:68], quantile90minmax)))
chip.denovo.quantstan$gene<-chip.denovo.quantstan[,1]
rownames(chip.denovo.quantstan)<-chip.denovo.quantstan$gene


boxplot(chip.denovo.quantstan[,grep("K27",names(chip.denovo.quantstan))],las=2)
boxplot(chip.denovo.quantstan[,grep("H1",names(chip.denovo.quantstan))],las=2)
boxplot(chip.denovo.quantstan[,grep("K9",names(chip.denovo.quantstan))],las=2)

boxplot(chip.denovo.minmax[,grep("K27",names(chip.denovo.minmax))],las=2)
boxplot(chip.denovo.minmax[,grep("H1",names(chip.denovo.minmax))],las=2)
boxplot(chip.denovo.minmax[,grep("K9",names(chip.denovo.minmax))],las=2)

plot(chip.denovo.minmax$H1.1741,chip.denovo.minmax$H1.5856)
cor(chip.denovo.minmax$H1.1741,chip.denovo.minmax$H1.5856)
cor(chip.denovo.log2$H1.1741,chip.denovo.log2$H1.5856)
cor(chip.denovo.minmax$K27.1741,chip.denovo.minmax$K27.5856)
cor(chip.denovo.log2$K27.1741,chip.denovo.log2$K27.5856)




########histone variation


chip.denovo.quantstan$sd.hist1<-apply(chip.denovo.quantstan[,grep("H1",names(chip.denovo.quantstan))],1,sd)
chip.denovo.quantstan$sd.key4<-apply(chip.denovo.quantstan[,grep("K4",names(chip.denovo.quantstan))],1,sd)
chip.denovo.quantstan$sd.key9<-apply(chip.denovo.quantstan[,grep("K9",names(chip.denovo.quantstan))],1,sd)
chip.denovo.quantstan$sd.key27<-apply(chip.denovo.quantstan[,grep("K27",names(chip.denovo.quantstan))],1,sd)
chip.denovo.quantstan$sd.key36<-apply(chip.denovo.quantstan[,grep("K36",names(chip.denovo.quantstan))],1,sd)


chip.denovo.quantstan$mean.hist1<-apply(chip.denovo.quantstan[,grep("H1",names(chip.denovo.quantstan))],1,mean)
chip.denovo.quantstan$mean.key4<-apply(chip.denovo.quantstan[,grep("K4",names(chip.denovo.quantstan))],1,mean)
chip.denovo.quantstan$mean.key9<-apply(chip.denovo.quantstan[,grep("K9",names(chip.denovo.quantstan))],1,mean)
chip.denovo.quantstan$mean.key27<-apply(chip.denovo.quantstan[,grep("K27",names(chip.denovo.quantstan))],1,mean)
chip.denovo.quantstan$mean.key36<-apply(chip.denovo.quantstan[,grep("K36",names(chip.denovo.quantstan))],1,mean)

chip.denovo.quantstan$min.hist1<-apply(chip.denovo.quantstan[,grep("H1",names(chip.denovo.quantstan))],1,min)
chip.denovo.quantstan$min.key4<-apply(chip.denovo.quantstan[,grep("K4",names(chip.denovo.quantstan))],1,min)
chip.denovo.quantstan$min.key9<-apply(chip.denovo.quantstan[,grep("K9",names(chip.denovo.quantstan))],1,min)
chip.denovo.quantstan$min.key27<-apply(chip.denovo.quantstan[,grep("K27",names(chip.denovo.quantstan))],1,min)
chip.denovo.quantstan$min.key36<-apply(chip.denovo.quantstan[,grep("K36",names(chip.denovo.quantstan))],1,min)

chip.denovo.quantstan$max.hist1<-apply(chip.denovo.quantstan[,grep("H1",names(chip.denovo.quantstan))],1,max)
chip.denovo.quantstan$max.key4<-apply(chip.denovo.quantstan[,grep("K4",names(chip.denovo.quantstan))],1,max)
chip.denovo.quantstan$max.key9<-apply(chip.denovo.quantstan[,grep("K9",names(chip.denovo.quantstan))],1,max)
chip.denovo.quantstan$max.key27<-apply(chip.denovo.quantstan[,grep("K27",names(chip.denovo.quantstan))],1,max)
chip.denovo.quantstan$max.key36<-apply(chip.denovo.quantstan[,grep("K36",names(chip.denovo.quantstan))],1,max)



chip.denovo.quantstan$variance.hist1<-chip.denovo.quantstan$sd.hist1/chip.denovo.quantstan$mean.hist1
chip.denovo.quantstan$variance.key4<-chip.denovo.quantstan$sd.key4/chip.denovo.quantstan$mean.key4
chip.denovo.quantstan$variance.key9<-chip.denovo.quantstan$sd.key9/chip.denovo.quantstan$mean.key9
chip.denovo.quantstan$variance.key27<-chip.denovo.quantstan$sd.key27/chip.denovo.quantstan$mean.key27
chip.denovo.quantstan$variance.key36<-chip.denovo.quantstan$sd.key36/chip.denovo.quantstan$mean.key36



rownames(chip.denovo.log2_all_samples)<-chip.denovo.log2_all_samples$gene


# export processed tables for GEO submission 

write.table(chip.denovo.log2,"03_Projects/ChIP-seq_data/2021/coverage/chip.5marks.ourannotation.log2_input_normalized_coverage.genebody.bed",quote = F,row.names = F, col.names = T,sep="\t")


write.table(chip.denovo_TSS.log2,"03_Projects/ChIP-seq_data/2021/coverage/chip.5marks.ourannotation.log2_input_normalized_coverage.promoter.bed",quote = F,row.names = F, col.names = T,sep="\t")


write.table(chip.denovo.quantstan,"03_Projects/ChIP-seq_data/2021/coverage/chip.5marks.ourannotation.log2_input_normalized_coverage.80_20_quantiles_normalized.genebody.bed",quote = F,row.names = F, col.names = T,sep="\t")
write.table(chip.denovo_TSS.quantstan,"03_Projects/ChIP-seq_data/2021/coverage/chip.5marks.ourannotation.log2_input_normalized_coverage.80_20_quantiles_normalized.promoter.bed",quote = F,row.names = F, col.names = T,sep="\t")


write.table(chip.araport.log2,"03_Projects/ChIP-seq_data/2021/coverage/chip.5marks.Araport11_PCs_TEs.log2_input_normalized_coverage.genebody.bed",quote = F,row.names = F, col.names = T,sep="\t")

write.table(chip.araport.quantstan,"03_Projects/ChIP-seq_data/2021/coverage/chip.5marks.Araport11_PCs_TEs.log2_input_normalized_coverage.80_20_quantiles_normalized.genebody.bed",quote = F,row.names = F, col.names = T,sep="\t")




