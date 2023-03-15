setwd("/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/")
setwd("Z:/01_POSTDOC/")





###########################
#upload TPM tables for de novo genes - 4 datasets: 1.1001G, 2. 1001Gnew, 3.ERACAPS, 4.Cortijo
############################

####
#denovo 20211013 1001G 461 accessions
####
denovo2021.TPMs.genes.1001G <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/denovo_Oct2021.TPMs.genes.1001G.bed")

denovo2021.TPMs.genes.1001G$mean<-apply(denovo2021.TPMs.genes.1001G[,2:462],1,mean)
denovo2021.TPMs.genes.1001G$max<-apply(denovo2021.TPMs.genes.1001G[,2:462],1,max)
denovo2021.TPMs.genes.1001G$variance<-apply(denovo2021.TPMs.genes.1001G[,2:462],1,sd)/denovo2021.TPMs.genes.1001G$mean
denovo2021.TPMs.genes.1001G$Nacc_where_expressed05<-apply(denovo2021.TPMs.genes.1001G[,2:462], 1, function(i) sum(i > 0.5))
denovo2021.TPMs.genes.1001G$Nacc_where_expressed2<-apply(denovo2021.TPMs.genes.1001G[,2:462], 1, function(i) sum(i > 2))
rownames(denovo2021.TPMs.genes.1001G)<-denovo2021.TPMs.genes.1001G$gene
denovo2021.TPMs.genes.1001G$gene_type<-"other"
denovo2021.TPMs.genes.1001G$gene_type[ denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene]<-"pc"
denovo2021.TPMs.genes.1001G$gene_type[ denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene]<-"as"
denovo2021.TPMs.genes.1001G$gene_type[ denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene]<-"linc"
denovo2021.TPMs.genes.1001G$gene_type[ denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.AS_to_TE.loci$gene]<-"as_to_te"
denovo2021.TPMs.genes.1001G$gene_type[ denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.AS_to_pseudo.loci$gene]<-"as_to_pseudo"
denovo2021.TPMs.genes.1001G$gene_type[ denovo2021.TPMs.genes.1001G$gene %in% TE_genes.loci$gene]<-"te_gene"
denovo2021.TPMs.genes.1001G$gene_type[ denovo2021.TPMs.genes.1001G$gene %in% TE_frags.transcripts$gene]<-"te_frags"
denovo2021.TPMs.genes.1001G$gene_type[ denovo2021.TPMs.genes.1001G$gene %in% denovo_pseudogene.transcripts$gene]<-"pseudogene"


pc_1001<-denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene %in% denovoPC.loci$gene,]
as_1001<-denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.antisense.loci$gene,]
linc_1001<-denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene %in% lncRNAs.intergenic.loci$gene,]
te_1001<-denovo2021.TPMs.genes.1001G[denovo2021.TPMs.genes.1001G$gene %in% TE_genes.loci$gene,]




#number of genes expressed 
a<-as.data.frame(t(denovo2021.TPMs.genes.1001G[,2:462]))

a$N_PC<-apply(a[,colnames(a) %in% denovoPC.loci$gene], 1, function(i) sum(i > 0.5))
a$N_AS<-apply(a[,names(a )%in% lncRNAs.antisense.loci$gene], 1, function(i) sum(i > 0.5))
a$N_TE<-apply(a[,names(a) %in% TE_genes.loci$gene], 1, function(i) sum(i > 0.5))
a$N_linc<-apply(a[,names(a) %in% lncRNAs.intergenic.loci$gene], 1, function(i) sum(i > 0.5))
ab<-a[,38155:38158]




##############################
# denovo - 1001G new dataset
#############################



denovo2021.TPMs.genes.1001Gnew <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/denovo_Oct2021.TPMs.genes.1001new.bed")
denovo2021.TPMs.genes.1001Gnew<-denovo2021.TPMs.genes.1001Gnew[,1:87]
denovo2021.TPMs.genes.1001Gnew$me_an<-apply(denovo2021.TPMs.genes.1001Gnew[,2:87],1,mean)
denovo2021.TPMs.genes.1001Gnew$ma_x<-apply(denovo2021.TPMs.genes.1001Gnew[,2:87],1,max)
denovo2021.TPMs.genes.1001Gnew$va_riance<-apply(denovo2021.TPMs.genes.1001Gnew[,2:87],1,sd)/denovo2021.TPMs.genes.1001Gnew$me_an

rownames(denovo2021.TPMs.genes.1001Gnew)<-denovo2021.TPMs.genes.1001Gnew$gene

a<-denovo2021.TPMs.genes.1001Gnew
a$mean.1741<-apply(a[,grep("1741",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.1741<-apply(a[,grep("1741",names(a[,1:87]))],1,sd)/apply(a[,grep("1741",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.1741<-apply(a[,grep("1741",names(a[,1:87]))],1,sd)

a$mean.4807<-apply(a[,grep("4807",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.4807<-apply(a[,grep("4807",names(a[,1:87]))],1,sd)/apply(a[,grep("4807",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.4807<-apply(a[,grep("4807",names(a[,1:87]))],1,sd)

a$mean.5210<-apply(a[,grep("5210",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.5210<-apply(a[,grep("5210",names(a[,1:87]))],1,sd)/apply(a[,grep("5210",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.5210<-apply(a[,grep("5210",names(a[,1:87]))],1,sd)

a$mean.5772<-apply(a[,grep("5772",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.5772<-apply(a[,grep("5772",names(a[,1:87]))],1,sd)/apply(a[,grep("5772",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.5772<-apply(a[,grep("5772",names(a[,1:87]))],1,sd)

a$mean.5784<-apply(a[,grep("5784",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.5784<-apply(a[,grep("5784",names(a[,1:87]))],1,sd)/apply(a[,grep("5784",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.5784<-apply(a[,grep("5784",names(a[,1:87]))],1,sd)

a$mean.5856<-apply(a[,grep("5856",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.5856<-apply(a[,grep("5856",names(a[,1:87]))],1,sd)/apply(a[,grep("5856",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.5856<-apply(a[,grep("5856",names(a[,1:87]))],1,sd)

a$mean.6021<-apply(a[,grep("6021",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6021<-apply(a[,grep("6021",names(a[,1:87]))],1,sd)/apply(a[,grep("6021",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6021<-apply(a[,grep("6021",names(a[,1:87]))],1,sd)

a$mean.6220<-apply(a[,grep("6220",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6220<-apply(a[,grep("6220",names(a[,1:87]))],1,sd)/apply(a[,grep("6220",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6220<-apply(a[,grep("6220",names(a[,1:87]))],1,sd)

a$mean.6909<-apply(a[,grep("6909",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6909<-apply(a[,grep("6909",names(a[,1:87]))],1,sd)/apply(a[,grep("6909",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6909<-apply(a[,grep("6909",names(a[,1:87]))],1,sd)

a$mean.6911<-apply(a[,grep("6911",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6911<-apply(a[,grep("6911",names(a[,1:87]))],1,sd)/apply(a[,grep("6911",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6911<-apply(a[,grep("6911",names(a[,1:87]))],1,sd)

a$mean.6966<-apply(a[,grep("6966",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6966<-apply(a[,grep("6966",names(a[,1:87]))],1,sd)/apply(a[,grep("6966",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6966<-apply(a[,grep("6966",names(a[,1:87]))],1,sd)

a$mean.8244<-apply(a[,grep("8244",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.8244<-apply(a[,grep("8244",names(a[,1:87]))],1,sd)/apply(a[,grep("8244",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.8244<-apply(a[,grep("8244",names(a[,1:87]))],1,sd)

a$mean.8366<-apply(a[,grep("8366",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.8366<-apply(a[,grep("8366",names(a[,1:87]))],1,sd)/apply(a[,grep("8366",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.8366<-apply(a[,grep("8366",names(a[,1:87]))],1,sd)

a$mean.9518<-apply(a[,grep("9518",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.9518<-apply(a[,grep("9518",names(a[,1:87]))],1,sd)/apply(a[,grep("9518",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.9518<-apply(a[,grep("9518",names(a[,1:87]))],1,sd)

a$mean.9588<-apply(a[,grep("9588",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.9588<-apply(a[,grep("9588",names(a[,1:87]))],1,sd)/apply(a[,grep("9588",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.9588<-apply(a[,grep("9588",names(a[,1:87]))],1,sd)

a$mean.9888<-apply(a[,grep("9888",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.9888<-apply(a[,grep("9888",names(a[,1:87]))],1,sd)/apply(a[,grep("9888",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.9888<-apply(a[,grep("9888",names(a[,1:87]))],1,sd)

a$mean.9905<-apply(a[,grep("9905",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.9905<-apply(a[,grep("9905",names(a[,1:87]))],1,sd)/apply(a[,grep("9905",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.9905<-apply(a[,grep("9905",names(a[,1:87]))],1,sd)

a$mean.10012<-apply(a[,grep("10012",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.10012<-apply(a[,grep("10012",names(a[,1:87]))],1,sd)/apply(a[,grep("10012",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.10012<-apply(a[,grep("10012",names(a[,1:87]))],1,sd)

a$mean.1254<-apply(a[,grep("1254",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.1254<-apply(a[,grep("1254",names(a[,1:87]))],1,sd)/apply(a[,grep("1254",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.1254<-apply(a[,grep("1254",names(a[,1:87]))],1,sd)

a$mean.6024<-apply(a[,grep("6024",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6024<-apply(a[,grep("6024",names(a[,1:87]))],1,sd)/apply(a[,grep("6024",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6024<-apply(a[,grep("6024",names(a[,1:87]))],1,sd)

a$mean.6069<-apply(a[,grep("6069",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6069<-apply(a[,grep("6069",names(a[,1:87]))],1,sd)/apply(a[,grep("6069",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6069<-apply(a[,grep("6069",names(a[,1:87]))],1,sd)

a$mean.6076<-apply(a[,grep("6076",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6076<-apply(a[,grep("6076",names(a[,1:87]))],1,sd)/apply(a[,grep("6076",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6076<-apply(a[,grep("6076",names(a[,1:87]))],1,sd)

a$mean.6184<-apply(a[,grep("6184",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6184<-apply(a[,grep("6184",names(a[,1:87]))],1,sd)/apply(a[,grep("6184",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6184<-apply(a[,grep("6184",names(a[,1:87]))],1,sd)

a$mean.6189<-apply(a[,grep("6189",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6189<-apply(a[,grep("6189",names(a[,1:87]))],1,sd)/apply(a[,grep("6189",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6189<-apply(a[,grep("6189",names(a[,1:87]))],1,sd)

a$mean.6244<-apply(a[,grep("6244",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6244<-apply(a[,grep("6244",names(a[,1:87]))],1,sd)/apply(a[,grep("6244",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6244<-apply(a[,grep("6244",names(a[,1:87]))],1,sd)

a$mean.9057<-apply(a[,grep("9057",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.9057<-apply(a[,grep("9057",names(a[,1:87]))],1,sd)/apply(a[,grep("9057",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.9057<-apply(a[,grep("9057",names(a[,1:87]))],1,sd)

a$mean.9412<-apply(a[,grep("9412",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.9412<-apply(a[,grep("9412",names(a[,1:87]))],1,sd)/apply(a[,grep("9412",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.9412<-apply(a[,grep("9412",names(a[,1:87]))],1,sd)

a$mean.9470<-apply(a[,grep("9412",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.9470<-apply(a[,grep("9470",names(a[,1:87]))],1,sd,na.rm = TRUE)/apply(a[,grep("9412",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.9470<-apply(a[,grep("9470",names(a[,1:87]))],1,sd,na.rm = TRUE)

a$mean_intravariance<-apply(a[,grep("var.",names(a))],1,mean,na.rm = TRUE)
a$mean_intra_sd<-apply(a[,grep("sd.",names(a))],1,mean,na.rm = TRUE)

a$mean_of_means<-apply(a[,grep("mean.",names(a))],1,mean,na.rm = TRUE)
a$variance_of_means<-apply(a[,grep("mean.",names(a))],1,sd)/a$mean_of_means
a$sd_of_means<-apply(a[,grep("mean.",names(a))],1,sd)

denovo2021.TPMs.genes.1001Gnew<-a

rm(a)

denovo2021.TPMs.genes.1001Gnew$gene_type<-"other"
denovo2021.TPMs.genes.1001Gnew$gene_type[ denovo2021.TPMs.genes.1001Gnew$gene %in% denovoPC.loci$gene]<-"pc"
denovo2021.TPMs.genes.1001Gnew$gene_type[ denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.antisense.loci$gene]<-"as"
denovo2021.TPMs.genes.1001Gnew$gene_type[ denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.intergenic.loci$gene]<-"linc"
denovo2021.TPMs.genes.1001Gnew$gene_type[ denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.AS_to_TE.loci$gene]<-"as_to_te"
denovo2021.TPMs.genes.1001Gnew$gene_type[ denovo2021.TPMs.genes.1001Gnew$gene %in% lncRNAs.AS_to_pseudo.loci$gene]<-"as_to_pseudo"
denovo2021.TPMs.genes.1001Gnew$gene_type[ denovo2021.TPMs.genes.1001Gnew$gene %in% TE_genes.loci$gene]<-"te_gene"
denovo2021.TPMs.genes.1001Gnew$gene_type[ denovo2021.TPMs.genes.1001Gnew$gene %in% TE_frags.transcripts$gene]<-"te_frags"
denovo2021.TPMs.genes.1001Gnew$gene_type[ denovo2021.TPMs.genes.1001Gnew$gene %in% denovo_pseudogene.transcripts$gene]<-"pseudogene"




#variability between random 3 samples

rand1<-denovo2021.TPMs.genes.1001Gnew[,sample(2:87,3)]
denovo2021.TPMs.genes.1001Gnew$var_3random_1<-apply(rand1,1,sd)/apply(rand1,1,mean)
rand2<-denovo2021.TPMs.genes.1001Gnew[,sample(2:87,3)]
denovo2021.TPMs.genes.1001Gnew$var_3random_2<-apply(rand2,1,sd)/apply(rand2,1,mean)
rand3<-denovo2021.TPMs.genes.1001Gnew[,sample(2:87,3)]
denovo2021.TPMs.genes.1001Gnew$var_3random_3<-apply(rand3,1,sd)/apply(rand3,1,mean)
rand4<-denovo2021.TPMs.genes.1001Gnew[,sample(2:87,3)]
denovo2021.TPMs.genes.1001Gnew$var_3random_4<-apply(rand4,1,sd)/apply(rand4,1,mean)
rand5<-denovo2021.TPMs.genes.1001Gnew[,sample(2:87,3)]
denovo2021.TPMs.genes.1001Gnew$var_3random_5<-apply(rand5,1,sd)/apply(rand5,1,mean)
rand6<-denovo2021.TPMs.genes.1001Gnew[,sample(2:87,3)]
denovo2021.TPMs.genes.1001Gnew$var_3random_6<-apply(rand6,1,sd)/apply(rand6,1,mean)
rand7<-denovo2021.TPMs.genes.1001Gnew[,sample(2:87,3)]
denovo2021.TPMs.genes.1001Gnew$var_3random_7<-apply(rand7,1,sd)/apply(rand7,1,mean)
rand8<-denovo2021.TPMs.genes.1001Gnew[,sample(2:87,3)]
denovo2021.TPMs.genes.1001Gnew$var_3random_8<-apply(rand8,1,sd)/apply(rand8,1,mean)
denovo2021.TPMs.genes.1001Gnew$var_3random_mean<-apply(denovo2021.TPMs.genes.1001Gnew[,c("var_3random_1","var_3random_2","var_3random_3","var_3random_4","var_3random_5","var_3random_6","var_3random_7","var_3random_8")],1,mean)




####################

#import ERACAPS 
####################
denovo2021.TPMs.genes.ERACAPS <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/denovo_Oct2021.TPMs.genes.ERACAPS.bed")
#exclude P.9638.batch1
denovo2021.TPMs.genes.ERACAPS<-denovo2021.TPMs.genes.ERACAPS[,!(names(denovo2021.TPMs.genes.ERACAPS)=="P.9638.batch1")]
rownames(denovo2021.TPMs.genes.ERACAPS)<-denovo2021.TPMs.genes.ERACAPS$gene

denovo2021.TPMs.genes.ERACAPS$R.10015<-apply(denovo2021.TPMs.genes.ERACAPS[,c("R.10015","R.Sha")],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$S.10015<-apply(denovo2021.TPMs.genes.ERACAPS[,c("S.10015","S.Sha")],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS<-denovo2021.TPMs.genes.ERACAPS[,!(names(denovo2021.TPMs.genes.ERACAPS) %in% c("S.Sha","R.Sha","F.9543.batch1","P.9543.batch1","S.9543","R.9543","F.9638.batch1","P.9638.batch1","R.9638","S.9638"))]

denovo2021.TPMs.genes.ERACAPS$S.22001<-denovo2021.TPMs.genes.ERACAPS$S.85.3
denovo2021.TPMs.genes.ERACAPS$S.22002<-denovo2021.TPMs.genes.ERACAPS$S.35.1
denovo2021.TPMs.genes.ERACAPS$S.22003<-denovo2021.TPMs.genes.ERACAPS$S.Taz.0
denovo2021.TPMs.genes.ERACAPS$S.22004<-denovo2021.TPMs.genes.ERACAPS$S.Elh.2
denovo2021.TPMs.genes.ERACAPS$S.22005<-denovo2021.TPMs.genes.ERACAPS$S.R1
denovo2021.TPMs.genes.ERACAPS$S.22006<-denovo2021.TPMs.genes.ERACAPS$S.A1
denovo2021.TPMs.genes.ERACAPS$S.22007<-denovo2021.TPMs.genes.ERACAPS$S.ET.86.4

denovo2021.TPMs.genes.ERACAPS$R.22001<-denovo2021.TPMs.genes.ERACAPS$R.85.3
denovo2021.TPMs.genes.ERACAPS$R.22002<-denovo2021.TPMs.genes.ERACAPS$R.35.1
denovo2021.TPMs.genes.ERACAPS$R.22003<-denovo2021.TPMs.genes.ERACAPS$R.Taz.0
denovo2021.TPMs.genes.ERACAPS$R.22004<-denovo2021.TPMs.genes.ERACAPS$R.Elh.2
denovo2021.TPMs.genes.ERACAPS$R.22005<-denovo2021.TPMs.genes.ERACAPS$R.R1
denovo2021.TPMs.genes.ERACAPS$R.22006<-denovo2021.TPMs.genes.ERACAPS$R.A1
denovo2021.TPMs.genes.ERACAPS$R.22007<-denovo2021.TPMs.genes.ERACAPS$R.ET.86.4

denovo2021.TPMs.genes.ERACAPS$F.10002<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.10002.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$F.10015<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.10015.|F.Sha",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$F.10024<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.10024.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$F.1741<-denovo2021.TPMs.genes.ERACAPS$F.1741.batch1
denovo2021.TPMs.genes.ERACAPS$F.22002<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.35.1.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$F.6024<-denovo2021.TPMs.genes.ERACAPS$F.6024.batch1
denovo2021.TPMs.genes.ERACAPS$F.6244<-denovo2021.TPMs.genes.ERACAPS$F.6244.batch1
denovo2021.TPMs.genes.ERACAPS$F.6909<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.6909.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$F.6966<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.6966.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$F.8236<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.8236.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$F.22001<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.85.3.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$F.9075<-denovo2021.TPMs.genes.ERACAPS$F.9075.batch1

#denovo2021.TPMs.genes.ERACAPS$F.9543<-denovo2021.TPMs.genes.ERACAPS$F.9543.batch1 #! seed mixup 
#denovo2021.TPMs.genes.ERACAPS$F.9638<-denovo2021.TPMs.genes.ERACAPS$F.9638.batch1 #! seed mixup 

denovo2021.TPMs.genes.ERACAPS$F.9537<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.9537.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$F.9764<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.9764.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$F.9728<-denovo2021.TPMs.genes.ERACAPS$F.9728.batch1 
denovo2021.TPMs.genes.ERACAPS$F.9888<-denovo2021.TPMs.genes.ERACAPS$F.9888.batch1 
denovo2021.TPMs.genes.ERACAPS$F.9905<-denovo2021.TPMs.genes.ERACAPS$F.9905.batch1 

denovo2021.TPMs.genes.ERACAPS$F.9981<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.9981.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$F.22006<-denovo2021.TPMs.genes.ERACAPS$F.A1.batch1 
denovo2021.TPMs.genes.ERACAPS$F.22004<-denovo2021.TPMs.genes.ERACAPS$F.Elh.2.batch1 
denovo2021.TPMs.genes.ERACAPS$F.22003<-denovo2021.TPMs.genes.ERACAPS$F.Taz.0.batch1 

denovo2021.TPMs.genes.ERACAPS$F.22007<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.ET.86.4.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$F.22005<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.R1.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)


denovo2021.TPMs.genes.ERACAPS$P.10002<-denovo2021.TPMs.genes.ERACAPS$P.10002.batch2 
denovo2021.TPMs.genes.ERACAPS$P.10015<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("P.10015.batch2|P.Sha.batch2",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$P.10024<-denovo2021.TPMs.genes.ERACAPS$P.10024.batch2 
denovo2021.TPMs.genes.ERACAPS$P.1741<-denovo2021.TPMs.genes.ERACAPS$P.1741.batch1 
denovo2021.TPMs.genes.ERACAPS$P.22002<-denovo2021.TPMs.genes.ERACAPS$P.35.1.batch2 
denovo2021.TPMs.genes.ERACAPS$P.6024<-denovo2021.TPMs.genes.ERACAPS$P.6024.batch1 
denovo2021.TPMs.genes.ERACAPS$P.6244<-denovo2021.TPMs.genes.ERACAPS$P.6244.batch1 
denovo2021.TPMs.genes.ERACAPS$P.6909<-denovo2021.TPMs.genes.ERACAPS$P.6909.batch2 

denovo2021.TPMs.genes.ERACAPS$P.6966<-denovo2021.TPMs.genes.ERACAPS$P.6966.batch2 

denovo2021.TPMs.genes.ERACAPS$P.8236<-denovo2021.TPMs.genes.ERACAPS$P.8236.batch2 

denovo2021.TPMs.genes.ERACAPS$P.22001<-denovo2021.TPMs.genes.ERACAPS$P.85.3.batch2 

denovo2021.TPMs.genes.ERACAPS$P.9075<-denovo2021.TPMs.genes.ERACAPS$P.9075.batch1 

denovo2021.TPMs.genes.ERACAPS$P.9537<-denovo2021.TPMs.genes.ERACAPS$P.9537.batch2 
# denovo2021.TPMs.genes.ERACAPS$P.9543<-denovo2021.TPMs.genes.ERACAPS$P.9543.batch1 #! seed mixup 

denovo2021.TPMs.genes.ERACAPS$P.9728<-denovo2021.TPMs.genes.ERACAPS$P.9728.batch1 

denovo2021.TPMs.genes.ERACAPS$P.9888<-denovo2021.TPMs.genes.ERACAPS$P.9888.batch1 
denovo2021.TPMs.genes.ERACAPS$P.9905<-denovo2021.TPMs.genes.ERACAPS$P.9905.batch1 
denovo2021.TPMs.genes.ERACAPS$P.9981<-denovo2021.TPMs.genes.ERACAPS$P.9981.batch2 
denovo2021.TPMs.genes.ERACAPS$P.22006<-denovo2021.TPMs.genes.ERACAPS$P.A1.batch1 
denovo2021.TPMs.genes.ERACAPS$P.22004<-denovo2021.TPMs.genes.ERACAPS$P.Elh.2.batch1 
denovo2021.TPMs.genes.ERACAPS$P.9888<-denovo2021.TPMs.genes.ERACAPS$P.9888.batch1 
denovo2021.TPMs.genes.ERACAPS$P.22003<-denovo2021.TPMs.genes.ERACAPS$P.Taz.0.batch1 
denovo2021.TPMs.genes.ERACAPS$P.9764<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("P.9764.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$P.22007<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("P.ET.86.4.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$P.22005<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("P.R1.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)


denovo2021.TPMs.genes.ERACAPS<-denovo2021.TPMs.genes.ERACAPS[,c(1,67:176)]
denovo2021.TPMs.genes.ERACAPS<-denovo2021.TPMs.genes.ERACAPS[,grep("^(?!.*85.3)^(?!.*35.1)^(?!.*Taz)^(?!.*Elh)^(?!.*R1)^(?!.*A1)^(?!.*ET.86)",names(denovo2021.TPMs.genes.ERACAPS),perl = T)]


denovo2021.TPMs.genes.ERACAPS$mean<-apply(denovo2021.TPMs.genes.ERACAPS[,2:97],1,mean)
denovo2021.TPMs.genes.ERACAPS$max<-apply(denovo2021.TPMs.genes.ERACAPS[,2:97],1,max)
denovo2021.TPMs.genes.ERACAPS$variance<-apply(denovo2021.TPMs.genes.ERACAPS[,2:97],1,sd)/denovo2021.TPMs.genes.ERACAPS$mean
denovo2021.TPMs.genes.ERACAPS$sd<-apply(denovo2021.TPMs.genes.ERACAPS[,2:97],1,sd)



denovo2021.TPMs.genes.ERACAPS$mean.flowers<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$max.flowers<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.",names(denovo2021.TPMs.genes.ERACAPS))],1,max,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$var.flowers<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.",names(denovo2021.TPMs.genes.ERACAPS))],1,sd)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$sd.flowers<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.",names(denovo2021.TPMs.genes.ERACAPS))],1,sd)

denovo2021.TPMs.genes.ERACAPS$mean.pollen<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("P.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$max.pollen<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("P.",names(denovo2021.TPMs.genes.ERACAPS))],1,max,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$var.pollen<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("P.",names(denovo2021.TPMs.genes.ERACAPS))],1,sd)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("P.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$sd.pollen<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("P.",names(denovo2021.TPMs.genes.ERACAPS))],1,sd)

denovo2021.TPMs.genes.ERACAPS$mean.seedl<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("^(?!.*Sha).*S",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)

denovo2021.TPMs.genes.ERACAPS$max.seedl<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("^(?!.*Sha).*S",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,max,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$var.seedl<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("^(?!.*Sha).*S",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("S.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$sd.seedl<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("^(?!.*Sha).*S",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd)

denovo2021.TPMs.genes.ERACAPS$mean.rosette<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("^(?!.*S.R1)^(?!.*F.R1)^(?!.*P.R1).*R",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$max.rosette<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("^(?!.*S.R1)^(?!.*F.R1)^(?!.*P.R1).*R",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,max,na.rm = TRUE)

denovo2021.TPMs.genes.ERACAPS$var.rosette<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("^(?!.*S.R1)^(?!.*F.R1)^(?!.*P.R1).*R",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("R.",names(denovo2021.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$sd.rosette<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("^(?!.*S.R1)^(?!.*F.R1)^(?!.*P.R1).*R",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd)


denovo2021.TPMs.genes.ERACAPS$rosette.expr.frequency<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("^(?!.*S.R1)^(?!.*F.R1)^(?!.*P.R1).*R",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,function(i) sum(i > 0.5))/length(as.vector(t(denovo2021.TPMs.genes.ERACAPS[1,grep("^(?!.*S.R1)^(?!.*F.R1)^(?!.*P.R1).*R",names(denovo2021.TPMs.genes.ERACAPS),perl = T)])))  

denovo2021.TPMs.genes.ERACAPS$seedl.expr.frequency<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("^(?!.*Sha).*S",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,function(i) sum(i > 0.5))/length(as.vector(t(denovo2021.TPMs.genes.ERACAPS[1,grep("^(?!.*Sha).*S",names(denovo2021.TPMs.genes.ERACAPS),perl = T)])))  

denovo2021.TPMs.genes.ERACAPS$flower.expr.frequency<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("F.",names(denovo2021.TPMs.genes.ERACAPS))],1,function(i) sum(i > 0.5))/length(as.vector(t(denovo2021.TPMs.genes.ERACAPS[1,grep("F.",names(denovo2021.TPMs.genes.ERACAPS))])))  

denovo2021.TPMs.genes.ERACAPS$pollen.expr.frequency<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("P.",names(denovo2021.TPMs.genes.ERACAPS))],1,function(i) sum(i > 0.5))/length(as.vector(t(denovo2021.TPMs.genes.ERACAPS[1,grep("P.",names(denovo2021.TPMs.genes.ERACAPS))])))  



denovo2021.TPMs.genes.ERACAPS$gene_type<-"other"
denovo2021.TPMs.genes.ERACAPS$gene_type[ denovo2021.TPMs.genes.ERACAPS$gene %in% denovoPC.loci$gene]<-"pc"
denovo2021.TPMs.genes.ERACAPS$gene_type[ denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.antisense.loci$gene]<-"as"
denovo2021.TPMs.genes.ERACAPS$gene_type[ denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.intergenic.loci$gene]<-"linc"
denovo2021.TPMs.genes.ERACAPS$gene_type[ denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.AS_to_TE.loci$gene]<-"as_to_te"
denovo2021.TPMs.genes.ERACAPS$gene_type[ denovo2021.TPMs.genes.ERACAPS$gene %in% lncRNAs.AS_to_pseudo.loci$gene]<-"as_to_pseudo"
denovo2021.TPMs.genes.ERACAPS$gene_type[ denovo2021.TPMs.genes.ERACAPS$gene %in% TE_genes.loci$gene]<-"te_gene"
denovo2021.TPMs.genes.ERACAPS$gene_type[ denovo2021.TPMs.genes.ERACAPS$gene %in% TE_frags.transcripts$gene]<-"te_frags"
denovo2021.TPMs.genes.ERACAPS$gene_type[ denovo2021.TPMs.genes.ERACAPS$gene %in% denovo_pseudogene.transcripts$gene]<-"pseudogene"




#add variation between tissues :
denovo2021.TPMs.genes.ERACAPS$tissue_variance.10002<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("10002",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("10002",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.10015<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("10015",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("10015",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.10024<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("10024",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("10024",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.1741<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("1741",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("1741",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.6024<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("6024",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("6024",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.6244<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("6244",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("6244",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.6909<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("6909",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("6909",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.6966<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("6966",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("6966",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.8236<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("8236",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("8236",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.9075<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("9075",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("9075",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.9537<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("9537",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("9537",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.9728<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("9728",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("9728",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.9764<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("9764",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("9764",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.9888<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("9888",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("9888",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.9905<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("9905",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("9905",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.9981<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("9981",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("9981",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.22001<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("22001",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("22001",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.22002<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("22002",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("22002",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.22003<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("22003",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("22003",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.22004<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("22004",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("22004",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.22005<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("22005",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("22005",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.22006<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("22006",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("22006",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$tissue_variance.22007<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("22007",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.ERACAPS[,grep("22007",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.ERACAPS$avertissuevar<-apply(denovo2021.TPMs.genes.ERACAPS[,grep("tissue_variance",names(denovo2021.TPMs.genes.ERACAPS),perl = T)],1,mean,na.rm = TRUE)















#denovo2021.TPMs.genes.ERACAPS$Nacc_05_F<-apply(denovo2021.TPMs.genes.ERACAPS[,2:462], 1, function(i) sum(i > 0.5))
# denovo2021.TPMs.genes.ERACAPS$Nacc_where_expressed2<-apply(denovo2021.TPMs.genes.1001G[,2:462], 1, function(i) sum(i > 2))


## import Cortijo 
denovo2021.TPMs.genes.Cortijo <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/denovo_Oct2021.TPMs.genes.Cortijo.bed")


denovo2021.TPMs.genes.Cortijo$mean<-apply(denovo2021.TPMs.genes.Cortijo[,2:169],1,mean)
denovo2021.TPMs.genes.Cortijo$max<-apply(denovo2021.TPMs.genes.Cortijo[,2:169],1,max)
denovo2021.TPMs.genes.Cortijo$variance<-apply(denovo2021.TPMs.genes.Cortijo[,2:169],1,sd)/denovo2021.TPMs.genes.Cortijo$mean
denovo2021.TPMs.genes.Cortijo$sd<-apply(denovo2021.TPMs.genes.Cortijo[,2:169],1,sd)

denovo2021.TPMs.genes.Cortijo$mean.ZT2<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT2",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$noise.ZT2<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT2",names(denovo2021.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT2",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$max.ZT2<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT2",names(denovo2021.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

denovo2021.TPMs.genes.Cortijo$mean.ZT4<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT4",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$noise.ZT4<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT4",names(denovo2021.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT4",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$max.ZT4<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT4",names(denovo2021.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

denovo2021.TPMs.genes.Cortijo$mean.ZT6<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT6",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$noise.ZT6<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT6",names(denovo2021.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT6",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$max.ZT6<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT6",names(denovo2021.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

denovo2021.TPMs.genes.Cortijo$mean.ZT8<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT8",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$noise.ZT8<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT8",names(denovo2021.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT8",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$max.ZT8<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT8",names(denovo2021.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

denovo2021.TPMs.genes.Cortijo$mean.ZT10<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT10",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$noise.ZT10<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT10",names(denovo2021.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT10",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$max.ZT10<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT10",names(denovo2021.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

denovo2021.TPMs.genes.Cortijo$mean.ZT12<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT12",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$noise.ZT12<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT12",names(denovo2021.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT12",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$max.ZT12<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT12",names(denovo2021.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

denovo2021.TPMs.genes.Cortijo$mean.ZT14<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT14",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$noise.ZT14<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT14",names(denovo2021.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT14",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$max.ZT14<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT14",names(denovo2021.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

denovo2021.TPMs.genes.Cortijo$mean.ZT16<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT16",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$noise.ZT16<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT16",names(denovo2021.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT16",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$max.ZT16<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT16",names(denovo2021.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

denovo2021.TPMs.genes.Cortijo$mean.ZT18<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT18",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$noise.ZT18<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT18",names(denovo2021.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT18",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$max.ZT18<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT18",names(denovo2021.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

denovo2021.TPMs.genes.Cortijo$mean.ZT20<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT20",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$noise.ZT20<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT20",names(denovo2021.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT20",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$max.ZT20<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT20",names(denovo2021.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

denovo2021.TPMs.genes.Cortijo$mean.ZT22<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT22",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$noise.ZT22<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT22",names(denovo2021.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT22",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$max.ZT22<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT22",names(denovo2021.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

denovo2021.TPMs.genes.Cortijo$mean.ZT24<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT24",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$noise.ZT24<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT24",names(denovo2021.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT24",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$max.ZT24<-apply(denovo2021.TPMs.genes.Cortijo[,grep("ZT24",names(denovo2021.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

denovo2021.TPMs.genes.Cortijo$noise_average<-apply(denovo2021.TPMs.genes.Cortijo[,grep("noise",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
denovo2021.TPMs.genes.Cortijo$variance_circadian<-apply(denovo2021.TPMs.genes.Cortijo[,grep("mean.ZT",names(denovo2021.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(denovo2021.TPMs.genes.Cortijo[,grep("mean.ZT",names(denovo2021.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)



denovo2021.TPMs.genes.Cortijo$gene_type<-"other"
denovo2021.TPMs.genes.Cortijo$gene_type[ denovo2021.TPMs.genes.Cortijo$gene %in% denovoPC.loci$gene]<-"pc"
denovo2021.TPMs.genes.Cortijo$gene_type[ denovo2021.TPMs.genes.Cortijo$gene %in% lncRNAs.antisense.loci$gene]<-"as"
denovo2021.TPMs.genes.Cortijo$gene_type[ denovo2021.TPMs.genes.Cortijo$gene %in% lncRNAs.intergenic.loci$gene]<-"linc"
denovo2021.TPMs.genes.Cortijo$gene_type[ denovo2021.TPMs.genes.Cortijo$gene %in% lncRNAs.AS_to_TE.loci$gene]<-"as_to_te"
denovo2021.TPMs.genes.Cortijo$gene_type[ denovo2021.TPMs.genes.Cortijo$gene %in% lncRNAs.AS_to_pseudo.loci$gene]<-"as_to_pseudo"
denovo2021.TPMs.genes.Cortijo$gene_type[ denovo2021.TPMs.genes.Cortijo$gene %in% TE_genes.loci$gene]<-"te_gene"
denovo2021.TPMs.genes.Cortijo$gene_type[ denovo2021.TPMs.genes.Cortijo$gene %in% TE_frags.transcripts$gene]<-"te_frags"
denovo2021.TPMs.genes.Cortijo$gene_type[ denovo2021.TPMs.genes.Cortijo$gene %in% denovo_pseudogene.transcripts$gene]<-"pseudogene"







pc_cortijo<-denovo2021.TPMs.genes.Cortijo[denovo2021.TPMs.genes.Cortijo$gene %in% denovoPC.loci$gene,]
as_cortijo<-denovo2021.TPMs.genes.Cortijo[denovo2021.TPMs.genes.Cortijo$gene %in% lncRNAs.antisense.loci$gene,]
linc_cortijo<-denovo2021.TPMs.genes.Cortijo[denovo2021.TPMs.genes.Cortijo$gene %in% lncRNAs.intergenic.loci$gene,]
as_te_cortijo<-denovo2021.TPMs.genes.Cortijo[denovo2021.TPMs.genes.Cortijo$gene %in% lncRNAs.AS_to_TE.loci$gene,]
te_gene_cortijo<-denovo2021.TPMs.genes.Cortijo[denovo2021.TPMs.genes.Cortijo$gene %in% TE_genes.transcripts$gene,]
te_frag_cortijo<-denovo2021.TPMs.genes.Cortijo[denovo2021.TPMs.genes.Cortijo$gene %in% TE_frags.transcripts$gene,]

boxplot(pc$noise_average[pc$max>1],as$noise_average[as$max>1],linc$noise_average[linc$max>1],
        as_te$noise_average[as_te$max>1],te_gene$noise_average[te_gene$max>1],te_frag$noise_average[te_frag$max>1])

boxplot(pc$variance_circadian[pc$max>1&pc$max<3],as$variance_circadian[as$max>1&as$max<3],linc$variance_circadian[linc$max>1&linc$max<3],
        as_te$variance_circadian[as_te$max>1&as_te$max<3],te_gene$variance_circadian[te_gene$max>1&te_gene$max<3],te_frag$variance_circadian[te_frag$max>1&te_frag$max<3])
length( as_te$variance_circadian[as_te$max>1&as_te$max<3])
length( te_frag$variance_circadian[te_frag$max>1&te_frag$max<3])
length(te_gene$variance_circadian[te_gene$max>1&te_gene$max<3])

boxplot(pc$variance_circadian[pc$max>1&pc$max<3],as$variance_circadian[as$max>1&as$max<3],linc$variance_circadian[linc$max>1&linc$max<3],
        te_gene$variance_circadian[te_gene$max>1&te_gene$max<3])
boxplot(pc$noise_average[pc$max>1&pc$max<3],as$noise_average[as$max>1&as$max<3],linc$noise_average[linc$max>1&linc$max<3],
        te_gene$noise_average[te_gene$max>1&te_gene$max<3])










#import Araport11


Araport11.TPMs.genes.1001G <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/Araport11.TPMs.genes.1001G.bed")

Araport11.TPMs.genes.1001G$mean<-apply(Araport11.TPMs.genes.1001G[,2:462],1,mean)
Araport11.TPMs.genes.1001G$max<-apply(Araport11.TPMs.genes.1001G[,2:462],1,max)
Araport11.TPMs.genes.1001G$variance<-apply(Araport11.TPMs.genes.1001G[,2:462],1,sd)/Araport11.TPMs.genes.1001G$mean
Araport11.TPMs.genes.1001G$Nacc_where_expressed05<-apply(Araport11.TPMs.genes.1001G[,2:462], 1, function(i) sum(i > 0.5))
Araport11.TPMs.genes.1001G$Nacc_where_expressed2<-apply(Araport11.TPMs.genes.1001G[,2:462], 1, function(i) sum(i > 2))
rownames(Araport11.TPMs.genes.1001G)<-Araport11.TPMs.genes.1001G$gene
#rm(Araport11.TPMs.genes.1001G)


Araport11_withTEs.TPMs.genes.1001G <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/Araport11_withTEs.TPMs.genes.1001G.bed")

Araport11_withTEs.TPMs.genes.1001G$mean<-apply(Araport11_withTEs.TPMs.genes.1001G[,2:462],1,mean)
Araport11_withTEs.TPMs.genes.1001G$max<-apply(Araport11_withTEs.TPMs.genes.1001G[,2:462],1,max)
Araport11_withTEs.TPMs.genes.1001G$variance<-apply(Araport11_withTEs.TPMs.genes.1001G[,2:462],1,sd)/Araport11_withTEs.TPMs.genes.1001G$mean
Araport11_withTEs.TPMs.genes.1001G$Nacc_where_expressed05<-apply(Araport11_withTEs.TPMs.genes.1001G[,2:462], 1, function(i) sum(i > 0.5))
Araport11_withTEs.TPMs.genes.1001G$Nacc_where_expressed2<-apply(Araport11_withTEs.TPMs.genes.1001G[,2:462], 1, function(i) sum(i > 2))
rownames(Araport11_withTEs.TPMs.genes.1001G)<-Araport11_withTEs.TPMs.genes.1001G$gene

Araport11_withTEs.TPMs.genes.1001G$gene_type<-"other"
Araport11_withTEs.TPMs.genes.1001G$gene_type[Araport11_withTEs.TPMs.genes.1001G$gene %in% Araport11_protein_coding.201606.genes$gene]<-"PC_gene"
Araport11_withTEs.TPMs.genes.1001G$gene_type[Araport11_withTEs.TPMs.genes.1001G$gene %in% Araport11_transposable_element_gene$gene]<-"TE_gene"
Araport11_withTEs.TPMs.genes.1001G$gene_type[Araport11_withTEs.TPMs.genes.1001G$gene %in% Araport11_transposable_element_f$gene]<-"TE_fragment"
Araport11_withTEs.TPMs.genes.1001G$gene_type[Araport11_withTEs.TPMs.genes.1001G$gene %in% Araport11_non_coding$gene]<-"NC"
Araport11_withTEs.TPMs.genes.1001G$gene_type[Araport11_withTEs.TPMs.genes.1001G$gene %in% Araport11_novel_transcribed_region$gene]<-"NTR"








####################

#import ERACAPS 
####################
Araport11.TPMs.genes.ERACAPS <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/Araport11.TPMs.genes.ERACAPS.bed")
#exclude P.9638.batch1
Araport11.TPMs.genes.ERACAPS<-Araport11.TPMs.genes.ERACAPS[,!(names(Araport11.TPMs.genes.ERACAPS)=="P.9638.batch1")]
rownames(Araport11.TPMs.genes.ERACAPS)<-Araport11.TPMs.genes.ERACAPS$gene

Araport11.TPMs.genes.ERACAPS$R.10015<-apply(Araport11.TPMs.genes.ERACAPS[,c("R.10015","R.Sha")],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$S.10015<-apply(Araport11.TPMs.genes.ERACAPS[,c("S.10015","S.Sha")],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS<-Araport11.TPMs.genes.ERACAPS[,!(names(Araport11.TPMs.genes.ERACAPS) %in% c("S.Sha","R.Sha","F.9543.batch1","P.9543.batch1","S.9543","R.9543","F.9638.batch1","P.9638.batch1","R.9638","S.9638"))]


Araport11.TPMs.genes.ERACAPS$S.22001<-Araport11.TPMs.genes.ERACAPS$S.85.3
Araport11.TPMs.genes.ERACAPS$S.22002<-Araport11.TPMs.genes.ERACAPS$S.35.1
Araport11.TPMs.genes.ERACAPS$S.22003<-Araport11.TPMs.genes.ERACAPS$S.Taz.0
Araport11.TPMs.genes.ERACAPS$S.22004<-Araport11.TPMs.genes.ERACAPS$S.Elh.2
Araport11.TPMs.genes.ERACAPS$S.22005<-Araport11.TPMs.genes.ERACAPS$S.R1
Araport11.TPMs.genes.ERACAPS$S.22006<-Araport11.TPMs.genes.ERACAPS$S.A1
Araport11.TPMs.genes.ERACAPS$S.22007<-Araport11.TPMs.genes.ERACAPS$S.ET.86.4

Araport11.TPMs.genes.ERACAPS$R.22001<-Araport11.TPMs.genes.ERACAPS$R.85.3
Araport11.TPMs.genes.ERACAPS$R.22002<-Araport11.TPMs.genes.ERACAPS$R.35.1
Araport11.TPMs.genes.ERACAPS$R.22003<-Araport11.TPMs.genes.ERACAPS$R.Taz.0
Araport11.TPMs.genes.ERACAPS$R.22004<-Araport11.TPMs.genes.ERACAPS$R.Elh.2
Araport11.TPMs.genes.ERACAPS$R.22005<-Araport11.TPMs.genes.ERACAPS$R.R1
Araport11.TPMs.genes.ERACAPS$R.22006<-Araport11.TPMs.genes.ERACAPS$R.A1
Araport11.TPMs.genes.ERACAPS$R.22007<-Araport11.TPMs.genes.ERACAPS$R.ET.86.4

Araport11.TPMs.genes.ERACAPS$F.10002<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.10002.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$F.10015<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.10015.|F.Sha",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$F.10024<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.10024.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$F.1741<-Araport11.TPMs.genes.ERACAPS$F.1741.batch1
Araport11.TPMs.genes.ERACAPS$F.22002<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.35.1.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$F.6024<-Araport11.TPMs.genes.ERACAPS$F.6024.batch1
Araport11.TPMs.genes.ERACAPS$F.6244<-Araport11.TPMs.genes.ERACAPS$F.6244.batch1
Araport11.TPMs.genes.ERACAPS$F.6909<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.6909.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$F.6966<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.6966.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$F.8236<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.8236.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$F.22001<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.85.3.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$F.9075<-Araport11.TPMs.genes.ERACAPS$F.9075.batch1

#Araport11.TPMs.genes.ERACAPS$F.9543<-Araport11.TPMs.genes.ERACAPS$F.9543.batch1 #! seed mixup 
#Araport11.TPMs.genes.ERACAPS$F.9638<-Araport11.TPMs.genes.ERACAPS$F.9638.batch1 #! seed mixup 

Araport11.TPMs.genes.ERACAPS$F.9537<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.9537.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$F.9764<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.9764.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$F.9728<-Araport11.TPMs.genes.ERACAPS$F.9728.batch1 
Araport11.TPMs.genes.ERACAPS$F.9888<-Araport11.TPMs.genes.ERACAPS$F.9888.batch1 
Araport11.TPMs.genes.ERACAPS$F.9905<-Araport11.TPMs.genes.ERACAPS$F.9905.batch1 

Araport11.TPMs.genes.ERACAPS$F.9981<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.9981.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$F.22006<-Araport11.TPMs.genes.ERACAPS$F.A1.batch1 
Araport11.TPMs.genes.ERACAPS$F.22004<-Araport11.TPMs.genes.ERACAPS$F.Elh.2.batch1 
Araport11.TPMs.genes.ERACAPS$F.22003<-Araport11.TPMs.genes.ERACAPS$F.Taz.0.batch1 

Araport11.TPMs.genes.ERACAPS$F.22007<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.ET.86.4.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$F.22005<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.R1.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)


Araport11.TPMs.genes.ERACAPS$P.10002<-Araport11.TPMs.genes.ERACAPS$P.10002.batch2 
Araport11.TPMs.genes.ERACAPS$P.10015<-apply(Araport11.TPMs.genes.ERACAPS[,grep("P.10015.batch2|P.Sha.batch2",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$P.10024<-Araport11.TPMs.genes.ERACAPS$P.10024.batch2 
Araport11.TPMs.genes.ERACAPS$P.1741<-Araport11.TPMs.genes.ERACAPS$P.1741.batch1 
Araport11.TPMs.genes.ERACAPS$P.22002<-Araport11.TPMs.genes.ERACAPS$P.35.1.batch2 
Araport11.TPMs.genes.ERACAPS$P.6024<-Araport11.TPMs.genes.ERACAPS$P.6024.batch1 
Araport11.TPMs.genes.ERACAPS$P.6244<-Araport11.TPMs.genes.ERACAPS$P.6244.batch1 
Araport11.TPMs.genes.ERACAPS$P.6909<-Araport11.TPMs.genes.ERACAPS$P.6909.batch2 

Araport11.TPMs.genes.ERACAPS$P.6966<-Araport11.TPMs.genes.ERACAPS$P.6966.batch2 

Araport11.TPMs.genes.ERACAPS$P.8236<-Araport11.TPMs.genes.ERACAPS$P.8236.batch2 

Araport11.TPMs.genes.ERACAPS$P.22001<-Araport11.TPMs.genes.ERACAPS$P.85.3.batch2 

Araport11.TPMs.genes.ERACAPS$P.9075<-Araport11.TPMs.genes.ERACAPS$P.9075.batch1 

Araport11.TPMs.genes.ERACAPS$P.9537<-Araport11.TPMs.genes.ERACAPS$P.9537.batch2 
# Araport11.TPMs.genes.ERACAPS$P.9543<-Araport11.TPMs.genes.ERACAPS$P.9543.batch1 #! seed mixup 

Araport11.TPMs.genes.ERACAPS$P.9728<-Araport11.TPMs.genes.ERACAPS$P.9728.batch1 

Araport11.TPMs.genes.ERACAPS$P.9888<-Araport11.TPMs.genes.ERACAPS$P.9888.batch1 
Araport11.TPMs.genes.ERACAPS$P.9905<-Araport11.TPMs.genes.ERACAPS$P.9905.batch1 
Araport11.TPMs.genes.ERACAPS$P.9981<-Araport11.TPMs.genes.ERACAPS$P.9981.batch2 
Araport11.TPMs.genes.ERACAPS$P.22006<-Araport11.TPMs.genes.ERACAPS$P.A1.batch1 
Araport11.TPMs.genes.ERACAPS$P.22004<-Araport11.TPMs.genes.ERACAPS$P.Elh.2.batch1 
Araport11.TPMs.genes.ERACAPS$P.9888<-Araport11.TPMs.genes.ERACAPS$P.9888.batch1 
Araport11.TPMs.genes.ERACAPS$P.22003<-Araport11.TPMs.genes.ERACAPS$P.Taz.0.batch1 
Araport11.TPMs.genes.ERACAPS$P.9764<-apply(Araport11.TPMs.genes.ERACAPS[,grep("P.9764.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$P.22007<-apply(Araport11.TPMs.genes.ERACAPS[,grep("P.ET.86.4.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$P.22005<-apply(Araport11.TPMs.genes.ERACAPS[,grep("P.R1.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)



Araport11.TPMs.genes.ERACAPS<-Araport11.TPMs.genes.ERACAPS[,c(1,67:176)]
Araport11.TPMs.genes.ERACAPS<-Araport11.TPMs.genes.ERACAPS[,grep("^(?!.*85.3)^(?!.*35.1)^(?!.*Taz)^(?!.*Elh)^(?!.*R1)^(?!.*A1)^(?!.*ET.86)",names(Araport11.TPMs.genes.ERACAPS),perl = T)]




#Araport11.TPMs.genes.ERACAPS<-Araport11.TPMs.genes.ERACAPS[,c(1,67:162)]


Araport11.TPMs.genes.ERACAPS$mean<-apply(Araport11.TPMs.genes.ERACAPS[,2:97],1,mean)
Araport11.TPMs.genes.ERACAPS$max<-apply(Araport11.TPMs.genes.ERACAPS[,2:97],1,max)
Araport11.TPMs.genes.ERACAPS$variance<-apply(Araport11.TPMs.genes.ERACAPS[,2:97],1,sd)/Araport11.TPMs.genes.ERACAPS$mean
Araport11.TPMs.genes.ERACAPS$sd<-apply(Araport11.TPMs.genes.ERACAPS[,2:97],1,sd)



Araport11.TPMs.genes.ERACAPS$mean.flowers<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$max.flowers<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.",names(Araport11.TPMs.genes.ERACAPS))],1,max,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$var.flowers<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.",names(Araport11.TPMs.genes.ERACAPS))],1,sd)/apply(Araport11.TPMs.genes.ERACAPS[,grep("F.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$sd.flowers<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.",names(Araport11.TPMs.genes.ERACAPS))],1,sd)

Araport11.TPMs.genes.ERACAPS$mean.pollen<-apply(Araport11.TPMs.genes.ERACAPS[,grep("P.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$max.pollen<-apply(Araport11.TPMs.genes.ERACAPS[,grep("P.",names(Araport11.TPMs.genes.ERACAPS))],1,max,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$var.pollen<-apply(Araport11.TPMs.genes.ERACAPS[,grep("P.",names(Araport11.TPMs.genes.ERACAPS))],1,sd)/apply(Araport11.TPMs.genes.ERACAPS[,grep("P.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$sd.pollen<-apply(Araport11.TPMs.genes.ERACAPS[,grep("P.",names(Araport11.TPMs.genes.ERACAPS))],1,sd)

Araport11.TPMs.genes.ERACAPS$mean.seedl<-apply(Araport11.TPMs.genes.ERACAPS[,grep("S",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)

Araport11.TPMs.genes.ERACAPS$max.seedl<-apply(Araport11.TPMs.genes.ERACAPS[,grep("S.",names(Araport11.TPMs.genes.ERACAPS))],1,max,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$var.seedl<-apply(Araport11.TPMs.genes.ERACAPS[,grep("S.",names(Araport11.TPMs.genes.ERACAPS))],1,sd)/apply(Araport11.TPMs.genes.ERACAPS[,grep("S.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$sd.seedl<-apply(Araport11.TPMs.genes.ERACAPS[,grep("S.",names(Araport11.TPMs.genes.ERACAPS))],1,sd)

Araport11.TPMs.genes.ERACAPS$mean.rosette<-apply(Araport11.TPMs.genes.ERACAPS[,grep("R.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$max.rosette<-apply(Araport11.TPMs.genes.ERACAPS[,grep("R.",names(Araport11.TPMs.genes.ERACAPS))],1,max,na.rm = TRUE)

Araport11.TPMs.genes.ERACAPS$var.rosette<-apply(Araport11.TPMs.genes.ERACAPS[,grep("R.",names(Araport11.TPMs.genes.ERACAPS))],1,sd)/apply(Araport11.TPMs.genes.ERACAPS[,grep("R.",names(Araport11.TPMs.genes.ERACAPS))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.ERACAPS$sd.rosette<-apply(Araport11.TPMs.genes.ERACAPS[,grep("R.",names(Araport11.TPMs.genes.ERACAPS))],1,sd)


Araport11.TPMs.genes.ERACAPS$rosette.expr.frequency<-apply(Araport11.TPMs.genes.ERACAPS[,grep("R.",names(Araport11.TPMs.genes.ERACAPS))],1,function(i) sum(i > 0.5))/length(as.vector(t(Araport11.TPMs.genes.ERACAPS[1,grep("R.",names(Araport11.TPMs.genes.ERACAPS))])))  

Araport11.TPMs.genes.ERACAPS$seedl.expr.frequency<-apply(Araport11.TPMs.genes.ERACAPS[,grep("S.",names(Araport11.TPMs.genes.ERACAPS))],1,function(i) sum(i > 0.5))/length(as.vector(t(Araport11.TPMs.genes.ERACAPS[1,grep("S.",names(Araport11.TPMs.genes.ERACAPS))])))  

Araport11.TPMs.genes.ERACAPS$flower.expr.frequency<-apply(Araport11.TPMs.genes.ERACAPS[,grep("F.",names(Araport11.TPMs.genes.ERACAPS))],1,function(i) sum(i > 0.5))/length(as.vector(t(Araport11.TPMs.genes.ERACAPS[1,grep("F.",names(Araport11.TPMs.genes.ERACAPS))])))  

Araport11.TPMs.genes.ERACAPS$pollen.expr.frequency<-apply(Araport11.TPMs.genes.ERACAPS[,grep("P.",names(Araport11.TPMs.genes.ERACAPS))],1,function(i) sum(i > 0.5))/length(as.vector(t(Araport11.TPMs.genes.ERACAPS[1,grep("P.",names(Araport11.TPMs.genes.ERACAPS))])))  







## import 1001G new Araport 11 
Araport11.TPMs.genes.1001Gnew <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/Araport11.TPMs.genes.1001new.bed")
Araport11.TPMs.genes.1001Gnew<-Araport11.TPMs.genes.1001Gnew[,1:87]
Araport11.TPMs.genes.1001Gnew$me_an<-apply(Araport11.TPMs.genes.1001Gnew[,2:87],1,mean)
Araport11.TPMs.genes.1001Gnew$ma_x<-apply(Araport11.TPMs.genes.1001Gnew[,2:87],1,max)
Araport11.TPMs.genes.1001Gnew$va_riance<-apply(Araport11.TPMs.genes.1001Gnew[,2:87],1,sd)/Araport11.TPMs.genes.1001Gnew$me_an
a<-Araport11.TPMs.genes.1001Gnew
a$mean.1741<-apply(a[,grep("1741",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.1741<-apply(a[,grep("1741",names(a[,1:87]))],1,sd)/apply(a[,grep("1741",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.1741<-apply(a[,grep("1741",names(a[,1:87]))],1,sd)

a$mean.4807<-apply(a[,grep("4807",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.4807<-apply(a[,grep("4807",names(a[,1:87]))],1,sd)/apply(a[,grep("4807",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.4807<-apply(a[,grep("4807",names(a[,1:87]))],1,sd)

a$mean.5210<-apply(a[,grep("5210",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.5210<-apply(a[,grep("5210",names(a[,1:87]))],1,sd)/apply(a[,grep("5210",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.5210<-apply(a[,grep("5210",names(a[,1:87]))],1,sd)

a$mean.5772<-apply(a[,grep("5772",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.5772<-apply(a[,grep("5772",names(a[,1:87]))],1,sd)/apply(a[,grep("5772",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.5772<-apply(a[,grep("5772",names(a[,1:87]))],1,sd)

a$mean.5784<-apply(a[,grep("5784",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.5784<-apply(a[,grep("5784",names(a[,1:87]))],1,sd)/apply(a[,grep("5784",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.5784<-apply(a[,grep("5784",names(a[,1:87]))],1,sd)

a$mean.5856<-apply(a[,grep("5856",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.5856<-apply(a[,grep("5856",names(a[,1:87]))],1,sd)/apply(a[,grep("5856",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.5856<-apply(a[,grep("5856",names(a[,1:87]))],1,sd)

a$mean.6021<-apply(a[,grep("6021",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6021<-apply(a[,grep("6021",names(a[,1:87]))],1,sd)/apply(a[,grep("6021",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6021<-apply(a[,grep("6021",names(a[,1:87]))],1,sd)

a$mean.6220<-apply(a[,grep("6220",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6220<-apply(a[,grep("6220",names(a[,1:87]))],1,sd)/apply(a[,grep("6220",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6220<-apply(a[,grep("6220",names(a[,1:87]))],1,sd)

a$mean.6909<-apply(a[,grep("6909",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6909<-apply(a[,grep("6909",names(a[,1:87]))],1,sd)/apply(a[,grep("6909",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6909<-apply(a[,grep("6909",names(a[,1:87]))],1,sd)

a$mean.6911<-apply(a[,grep("6911",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6911<-apply(a[,grep("6911",names(a[,1:87]))],1,sd)/apply(a[,grep("6911",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6911<-apply(a[,grep("6911",names(a[,1:87]))],1,sd)

a$mean.6966<-apply(a[,grep("6966",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6966<-apply(a[,grep("6966",names(a[,1:87]))],1,sd)/apply(a[,grep("6966",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6966<-apply(a[,grep("6966",names(a[,1:87]))],1,sd)

a$mean.8244<-apply(a[,grep("8244",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.8244<-apply(a[,grep("8244",names(a[,1:87]))],1,sd)/apply(a[,grep("8244",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.8244<-apply(a[,grep("8244",names(a[,1:87]))],1,sd)

a$mean.8366<-apply(a[,grep("8366",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.8366<-apply(a[,grep("8366",names(a[,1:87]))],1,sd)/apply(a[,grep("8366",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.8366<-apply(a[,grep("8366",names(a[,1:87]))],1,sd)

a$mean.9518<-apply(a[,grep("9518",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.9518<-apply(a[,grep("9518",names(a[,1:87]))],1,sd)/apply(a[,grep("9518",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.9518<-apply(a[,grep("9518",names(a[,1:87]))],1,sd)

a$mean.9588<-apply(a[,grep("9588",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.9588<-apply(a[,grep("9588",names(a[,1:87]))],1,sd)/apply(a[,grep("9588",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.9588<-apply(a[,grep("9588",names(a[,1:87]))],1,sd)

a$mean.9888<-apply(a[,grep("9888",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.9888<-apply(a[,grep("9888",names(a[,1:87]))],1,sd)/apply(a[,grep("9888",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.9888<-apply(a[,grep("9888",names(a[,1:87]))],1,sd)

a$mean.9905<-apply(a[,grep("9905",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.9905<-apply(a[,grep("9905",names(a[,1:87]))],1,sd)/apply(a[,grep("9905",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.9905<-apply(a[,grep("9905",names(a[,1:87]))],1,sd)

a$mean.10012<-apply(a[,grep("10012",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.10012<-apply(a[,grep("10012",names(a[,1:87]))],1,sd)/apply(a[,grep("10012",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.10012<-apply(a[,grep("10012",names(a[,1:87]))],1,sd)

a$mean.1254<-apply(a[,grep("1254",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.1254<-apply(a[,grep("1254",names(a[,1:87]))],1,sd)/apply(a[,grep("1254",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.1254<-apply(a[,grep("1254",names(a[,1:87]))],1,sd)

a$mean.6024<-apply(a[,grep("6024",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6024<-apply(a[,grep("6024",names(a[,1:87]))],1,sd)/apply(a[,grep("6024",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6024<-apply(a[,grep("6024",names(a[,1:87]))],1,sd)

a$mean.6069<-apply(a[,grep("6069",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6069<-apply(a[,grep("6069",names(a[,1:87]))],1,sd)/apply(a[,grep("6069",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6069<-apply(a[,grep("6069",names(a[,1:87]))],1,sd)

a$mean.6076<-apply(a[,grep("6076",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6076<-apply(a[,grep("6076",names(a[,1:87]))],1,sd)/apply(a[,grep("6076",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6076<-apply(a[,grep("6076",names(a[,1:87]))],1,sd)

a$mean.6184<-apply(a[,grep("6184",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6184<-apply(a[,grep("6184",names(a[,1:87]))],1,sd)/apply(a[,grep("6184",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6184<-apply(a[,grep("6184",names(a[,1:87]))],1,sd)

a$mean.6189<-apply(a[,grep("6189",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6189<-apply(a[,grep("6189",names(a[,1:87]))],1,sd)/apply(a[,grep("6189",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6189<-apply(a[,grep("6189",names(a[,1:87]))],1,sd)

a$mean.6244<-apply(a[,grep("6244",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.6244<-apply(a[,grep("6244",names(a[,1:87]))],1,sd)/apply(a[,grep("6244",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.6244<-apply(a[,grep("6244",names(a[,1:87]))],1,sd)

a$mean.9057<-apply(a[,grep("9057",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.9057<-apply(a[,grep("9057",names(a[,1:87]))],1,sd)/apply(a[,grep("9057",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.9057<-apply(a[,grep("9057",names(a[,1:87]))],1,sd)

a$mean.9412<-apply(a[,grep("9412",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.9412<-apply(a[,grep("9412",names(a[,1:87]))],1,sd)/apply(a[,grep("9412",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.9412<-apply(a[,grep("9412",names(a[,1:87]))],1,sd)

a$mean.9470<-apply(a[,grep("9412",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$var.9470<-apply(a[,grep("9470",names(a[,1:87]))],1,sd,na.rm = TRUE)/apply(a[,grep("9412",names(a[,1:87]))],1,mean,na.rm = TRUE)
a$sd.9470<-apply(a[,grep("9470",names(a[,1:87]))],1,sd,na.rm = TRUE)

a$mean_intravariance<-apply(a[,grep("var.",names(a))],1,mean,na.rm = TRUE)
a$mean_intra_sd<-apply(a[,grep("sd.",names(a))],1,mean,na.rm = TRUE)

a$mean_of_means<-apply(a[,grep("mean.",names(a))],1,mean,na.rm = TRUE)
a$variance_of_means<-apply(a[,grep("mean.",names(a))],1,sd)/a$mean_of_means
a$sd_of_means<-apply(a[,grep("mean.",names(a))],1,sd)

Araport11.TPMs.genes.1001Gnew<-a

rm(a)



## import Cortijo 
Araport11.TPMs.genes.Cortijo <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/Araport11.TPMs.genes.Cortijo.bed")
rownames(Araport11.TPMs.genes.Cortijo)<-Araport11.TPMs.genes.Cortijo$gene
Araport11.TPMs.genes.Cortijo$mean<-apply(Araport11.TPMs.genes.Cortijo[,2:169],1,mean)
Araport11.TPMs.genes.Cortijo$max<-apply(Araport11.TPMs.genes.Cortijo[,2:169],1,max)
Araport11.TPMs.genes.Cortijo$variance<-apply(Araport11.TPMs.genes.Cortijo[,2:169],1,sd)/Araport11.TPMs.genes.Cortijo$mean
Araport11.TPMs.genes.Cortijo$sd<-apply(Araport11.TPMs.genes.Cortijo[,2:169],1,sd)

Araport11.TPMs.genes.Cortijo$mean.ZT2<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT2",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$noise.ZT2<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT2",names(Araport11.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(Araport11.TPMs.genes.Cortijo[,grep("ZT2",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$max.ZT2<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT2",names(Araport11.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

Araport11.TPMs.genes.Cortijo$mean.ZT4<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT4",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$noise.ZT4<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT4",names(Araport11.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(Araport11.TPMs.genes.Cortijo[,grep("ZT4",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$max.ZT4<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT4",names(Araport11.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

Araport11.TPMs.genes.Cortijo$mean.ZT6<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT6",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$noise.ZT6<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT6",names(Araport11.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(Araport11.TPMs.genes.Cortijo[,grep("ZT6",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$max.ZT6<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT6",names(Araport11.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

Araport11.TPMs.genes.Cortijo$mean.ZT8<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT8",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$noise.ZT8<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT8",names(Araport11.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(Araport11.TPMs.genes.Cortijo[,grep("ZT8",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$max.ZT8<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT8",names(Araport11.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

Araport11.TPMs.genes.Cortijo$mean.ZT10<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT10",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$noise.ZT10<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT10",names(Araport11.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(Araport11.TPMs.genes.Cortijo[,grep("ZT10",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$max.ZT10<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT10",names(Araport11.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

Araport11.TPMs.genes.Cortijo$mean.ZT12<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT12",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$noise.ZT12<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT12",names(Araport11.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(Araport11.TPMs.genes.Cortijo[,grep("ZT12",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$max.ZT12<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT12",names(Araport11.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

Araport11.TPMs.genes.Cortijo$mean.ZT14<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT14",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$noise.ZT14<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT14",names(Araport11.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(Araport11.TPMs.genes.Cortijo[,grep("ZT14",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$max.ZT14<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT14",names(Araport11.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

Araport11.TPMs.genes.Cortijo$mean.ZT16<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT16",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$noise.ZT16<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT16",names(Araport11.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(Araport11.TPMs.genes.Cortijo[,grep("ZT16",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$max.ZT16<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT16",names(Araport11.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

Araport11.TPMs.genes.Cortijo$mean.ZT18<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT18",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$noise.ZT18<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT18",names(Araport11.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(Araport11.TPMs.genes.Cortijo[,grep("ZT18",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$max.ZT18<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT18",names(Araport11.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

Araport11.TPMs.genes.Cortijo$mean.ZT20<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT20",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$noise.ZT20<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT20",names(Araport11.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(Araport11.TPMs.genes.Cortijo[,grep("ZT20",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$max.ZT20<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT20",names(Araport11.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

Araport11.TPMs.genes.Cortijo$mean.ZT22<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT22",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$noise.ZT22<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT22",names(Araport11.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(Araport11.TPMs.genes.Cortijo[,grep("ZT22",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$max.ZT22<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT22",names(Araport11.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

Araport11.TPMs.genes.Cortijo$mean.ZT24<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT24",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$noise.ZT24<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT24",names(Araport11.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(Araport11.TPMs.genes.Cortijo[,grep("ZT24",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$max.ZT24<-apply(Araport11.TPMs.genes.Cortijo[,grep("ZT24",names(Araport11.TPMs.genes.Cortijo))],1,max,na.rm = TRUE)

Araport11.TPMs.genes.Cortijo$noise_average<-apply(Araport11.TPMs.genes.Cortijo[,grep("noise",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)
Araport11.TPMs.genes.Cortijo$variance_circadian<-apply(Araport11.TPMs.genes.Cortijo[,grep("mean.ZT",names(Araport11.TPMs.genes.Cortijo))],1,sd,na.rm = TRUE)/apply(Araport11.TPMs.genes.Cortijo[,grep("mean.ZT",names(Araport11.TPMs.genes.Cortijo))],1,mean,na.rm = TRUE)



####################################

# Bhaggy - ddm1 knockout rosette Col0 


############################################
#upload bhagys data TPMs
#############################################
denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1 <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1.bed")
rownames(denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1)<-denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1$gene
Araport11.TPMs.genes.Bhagy_WT_ddm1 <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/Araport11.TPMs.genes.Bhagy_WT_ddm1.bed")
rownames(Araport11.TPMs.genes.Bhagy_WT_ddm1)<-Araport11.TPMs.genes.Bhagy_WT_ddm1$gene


denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1$WT<-apply(denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1[,5:7],1,mean)
denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1$ddm1<-apply(denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1[,2:4],1,mean)

Araport11.TPMs.genes.Bhagy_WT_ddm1$WT<-apply(Araport11.TPMs.genes.Bhagy_WT_ddm1[,5:7],1,mean)
Araport11.TPMs.genes.Bhagy_WT_ddm1$ddm1<-apply(Araport11.TPMs.genes.Bhagy_WT_ddm1[,2:4],1,mean)

linc_Bhagy<-denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1[denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1$gene %in% lncRNAs.intergenic.loci$gene,]
pc_Bhagy<-denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1[denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1$gene %in% denovoPC.loci$gene,]
as_Bhagy<-denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1[denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1$gene %in% lncRNAs.antisense.loci$gene,]
te_Bhagy<-denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1[denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1$gene %in% TE_genes.loci$gene,]



#####
# Vu's data - ddm1 and pol4 knockouts in Col-0 different cell types 
# 





############################################
#upload Vu's data TPMs
#############################################
Araport11.TPMs.genes.Vu.SEsamples <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/Araport11.TPMs.genes.Vu.SEsamples.sep_reps.bed")
denovo_Oct2021.TPMs.genes.Vu.SEsamples <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/denovo_Oct2021.TPMs.genes.Vu.SEsamples.sep_reps.bed")
Araport11.TPMs.genes.Vu.PEsamples <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/Araport11.TPMs.genes.Vu.PEsamples.sep_reps.bed")
denovo_Oct2021.TPMs.genes.Vu.PEsamples <- read.delim("03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation/denovo_Oct2021.TPMs.genes.Vu.PEsamples.sep_reps.bed")


###### exclude samples with <2mln reads 
#1.inflorescenceM_HS28D_ddm1.rep4
#2.inflorescenceM_M28D_21.rep3
#3. stem_M28D_poliv.rep4
#4. stem_M7D_21.rep1

excludesamples<-c("inflorescenceM_HS28D_ddm1.rep4","inflorescenceM_M28D_21.rep3","stem_M28D_poliv.rep4","stem_M7D_21.rep1")

Araport11.TPMs.genes.Vu.SEsamples<-Araport11.TPMs.genes.Vu.SEsamples[,!(names(Araport11.TPMs.genes.Vu.SEsamples) %in% excludesamples)]
denovo_Oct2021.TPMs.genes.Vu.SEsamples<-denovo_Oct2021.TPMs.genes.Vu.SEsamples[,!(names(denovo_Oct2021.TPMs.genes.Vu.SEsamples) %in% excludesamples)]
Araport11.TPMs.genes.Vu.PEsamples<-Araport11.TPMs.genes.Vu.PEsamples[,!(names(Araport11.TPMs.genes.Vu.PEsamples) %in% excludesamples)]
denovo_Oct2021.TPMs.genes.Vu.PEsamples<-denovo_Oct2021.TPMs.genes.Vu.PEsamples[,!(names(denovo_Oct2021.TPMs.genes.Vu.PEsamples) %in% excludesamples)]


Araport11.TPMs.genes.Vu.SEsamples$inflor_heat_28D_WT<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "inflorescenceM_HS28D_21",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$inflor_heat_28D_ddm1<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "inflorescenceM_HS28D_ddm1",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$inflor_heat_28D_pol4<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "inflorescenceM_HS28D_poliv",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$inflor_mock_28D_WT<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "inflorescenceM_M28D_21",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$inflor_mock_28D_ddm1<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "inflorescenceM_M28D_ddm1",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$inflor_mock_28D_pol4<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "inflorescenceM_M28D_poliv",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)


Araport11.TPMs.genes.Vu.SEsamples$stem_heat_28D_WT<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "stem_HS28D_21",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$stem_heat_28D_ddm1<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "stem_HS28D_ddm1",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$stem_heat_28D_pol4<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "stem_HS28D_poliv",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$stem_mock_28D_WT<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "stem_M28D_21",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$stem_mock_28D_ddm1<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "stem_M28D_ddm1",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$stem_mock_28D_pol4<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "stem_M28D_poliv",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)


Araport11.TPMs.genes.Vu.SEsamples$seedling_heat_7D_WT<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "seedling_HS7D_21",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$seedling_heat_7D_ddm1<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "seedling_HS7D_ddm1",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$seedling_heat_7D_pol4<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "seedling_HS7D_poliv",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$seedling_mock_7D_WT<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "seedling_M7D_21",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$seedling_mock_7D_ddm1<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "seedling_M7D_ddm1",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)
Araport11.TPMs.genes.Vu.SEsamples$seedling_mock_7D_pol4<-apply(Araport11.TPMs.genes.Vu.SEsamples[,grep(pattern = "seedling_M7D_poliv",x=names(Araport11.TPMs.genes.Vu.SEsamples))],1,mean)

Araport11.TPMs.genes.Vu.PEsamples$stem_heat_7D_WT<-apply(Araport11.TPMs.genes.Vu.PEsamples[,grep(pattern = "stem_HS7D_21",x=names(Araport11.TPMs.genes.Vu.PEsamples))],1,mean)
Araport11.TPMs.genes.Vu.PEsamples$stem_heat_7D_ddm1<-apply(Araport11.TPMs.genes.Vu.PEsamples[,grep(pattern = "stem_HS7D_ddm1",x=names(Araport11.TPMs.genes.Vu.PEsamples))],1,mean)
Araport11.TPMs.genes.Vu.PEsamples$stem_heat_7D_pol4<-apply(Araport11.TPMs.genes.Vu.PEsamples[,grep(pattern = "stem_HS7D_poliv",x=names(Araport11.TPMs.genes.Vu.PEsamples))],1,mean)
Araport11.TPMs.genes.Vu.PEsamples$stem_mock_7D_WT<-apply(Araport11.TPMs.genes.Vu.PEsamples[,grep(pattern = "stem_M7D_21",x=names(Araport11.TPMs.genes.Vu.PEsamples))],1,mean)
Araport11.TPMs.genes.Vu.PEsamples$stem_mock_7D_ddm1<-apply(Araport11.TPMs.genes.Vu.PEsamples[,grep(pattern = "stem_M7D_ddm1",x=names(Araport11.TPMs.genes.Vu.PEsamples))],1,mean)
Araport11.TPMs.genes.Vu.PEsamples$stem_mock_7D_pol4<-apply(Araport11.TPMs.genes.Vu.PEsamples[,grep(pattern = "stem_M7D_poliv",x=names(Araport11.TPMs.genes.Vu.PEsamples))],1,mean)
Araport11.TPMs.genes.Vu.PEsamples$nonstem_heat_7D_WT<-apply(Araport11.TPMs.genes.Vu.PEsamples[,grep(pattern = "non_stem_HS7D_21",x=names(Araport11.TPMs.genes.Vu.PEsamples))],1,mean)
Araport11.TPMs.genes.Vu.PEsamples$nonstem_mock_7D_WT<-apply(Araport11.TPMs.genes.Vu.PEsamples[,grep(pattern = "non_stem_M7D_21",x=names(Araport11.TPMs.genes.Vu.PEsamples))],1,mean)


Vu_Araport11_TPM<-cbind(Araport11.TPMs.genes.Vu.PEsamples[,c(1,61:68)],Araport11.TPMs.genes.Vu.SEsamples[,c(70:87)])

rownames(Vu_Araport11_TPM)<-Vu_Araport11_TPM$gene


denovo_Oct2021.TPMs.genes.Vu.SEsamples$inflor_heat_28D_WT<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "inflorescenceM_HS28D_21",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$inflor_heat_28D_ddm1<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "inflorescenceM_HS28D_ddm1",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$inflor_heat_28D_pol4<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "inflorescenceM_HS28D_poliv",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$inflor_mock_28D_WT<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "inflorescenceM_M28D_21",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$inflor_mock_28D_ddm1<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "inflorescenceM_M28D_ddm1",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$inflor_mock_28D_pol4<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "inflorescenceM_M28D_poliv",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)


denovo_Oct2021.TPMs.genes.Vu.SEsamples$stem_heat_28D_WT<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "stem_HS28D_21",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$stem_heat_28D_ddm1<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "stem_HS28D_ddm1",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$stem_heat_28D_pol4<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "stem_HS28D_poliv",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$stem_mock_28D_WT<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "stem_M28D_21",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$stem_mock_28D_ddm1<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "stem_M28D_ddm1",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$stem_mock_28D_pol4<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "stem_M28D_poliv",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)


denovo_Oct2021.TPMs.genes.Vu.SEsamples$seedling_heat_7D_WT<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "seedling_HS7D_21",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$seedling_heat_7D_ddm1<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "seedling_HS7D_ddm1",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$seedling_heat_7D_pol4<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "seedling_HS7D_poliv",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$seedling_mock_7D_WT<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "seedling_M7D_21",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$seedling_mock_7D_ddm1<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "seedling_M7D_ddm1",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.SEsamples$seedling_mock_7D_pol4<-apply(denovo_Oct2021.TPMs.genes.Vu.SEsamples[,grep(pattern = "seedling_M7D_poliv",x=names(denovo_Oct2021.TPMs.genes.Vu.SEsamples))],1,mean)

#inflorescenceM_HS28D_21
#inflorescenceM_HS28D_ddm1
#inflorescenceM_HS28D_poliv
#inflorescenceM_M28D_21
#inflorescenceM_M28D_ddm1
#inflorescenceM_M28D_poliv
#seedling_HS7D_21
#seedling_HS7D_ddm1
#seedling_HS7D_poliv
#seedling_M7D_21
#seedling_M7D_ddm1
#seedling_M7D_poliv
#stem_HS28D_21
#stem_HS28D_ddm1
#stem_HS28D_poliv
#stem_M28D_21
#stem_M28D_ddm1
#stem_M28D_poliv

#non_stem_HS7D_21
#non_stem_M7D_21
#stem_HS7D_21
#stem_HS7D_ddm1
#stem_HS7D_poliv
#stem_M7D_21
#stem_M7D_ddm1
#stem_M7D_poliv

denovo_Oct2021.TPMs.genes.Vu.PEsamples$stem_heat_7D_WT<-apply(denovo_Oct2021.TPMs.genes.Vu.PEsamples[,grep(pattern = "stem_HS7D_21",x=names(denovo_Oct2021.TPMs.genes.Vu.PEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.PEsamples$stem_heat_7D_ddm1<-apply(denovo_Oct2021.TPMs.genes.Vu.PEsamples[,grep(pattern = "stem_HS7D_ddm1",x=names(denovo_Oct2021.TPMs.genes.Vu.PEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.PEsamples$stem_heat_7D_pol4<-apply(denovo_Oct2021.TPMs.genes.Vu.PEsamples[,grep(pattern = "stem_HS7D_poliv",x=names(denovo_Oct2021.TPMs.genes.Vu.PEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.PEsamples$stem_mock_7D_WT<-apply(denovo_Oct2021.TPMs.genes.Vu.PEsamples[,grep(pattern = "stem_M7D_21",x=names(denovo_Oct2021.TPMs.genes.Vu.PEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.PEsamples$stem_mock_7D_ddm1<-apply(denovo_Oct2021.TPMs.genes.Vu.PEsamples[,grep(pattern = "stem_M7D_ddm1",x=names(denovo_Oct2021.TPMs.genes.Vu.PEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.PEsamples$stem_mock_7D_pol4<-apply(denovo_Oct2021.TPMs.genes.Vu.PEsamples[,grep(pattern = "stem_M7D_poliv",x=names(denovo_Oct2021.TPMs.genes.Vu.PEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.PEsamples$nonstem_heat_7D_WT<-apply(denovo_Oct2021.TPMs.genes.Vu.PEsamples[,grep(pattern = "non_stem_HS7D_21",x=names(denovo_Oct2021.TPMs.genes.Vu.PEsamples))],1,mean)
denovo_Oct2021.TPMs.genes.Vu.PEsamples$nonstem_mock_7D_WT<-apply(denovo_Oct2021.TPMs.genes.Vu.PEsamples[,grep(pattern = "non_stem_M7D_21",x=names(denovo_Oct2021.TPMs.genes.Vu.PEsamples))],1,mean)


Vu_denovo_TPM<-cbind(denovo_Oct2021.TPMs.genes.Vu.PEsamples[,c(1,61:68)],denovo_Oct2021.TPMs.genes.Vu.SEsamples[,c(70:87)])

rownames(Vu_denovo_TPM)<-Vu_denovo_TPM$gene

linc_Vu<-Vu_denovo_TPM[Vu_denovo_TPM$gene %in% lncRNAs.intergenic.loci$gene,]
pc_Vu<-Vu_denovo_TPM[Vu_denovo_TPM$gene %in% denovoPC.loci$gene,]
as_Vu<-Vu_denovo_TPM[Vu_denovo_TPM$gene %in% lncRNAs.antisense.loci$gene,]
te_Vu<-Vu_denovo_TPM[Vu_denovo_TPM$gene %in% TE_genes.loci$gene,]














# small RNA data 
sRNA.24nt.Ar11.RPM <- read.delim("03_Projects/RNA-seq_data/miRNA/flowers_14acc.24nt.Araport11.coverage.perbase_calc.normalized.bed")
sRNA.24nt.Ar11.RPM<-sRNA.24nt.Ar11.RPM[!duplicated(sRNA.24nt.Ar11.RPM$gene),]
sRNA.21nt.Ar11.RPM <- read.delim("03_Projects/RNA-seq_data/miRNA/flowers_14acc.21-22nt.Araport11.coverage.perbase_calc.normalized.bed")
sRNA.21nt.Ar11.RPM<-sRNA.21nt.Ar11.RPM[!duplicated(sRNA.21nt.Ar11.RPM$gene),]

sRNA.24nt.denovo2021.RPM <- read.delim("03_Projects/RNA-seq_data/miRNA/flowers_14acc.24nt.denovo_Oct2021.coverage.perbase_calc.normalized.bed")
sRNA.21nt.denovo2021.RPM <- read.delim("03_Projects/RNA-seq_data/miRNA/flowers_14acc.21-22nt.denovo_Oct2021.coverage.perbase_calc.normalized.bed")
rownames(sRNA.24nt.Ar11.RPM)<-sRNA.24nt.Ar11.RPM$gene
rownames(sRNA.21nt.Ar11.RPM)<-sRNA.21nt.Ar11.RPM$gene
rownames(sRNA.24nt.denovo2021.RPM)<-sRNA.24nt.denovo2021.RPM$gene
rownames(sRNA.21nt.denovo2021.RPM)<-sRNA.21nt.denovo2021.RPM$gene

sRNA.24nt.Ar11.RPM$mean<-apply(sRNA.24nt.Ar11.RPM[,7:20],1,mean)
sRNA.21nt.Ar11.RPM$mean<-apply(sRNA.21nt.Ar11.RPM[,7:20],1,mean)
sRNA.24nt.denovo2021.RPM$mean<-apply(sRNA.24nt.denovo2021.RPM[,7:20],1,mean)
sRNA.21nt.denovo2021.RPM$mean<-apply(sRNA.21nt.denovo2021.RPM[,7:20],1,mean)


sRNA.24nt.Ar11.RPM$max<-apply(sRNA.24nt.Ar11.RPM[,7:20],1,max)
sRNA.21nt.Ar11.RPM$max<-apply(sRNA.21nt.Ar11.RPM[,7:20],1,max)
sRNA.24nt.denovo2021.RPM$max<-apply(sRNA.24nt.denovo2021.RPM[,7:20],1,max)
sRNA.21nt.denovo2021.RPM$max<-apply(sRNA.21nt.denovo2021.RPM[,7:20],1,max)

sRNA.24nt.Ar11.RPM$variance<-apply(sRNA.24nt.Ar11.RPM[,7:20],1,sd)/apply(sRNA.24nt.Ar11.RPM[,7:20],1,mean)
sRNA.21nt.Ar11.RPM$variance<-apply(sRNA.21nt.Ar11.RPM[,7:20],1,sd)/apply(sRNA.21nt.Ar11.RPM[,7:20],1,mean)
sRNA.24nt.denovo2021.RPM$variance<-apply(sRNA.24nt.denovo2021.RPM[,7:20],1,sd)/apply(sRNA.24nt.denovo2021.RPM[,7:20],1,mean)
sRNA.21nt.denovo2021.RPM$variance<-apply(sRNA.21nt.denovo2021.RPM[,7:20],1,sd)/apply(sRNA.21nt.denovo2021.RPM[,7:20],1,mean)




sRNA.24nt.denovo2021.RPM.Ranj <- read.delim("03_Projects/Public_data/small_RNA_and_Methyl_Seq_Papareddy_et_al_2020/Ranj.sRNA.24nt.denovo_Oct2021.coverage.perbase_calc.normalized.bed")
a<-sRNA.24nt.denovo2021.RPM.Ranj
a$WT.eheart<-apply(a[,grep("WT.eheart",names(a))],1,mean)
a$nrpda3.eheart<-apply(a[,grep("nrpda3.eheart",names(a))],1,mean)

a$WT.leaf<-apply(a[,grep("WT.leaf",names(a))],1,mean)
a$nrpda3.leaf<-apply(a[,grep("nrpda3.leaf",names(a))],1,mean)
a$dcl234.leaf<-apply(a[,grep("dcl234.leaf",names(a))],1,mean)
a$Hon12.leaf<-apply(a[,grep("Hon12.leaf",names(a))],1,mean)
a$suvh456.leaf<-apply(a[,grep("suvh456.leaf",names(a))],1,mean)

a$nrpda3.fb<-apply(a[,grep("nrpda3.fb",names(a))],1,mean)
a$dcl234.fb<-a$dcl234.fb.r3
a$dcl234.gl<-apply(a[,grep("dcl234.gl",names(a))],1,mean)
sRNA.24nt.denovo2021.RPM.Ranj<-a

rownames(sRNA.24nt.denovo2021.RPM.Ranj)<-sRNA.24nt.denovo2021.RPM.Ranj$gene




sRNA.2122nt.denovo2021.RPM.Ranj <- read.delim("03_Projects/Public_data/small_RNA_and_Methyl_Seq_Papareddy_et_al_2020/Ranj.sRNA.21-22nt.denovo_Oct2021.coverage.perbase_calc.normalized.bed")
a<-sRNA.2122nt.denovo2021.RPM.Ranj
a$WT.eheart<-apply(a[,grep("WT.eheart",names(a))],1,mean)
a$nrpda3.eheart<-apply(a[,grep("nrpda3.eheart",names(a))],1,mean)

a$WT.leaf<-apply(a[,grep("WT.leaf",names(a))],1,mean)
a$nrpda3.leaf<-apply(a[,grep("nrpda3.leaf",names(a))],1,mean)
a$dcl234.leaf<-apply(a[,grep("dcl234.leaf",names(a))],1,mean)
a$Hon12.leaf<-apply(a[,grep("Hon12.leaf",names(a))],1,mean)
a$suvh456.leaf<-apply(a[,grep("suvh456.leaf",names(a))],1,mean)

a$nrpda3.fb<-apply(a[,grep("nrpda3.fb",names(a))],1,mean)
a$dcl234.fb<-a$dcl234.fb.r3
a$dcl234.gl<-apply(a[,grep("dcl234.gl",names(a))],1,mean)
sRNA.2122nt.denovo2021.RPM.Ranj<-a
rownames(sRNA.2122nt.denovo2021.RPM.Ranj)<-sRNA.2122nt.denovo2021.RPM.Ranj$gene



sRNA.2122nt.denovo2021.RPM.Ranj <- read.delim("03_Projects/Public_data/small_RNA_and_Methyl_Seq_Papareddy_et_al_2020/Ranj.sRNA.21-22nt.denovo_Oct2021.coverage.perbase_calc.normalized.bed")
a<-sRNA.2122nt.denovo2021.RPM.Ranj
a$WT.eheart<-apply(a[,grep("WT.eheart",names(a))],1,mean)
a$nrpda3.eheart<-apply(a[,grep("nrpda3.eheart",names(a))],1,mean)

a$WT.leaf<-apply(a[,grep("WT.leaf",names(a))],1,mean)
a$nrpda3.leaf<-apply(a[,grep("nrpda3.leaf",names(a))],1,mean)
a$dcl234.leaf<-apply(a[,grep("dcl234.leaf",names(a))],1,mean)
a$Hon12.leaf<-apply(a[,grep("Hon12.leaf",names(a))],1,mean)
a$suvh456.leaf<-apply(a[,grep("suvh456.leaf",names(a))],1,mean)

a$nrpda3.fb<-apply(a[,grep("nrpda3.fb",names(a))],1,mean)
a$dcl234.fb<-a$dcl234.fb.r3
a$dcl234.gl<-apply(a[,grep("dcl234.gl",names(a))],1,mean)
sRNA.2122nt.denovo2021.RPM.Ranj<-a
rownames(sRNA.2122nt.Araport11.RPM.Ranj)<-sRNA.2122nt.Araport11.RPM.Ranj$gene




sRNA.24nt.Araport11.RPM.Ranj <- read.delim("03_Projects/Public_data/small_RNA_and_Methyl_Seq_Papareddy_et_al_2020/Ranj.sRNA.24nt.Araport11.coverage.perbase_calc.normalized.bed")
a<-sRNA.24nt.Araport11.RPM.Ranj
a$WT.eheart<-apply(a[,grep("WT.eheart",names(a))],1,mean)
a$nrpda3.eheart<-apply(a[,grep("nrpda3.eheart",names(a))],1,mean)

a$WT.leaf<-apply(a[,grep("WT.leaf",names(a))],1,mean)
a$nrpda3.leaf<-apply(a[,grep("nrpda3.leaf",names(a))],1,mean)
a$dcl234.leaf<-apply(a[,grep("dcl234.leaf",names(a))],1,mean)
a$Hon12.leaf<-apply(a[,grep("Hon12.leaf",names(a))],1,mean)
a$suvh456.leaf<-apply(a[,grep("suvh456.leaf",names(a))],1,mean)

a$nrpda3.fb<-apply(a[,grep("nrpda3.fb",names(a))],1,mean)
a$dcl234.fb<-a$dcl234.fb.r3
a$dcl234.gl<-apply(a[,grep("dcl234.gl",names(a))],1,mean)
sRNA.24nt.Araport11.RPM.Ranj<-a




sRNA.2122nt.Araport11.RPM.Ranj <- read.delim("03_Projects/Public_data/small_RNA_and_Methyl_Seq_Papareddy_et_al_2020/Ranj.sRNA.21-22nt.Araport11.coverage.perbase_calc.normalized.bed")
a<-sRNA.2122nt.Araport11.RPM.Ranj
a$WT.eheart<-apply(a[,grep("WT.eheart",names(a))],1,mean)
a$nrpda3.eheart<-apply(a[,grep("nrpda3.eheart",names(a))],1,mean)

a$WT.leaf<-apply(a[,grep("WT.leaf",names(a))],1,mean)
a$nrpda3.leaf<-apply(a[,grep("nrpda3.leaf",names(a))],1,mean)
a$dcl234.leaf<-apply(a[,grep("dcl234.leaf",names(a))],1,mean)
a$Hon12.leaf<-apply(a[,grep("Hon12.leaf",names(a))],1,mean)
a$suvh456.leaf<-apply(a[,grep("suvh456.leaf",names(a))],1,mean)

a$nrpda3.fb<-apply(a[,grep("nrpda3.fb",names(a))],1,mean)
a$dcl234.fb<-a$dcl234.fb.r3
a$dcl234.gl<-apply(a[,grep("dcl234.gl",names(a))],1,mean)
sRNA.2122nt.Araport11.RPM.Ranj<-a
rownames(sRNA.2122nt.Araport11.RPM.Ranj)<-sRNA.2122nt.Araport11.RPM.Ranj$gene
