setwd("Z:/01_POSTDOC/")



#upload public annotations 

Araport11_protein_coding.201606.genes <- read.delim("04_writing/lncRNA variation paper/Annotation/Araport11_protein_coding.201606.genes.bed", header=FALSE)
names(Araport11_protein_coding.201606.genes)<-c("chr","start","end","gene","score","strand")


Araport11_protein_coding <- read.delim("04_writing/lncRNA variation paper/Annotation/Araport11_protein_coding.201606.bed", header=FALSE)
names(Araport11_protein_coding)<-c("chr","start","end","transcript","score","strand","start1","end1","rgb","exon_N","exon_length","exon_start")
Araport11_protein_coding$gene<-lapply(strsplit(as.character(Araport11_protein_coding$transcript),".", fixed = T), "[",1)


Araport11_non_coding<- read.delim("04_writing/lncRNA variation paper/Annotation/Araport11_non_coding.2016016.sorted.bed", header=FALSE)
names(Araport11_non_coding)<-c("chr","start","end","transcript","score","strand","start1","end1","rgb","exon_N","exon_length","exon_start")
Araport11_non_coding$gene<-lapply(strsplit(as.character(Araport11_non_coding$transcript),".", fixed = T), "[",1)


Araport11_novel_transcribed_region<- read.delim("04_writing/lncRNA variation paper/Annotation/Araport11_novel_transcribed_region.201606.bed", header=FALSE)
names(Araport11_novel_transcribed_region)<-c("chr","start","end","gene","score","strand","start1","end1","rgb","exon_N","exon_length","exon_start")


Araport11_transposable_element_gene<- read.delim("04_writing/lncRNA variation paper/Annotation/Araport11_transposable_element_gene.201606.bed", header=FALSE)
names(Araport11_transposable_element_gene)<-c("chr","start","end","transcript","score","strand","start1","end1","rgb","exon_N","exon_length","exon_start")
Araport11_transposable_element_gene$gene<-lapply(strsplit(as.character(Araport11_transposable_element_gene$transcript),".", fixed = T), "[",1)



Araport11_transposable_element_f<- read.delim("04_writing/lncRNA variation paper/Annotation/Araport11_TEs.transposon_fragments.bed", header=FALSE)
names(Araport11_transposable_element_f)<-c("chr","start","end","gene","score","strand")


#upload denovo annotations 




# novel and known lncRNAs 
denovoNC.loci.novel <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/denovoNC.loci.novel.bed", header=FALSE)
names(denovoNC.loci.novel)<-c("chr","start","end","gene","score","strand")
denovoNC.loci.novel.names<-as.data.frame(denovoNC.loci.novel[,4])

denovoNC.loci.known <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/denovoNC.loci.known.bed", header=FALSE)
names(denovoNC.loci.known)<-c("chr","start","end","gene","score","strand")
denovoNC.loci.known.names<-as.data.frame(denovoNC.loci.known[,4])

denovoNC.transcripts.novel <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/denovoNC.transcripts.novel.bed", header=FALSE)
names(denovoNC.transcripts.novel)<-c("chr","start","end","gene","score","strand","start1","end1","rgb","exon_N","exon_length","exon_start")
denovoNC.transcripts.known <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/denovoNC.transcripts.known.bed", header=FALSE)
names(denovoNC.transcripts.known)<-c("chr","start","end","gene","score","strand","start1","end1","rgb","exon_N","exon_length","exon_start")


#6 types of lncRNAs
#AS
lncRNAs.antisense.loci <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/lncRNAs.antisense.loci.bed", header=FALSE)
names(lncRNAs.antisense.loci)<-c("chr","start","end","gene","score","strand")
lncRNAs.antisense.transcripts <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/lncRNAs.antisense.transcripts.bed", header=FALSE)
names(lncRNAs.antisense.transcripts)<-c("chr","start","end","gene","score","strand","start1","end1","rgb","exon_N","exon_length","exon_start")
lncRNAs.antisense.loci$length<-lncRNAs.antisense.loci$end-lncRNAs.antisense.loci$start


#AS to TE genes
lncRNAs.AS_to_TE.loci <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/lncRNAs.AS_to_TE.loci.bed", header=FALSE)
names(lncRNAs.AS_to_TE.loci)<-c("chr","start","end","gene","score","strand")
lncRNAs.AS_to_TE.transcripts <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/lncRNAs.AS_to_TE.transcripts.bed", header=FALSE)
names(lncRNAs.AS_to_TE.transcripts)<-c("chr","start","end","gene","score","strand","start1","end1","rgb","exon_N","exon_length","exon_start")
lncRNAs.AS_to_TE.loci$length<-lncRNAs.AS_to_TE.loci$end-lncRNAs.AS_to_TE.loci$start



#AS to pseudogenes
lncRNAs.AS_to_pseudo.loci <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/lncRNAs.AS_to_pseudogenes.loci.bed", header=FALSE)
names(lncRNAs.AS_to_pseudo.loci)<-c("chr","start","end","gene","score","strand")
lncRNAs.AS_to_pseudo.transcripts <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/lncRNAs.AS_to_pseudogenes.transcripts.bed", header=FALSE)
names(lncRNAs.AS_to_pseudo.transcripts)<-c("chr","start","end","gene","score","strand","start1","end1","rgb","exon_N","exon_length","exon_start")


#lincRNAs
lncRNAs.intergenic.loci <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/lncRNAs.intergenic.loci.bed", header=FALSE)
names(lncRNAs.intergenic.loci)<-c("chr","start","end","gene","score","strand")
lncRNAs.intergenic.transcripts <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/lncRNAs.intergenic.transcripts.bed", header=FALSE)
names(lncRNAs.intergenic.transcripts)<-c("chr","start","end","gene","score","strand","start1","end1","rgb","exon_N","exon_length","exon_start")
lncRNAs.intergenic.loci$length<-lncRNAs.intergenic.loci$end-lncRNAs.intergenic.loci$start

#expressed TE genes
TE_genes.transcripts <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/transcripts.TE_genes.bed", header=FALSE)
names(TE_genes.transcripts)<-c("chr","start","end","transcript","score","strand","start1","end1","rgb","exon_N","exon_size","exon_loc")
TE_genes.transcripts$gene<-paste(lapply(strsplit(as.character(TE_genes.transcripts$transcript),".", fixed = T), "[",1),lapply(strsplit(as.character(TE_genes.transcripts$transcript),".", fixed = T), "[",2),sep = ".")
TE_genes.loci<-read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/loci.TE_genes.bed", header=FALSE)
names(TE_genes.loci)<-c("chr","start","end","gene","score","strand")

TE_genes.loci$length<-TE_genes.loci$end-TE_genes.loci$start

#expressed TE frags
TE_frags.transcripts <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/transcripts.TE_frags.bed", header=FALSE)
names(TE_frags.transcripts)<-c("chr","start","end","transcript","score","strand","start1","end1","rgb","exon_N","exon_size","exon_loc")
TE_frags.transcripts$gene<-paste(lapply(strsplit(as.character(TE_frags.transcripts$transcript),".", fixed = T), "[",1),lapply(strsplit(as.character(TE_frags.transcripts$transcript),".", fixed = T), "[",2),sep = ".")
TE_frags.transcripts$length<-TE_frags.transcripts$end-TE_frags.transcripts$start


#expressed pseudogenes 
denovo_pseudogene.transcripts <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/transcripts.f1.denovo_pseudogene.bed", header=FALSE)
names(denovo_pseudogene.transcripts)<-c("chr","start","end","transcript","score","strand","start1","end1","rgb","exon_N","exon_size","exon_loc")
denovo_pseudogene.transcripts$gene<-paste(lapply(strsplit(as.character(denovo_pseudogene.transcripts$transcript),".", fixed = T), "[",1),lapply(strsplit(as.character(denovo_pseudogene.transcripts$transcript),".", fixed = T), "[",2),sep = ".")


#PC genes
denovoPC.loci <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/denovoPC.loci.bed", header=FALSE)
names(denovoPC.loci)<-c("chr","start","end","gene","score","strand")
denovoPC.transcripts <- read.delim("03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/20211013_annotation/denovoPC.transcripts.bed", header=FALSE)
names(denovoPC.transcripts)<-c("chr","start","end","transcript","score","strand","start1","end1","rgb","exon_N","exon_size","exon_loc")
denovoPC.loci$length<-denovoPC.loci$end-denovoPC.loci$start



#upload Araport11 annotations 

Ar11_TE_genes<-read.delim("04_writing/lncRNA variation paper/Annotation/Araport11_transposable_element_gene.201606.bed", header=FALSE)
names(Ar11_TE_genes)<-c("chr","start","end","transcript","score","strand","start1","end1","rgb","exon_N","exon_size","exon_loc")
Ar11_TE_genes$gene<-lapply(strsplit(as.character(Ar11_TE_genes$transcript),".", fixed = T), "[",1)


Ar11_TE_fragments<-read.delim("04_writing/lncRNA variation paper/Annotation/Araport11_TEs.transposon_fragments.bed", header=FALSE)
names(Ar11_TE_fragments)<-c("chr","start","end","transcript","score","strand")

Ar11_TE_frag_types<-read.delim("04_writing/lncRNA variation paper/Annotation/TAIR10_Transposable_Elements_8superfamilies.txt", header=FALSE)
names(Ar11_TE_frag_types)<-c("TE","TE_family")

#TAIR10_TEs<-("03_Projects/TAIR10_Transposable_Elements.bed", header=FALSE)
#names(TAIR10_TEs)<-c("chr","start","end","transcript","score","strand")



#add distance to centromeres 
a<-lncRNAs.antisense.loci

a$dist_from_centromere<-0
a$dist_from_centromere[a$chr=="Chr1"]<-abs(a$start[a$chr=="Chr1"]-15000000)
a$dist_from_centromere[a$chr=="Chr2"]<-abs(a$start[a$chr=="Chr2"]-4700000)
a$dist_from_centromere[a$chr=="Chr3"]<-abs(a$start[a$chr=="Chr3"]-13000000)
a$dist_from_centromere[a$chr=="Chr4"]<-abs(a$start[a$chr=="Chr4"]-3800000)
a$dist_from_centromere[a$chr=="Chr5"]<-abs(a$start[a$chr=="Chr5"]-11900000)
lncRNAs.antisense.loci<-a

a<-lncRNAs.intergenic.loci
a$dist_from_centromere<-0
a$dist_from_centromere[a$chr=="Chr1"]<-abs(a$start[a$chr=="Chr1"]-15000000)
a$dist_from_centromere[a$chr=="Chr2"]<-abs(a$start[a$chr=="Chr2"]-4700000)
a$dist_from_centromere[a$chr=="Chr3"]<-abs(a$start[a$chr=="Chr3"]-13000000)
a$dist_from_centromere[a$chr=="Chr4"]<-abs(a$start[a$chr=="Chr4"]-3800000)
a$dist_from_centromere[a$chr=="Chr5"]<-abs(a$start[a$chr=="Chr5"]-11900000)
lncRNAs.intergenic.loci<-a


a<-denovoPC.loci
a$dist_from_centromere<-0
a$dist_from_centromere[a$chr=="Chr1"]<-abs(a$start[a$chr=="Chr1"]-15000000)
a$dist_from_centromere[a$chr=="Chr2"]<-abs(a$start[a$chr=="Chr2"]-4700000)
a$dist_from_centromere[a$chr=="Chr3"]<-abs(a$start[a$chr=="Chr3"]-13000000)
a$dist_from_centromere[a$chr=="Chr4"]<-abs(a$start[a$chr=="Chr4"]-3800000)
a$dist_from_centromere[a$chr=="Chr5"]<-abs(a$start[a$chr=="Chr5"]-11900000)
denovoPC.loci<-a

a<-TE_genes.loci
a$dist_from_centromere<-0
a$dist_from_centromere[a$chr=="Chr1"]<-abs(a$start[a$chr=="Chr1"]-15000000)
a$dist_from_centromere[a$chr=="Chr2"]<-abs(a$start[a$chr=="Chr2"]-4700000)
a$dist_from_centromere[a$chr=="Chr3"]<-abs(a$start[a$chr=="Chr3"]-13000000)
a$dist_from_centromere[a$chr=="Chr4"]<-abs(a$start[a$chr=="Chr4"]-3800000)
a$dist_from_centromere[a$chr=="Chr5"]<-abs(a$start[a$chr=="Chr5"]-11900000)
TE_genes.loci<-a