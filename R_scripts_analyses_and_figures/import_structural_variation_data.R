setwd("/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/")
setwd("Z:/01_POSTDOC/")

#copy number 

CN_as_27genomes<- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/CN_tables/denovo_AS_Oct2021_CN_eracaps_genomes.processed.bed")
CN_te_27genomes<- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/CN_tables/denovo_TEgenes_Oct2021_CN_eracaps_genomes.processed.bed")
CN_linc_27genomes<- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/CN_tables/denovo_lincRNAs_Oct2021_CN_eracaps_genomes.processed.bed")
CN_linc_27genomes_lessstrict<- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/CN_tables/denovo_lincRNAs_Oct2021_CN_eracaps_genomes.40_80.processed.bed")
CN_pc_27genomes<- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/CN_tables/denovo_PC_Oct2021_CN_eracaps_genomes.processed.bed")

CN_ar11pc_27genomes<- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/CN_tables/Ar11_PC_Oct2021_CN_eracaps_genomes.processed.bed")


CN_linc_27genomes$gene<-rownames(CN_linc_27genomes)
CN_linc_27genomes_lessstrict$gene<-rownames(CN_linc_27genomes_lessstrict)
CN_as_27genomes$gene<-rownames(CN_as_27genomes)
CN_te_27genomes$gene<-rownames(CN_te_27genomes)
CN_pc_27genomes$gene<-rownames(CN_pc_27genomes)
CN_ar11pc_27genomes$gene<-rownames(CN_ar11pc_27genomes)





#upload genetic var stats 

linc_Oct2021_genetic_var_stats <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/genestats_tables/LINC.genetic_var_stats.txt")
linc_Oct2021_genetic_var_stats$snp_per_kb<-linc_Oct2021_genetic_var_stats$nr_of_snps_in_locus*1000/linc_Oct2021_genetic_var_stats$length

linc_TSS_Oct2021_genetic_var_stats <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/genestats_tables/LINC.TSS_plusminus200bp.genetic_var_stats.txt")
linc_TSS_Oct2021_genetic_var_stats$snp_per_kb<-linc_TSS_Oct2021_genetic_var_stats$nr_of_snps_in_locus*1000/linc_TSS_Oct2021_genetic_var_stats$length

linc_TES_Oct2021_genetic_var_stats <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/genestats_tables/LINC.TES_plusminus200bp.genetic_var_stats.txt")
linc_TES_Oct2021_genetic_var_stats$snp_per_kb<-linc_TES_Oct2021_genetic_var_stats$nr_of_snps_in_locus*1000/linc_TES_Oct2021_genetic_var_stats$length

te_Oct2021_genetic_var_stats <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/genestats_tables/TEgene.genetic_var_stats.txt")
te_Oct2021_genetic_var_stats$snp_per_kb<-te_Oct2021_genetic_var_stats$nr_of_snps_in_locus*1000/te_Oct2021_genetic_var_stats$length

te_TSS_Oct2021_genetic_var_stats <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/genestats_tables/TEgene.TSS_plusminus200bp.genetic_var_stats.txt")
te_TSS_Oct2021_genetic_var_stats$snp_per_kb<-te_TSS_Oct2021_genetic_var_stats$nr_of_snps_in_locus*1000/te_TSS_Oct2021_genetic_var_stats$length

te_TES_Oct2021_genetic_var_stats <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/genestats_tables/TEgene.TES_plusminus200bp.genetic_var_stats.txt")
te_TES_Oct2021_genetic_var_stats$snp_per_kb<-te_TES_Oct2021_genetic_var_stats$nr_of_snps_in_locus*1000/te_TES_Oct2021_genetic_var_stats$length


as_Oct2021_genetic_var_stats <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/genestats_tables/AS.genetic_var_stats.txt")
as_Oct2021_genetic_var_stats$snp_per_kb<-as.numeric(as_Oct2021_genetic_var_stats$nr_of_snps_in_locus)*1000/as.numeric(as_Oct2021_genetic_var_stats$length)

as_TSS_Oct2021_genetic_var_stats <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/genestats_tables/AS.TSS_plusminus200bp.genetic_var_stats.txt")
as_TSS_Oct2021_genetic_var_stats$snp_per_kb<-as_TSS_Oct2021_genetic_var_stats$nr_of_snps_in_locus*1000/as_TSS_Oct2021_genetic_var_stats$length

as_TES_Oct2021_genetic_var_stats <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/genestats_tables/AS.TES_plusminus200bp.genetic_var_stats.txt")
as_TES_Oct2021_genetic_var_stats$snp_per_kb<-as_TES_Oct2021_genetic_var_stats$nr_of_snps_in_locus*1000/as_TES_Oct2021_genetic_var_stats$length

pc_Oct2021_genetic_var_stats <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/genestats_tables/PC.genetic_var_stats.txt")
pc_Oct2021_genetic_var_stats$snp_per_kb<-as.numeric(pc_Oct2021_genetic_var_stats$nr_of_snps_in_locus)*1000/as.numeric(pc_Oct2021_genetic_var_stats$length)

pc_TSS_Oct2021_genetic_var_stats <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/genestats_tables/PC.TSS_plusminus200bp.genetic_var_stats.txt")
pc_TSS_Oct2021_genetic_var_stats$snp_per_kb<-as.numeric(pc_TSS_Oct2021_genetic_var_stats$nr_of_snps_in_locus)*1000/as.numeric(pc_TSS_Oct2021_genetic_var_stats$length)

pc_TES_Oct2021_genetic_var_stats <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/genestats_tables/PC.TES_plusminus200bp.genetic_var_stats.txt")
pc_TES_Oct2021_genetic_var_stats$snp_per_kb<-as.numeric(pc_TES_Oct2021_genetic_var_stats$nr_of_snps_in_locus)*1000/as.numeric(pc_TES_Oct2021_genetic_var_stats$length)



rownames(linc_Oct2021_genetic_var_stats)<-linc_Oct2021_genetic_var_stats$gene
rownames(linc_TES_Oct2021_genetic_var_stats)<-linc_TES_Oct2021_genetic_var_stats$gene
rownames(linc_TSS_Oct2021_genetic_var_stats)<-linc_TSS_Oct2021_genetic_var_stats$gene

rownames(as_Oct2021_genetic_var_stats)<-as_Oct2021_genetic_var_stats$gene
rownames(as_TES_Oct2021_genetic_var_stats)<-as_TES_Oct2021_genetic_var_stats$gene
rownames(as_TSS_Oct2021_genetic_var_stats)<-as_TSS_Oct2021_genetic_var_stats$gene

rownames(te_Oct2021_genetic_var_stats)<-te_Oct2021_genetic_var_stats$gene
rownames(te_TES_Oct2021_genetic_var_stats)<-te_TES_Oct2021_genetic_var_stats$gene
rownames(te_TSS_Oct2021_genetic_var_stats)<-te_TSS_Oct2021_genetic_var_stats$gene

rownames(pc_Oct2021_genetic_var_stats)<-pc_Oct2021_genetic_var_stats$gene
rownames(pc_TES_Oct2021_genetic_var_stats)<-pc_TES_Oct2021_genetic_var_stats$gene
rownames(pc_TSS_Oct2021_genetic_var_stats)<-pc_TSS_Oct2021_genetic_var_stats$gene





length( linc_Oct2021_genetic_var_stats$nr_lines_with_besthit_on_correct_chr)
#2238




#upload TE insert data 
lincs.27genomes.TEsizes.forward <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/aug2022_TE_tables/LINC.TEinserts.TEsizes.F.bed")
lincs.27genomes.TEsizes.reverse <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/aug2022_TE_tables/LINC.TEinserts.TEsizes.R.bed")
lincs.27genomes.TEtypes.forward <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/aug2022_TE_tables/LINC.TEinserts.TEtypes.F.bed")
lincs.27genomes.TEtypes.reverse <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/aug2022_TE_tables/LINC.TEinserts.TEtypes.R.bed")
lincs.27genomes.TEtypes <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/aug2022_TE_tables/LINC.TEinserts.TEtypes.bed")
lincs.27genomes.TEsizes <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/aug2022_TE_tables/LINC.TEinserts.TEsizes.bed")


lincs.27genomes.TEsizes.forward[lincs.27genomes.TEsizes.forward==""]<-0
lincs.27genomes.TEsizes.reverse[lincs.27genomes.TEsizes.reverse==""]<-0
lincs.27genomes.TEtypes.forward[lincs.27genomes.TEtypes.forward==""]<-"no_TE"
lincs.27genomes.TEtypes.reverse[lincs.27genomes.TEtypes.reverse==""]<-"no_TE"
rownames(lincs.27genomes.TEsizes.forward)<-lincs.27genomes.TEsizes.forward$gene
rownames(lincs.27genomes.TEsizes.reverse)<-lincs.27genomes.TEsizes.reverse$gene
rownames(lincs.27genomes.TEtypes.forward)<-lincs.27genomes.TEtypes.forward$gene
rownames(lincs.27genomes.TEtypes.reverse)<-lincs.27genomes.TEtypes.reverse$gene



# PC 

pc.27genomes.TEsizes.forward <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/denovoPC_Oct2021_each_acc_TEinserts.forward.TEsizes.bed")
pc.27genomes.TEsizes.reverse <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/denovoPC_Oct2021_each_acc_TEinserts.reverse.TEsizes.bed")
pc.27genomes.TEtypes.forward <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/denovoPC_Oct2021_each_acc_TEinserts.forward.TEtypes.bed")
pc.27genomes.TEtypes.reverse <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/denovoPC_Oct2021_each_acc_TEinserts.reverse.TEtypes.bed")

rownames(pc.27genomes.TEsizes.forward)<-pc.27genomes.TEsizes.forward$gene
rownames(pc.27genomes.TEsizes.reverse)<-pc.27genomes.TEsizes.reverse$gene
rownames(pc.27genomes.TEtypes.forward)<-pc.27genomes.TEtypes.forward$gene
rownames(pc.27genomes.TEtypes.reverse)<-pc.27genomes.TEtypes.reverse$gene


pc.27genomes.TEsizes.forward[pc.27genomes.TEsizes.forward==""]<-0
pc.27genomes.TEsizes.reverse[pc.27genomes.TEsizes.reverse==""]<-0
pc.27genomes.TEtypes.forward[pc.27genomes.TEtypes.forward==""]<-"no_TE"
pc.27genomes.TEtypes.reverse[pc.27genomes.TEtypes.reverse==""]<-"no_TE"



pcTSS.27genomes.TEsizes.forward <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/denovoPC_TSS_Oct2021_each_acc_TEinserts.forward.TEsizes.bed")
pcTSS.27genomes.TEsizes.reverse <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/denovoPC_TSS_Oct2021_each_acc_TEinserts.reverse.TEsizes.bed")
pcTSS.27genomes.TEtypes.forward <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/denovoPC_TSS_Oct2021_each_acc_TEinserts.forward.TEtypes.bed")
pcTSS.27genomes.TEtypes.reverse <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/denovoPC_TSS_Oct2021_each_acc_TEinserts.reverse.TEtypes.bed")



rownames(pcTSS.27genomes.TEsizes.forward)<-pcTSS.27genomes.TEsizes.forward$gene
rownames(pcTSS.27genomes.TEsizes.reverse)<-pcTSS.27genomes.TEsizes.reverse$gene
rownames(pcTSS.27genomes.TEtypes.forward)<-pcTSS.27genomes.TEtypes.forward$gene
rownames(pcTSS.27genomes.TEtypes.reverse)<-pcTSS.27genomes.TEtypes.reverse$gene



pcTSS.27genomes.TEsizes.forward[pcTSS.27genomes.TEsizes.forward==""]<-0
pcTSS.27genomes.TEsizes.reverse[pcTSS.27genomes.TEsizes.reverse==""]<-0
pcTSS.27genomes.TEtypes.forward[pcTSS.27genomes.TEtypes.forward==""]<-"no_TE"
pcTSS.27genomes.TEtypes.reverse[pcTSS.27genomes.TEtypes.reverse==""]<-"no_TE"


################################## 
#upload distance to the centromere in TAIR10 and 27 genomes
LINC.distance_from_centromere <- read.delim("03_Projects/2018_lncRNA_variation_paper/03_genetic_variation/blast_presence_absence_TE_insertions/LINC.distance_from_centromere.TAIR10_27genomes.processed.bed")

