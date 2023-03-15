#Creating cumulative transcriptome annotation 



#de novo transcriptome annotation
################################
#### only change this: 
export name_assembly=20211013_annotation

#lists of assemblies for the first step of merging
pollen=$folder_identification_pipeline/assemblies_merge_pollen_eracaps_may.txt
flower=$folder_identification_pipeline/assemblies_merge_flower_eracaps_may.txt
seedling_eracaps=$folder_identification_pipeline/assemblies_merge_seedling_eracaps_may.txt
seedling_cortijo=$folder_identification_pipeline/assemblies_merge_seedling_cortijo_may.txt
rosette_eracaps=$folder_identification_pipeline/assemblies_merge_rosette_eracaps_may.txt
rosette_1001g=$folder_identification_pipeline/assemblies_merge_rosette_1001G_may.txt
rosette_1001gnew=$folder_identification_pipeline/assemblies_merge_rosette_1001Gnew_may.txt


#merging 7 different assemblies
################################

#discard ChrM and ChrC 

#the folder where the output annotation will be :
export  working_folder=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/$name_assembly
mkdir $working_folder
cd $working_folder
mkdir /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/01_Identification_pipeline/$name_assembly
export R_folder=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/short_Rscripts
genome=/groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/TAIR10_all.fa
chr_sizes=/groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/chr_length_wo_ChrM_ChrC.txt

export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation

ml stringtie/2.1.5-gcccore-7.3.0
ml cufflinks/2.2.1-foss-2018b
ml python/2.7.15-foss-2018b
ml bedtools/2.27.1-foss-2018b
module load r/3.5.1-foss-2018b


#merge with cuffmerge (without reference!) - step 1 - separate tissues/datasets
samples=pollen
cuffmerge --keep-tmp -s $genome --min-isoform-fraction 0 -p 8 -o $working_folder/$samples ${!samples}

samples=flower
cuffmerge --keep-tmp -s $genome --min-isoform-fraction 0 -p 8 -o $working_folder/$samples ${!samples}

samples=seedling_eracaps
cuffmerge --keep-tmp -s $genome --min-isoform-fraction 0 -p 8 -o $working_folder/$samples ${!samples}

samples=seedling_cortijo
cuffmerge --keep-tmp -s $genome --min-isoform-fraction 0 -p 8 -o $working_folder/$samples ${!samples}

samples=rosette_eracaps
cuffmerge --keep-tmp -s $genome --min-isoform-fraction 0 -p 8 -o $working_folder/$samples ${!samples}

samples=rosette_1001g
cuffmerge --keep-tmp -s $genome --min-isoform-fraction 0 -p 8 -o $working_folder/$samples ${!samples}

samples=rosette_1001gnew
cuffmerge --keep-tmp -s $genome --min-isoform-fraction 0 -p 8 -o $working_folder/$samples ${!samples}


samples=pollen
cat  $working_folder/$samples/transcripts.gtf |  awk -v OFS="\t" '!($7=="."){print $0}' | grep -v chloroplast | grep -v  mitochondria > $working_folder/$samples/transcripts.nodots.gtf

samples=flower
cat  $working_folder/$samples/transcripts.gtf |  awk -v OFS="\t" '!($7=="."){print $0}' | grep -v chloroplast | grep -v  mitochondria> $working_folder/$samples/transcripts.nodots.gtf

samples=rosette_1001g
cat  $working_folder/$samples/transcripts.gtf |  awk -v OFS="\t" '!($7=="."){print $0}' | grep -v chloroplast | grep -v  mitochondria> $working_folder/$samples/transcripts.nodots.gtf

samples=rosette_1001gnew
cat  $working_folder/$samples/transcripts.gtf |  awk -v OFS="\t" '!($7=="."){print $0}'| grep -v chloroplast | grep -v  mitochondria > $working_folder/$samples/transcripts.nodots.gtf

samples=rosette_eracaps
cat  $working_folder/$samples/transcripts.gtf |  awk -v OFS="\t" '!($7=="."){print $0}'| grep -v chloroplast | grep -v  mitochondria > $working_folder/$samples/transcripts.nodots.gtf

samples=seedling_cortijo
cat  $working_folder/$samples/transcripts.gtf |  awk -v OFS="\t" '!($7=="."){print $0}' | grep -v chloroplast | grep -v  mitochondria> $working_folder/$samples/transcripts.nodots.gtf
 
samples=seedling_eracaps
cat  $working_folder/$samples/transcripts.gtf |  awk -v OFS="\t" '!($7=="."){print $0}' | grep -v chloroplast | grep -v  mitochondria> $working_folder/$samples/transcripts.nodots.gtf


#create list of single-tissue assemblies for the 2nd step of Cuffmerge
cd $working_folder
ls $working_folder/*/transcripts.nodots.gtf | grep -v pools> $working_folder/CFM_merged_assemblies_list.txt
ls $working_folder/*/transcripts.nodots.gtf > $working_folder/CFM_merged_assemblies_list_with_pools.txt


#######################################################
#####second step of merging      ######################
#######################################################
#two step merging allows avoiding (many, but not all) chimeric trascripts 

cuffmerge --keep-tmp -s $genome --min-isoform-fraction 0 -p 8 -o $working_folder $working_folder/CFM_merged_assemblies_list.txt

#convert to bed 
file=transcripts
cat  $working_folder/$file.gtf |  awk -v OFS="\t" '!($7=="."){print $0}' > $working_folder/$file.nodots.gtf
/groups/nordborg/projects/cegs/alexandra/software/gtfToGenePred $working_folder/$file.nodots.gtf $working_folder/$file.genepred
/groups/nordborg/projects/cegs/alexandra/software/genePredToBed $working_folder/$file.genepred $working_folder/$file.1.bed
cat $working_folder/$file.1.bed | awk -v OFS="\t"  '{$7=$2;print $0}'| sed -e 's/chloroplast/ChrC/g'| sed -e 's/mitochondria/ChrM/g'  > $working_folder/$file.bed  
cat $working_folder/$file.bed  | grep Chr4 > /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/01_Identification_pipeline/$name_assembly/CFM_2step_CFM.chr4.bed 


echo "transcripts.bed transcript number">  $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.bed >>  $working_folder/linenumbers.txt
echo "transcripts.bed non-overlapping loci number">>  $working_folder/linenumbers.txt
cat $working_folder/transcripts.bed  | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l  >> $working_folder/linenumbers.txt
  
 
wc -l  $working_folder/transcripts.bed
  
#115124 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.bed

#loci  
 
sortBed -i $working_folder/transcripts.bed | mergeBed -s -i stdin| wc -l 
# 41635

  
# make gene positions annotation 
cat $working_folder/genes.fpkm_tracking |grep -v tracking_id |  sed -e 's/chloroplast/ChrC/g'|  sed -e 's/mitochondria/ChrM/g' | awk -v OFS="\t"  '{split($7,s,":"); split(s[2],a,"-"); print $1, s[1], a[1],a[2]}' >  $working_folder/genes.txt 
cat $working_folder/transcripts.bed | awk -v OFS="\t"  '{split($4,s,"."); print s[1]"."s[2],$6}' |uniq > $working_folder/genes_strand.txt 


#add strand 
export genes_strand_location=$working_folder/genes_strand.txt 
export genes_location=$working_folder/genes.txt


module load r/3.5.1-foss-2018b
Rscript $R_folder/filtering_Rscript_1.r 
################

cat $working_folder/genes.txt | awk -v OFS="\t"  '{print $3,$4,$5,$1,0,$2}'| sort -k1,1 -k2,2g > $working_folder/genes.bed


##########################################################################
####   1. artefact filtering 							##################
##########################################################################
 
### 1. remove transcripts with mRNA length <200nt, remove single-exon transcripts <400nt 
cat $working_folder/transcripts.bed | \
awk '{split($11,s,","); total=0; for(x in s){ total+=s[x];} print total,$0}'  | awk -v OFS="\t" '$1>=200 {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13 }' \
  | awk -v OFS="\t" '!(($3-$2)<=400 && $10==1) {print $0 }' > $working_folder/transcripts.f1.bed 

echo "number of transcripts after removing transcripts <200nt and single-exon transcripts <400nt " >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.f1.bed >> $working_folder/linenumbers.txt
echo "number of loci" >> $working_folder/linenumbers.txt
cat $working_folder/transcripts.f1.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l  >> $working_folder/linenumbers.txt

##########################################################################
### 2. Additional antisense filter 						##################
##########################################################################


# a. Antisense leakage single exon
#remove transcripts that 
#			1.are single exon , 
#			2. overlap another transcript's exon >98% of their length (fully included), 
#			3. that overlap multiexon genes (because singleexon genes often have a fully overlapping AS RNA), 
#			4. and  that are short - <0.6kb
#super strict - don't want artifacts

intersectBed -split -f 0.98 -S -wo -a $working_folder/transcripts.f1.bed -b $working_folder/transcripts.f1.bed | awk '($10==1 && $22>1 &&($3-$2)<600){print $4 }'| uniq > $working_folder/potential_AS_artifacts.txt

#############
# b.  spliced leakage mirror transcripts (pieces of transcripts - usually 2 exons)
# exclude transcripts that have very similar introns on the antisense strand to PCgenes 
#create annotation of introns 
ml bedparse/0.2.3-foss-2018b-python-3.6.6
bedparse introns $working_folder/transcripts.f1.bed > $working_folder/transcripts.f1.introns.bed
bed12ToBed6 -i $working_folder/transcripts.f1.introns.bed >$working_folder/transcripts.f1.introns.bed6

bedparse introns $annotationfolder/Araport11_protein_coding.201606.bed > $annotationfolder/Araport11_protein_coding.201606.introns.bed
bed12ToBed6 -i $annotationfolder/Araport11_protein_coding.201606.introns.bed > $annotationfolder/Araport11_protein_coding.201606.introns.bed6


#intersectBed -u -r -f 0.95 -S -a $annotationfolder/Araport11_protein_coding.201606.introns.bed6 -b $annotationfolder/Araport11_protein_coding.201606.introns.bed6   | awk -v OFS="\t" '{print $4}' | uniq > $annotationfolder/Araport11_protein_coding.201606.genes_with_mirror_introns_AS_genes.names.txt

cat $annotationfolder/Araport11_protein_coding.201606.introns.bed6  | grep -v -w -f $annotationfolder/Araport11_protein_coding.201606.genes_with_mirror_introns_AS_genes.names.txt | intersectBed -u -r -f 0.95 -S -a $working_folder/transcripts.f1.introns.bed6 -b stdin  | awk -v OFS="\t" '{print $4}' | uniq > $working_folder/mirror_transcripts_spliced_artefacts_names.txt

# combine all transcript names to exclude
cat $working_folder/mirror_transcripts_spliced_artefacts_names.txt $working_folder/potential_AS_artifacts.txt |  uniq > $working_folder/potential_artifacts.txt
 
 export transcripts_tofiltert_loc=$working_folder/potential_artifacts.txt
export transcripts_loc=$working_folder/transcripts.f1.bed  
export transcripts_filtered=$working_folder/transcripts.f1A.bed  

# R script to remove artefacts 
Rscript $R_folder/filtering_Rscript_1A.r  

#write down number of transcripts/loci
echo "number of transcripts after (stringently) removing potential AS and chimeric atrifacts transcripts " >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.f1A.bed >> $working_folder/linenumbers.txt
echo "number of loci after (stringently) removing potential AS and chimeric atrifacts transcripts" >> $working_folder/linenumbers.txt
cat $working_folder/transcripts.f1A.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l  >> $working_folder/linenumbers.txt

wc -l $working_folder/transcripts.f1A.bed
#110947 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.f1A.bed


##########################################################################
# 3. splitting annotation into PC and non-coding part 			##########
##########################################################################

#use both Araport and TAIR10 for PC annotation 
export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation
#pc from TAIR10 
export PC_tair=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.protein_coding.gtf
#pc annotation from Araport 
export PC_araport=$annotationfolder/Araport11_protein_coding.201606.bed
#take pseudogene annotation from araport 
export pseudogenes=$annotationfolder/Araport11_pseudogene.201606.bed
#take transposon gene annotation from araport 
export TE_genes=$annotationfolder/Araport11_transposable_element_gene.201606.bed
#TE annotation 
export TE_fragments=$annotationfolder/Araport11_TEs.transposon_fragments.bed 
export TE_elements=$annotationfolder/Araport11_TEs.transposable_elements.bed
#araport nc RNAs 
export AraportNC=$annotationfolder/Araport11_non_coding.2016016.sorted.bed
export AraportNTR=$annotationfolder/Araport11_novel_transcribed_region.201606.bed

#TAIR10
protein_coding=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.protein_coding.gtf
lncRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.lncRNA.gtf
atlncRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.atlncRNA.gtf
miRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.miRNA.gtf
pre_miRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.pre_miRNA.gtf
tRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.tRNA.gtf
atRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.atRNA.gtf
snoRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.snoRNA.gtf
otherRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.otherRNA.gtf
snRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.snRNA.gtf
SRP_RNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.SRP_RNA.gtf
nontranslating_CDS=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.nontranslating_CDS.gtf
rRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.rRNA.gtf
ncRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.ncRNA.gtf
RNase_MRP_RNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.RNase_MRP_RNA.gtf
sense_intronic=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.sense_intronic.gtf
antisense_RNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.antisense_RNA.gtf

###############################
#make denovo PC annotation: 
###############################

#(there will be many chimeric loci!! so do not take this PC annotation too seriously, use Araport)
f1=$working_folder/transcripts.f1A.bed
intersectBed -u -split -s -a $f1 -b $PC_tair $PC_araport | intersectBed -v -split -s -a stdin -b $tRNA $rRNA $snRNA $snoRNA|sortBed -i stdin > $working_folder/transcripts.f1.denovoPC.with_chimeras.bed

#will need this fuller denovoPC annotation later for excluding PC gene extensions exons 
 
wc -l $working_folder/transcripts.f1.denovoPC.with_chimeras.bed
# 85123 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20210513_annotation/transcripts.f1.denovoPC.bed


# Remove chimeric transcript assemblies merging 2 or more genes
# remove chimeric PC-PC PC-pseudogene and pseudo-pseudo transcripts

export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation

#sortBed -i $annotationfolder/Araport11_pseudogene.201606.bed | mergeBed -s -c 4,5,6 -o first -i stdin > $annotationfolder/Araport11_pseudogene.201606.mergebed.loci.bed

#intersectBed -split -s -v -a $annotationfolder/Araport11_pseudogene.201606.mergebed.loci.bed -b /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.mergeBed.loci.bed |awk -v OFS="\t" '{split($4,a,".");print $1,$2,$3,a[1],$5,$6}'  >  $annotationfolder/Araport11_pseudogene.201606.mergebed.loci.nooverlap_with_Ar11PC.bed

intersectBed -u -s -split -a $annotationfolder/Araport11_pseudogene.201606.mergebed.loci.nooverlap_with_Ar11PC.bed -b /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes_nooverlapping_genes.bed | wc -l
#0
cat $annotationfolder/Araport11_pseudogene.201606.mergebed.loci.nooverlap_with_Ar11PC.bed  /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes_nooverlapping_genes.bed | sortBed -i stdin > pseudo_and_PC_for_chimeras_filt.bed

intersectBed -s -f 0.3 -e -F 0.3 -wo -a $working_folder/transcripts.f1.denovoPC.with_chimeras.bed  -b pseudo_and_PC_for_chimeras_filt.bed | awk -v OFS="\t" '{print $4,$16}' | uniq|awk -v OFS="\t" '{print $1}'| uniq -d >  $working_folder/chimeric_transcripts_artefacts_names.txt

export transcripts_tofiltert_loc=$working_folder/chimeric_transcripts_artefacts_names.txt
export transcripts_loc=$working_folder/transcripts.f1.denovoPC.with_chimeras.bed 
export transcripts_filtered=$working_folder/transcripts.f1.denovoPC.without_chimeras.bed 

# R script to remove artefacts 
Rscript $R_folder/filtering_Rscript_1A.r  

wc -l $working_folder/transcripts.f1.denovoPC.without_chimeras.bed

#77118 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.f1.denovoPC.without_chimeras.bed


#redefine PC loci 

head $working_folder/transcripts.f1.denovoPC.without_chimeras.bed
Chr1    3630    5899    CUFF.2.1        0       +       3630    5899    0       6       283,281,120,390,153,461,        0,365,855,1075,1543,1808,
Chr1    3635    5899    CUFF.2.2        0       +       3635    5899    0       5       278,610,390,153,461,    0,360,1070,1538,1803,
Chr1    3642    5899    CUFF.2.3        0       +       3642    5899    0       6       271,276,120,390,153,461,        0,358,843,1063,1531,1796,
Chr1    3643    5899    CUFF.2.4        0       +       3643    5899    0       6       270,281,100,390,153,461,        0,352,862,1062,1530,1795,
Chr1    3691    5899    CUFF.2.5        0       +       3691    5899    0       5       222,281,390,153,461,    0,304,1014,1482,1747,


#PC
#Cuffmerge loci names in the $working_folder/transcripts.f1.denovoPC.bed
cat $working_folder/transcripts.f1.denovoPC.without_chimeras.bed | awk -v OFS="\t"  '{print $4,$5}' | awk -v OFS="\t"  '{split($1,s,".");   print s[1]"."s[2]}' |sort | uniq  | wc -l
# 23098


cat $working_folder/transcripts.f1.denovoPC.without_chimeras.bed | awk -v OFS="\t"  '{print $4,$5}' | awk -v OFS="\t"  '{split($1,s,".");   print s[1]"."s[2]}' |sort -V -k1,1 -k2,2 | uniq  > $working_folder/transcripts.f1.denovoPC.cuff_locinames.bed
wc -l $working_folder/transcripts.f1.denovoPC.cuff_locinames.bed
#23098 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.f1.denovoPC.cuff_locinames.bed


rm $working_folder/loci_table
touch 
count=1
while read line
do
locus_name=$line
cat  $working_folder/transcripts.f1.denovoPC.without_chimeras.bed  |awk -v OFS="\t"  '{split($4,s,".");   print $0,s[1]"."s[2]}' | grep -w $line | awk -v OFS="\t"  '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > $working_folder/tmp_1
mergeBed  -c 4,5,6 -o distinct,first,first -s -i $working_folder/tmp_1 | awk -v locus_number="$count" '{if (NR==1) {print $0,"CUFF_PC."locus_number} else {print $0,"CUFF_PC."locus_number"-"NR}}'>> $working_folder/loci_table
mergeBed  -c 4,5,6 -o distinct,first,first -s -i $working_folder/tmp_1 | awk -v locus_number="$count" '{if (NR==1) {print $0,"CUFF_PC."locus_number} else {print $0,"CUFF_PC."locus_number"-"NR}}'
count=$((count+1))
done < $working_folder/transcripts.f1.denovoPC.cuff_locinames.bed

wc -l loci_table
#23676 loci_table



cat $working_folder/loci_table > $working_folder/transcripts.f1.denovoPC.loci_old-new_names.bed

cat $working_folder/loci_table |awk -v OFS="\t"  '{ print $1,$2,$3,$7,$5,$6}'> $working_folder/denovoPC.loci.bed 
cat $working_folder/denovoPC.loci.bed  > inter 
sortBed -i inter > $working_folder/denovoPC.loci.bed 
 head $working_folder/denovoPC.loci.bed
Chr1    3630    5899    CUFF_PC.2       0       +
Chr1    6787    9130    CUFF_PC.8       0       -
Chr1    9203    13714   CUFF_PC.1       0       -
Chr1    23120   31252   CUFF_PC.22      0       +
Chr1    31169   33171   CUFF_PC.23      0       -
Chr1    33378   37871   CUFF_PC.9       0       -
Chr1    38432   41017   CUFF_PC.4       0       -
Chr1    44969   47059   CUFF_PC.3       0       -
Chr1    47233   49392   CUFF_PC.5       0       -


cat $working_folder/loci_table  |  awk  -v OFS="\t"  'split($4,s,","){ for(i in s) {print s[i],$7,$7"."i  }}' | sort -V -k1,1 -k2,2 |head
CUFF.1.1        CUFF_PC.1       CUFF_PC.1.1
CUFF.2.1        CUFF_PC.2       CUFF_PC.2.1
CUFF.2.2        CUFF_PC.2       CUFF_PC.2.2
CUFF.2.3        CUFF_PC.2       CUFF_PC.2.3
CUFF.2.4        CUFF_PC.2       CUFF_PC.2.4
CUFF.2.5        CUFF_PC.2       CUFF_PC.2.5


cat $working_folder/loci_table  |  awk  -v OFS="\t"  'split($4,s,","){ for(i in s) {print s[i],$7,$7"."i  }}' > $working_folder/denovoPC.loci_transcripts_names_old_new.bed 
wc -l  $working_folder/denovoPC.loci_transcripts_names_old_new.bed
#77118 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/denovoPC.loci_transcripts_names_old_new.bed



# make annotation with new names 

sort -k 4b,4  $working_folder/transcripts.f1.denovoPC.without_chimeras.bed  > $working_folder/old_names_annotation_sorted
sort -k 1b,1 $working_folder/denovoPC.loci_transcripts_names_old_new.bed  >  $working_folder/new_old_names_sorted
join -1 1 -2 4 $working_folder/new_old_names_sorted $working_folder/old_names_annotation_sorted |  awk -v OFS="\t" '{print $2,$3,$4,$5,$6,$1,$7,$8,$9,$10,$11,$12,$13,$14}' | sort -V -k1,1 -k2,2 > $working_folder/denovoPC.transcripts.old_names_new_names.bed

 head $working_folder/denovoPC.transcripts.old_names_new_names.bed
CUFF_PC.1       CUFF_PC.1.1     Chr1    9203    13714   CUFF.1.1        0       -       9203    13714   0       3       3151,750,380,   0,3220,4131,
CUFF_PC.2       CUFF_PC.2.1     Chr1    3630    5899    CUFF.2.1        0       +       3630    5899    0       6       283,281,120,390,153,461,        0,365,855,1075,1543,1808,
CUFF_PC.2       CUFF_PC.2.2     Chr1    3635    5899    CUFF.2.2        0       +       3635    5899    0       5       278,610,390,153,461,    0,360,1070,1538,1803,
CUFF_PC.2       CUFF_PC.2.3     Chr1    3642    5899    CUFF.2.3        0       +       3642    5899    0       6       271,276,120,390,153,461,        0,358,843,1063,1531,1796,
CUFF_PC.2       CUFF_PC.2.4     Chr1    3643    5899    CUFF.2.4        0       +       3643    5899    0       6       270,281,100,390,153,461,        0,352,862,1062,1530,1795,


cat $working_folder/denovoPC.transcripts.old_names_new_names.bed |  awk -v OFS="\t" '{print $3,$4,$5,$2,$7,$8,$9,$10,$11,$12,$13,$14}' | sortBed -i stdin > $working_folder/denovoPC.transcripts.bed

head  $working_folder/denovoPC.transcripts.bed                               
Chr1    3630    5899    CUFF_PC.2.1     0       +       3630    5899    0       6       283,281,120,390,153,461,        0,365,855,1075,1543,1808,
Chr1    3635    5899    CUFF_PC.2.2     0       +       3635    5899    0       5       278,610,390,153,461,    0,360,1070,1538,1803,
Chr1    3642    5899    CUFF_PC.2.3     0       +       3642    5899    0       6       271,276,120,390,153,461,        0,358,843,1063,1531,1796,
Chr1    3643    5899    CUFF_PC.2.4     0       +       3643    5899    0       6       270,281,100,390,153,461,        0,352,862,1062,1530,1795,
Chr1    3691    5899    CUFF_PC.2.5     0       +       3691    5899    0       5       222,281,390,153,461,    0,304,1014,1482,1747,
Chr1    6787    8710    CUFF_PC.8.1     0       -       6787    8710    0       8       282,76,67,86,74,355,48,117,     0,369,596,776,974,1154,1629,1806,
Chr1    6787    8737    CUFF_PC.8.2     0       -       6787    8737    0       9       282,76,67,86,74,46,36,48,167,   0,369,596,776,974,1154,1502,1629,1783,

 wc -l $working_folder/denovoPC.transcripts.bed
#77118 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/denovoPC.transcripts.bed

echo "number of denovo protein-coding transcripts" >> $working_folder/linenumbers.txt
wc -l $working_folder/denovoPC.transcripts.bed >> $working_folder/linenumbers.txt

echo "number of denovo protein-coding loci" >> $working_folder/linenumbers.txt
wc -l $working_folder/denovoPC.loci.bed  >> $working_folder/linenumbers.txt

 wc -l $working_folder/denovoPC.loci.bed
#23676 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/denovoPC.loci.bed

denovoPC=$working_folder/denovoPC.transcripts.bed


###############################
# Define pseudogenes 
###############################

#will not allow PC and pseudogenes to share exons in my annotation although there are 45 transcripts in araport pseudogene annotation that have a sense overlap with PC genes 
 #intersectBed -s  -u -a $pseudogenes -b  $PC_araport $PC_tair | wc -l                               
# 46
#do not allow overlap with trna,snRNA snoRNA and rRNA 

intersectBed -u -split -s -a $f1 -b $pseudogenes | intersectBed -v  -s -a stdin -b $denovoPC $PC_araport $PC_tair | intersectBed -v -split -s -a stdin -b $tRNA $rRNA $snRNA $snoRNA| sortBed -i stdin  > $working_folder/transcripts.f1.denovo_pseudogene.bed
wc -l $working_folder/transcripts.f1.denovo_pseudogene.bed
#1412 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.f1.denovo_pseudogene.bed


echo "number of denovo pseudogene transcripts" >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.f1.denovo_pseudogene.bed >> $working_folder/linenumbers.txt

echo "number of denovo pseudogene loci" >> $working_folder/linenumbers.txt
cat $working_folder/transcripts.f1.denovo_pseudogene.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l >> $working_folder/linenumbers.txt

intersectBed -u -s -a $f1 -b $PC_tair $PC_araport | intersectBed -v -s -split -a stdin -b $PC_tair $PC_araport|intersectBed -v -split -s -a stdin -b $tRNA $rRNA $snRNA $snoRNA | sortBed -i stdin > $working_folder/transcripts.f1.senseoverlapping.bed
wc -l $working_folder/transcripts.f1.senseoverlapping.bed
#187 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.f1.senseoverlapping.bed

echo "number of denovo senseoverlapping transcripts" >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.f1.senseoverlapping.bed >> $working_folder/linenumbers.txt

echo "number of denovo senseoverlapping loci" >> $working_folder/linenumbers.txt
cat $working_folder/transcripts.f1.senseoverlapping.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l >> $working_folder/linenumbers.txt


######################
#TEs				##
######################

#####################
#a. TE genes
#####################

#define TE genes - transcripts that have an exonic overlap with an Araport 11 TE gene
intersectBed -v -s -a $f1 -b $PC_tair $PC_araport $pseudogenes|   intersectBed -v -split -s -a stdin -b  $denovoPC_with_chimeras|intersectBed -v -split -s -a stdin -b $tRNA $rRNA $snRNA $snoRNA |   intersectBed -u -s  -split -a stdin -b  $TE_genes | sortBed -i stdin >  $working_folder/transcripts.TE_genes.bed
wc -l  $working_folder/transcripts.TE_genes.bed
#3777 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.TE_genes.bed


#####################
#b. expressed TE fragments 
#####################
#Define TE fragment transcripts (arbitrary separation from lncRNAs - threshold 60% exonic sense overlap)

intersectBed -v -s -a $f1 -b $PC_tair $PC_araport $pseudogenes|   intersectBed -v -split -s -a stdin -b  $denovoPC_with_chimeras| intersectBed -v -split -s -a stdin -b $tRNA $rRNA $snRNA $snoRNA | intersectBed -v -s -split  -a stdin -b  $TE_genes | sortBed -i stdin  >  $working_folder/transcripts.TEs.prelim.bed
#19087 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.TEs.prelim.bed


# definition of lncRNAs : 
# join -1 4 -2 1 $working_folder/tab3_sorted $working_folder/tab4_sorted |  awk -v OFS="\t" '($13<0.6 &&($13<0.5 || $15<0.8)){print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12}'|sortBed -i stdin > $working_folder/transcripts.denovoNC.f2.bed

bed12ToBed6 -i $working_folder/transcripts.TEs.prelim.bed > $working_folder/transcripts.TEs.prelim.bed6

coverageBed -split -s -a  $working_folder/transcripts.TEs.prelim.bed6 -b $TE_fragments | awk -v OFS="\t" '{arr[$4]+=$8} END {for (i in arr) {print i,arr[i]}}' | sort  -k1,1n > $working_folder/exon_length_covered

coverageBed -split -s -a  $working_folder/transcripts.TEs.prelim.bed6 -b $TE_fragments | awk -v OFS="\t" '{arr[$4]+=$9} END {for (i in arr) {print i,arr[i]}}' | sort -k1,1n > $working_folder/exon_length

sort -k 1b,1  $working_folder/exon_length_covered  > $working_folder/exon_length_covered_sorted
sort -k 1b,1 $working_folder/exon_length  >  $working_folder/exon_length_sorted
join -1 1 -2 1 $working_folder/exon_length_covered_sorted $working_folder/exon_length_sorted |  awk -v OFS="\t" '{print $1, $2/$3}'  > $working_folder/transcripts.TEs.prelim.exon_coverage_by.TE_frags.bed

wc -l  $working_folder/transcripts.TEs.prelim.exon_coverage_by.TE_frags.bed
#19087 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.TEs.prelim.exon_coverage_by.TE_frags.bed

coverageBed -s -a  $working_folder/transcripts.TEs.prelim.bed -b $TE_fragments | awk -v OFS="\t" '{print $4, $15,$16}'>  $working_folder/transcripts.TEs.prelim.locus_coverage_by.TE_frags.bed

sort -k 4b,4  $working_folder/transcripts.TEs.prelim.bed   > $working_folder/tab1_sorted
sort -k 1b,1 $working_folder/transcripts.TEs.prelim.exon_coverage_by.TE_frags.bed   >  $working_folder/tab2_sorted
join -1 4 -2 1 $working_folder/tab1_sorted $working_folder/tab2_sorted |  awk -v OFS="\t" '{print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13}'  | sort -k 4b,4 > $working_folder/tab3_sorted
sort -k 1b,1 $working_folder/transcripts.TEs.prelim.locus_coverage_by.TE_frags.bed > $working_folder/tab4_sorted

join -1 4 -2 1 $working_folder/tab3_sorted $working_folder/tab4_sorted |  awk -v OFS="\t" '($13>=0.6 ||($13>=0.5 && $15>=0.8)){print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12}'|sortBed -i stdin > $working_folder/transcripts.TEs.bed

wc -l  $working_folder/transcripts.TEs.bed 

wc -l  $working_folder/transcripts.TEs.bed
#882 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.TEs.bed


echo "number of expressed TEs genes' transcripts expressed" >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.TE_genes.bed >> $working_folder/linenumbers.txt

echo "number of expressed TEs genes' loci" >> $working_folder/linenumbers.txt
cat $working_folder/transcripts.TE_genes.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l >> $working_folder/linenumbers.txt
#2094 loci

echo "number of expressed TEs transcripts expressed" >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.TEs.bed >> $working_folder/linenumbers.txt

echo "number of expressed TEs loci" >> $working_folder/linenumbers.txt
cat $working_folder/transcripts.TEs.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l >> $working_folder/linenumbers.txt
#692 loci


# how many TE frag transcripts overlap TE gene transcripts? 
intersectBed -s -u -a $working_folder/transcripts.TEs.bed -b $working_folder/transcripts.TE_genes.bed | wc -l 
#98 transcripts (11%)
intersectBed -s -u -a $working_folder/transcripts.TEs.bed -b $working_folder/transcripts.TE_genes.bed |  sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l
# 60 loci (9%)
intersectBed -split -s -u -a $working_folder/transcripts.TEs.bed -b $working_folder/transcripts.TE_genes.bed | wc -l 
#75 (11%)

#remove those ambiguous TE transcripts to avoid confusions 

intersectBed -s -v -a $working_folder/transcripts.TEs.bed -b $working_folder/transcripts.TE_genes.bed > $working_folder/transcripts.TE_frags.bed
wc -l $working_folder/transcripts.TE_frags.bed
#784 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.TE_frags.bed




#################################
# 4. initial lncRNA file: 
#################################
# exclude tr-ts that have any sense overlap with 
# 1. PC genes (TAIR and Araport11) 
# 2. exons of potential extensions of PC genes (de novo assembled PC genes (with chimeras!)) 
# 3. pseudogenes (Araport11) 
# 4. TE genes from  Araport11
# 5. no more than 60% sense exonic overlap with TE fragments (and no more than 80% gene body overlap with TE frags )
denovoPC_with_chimeras=$working_folder/transcripts.f1.denovoPC.with_chimeras.bed
f1=$working_folder/transcripts.f1A.bed
export TE_genes=$annotationfolder/Araport11_transposable_element_gene.201606.bed
export TE_fragments=$annotationfolder/Araport11_TEs.transposon_fragments.bed 
export TE_elements=$annotationfolder/Araport11_TEs.transposable_elements.bed


intersectBed -v -s -a $f1 -b $PC_tair  $PC_araport  $pseudogenes |  intersectBed -v -split -s -a stdin -b  $denovoPC_with_chimeras|  intersectBed -v -s -a stdin -b $TE_genes  > $working_folder/transcripts.denovoNC.f1b.bed
wc -l $working_folder/transcripts.denovoNC.f1b.bed
# 19084 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.denovoNC.f1b.bed
cat $working_folder/transcripts.denovoNC.f1b.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l 
#12562

# check exon overlap with TE fragments 
 
bed12ToBed6 -i $working_folder/transcripts.denovoNC.f1b.bed > $working_folder/transcripts.denovoNC.f1b.bed6

coverageBed -split -s -a  $working_folder/transcripts.denovoNC.f1b.bed6 -b $TE_fragments | awk -v OFS="\t" '{arr[$4]+=$8} END {for (i in arr) {print i,arr[i]}}' | sort -k1> $working_folder/exon_length_covered

coverageBed -split -s -a  $working_folder/transcripts.denovoNC.f1b.bed6 -b $TE_fragments | awk -v OFS="\t" '{arr[$4]+=$9} END {for (i in arr) {print i,arr[i]}}' | sort -k1> $working_folder/exon_length

sort -k 1b,1  $working_folder/exon_length_covered  > $working_folder/exon_length_covered_sorted
sort -k 1b,1 $working_folder/exon_length  >  $working_folder/exon_length_sorted
join -1 1 -2 1 $working_folder/exon_length_covered_sorted $working_folder/exon_length_sorted |  awk -v OFS="\t" '{print $1, $2/$3}'  > $working_folder/transcripts.denovoNC.f1b.exon_coverage_by.TE_frags.bed
wc -l $working_folder/transcripts.denovoNC.f1b.exon_coverage_by.TE_frags.bed
#19084 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.denovoNC.f1b.exon_coverage_by.TE_frags.bed

coverageBed -s -a  $working_folder/transcripts.denovoNC.f1b.bed -b $TE_fragments | awk -v OFS="\t" '{print $4, $15,$16}'>  $working_folder/transcripts.denovoNC.f1b.locus_coverage_by.TE_frags.bed


sort -k 4b,4  $working_folder/transcripts.denovoNC.f1b.bed  > $working_folder/tab1_sorted
sort -k 1b,1 $working_folder/transcripts.denovoNC.f1b.exon_coverage_by.TE_frags.bed   >  $working_folder/tab2_sorted
join -1 4 -2 1 $working_folder/tab1_sorted $working_folder/tab2_sorted |  awk -v OFS="\t" '{print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13}'  | sort -k 4b,4 > $working_folder/tab3_sorted
sort -k 1b,1 $working_folder/transcripts.denovoNC.f1b.locus_coverage_by.TE_frags.bed > $working_folder/tab4_sorted

join -1 4 -2 1 $working_folder/tab3_sorted $working_folder/tab4_sorted |  awk -v OFS="\t" '($13<0.6 &&($13<0.5 || $15<0.8)){print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12}'|sortBed -i stdin > $working_folder/transcripts.denovoNC.f2.bed

wc -l $working_folder/transcripts.denovoNC.f2.bed
#18258 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.denovoNC.f2.bed

cat $working_folder/transcripts.denovoNC.f2.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l
#11928 loci


echo "number of preliminary lncRNA transcripts" >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.denovoNC.f2.bed  >> $working_folder/linenumbers.txt
echo "number of preliminary lncRNA loci" >> $working_folder/linenumbers.txt
cat $working_folder/transcripts.denovoNC.f2.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l >> $working_folder/linenumbers.txt


#####
#5. remove transcripts with protein coding potential 
#####

#run CPC 
#make fasta file 
bedtools getfasta -name -s -split -fi /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/TAIR10_all.fa -bed $working_folder/transcripts.denovoNC.f2.bed -fo $working_folder/transcripts.denovoNC.f2.fa
bedtools getfasta -name -s -split -fi /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/TAIR10_all.fa -bed $working_folder/denovoPC.transcripts.bed -fo $working_folder/denovoPC.transcripts.fa
bedtools getfasta -name -s -split -fi /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/TAIR10_all.fa -bed  $working_folder/transcripts.TE_frags.bed -fo  $working_folder/transcripts.TE_frags.fa
bedtools getfasta -name -s -split -fi /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/TAIR10_all.fa -bed  $working_folder/transcripts.TE_genes.bed -fo  $working_folder/transcripts.TE_genes.fa
bedtools getfasta -name -s -split -fi /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/TAIR10_all.fa -bed  $working_folder/transcripts.f1.denovo_pseudogene.bed -fo  $working_folder/transcripts.f1.denovo_pseudogene.fa

#araport nc and araport pc
#bedtools getfasta -name -s -split -fi /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/TAIR10_all.fa -bed $annotationfolder/Araport11_protein_coding.201606.bed -fo $annotationfolder/Araport11_protein_coding.201606.fa 

#bedtools getfasta -name -s -split -fi /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/TAIR10_all.fa -bed $annotationfolder/Araport11_non_coding.2016016.sorted.bed -fo $annotationfolder/Araport11_non_coding.2016016.fa 


ml biopython/1.72-foss-2018b-python-2.7.15


#!!!! CPC2 needs old Python! can be a conflict of versions - be careful!
#NC
/groups/nordborg/projects/cegs/alexandra/software/CPC2/CPC2-beta/bin/CPC2.py -i $working_folder/transcripts.denovoNC.f2.fa -o $working_folder/transcripts.denovoNC.f2.CPC_results.txt

cat $working_folder/transcripts.denovoNC.f2.CPC_results.txt |  awk -v OFS="\t" '($8=="noncoding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/transcripts.denovoNC.f2.CPC_results.noncoding.txt

cat $working_folder/transcripts.denovoNC.f2.CPC_results.txt |  awk -v OFS="\t" '($8=="coding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/transcripts.denovoNC.f2.CPC_results.coding.txt

#PC
/groups/nordborg/projects/cegs/alexandra/software/CPC2/CPC2-beta/bin/CPC2.py -i $working_folder/denovoPC.transcripts.fa -o $working_folder/denovoPC.transcripts.CPC_results.txt

cat $working_folder/denovoPC.transcripts.CPC_results.txt |  awk -v OFS="\t" '($8=="noncoding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/denovoPC.transcripts.CPC_results.noncoding.txt

cat $working_folder/denovoPC.transcripts.CPC_results.txt |  awk -v OFS="\t" '($8=="coding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/denovoPC.transcripts.CPC_results.coding.txt

#TE genes

/groups/nordborg/projects/cegs/alexandra/software/CPC2/CPC2-beta/bin/CPC2.py -i $working_folder/transcripts.TE_genes.fa -o $working_folder/transcripts.TE_genes.CPC_results.txt

cat $working_folder/transcripts.TE_genes.CPC_results.txt |  awk -v OFS="\t" '($8=="noncoding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/transcripts.TE_genes.CPC_results.noncoding.txt

cat $working_folder/transcripts.TE_genes.CPC_results.txt |  awk -v OFS="\t" '($8=="coding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/transcripts.TE_genes.CPC_results.coding.txt


#TE fragments 

/groups/nordborg/projects/cegs/alexandra/software/CPC2/CPC2-beta/bin/CPC2.py -i $working_folder/transcripts.TE_frags.fa -o $working_folder/transcripts.TE_frags.CPC_results.txt

cat $working_folder/transcripts.TE_frags.CPC_results.txt |  awk -v OFS="\t" '($8=="noncoding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/transcripts.TE_frags.CPC_results.noncoding.txt

cat $working_folder/transcripts.TE_frags.CPC_results.txt |  awk -v OFS="\t" '($8=="coding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/transcripts.TE_frags.CPC_results.coding.txt


#Pseudogenes

/groups/nordborg/projects/cegs/alexandra/software/CPC2/CPC2-beta/bin/CPC2.py -i $working_folder/transcripts.f1.denovo_pseudogene.fa -o $working_folder/transcripts.f1.denovo_pseudogene.CPC_results.txt

cat $working_folder/transcripts.f1.denovo_pseudogene.CPC_results.txt |  awk -v OFS="\t" '($8=="noncoding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/transcripts.f1.denovo_pseudogene.CPC_results.noncoding.txt

cat $working_folder/transcripts.f1.denovo_pseudogene.CPC_results.txt |  awk -v OFS="\t" '($8=="coding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/transcripts.f1.denovo_pseudogene.CPC_results.coding.txt


#Araport PC

/groups/nordborg/projects/cegs/alexandra/software/CPC2/CPC2-beta/bin/CPC2.py -i $annotationfolder/Araport11_protein_coding.201606.fa  -o $working_folder/Araport11_protein_coding.201606.CPC_results.txt

cat $working_folder/Araport11_protein_coding.201606.CPC_results.txt |  awk -v OFS="\t" '($8=="noncoding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/Araport11_protein_coding.201606.CPC_results.noncoding.txt

cat $working_folder/Araport11_protein_coding.201606.CPC_results.txt |  awk -v OFS="\t" '($8=="coding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/Araport11_protein_coding.201606.CPC_results.coding.txt



#Araport NC

/groups/nordborg/projects/cegs/alexandra/software/CPC2/CPC2-beta/bin/CPC2.py -i $annotationfolder/Araport11_non_coding.2016016.fa  -o $working_folder/Araport11_non_coding.2016016.CPC_results.txt

cat $working_folder/Araport11_non_coding.2016016.CPC_results.txt |  awk -v OFS="\t" '($8=="noncoding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/Araport11_non_coding.2016016.CPC_results.noncoding.txt

cat $working_folder/Araport11_non_coding.2016016.CPC_results.txt |  awk -v OFS="\t" '($8=="coding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/Araport11_non_coding.2016016.CPC_results.coding.txt




wc -l *CPC_results*
      109 Araport11_non_coding.2016016.CPC_results.coding.txt
    5470 Araport11_non_coding.2016016.CPC_results.noncoding.txt
    5580 Araport11_non_coding.2016016.CPC_results.txt
   44549 Araport11_protein_coding.201606.CPC_results.coding.txt
    3600 Araport11_protein_coding.201606.CPC_results.noncoding.txt
   48150 Araport11_protein_coding.201606.CPC_results.txt
   68362 denovoPC.transcripts.CPC_results.coding.txt
    8756 denovoPC.transcripts.CPC_results.noncoding.txt
   77119 denovoPC.transcripts.CPC_results.txt
    1247 transcripts.denovoNC.f2.CPC_results.coding.txt
   17011 transcripts.denovoNC.f2.CPC_results.noncoding.txt
   18259 transcripts.denovoNC.f2.CPC_results.txt
     559 transcripts.f1.denovo_pseudogene.CPC_results.coding.txt
     853 transcripts.f1.denovo_pseudogene.CPC_results.noncoding.txt
    1413 transcripts.f1.denovo_pseudogene.CPC_results.txt
      48 transcripts.TE_frags.CPC_results.coding.txt
     736 transcripts.TE_frags.CPC_results.noncoding.txt
     785 transcripts.TE_frags.CPC_results.txt
    2123 transcripts.TE_genes.CPC_results.coding.txt
    1654 transcripts.TE_genes.CPC_results.noncoding.txt
    3778 transcripts.TE_genes.CPC_results.txt
      37 transcripts.TEs.CPC_results.coding.txt
     590 transcripts.TEs.CPC_results.noncoding.txt
     628 transcripts.TEs.CPC_results.txt
  311416 total



# !!!!!!!!!!!!!!!!!filter not just the potentially coding transcript, but the whole locus ( all the transcripts overlapping the transcript with high PC potential)



export CPC_loc=$working_folder/transcripts.denovoNC.f2.CPC_results.noncoding.txt
export CPC_coding_loc=$working_folder/transcripts.denovoNC.f2.CPC_results.coding.txt

export transcripts_loc=$working_folder/transcripts.denovoNC.f2.bed
export transcripts_filt_loc=$working_folder/transcripts.denovoNC.f3.bed
export transcripts_coding_loc=$working_folder/transcripts.denovoNC.f2.potentially_coding.bed

Rscript $R_folder/filtering_Rscript_2.r

#exonic overlap vs any overlap - very tiny difference - just few transcripts (~10)

intersectBed -u -split -s -a $working_folder/transcripts.denovoNC.f3.bed -b $working_folder/transcripts.denovoNC.f2.potentially_coding.bed | wc -l
#842 
intersectBed -v -s -a $working_folder/transcripts.denovoNC.f3.bed -b $working_folder/transcripts.denovoNC.f2.potentially_coding.bed  > $working_folder/transcripts.denovoNC.f3a.bed
 wc -l $working_folder/transcripts.denovoNC.f3.bed
 #17011 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20210910_annotation/transcripts.denovoNC.f3.bed

wc -l $working_folder/transcripts.denovoNC.f3a.bed
#16140 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20210910_annotation/transcripts.denovoNC.f3a.bed


#cp $working_folder/transcripts.denovoNC.f2.potentially_coding.bed /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/ 
#cp $working_folder/transcripts.denovoNC.f3.bed /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/
#cp $working_folder/transcripts.denovoNC.f2.bed  /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/


#check for overlaps between categories 

 intersectBed  -u  -s -a $working_folder/transcripts.denovoNC.f3a.bed -b $working_folder/transcripts.TE_genes.bed | wc -l 
#99
 intersectBed -split -u  -s -a $working_folder/transcripts.denovoNC.f3a.bed -b $working_folder/transcripts.TE_genes.bed | wc -l 
#85

intersectBed -u  -s -a $working_folder/transcripts.denovoNC.f3a.bed -b $working_folder/transcripts.TE_frags.bed | sortBed -i stdin| mergeBed -s -i stdin | wc -l 
#28
intersectBed -split -u  -s -a $working_folder/transcripts.denovoNC.f3a.bed -b $working_folder/transcripts.TE_frags.bed | sortBed -i stdin| mergeBed -s -i stdin | wc -l 
#21

intersectBed -u -s -a $working_folder/transcripts.denovoNC.f3a.bed -b $working_folder/transcripts.f1.denovo_pseudogene.bed  | wc -l 
#25
intersectBed -u -s -a $working_folder/transcripts.denovoNC.f3a.bed -b $working_folder/denovoPC.transcripts.bed  | wc -l 
#191
intersectBed -u -split -s -a $working_folder/transcripts.denovoNC.f3a.bed -b $working_folder/denovoPC.transcripts.bed  | wc -l 
#0


intersectBed -u -split -s -a $working_folder/transcripts.denovoNC.f3a.bed -b $working_folder/transcripts.f1.denovo_pseudogene.bed  | wc -l 
#21
intersectBed -u -split -s -a $working_folder/transcripts.denovoNC.f3a.bed -b $working_folder/denovoPC.transcripts.bed  | wc -l 
#0



#remove all ambiguous NC rnas into a separate file 
#remove transcripts with strand exonic overlap with de novo TE_frags and de novo pseudogenes
intersectBed -v  -split -s -a $working_folder/transcripts.denovoNC.f3a.bed -b $working_folder/transcripts.TE_frags.bed >$working_folder/transcripts.denovoNC.f3b.bed
intersectBed -v  -split -s -a $working_folder/transcripts.denovoNC.f3b.bed -b $working_folder/transcripts.TE_genes.bed >$working_folder/transcripts.denovoNC.f3b1.bed
intersectBed -v -split -s -a $working_folder/transcripts.denovoNC.f3b1.bed -b $working_folder/transcripts.f1.denovo_pseudogene.bed >$working_folder/transcripts.denovoNC.f3c.bed 

wc -l $working_folder/transcripts.denovoNC.f3a.bed 
#16140 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.denovoNC.f3a.bed
wc -l $working_folder/transcripts.denovoNC.f3b.bed 
#16095
wc -l $working_folder/transcripts.denovoNC.f3b1.bed 
#16010
wc -l $working_folder/transcripts.denovoNC.f3c.bed 
#15989
sort -k1,1 -k2,2g $working_folder/transcripts.denovoNC.f3c.bed    |  mergeBed -s -i  stdin| wc -l 
#11156

#############

echo "number of preliminary lncRNA (CPC filtered) transcripts" >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.denovoNC.f3a.bed  >> $working_folder/linenumbers.txt

echo "number of preliminary lncRNA (CPC filtered) non-overlapping loci " >> $working_folder/linenumbers.txt
sort -k1,1 -k2,2g $working_folder/transcripts.denovoNC.f3a.bed  |  mergeBed -s -i  stdin| wc -l >> $working_folder/linenumbers.txt


echo "number of potentially coding lncRNA (CPC) transcripts" >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.denovoNC.f2.potentially_coding.bed >> $working_folder/linenumbers.txt

echo "number of potentially coding lncRNA (CPC) non-overlapping loci " >> $working_folder/linenumbers.txt
sort -k1,1 -k2,2g $working_folder/transcripts.denovoNC.f2.potentially_coding.bed |  mergeBed -s -i  stdin| wc -l >> $working_folder/linenumbers.txt


echo "number of preliminary lncRNA transcripts (no TE/pseudogene exonic overlap)" >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.denovoNC.f3c.bed     >> $working_folder/linenumbers.txt

echo "number of preliminary lncRNA loci (no transcripts with TE/pseudogene exonic overlap)" >> $working_folder/linenumbers.txt
sort -k1,1 -k2,2g $working_folder/transcripts.denovoNC.f3c.bed    |  mergeBed -s -i  stdin| wc -l >> $working_folder/linenumbers.txt

 wc -l $working_folder/transcripts.denovoNC.f3c.bed
#15989 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.denovoNC.f3c.bed

#############################################################
#filter out different types of RNA - rRNAs, tRNAs, snoRNAs  etc 
	
#TAIR10
protein_coding=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.protein_coding.gtf
lncRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.lncRNA.gtf
atlncRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.atlncRNA.gtf
miRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.miRNA.gtf
pre_miRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.pre_miRNA.gtf
tRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.tRNA.gtf
atRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.atRNA.gtf
snoRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.snoRNA.gtf
otherRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.otherRNA.gtf
snRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.snRNA.gtf
SRP_RNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.SRP_RNA.gtf
nontranslating_CDS=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.nontranslating_CDS.gtf
rRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.rRNA.gtf
ncRNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.ncRNA.gtf
RNase_MRP_RNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.RNase_MRP_RNA.gtf
sense_intronic=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.sense_intronic.gtf
antisense_RNA=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.antisense_RNA.gtf


#filter out ribosomal RNA and tRNA precursors 
intersectBed -s -split -u -a $working_folder/transcripts.bed -b $rRNA | sortBed -i stdin > $working_folder/transcripts.rRNA.bed

echo "number of de novo transcripts with exonic overlap to an rRNA" >> $working_folder/linenumbers.txt
wc -l  $working_folder/transcripts.rRNA.bed >> $working_folder/linenumbers.txt
sort -k1,1 -k2,2g $working_folder/transcripts.rRNA.bed |  mergeBed -s -i  stdin| wc -l 
#9

#tRNAs
intersectBed -s -split -u -a $working_folder/transcripts.bed -b $tRNA | sortBed -i stdin > $working_folder/transcripts.tRNA.bed 

echo "number of de novo transcripts with exonic overlap to a tRNA" >> $working_folder/linenumbers.txt
wc -l  $working_folder/transcripts.tRNA.bed >> $working_folder/linenumbers.txt
sort -k1,1 -k2,2g $working_folder/transcripts.tRNA.bed |  mergeBed -s -i  stdin| wc -l 
#46

# exclude ribosomal RNA and ribosomal RNA extensions 
intersectBed -s -v -a $working_folder/transcripts.denovoNC.f3c.bed -b $working_folder/transcripts.rRNA.bed $working_folder/transcripts.tRNA.bed   > $working_folder/transcripts.denovoNC.f4.bed

echo "number of de novo lncRNA transcripts" >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.denovoNC.f4.bed >> $working_folder/linenumbers.txt
echo "number of de novo lncRNA loci " >> $working_folder/linenumbers.txt
sort -k1,1 -k2,2g $working_folder/transcripts.denovoNC.f4.bed  |  mergeBed -s -i  stdin | wc -l >> $working_folder/linenumbers.txt
#11136 loci 

###############################################################
###############################################################

#rename loci and transcripts for NC 


#rename NC loci 

#
#redefine loci 
head $working_folder/transcripts.denovoNC.f4.bed
Chr1    5814    6640    CUFF.3.1        0       -       5814    6640    0       2       449,204,        0,622,
Chr1    33820   36039   CUFF.13.1       0       +       33820   36039   0       1       2219,   0,
Chr1    45074   46980   CUFF.5.1        0       +       45074   46980   0       1       1906,   0,
Chr1    51977   52743   CUFF.10.1       0       -       51977   52743   0       1       766,    0,
Chr1    63146   64476   CUFF.29.1       0       +       63146   64476   0       1       1330,   0,
Chr1    67285   67782   CUFF.30.1       0       +       67285   67782   0       1       497,    0,
Chr1    69606   71157   CUFF.15.1       0       +       69606   71157   0       2       1144,333,       0,1218,

#NC
#Cuffmerge loci names in the $working_folder/transcripts.f1.denovoPC.bed
cat $working_folder/transcripts.denovoNC.f4.bed | awk -v OFS="\t"  '{print $4,$5}' | awk -v OFS="\t"  '{split($1,s,".");   print s[1]"."s[2]}' |sort | uniq  | wc -l
#11443



cat $working_folder/transcripts.denovoNC.f4.bed | awk -v OFS="\t"  '{print $4,$5}' | awk -v OFS="\t"  '{split($1,s,".");   print s[1]"."s[2]}' |sort -V -k1,1 -k2,2 | uniq  > $working_folder/transcripts.denovoNC.f4.cuff_locinames.bed
wc -l $working_folder/transcripts.denovoNC.f4.cuff_locinames.bed
 wc -l $working_folder/transcripts.denovoNC.f4.cuff_locinames.bed
#11443 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/transcripts.denovoNC.f4.cuff_locinames.bed



 head /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20210513_annotation/transcripts.denovoNC.f4.cuff_locinames.bed
CUFF.3
CUFF.6
CUFF.10
CUFF.13
CUFF.15
CUFF.18
CUFF.20
CUFF.25
CUFF.29
CUFF.30



rm $working_folder/loci_table
touch 
count=1
while read line
do
locus_name=$line
cat  $working_folder/transcripts.denovoNC.f4.bed   |awk -v OFS="\t"  '{split($4,s,".");   print $0,s[1]"."s[2]}' | grep -w $line | awk -v OFS="\t"  '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > $working_folder/tmp_1
mergeBed  -c 4,5,6 -o distinct,first,first -s -i $working_folder/tmp_1 | awk -v locus_number="$count" '{if (NR==1) {print $0,"CUFF_NC."locus_number} else {print $0,"CUFF_NC."locus_number"-"NR}}'>> $working_folder/loci_table
mergeBed  -c 4,5,6 -o distinct,first,first -s -i $working_folder/tmp_1 | awk -v locus_number="$count" '{if (NR==1) {print $0,"CUFF_NC."locus_number} else {print $0,"CUFF_NC."locus_number"-"NR}}'
count=$((count+1))
done < $working_folder/transcripts.denovoNC.f4.cuff_locinames.bed

wc -l loci_table
# 11447 loci_table


head loci_table
Chr1    5814    6640    CUFF.3.1        0       - CUFF_NC.1
Chr1    45074   46980   CUFF.5.1        0       + CUFF_NC.2
Chr1    51977   52743   CUFF.10.1       0       - CUFF_NC.3
Chr1    33820   36039   CUFF.13.1       0       + CUFF_NC.4
Chr1    69606   71157   CUFF.15.1       0       + CUFF_NC.5
Chr1    86737   87289   CUFF.18.1       0       + CUFF_NC.6
Chr1    104712  105333  CUFF.20.1       0       + CUFF_NC.7
Chr1    73084   73716   CUFF.25.1       0       - CUFF_NC.8
Chr1    63146   64476   CUFF.29.1       0       + CUFF_NC.9
Chr1    67285   67782   CUFF.30.1       0       + CUFF_NC.10


cat $working_folder/loci_table > $working_folder/transcripts.denovoNC.loci_old-new_names.bed

cat $working_folder/loci_table |awk -v OFS="\t"  '{ print $1,$2,$3,$7,$5,$6}' | sortBed -i > $working_folder/denovoNC.loci.bed 
 head $working_folder/denovoNC.loci.bed
Chr1    5814    6640    CUFF_NC.1       0       -
Chr1    33820   36039   CUFF_NC.4       0       +
Chr1    45074   46980   CUFF_NC.2       0       +
Chr1    51977   52743   CUFF_NC.3       0       -
Chr1    63146   64476   CUFF_NC.9       0       +
Chr1    67285   67782   CUFF_NC.10      0       +
Chr1    69606   71157   CUFF_NC.5       0       +
Chr1    73084   73716   CUFF_NC.8       0       -
Chr1    75366   76624   CUFF_NC.12      0       -


cat $working_folder/loci_table  |  awk  -v OFS="\t"  'split($4,s,","){ for(i in s) {print s[i],$7,$7"."i  }}' | head 
CUFF.3.1        CUFF_NC.1       CUFF_NC.1.1
CUFF.5.1        CUFF_NC.2       CUFF_NC.2.1
CUFF.11.1       CUFF_NC.3       CUFF_NC.3.1
CUFF.13.1       CUFF_NC.4       CUFF_NC.4.1
CUFF.16.1       CUFF_NC.5       CUFF_NC.5.1
CUFF.18.1       CUFF_NC.6       CUFF_NC.6.1
CUFF.20.1       CUFF_NC.7       CUFF_NC.7.1
CUFF.23.1       CUFF_NC.8       CUFF_NC.8.1
CUFF.29.1       CUFF_NC.9       CUFF_NC.9.1
CUFF.30.1       CUFF_NC.10      CUFF_NC.10.1


cat $working_folder/loci_table  |  awk  -v OFS="\t"  'split($4,s,","){ for(i in s) {print s[i],$7,$7"."i  }}' > $working_folder/denovoNC.loci_transcripts_names_old_new.bed 


wc -l  $working_folder/denovoNC.loci_transcripts_names_old_new.bed
#15938 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/denovoNC.loci_transcripts_names_old_new.bed



# make annotation with new names 

sort -k 4b,4  $working_folder/transcripts.denovoNC.f4.bed  > $working_folder/old_names_annotation_sorted
sort -k 1b,1 $working_folder/denovoNC.loci_transcripts_names_old_new.bed  >  $working_folder/new_old_names_sorted
join -1 1 -2 4 $working_folder/new_old_names_sorted $working_folder/old_names_annotation_sorted |  awk -v OFS="\t" '{print $2,$3,$4,$5,$6,$1,$7,$8,$9,$10,$11,$12,$13,$14}' | sort -V -k1,1 -k2,2 > $working_folder/denovoNC.transcripts.old_names_new_names.bed

 head $working_folder/denovoNC.transcripts.old_names_new_names.bed
CUFF_NC.1       CUFF_NC.1.1     Chr1    5814    6640    CUFF.3.1        0       -       5814    6640    0       2  449,204, 0,622,
CUFF_NC.2       CUFF_NC.2.1     Chr1    45074   46980   CUFF.5.1        0       +       45074   46980   0       1  1906,    0,
CUFF_NC.3       CUFF_NC.3.1     Chr1    51977   52743   CUFF.10.1       0       -       51977   52743   0       1  766,     0,
CUFF_NC.4       CUFF_NC.4.1     Chr1    33820   36039   CUFF.13.1       0       +       33820   36039   0       1  2219,    0,
CUFF_NC.5       CUFF_NC.5.1     Chr1    69606   71157   CUFF.15.1       0       +       69606   71157   0       2  1144,333,        0,1218,
CUFF_NC.6       CUFF_NC.6.1     Chr1    86737   87289   CUFF.18.1       0       +       86737   87289   0       1  5
cat $working_folder/denovoNC.transcripts.old_names_new_names.bed |  awk -v OFS="\t" '{print $3,$4,$5,$2,$7,$8,$9,$10,$11,$12,$13,$14}'| sortBed -i stdin > $working_folder/denovoNC.transcripts.bed

head $working_folder/denovoNC.transcripts.bed
Chr1    5814    6640    CUFF_NC.1.1     0       -       5814    6640    0       2       449,204,        0,622,
Chr1    33820   36039   CUFF_NC.4.1     0       +       33820   36039   0       1       2219,   0,
Chr1    45074   46980   CUFF_NC.2.1     0       +       45074   46980   0       1       1906,   0,
Chr1    51977   52743   CUFF_NC.3.1     0       -       51977   52743   0       1       766,    0,
Chr1    63146   64476   CUFF_NC.9.1     0       +       63146   64476   0       1       1330,   0,
Chr1    67285   67782   CUFF_NC.10.1    0       +       67285   67782   0       1       497,    0,
Chr1    69606   71157   CUFF_NC.5.1     0       +       69606   71157   0       2       1144,333,       0,1218,
Chr1    73084   73716   CUFF_NC.8.1     0       -       73084   73716   0       2       90,219, 0,413,

wc -l $working_folder/denovoNC.transcripts.bed 
#15938 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/denovoNC.transcripts.bed

wc -l $working_folder/denovoNC.loci.bed 
#11447 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/denovoNC.loci.bed


###############################################################
# classify lncRNA loci by position
#operate on loci level, then make the transcript file
###############################################################

#classify lncRNAs
#1. miRNA precursors
#2. sno/snRNA precursors
#3. Antisense 
#4. Antisense to pseudogene
#5. Antisense to TE genes and TE_frags 
#6. Intergenic

#Position based classification 
	
#1. miRNA precursor 
intersectBed -s  -u -a $working_folder/denovoNC.loci.bed -b $miRNA | sortBed -i stdin > $working_folder/lncRNAs.miRNA_precursors.loci.bed
cat $working_folder/lncRNAs.miRNA_precursors.loci.bed | awk -v OFS="\t" '{print $4}' > $working_folder/lncRNAs.miRNA_precursors.loci_names.txt
 cat $working_folder/denovoNC.transcripts.bed | grep -w -f $working_folder/lncRNAs.miRNA_precursors.loci_names.txt > $working_folder/lncRNAs.miRNA_precursors.transcripts.bed

#2. snoRNA and snRNA precursors
intersectBed -s  -u -a $working_folder/denovoNC.loci.bed -b $snRNA $snoRNA | sortBed -i stdin > $working_folder/lncRNAs.sn_snoRNA_precursors.loci.bed
cat $working_folder/lncRNAs.sn_snoRNA_precursors.loci.bed | awk -v OFS="\t" '{print $4}' > $working_folder/lncRNAs.sn_snoRNA_precursors.loci_names.txt
cat $working_folder/denovoNC.transcripts.bed | grep -w -f $working_folder/lncRNAs.sn_snoRNA_precursors.loci_names.txt > $working_folder/lncRNAs.sn_snoRNA_precursors.transcripts.bed

 
#3. Antisense 
intersectBed -s  -v -a $working_folder/denovoNC.loci.bed -b $miRNA $snRNA $snoRNA |intersectBed -S -u -a stdin -b $PC_tair $PC_araport  | sortBed -i stdin > $working_folder/lncRNAs.antisense.loci.bed

cat $working_folder/lncRNAs.antisense.loci.bed | awk -v OFS="\t" '{print $4}' > $working_folder/lncRNAs.antisense.loci_names.txt
 
wc -l $working_folder/lncRNAs.antisense.loci_names.txt
#8195 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/lncRNAs.antisense.loci_names.txt


cat $working_folder/denovoNC.transcripts.bed | grep -w -f $working_folder/lncRNAs.antisense.loci_names.txt > $working_folder/lncRNAs.antisense.transcripts.bed
wc -l $working_folder/lncRNAs.antisense.transcripts.bed
#10934 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/lncRNAs.antisense.transcripts.bed


#	4. Antisense to pseudogene  (should not overlap with antisense (to PC) lncRNAs)
intersectBed -s  -v -a $working_folder/denovoNC.loci.bed -b $working_folder/lncRNAs.antisense.loci.bed $miRNA $snRNA $snoRNA | intersectBed -S -u -a stdin -b $pseudogenes   | sortBed -i stdin > $working_folder/lncRNAs.AS_to_pseudogenes.loci.bed
cat $working_folder/lncRNAs.AS_to_pseudogenes.loci.bed | awk -v OFS="\t" '{print $4}' > $working_folder/lncRNAs.AS_to_pseudogenes.loci_names.txt
cat $working_folder/denovoNC.transcripts.bed | grep -w -f $working_folder/lncRNAs.AS_to_pseudogenes.loci_names.txt > $working_folder/lncRNAs.AS_to_pseudogenes.transcripts.bed
wc -l  $working_folder/lncRNAs.AS_to_pseudogenes.transcripts.bed
#185 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/lncRNAs.AS_to_pseudogenes.transcripts.bed


# 5.Antisense to TE genes 
export TE_genes=$annotationfolder/Araport11_transposable_element_gene.201606.bed

intersectBed -s  -v -a $working_folder/denovoNC.loci.bed -b $working_folder/lncRNAs.antisense.loci.bed $working_folder/lncRNAs.AS_to_pseudogenes.loci.bed $miRNA $snRNA $snoRNA | intersectBed -S -u -a stdin -b $TE_genes | sortBed -i stdin > $working_folder/lncRNAs.AS_to_TE.loci.bed
cat $working_folder/lncRNAs.AS_to_TE.loci.bed | awk -v OFS="\t" '{print $4}' > $working_folder/lncRNAs.AS_to_TE.loci_names.txt
cat $working_folder/denovoNC.transcripts.bed | grep -w -f $working_folder/lncRNAs.AS_to_TE.loci_names.txt > $working_folder/lncRNAs.AS_to_TE.transcripts.bed

 wc -l $working_folder/lncRNAs.AS_to_TE.loci.bed
#630 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/lncRNAs.AS_to_TE.loci.bed

##!! not only annotated genes but also de novo PC (PC gene extensions)
# 6. Intergenic 
denovoPC_with_chimeras=$working_folder/transcripts.f1.denovoPC.with_chimeras.bed
export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation
export PC_tair=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.protein_coding.sorted.bed
export PC_araport=$annotationfolder/Araport11_protein_coding.201606.sorted.bed
export pseudogenes=$annotationfolder/Araport11_pseudogene.201606.sorted.bed

#other NC loci 
cat $working_folder/lncRNAs.AS_to_TE.loci.bed $working_folder/lncRNAs.antisense.loci.bed $working_folder/lncRNAs.AS_to_pseudogenes.loci.bed  $working_folder/lncRNAs.miRNA_precursors.loci.bed  $working_folder/lncRNAs.sn_snoRNA_precursors.loci.bed | awk -v OFS="\t" '{print $4}' > $working_folder/lncRNAs.non_linc.loci_names.txt

#lincRNAs
cat $working_folder/denovoNC.loci.bed | grep -w -v -f $working_folder/lncRNAs.non_linc.loci_names.txt | intersectBed -v -a stdin -b $PC_tair $PC_araport $pseudogenes  |intersectBed -split -v -a stdin -b $denovoPC_with_chimeras |sortBed -i stdin  > $working_folder/lncRNAs.intergenic.loci.preliminary.bed
wc -l $working_folder/lncRNAs.intergenic.loci.preliminary.bed
#2285 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/lncRNAs.intergenic.loci.preliminary.bed
 

intersectBed -v -s -a $f1 -b $PC_tair  $PC_araport  $pseudogenes |  intersectBed -v -split -s -a stdin -b  $denovoPC_with_chimeras|coverageBed -split

#remove transcripts that are too close (and downstream) to the PC genes and pseudogenes (might be trnascription read-through)

bedtools closest -a $working_folder/lncRNAs.intergenic.loci.preliminary.bed -b $PC_tair $PC_araport $pseudogenes -s -D b -iu| awk '($20<100){print $4}' | uniq > $working_folder/lincRNAs.potential_PC_runons.names.txt
wc -l $working_folder/lincRNAs.potential_PC_runons.names.txt
#39 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/lincRNAs.potential_PC_runons.names.txt




cat $working_folder/lncRNAs.intergenic.loci.preliminary.bed | grep -v -w -f $working_folder/lincRNAs.potential_PC_runons.names.txt >$working_folder/lncRNAs.intergenic.loci.bed 
 
 wc -l $working_folder/lncRNAs.intergenic.loci.bed
#2246 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/lncRNAs.intergenic.loci.bed

 
cat $working_folder/lncRNAs.intergenic.loci.bed | awk -v OFS="\t" '{print $4}' > $working_folder/lncRNAs.intergenic.loci_names.txt
cat $working_folder/denovoNC.transcripts.bed | grep -w -f $working_folder/lncRNAs.intergenic.loci_names.txt > $working_folder/lncRNAs.intergenic.transcripts.bed
wc -l $working_folder/lncRNAs.intergenic.transcripts.bed
 
# 3207 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/lncRNAs.intergenic.transcripts.bed

 
 cat lncRNAs.intergenic.loci.bed lncRNAs.antisense.loci.bed lncRNAs.AS_to_TE.loci.bed lncRNAs.AS_to_pseudogenes.loci.bed lncRNAs.sn_snoRNA_precursors.loci.bed lncRNAs.miRNA_precursors.loci.bed | sortBed -i > lncRNAs.6types.loci.bed

 cat lncRNAs.6types.loci.bed |awk -v OFS="\t" '{print $4}' | uniq | wc -l 
 #11295
wc -l lncRNAs.6types.loci.bed
#11295 lncRNAs.6types.loci.bed

wc -l denovoNC.loci.bed
#11565 denovoNC.loci.bed

 
#######################
 
#check overlap between categories 

#number of transcripts/loci after overlap between categories clean 
 
 
 
 wc -l lncRNAs*.loci.bed
    11295 lncRNAs.6types.loci.bed
  8195 lncRNAs.antisense.loci.bed
   121 lncRNAs.AS_to_pseudogenes.loci.bed
   630 lncRNAs.AS_to_TE.loci.bed
  2246 lncRNAs.intergenic.loci.bed
    65 lncRNAs.miRNA_precursors.loci.bed
    38 lncRNAs.sn_snoRNA_precursors.loci.bed



	 wc -l lncRNAs*.transcripts.bed
  15964 lncRNAs.6types.transcripts.bed
  10934 lncRNAs.antisense.transcripts.bed
    185 lncRNAs.AS_to_pseudogenes.transcripts.bed
   1107 lncRNAs.AS_to_TE.transcripts.bed
   3207 lncRNAs.intergenic.transcripts.bed
    250 lncRNAs.miRNA_precursors.transcripts.bed
     62 lncRNAs.sn_snoRNA_precursors.transcripts.bed

cat lncRNAs.antisense.transcripts.bed lncRNAs.AS_to_pseudogenes.transcripts.bed  lncRNAs.AS_to_TE.transcripts.bed lncRNAs.intergenic.transcripts.bed lncRNAs.miRNA_precursors.transcripts.bed  lncRNAs.sn_snoRNA_precursors.transcripts.bed | sortBed -i > lncRNAs.6types.transcripts.bed 
  
wc -l lncRNAs.6types.transcripts.bed 
#15745 lncRNAs.6types.transcripts.bed

wc -l denovoNC.transcripts.bed
#15938
 
 
#known lncRNAs (present in Araport11) 
intersectBed -s -u -a $working_folder/lncRNAs.6types.transcripts.bed  -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed $annotationfolder/Araport11_novel_transcribed_region.201606.bed > $working_folder/lncRNAs.6types.transcripts.known.bed 
intersectBed -s -u -a $working_folder/lncRNAs.6types.loci.bed -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed $annotationfolder/Araport11_novel_transcribed_region.201606.bed > $working_folder/lncRNAs.6types.loci.known.bed 
 
 
 #novel lncRNAs (not present in Araport11)
intersectBed -s -v -a $working_folder/lncRNAs.6types.transcripts.bed  -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed $annotationfolder/Araport11_novel_transcribed_region.201606.bed > $working_folder/lncRNAs.6types.transcripts.novel.bed 

intersectBed -s -v -a $working_folder/lncRNAs.6types.loci.bed -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed $annotationfolder/Araport11_novel_transcribed_region.201606.bed> $working_folder/lncRNAs.6types.loci.novel.bed 


wc -l *known*
980 lncRNAs.6types.loci.known.bed
2606 lncRNAs.6types.transcripts.known.bed

  
wc -l *novel*
10315 lncRNAs.6types.loci.novel.bed
13139 lncRNAs.6types.transcripts.novel.bed

 
 #lncRNAs 
 # all 
 #$working_folder/lncRNAs.6types.loci.bed
 #$working_folder/lncRNAs.6types.transcripts.bed 
 
 # types 
 
  # lncRNAs.antisense.loci.bed
   # lncRNAs.AS_to_pseudogenes.loci.bed
  #  lncRNAs.AS_to_TE.loci.bed
  # lncRNAs.intergenic.loci.bed
  #   lncRNAs.miRNA_precursors.loci.bed
  #   lncRNAs.sn_snoRNA_precursors.loci.bed
 
 #create SAF for feature counts 
 GeneID	Chr	Start	End	Strand
497097	chr1	3204563	3207049	-
497097	chr1	3411783	3411982	-
497097	chr1	3660633	3661579	-
...
cd $working_folder/

#make a list of transcripts 
1. lncRNAs.antisense
2. lncRNAs.AS_to_pseudogenes
3. lncRNAs.AS_to_TE
4. lncRNAs.intergenic
5. lncRNAs.miRNA_precursors
6. denovo_pseudogene
7. denovoPC
8. expressed_TEs
9.expressed TE genes

lncRNAs.antisense.transcripts.bed
lncRNAs.AS_to_pseudogenes.transcripts.bed
lncRNAs.AS_to_TE.transcripts.bed
lncRNAs.intergenic.transcripts.bed
lncRNAs.miRNA_precursors.transcripts.bed
lncRNAs.sn_snoRNA_precursors.transcripts.bed


transcripts.TE_genes.bed
transcripts.TE_frags.bed
transcripts.f1.denovo_pseudogene.bed

denovoPC.transcripts.bed


cat  lncRNAs.antisense.transcripts.bed lncRNAs.AS_to_pseudogenes.transcripts.bed lncRNAs.AS_to_TE.transcripts.bed lncRNAs.intergenic.transcripts.bed lncRNAs.miRNA_precursors.transcripts.bed lncRNAs.sn_snoRNA_precursors.transcripts.bed transcripts.TE_genes.bed transcripts.TE_frags.bed transcripts.f1.denovo_pseudogene.bed denovoPC.transcripts.bed  | sortBed -i stdin > 20211013_annotation.transcripts.bed
wc -l 20211013_annotation.transcripts.bed
#98836 20211013_annotation.transcripts.bed



cat 20211013_annotation.transcripts.bed| bed12ToBed6 -i stdin |  awk -v OFS="\t"  '{split($4,s,".");print s[1]"."s[2],$1,$2,$3,$6}'>20211013_annotation.transcripts.saf 


/groups/nordborg/projects/cegs/alexandra/software/bedToExons 20211013_annotation.transcripts.bed 20211013_annotation.transcripts.exons
cat 20211013_annotation.transcripts.bed | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'> 20211013_annotation.transcripts.6columns.bed
cat 20211013_annotation.transcripts.exons | awk -v OFS="\t" '{print $1,$2+1,$3,$4,$5,$6,"exon"}' > gtf
cat 20211013_annotation.transcripts.6columns.bed | awk -v OFS="\t" '{print $1,$2+1,$3,$4,$5,$6,"transcript"}' >> gtf

 cat gtf |  awk  -v OFS="\t"  '{split($4,s,"."); print $0,s[1]"."s[2]}'| awk  -v OFS="\t"  '{print $1,"denovo_annotation",$7,$2,$3,0, $6,".", "gene_id \""$8"\"; transcript_id \""$4"\";"}' >20211013_annotation.transcripts.gtf

 # export loci and transcript numbers
 PC_loci=`cat $working_folder/denovoPC.loci.bed| wc -l`
 PC_tr=`cat $working_folder/denovoPC.transcripts.bed| wc -l`
  
 AS_loci=`cat $working_folder/lncRNAs.antisense.loci.bed| wc -l`
 AS_tr=`cat $working_folder/lncRNAs.antisense.transcripts.bed| wc -l`

 linc_loci=`cat $working_folder/lncRNAs.intergenic.loci.bed| wc -l`
 linc_tr=`cat  $working_folder/lncRNAs.intergenic.transcripts.bed| wc -l`
 
 AS_to_TE_loci=`cat $working_folder/lncRNAs.AS_to_TE.loci.bed| wc -l`
 AS_to_TE_tr=`cat $working_folder/lncRNAs.AS_to_TE.transcripts.bed| wc -l`
 
 AS_to_pseudo_loci=`cat $working_folder/lncRNAs.AS_to_pseudogenes.loci.bed| wc -l`
 AS_to_pseudo_tr=`cat $working_folder/lncRNAs.AS_to_pseudogenes.transcripts.bed| wc -l`

 
 TEgene_loci=`cat $working_folder/transcripts.TE_genes.bed|   mergeBed -s -i  stdin|wc -l`
 TEgene_tr=`cat  $working_folder/transcripts.TE_genes.bed| wc -l`

 TEs_tr=`cat $working_folder/transcripts.TE_frags.bed| wc -l`
 TEs_loci=`cat $working_folder/transcripts.TE_frags.bed |  mergeBed -s -i  stdin| wc -l`

 pot_cod_tr=`cat $working_folder/transcripts.denovoNC.f2.potentially_coding.bed| wc -l`
 pot_cod_loci=`sort -k1,1 -k2,2g $working_folder/transcripts.denovoNC.f2.potentially_coding.bed |  mergeBed -s -i  stdin| wc -l `

 miRNA_loci=`cat  $working_folder/lncRNAs.miRNA_precursors.loci.bed| wc -l`
 miRNA_tr=`cat  $working_folder/lncRNAs.miRNA_precursors.transcripts.bed| wc -l`
 
 snoRNA_loci=`cat  $working_folder/lncRNAs.sn_snoRNA_precursors.loci.bed | wc -l`
 snoRNA_tr=`cat  $working_folder/lncRNAs.sn_snoRNA_precursors.transcripts.bed | wc -l`
 
 pseudo_loci=`cat  $working_folder/transcripts.f1.denovo_pseudogene.bed| mergeBed -s -i  stdin|wc -l`
 pseudo_tr=`cat  $working_folder/transcripts.f1.denovo_pseudogene.bed| wc -l`
 
 lncRNA_6types_loci=`cat  $working_folder/lncRNAs.6types.loci.bed| wc -l`
 lncRNA_6types_tr=`cat  $working_folder/lncRNAs.6types.transcripts.bed| wc -l`
 
  rRNA_loci=`cat  $working_folder/transcripts.rRNA.bed| mergeBed -s -i  stdin| wc -l`
 rRNA_tr=`cat  $working_folder/transcripts.rRNA.bed| wc -l`
 
 
  tRNA_loci=`cat  $working_folder/transcripts.tRNA.bed| mergeBed -s -i  stdin| wc -l`
 tRNA_tr=`cat  $working_folder/transcripts.tRNA.bed| wc -l`
 
 
 intronic_loci=`cat  $working_folder/transcripts.f1.senseoverlapping.bed| mergeBed -s -i  stdin| wc -l`
 intronic_tr=`cat  $working_folder/transcripts.f1.senseoverlapping.bed| wc -l`

 echo -e "PC_loci\t"  "PC_tr\t"   "pseudo_loci\t" "pseudo_tr\t" "pot_cod_loci\t"  "pot_cod_tr\t" \
 "lncRNA_6types_loci\t"  "lncRNA_6types_tr\t"  "AS_loci\t"  "AS_tr\t"  "linc_loci\t"  "linc_tr\t" "AS_to_TE_loci\t"  "AS_to_TE_tr\t" \
 "rRNA_loci\t"  "rRNA_tr\t" "tRNA_loci\t" "tRNA_tr\t" "intronic_loci\t" "intronic_tr\t" \
 "AS_to_pseudo_loci\t"  "AS_to_pseudo_loci_tr\t"  "TEgene_loci\t"  "TEgene_tr\t"  "TEs_loci\t"  "TEs_tr\t" \
 "miRNA_loci\t"  "miRNA_tr\t" "snoRNA_loci\t"  "snoRNA_tr\t"  > $working_folder/$name_assembly.gene_numbers.txt

echo -e  $PC_loci"\t" $PC_tr"\t" $pseudo_loci"\t" $pseudo_tr"\t" $pot_cod_loci"\t" $pot_cod_tr"\t"\
$lncRNA_6types_loci"\t" $lncRNA_6types_tr"\t" $AS_loci"\t" $AS_tr"\t"  $linc_loci"\t" $linc_tr"\t"  $AS_to_TE_loci"\t" $AS_to_TE_tr"\t"\
$rRNA_loci"\t" $rRNA_tr"\t" $tRNA_loci"\t" $tRNA_tr"\t" $intronic_loci"\t" $intronic_tr"\t"\
$AS_to_pseudo_loci"\t" $AS_to_pseudo_tr"\t" $TEgene_loci"\t" $TEgene_tr"\t" $TEs_loci"\t" $TEs_tr"\t"\
$miRNA_loci"\t" $miRNA_tr"\t" $snoRNA_loci"\t" $snoRNA_tr"\t">> $working_folder/$name_assembly.gene_numbers.txt


 
 cp -R $working_folder/  /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/ 
