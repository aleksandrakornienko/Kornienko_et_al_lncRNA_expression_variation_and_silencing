#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=08:00:00
#SBATCH --array=1-96
#SBATCH --output=/groups/nordborg/user/aleksandra.kornienko/analyses/logs/array_job_slurm_%A_%a.out

#lists of assemblies for the first step of merging - shuffling done on 20221109
#folder_ident_pipeline=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/01_Identification_pipeline/2021
#seedling_cortijo=$folder_ident_pipeline/assemblies_merge_seedling_cortijo_may.txt
#dir=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/saturation_curve/cortijo 

#for i in {1..12}
#do
#for j in {1..8}
#do
#name_assembly=sat_curve.$i.pools.cortijo.rep.$j
#echo $name_assembly
#mkdir $dir/$name_assembly
#shuf $seedling_cortijo | head -$i > $dir/$name_assembly/$name_assembly.assemblies_for_merge.txt
#done 
#done 

export numberofline=$SLURM_ARRAY_TASK_ID

#####

samples=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/02_Identification_Saturation_Curve/2021/cortijo.samples.sat.curve.txt

export name_assembly=`sed -n $numberofline,"$numberofline"p $samples`
echo $name_assembly

################################

#the folder where the output annotation will be :
export  working_folder=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/saturation_curve/cortijo/$name_assembly
#mkdir $working_folder
cd $working_folder


export R_folder=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/short_Rscripts
genome=/groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/TAIR10_all.fa
chr_sizes=/groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/chr_length_wo_ChrM_ChrC.txt

export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation


ml bedtools/2.27.1-foss-2018b
module load r/3.5.1-foss-2018b

##########
ml cufflinks/2.2.1-foss-2018b
ml python/2.7.15-foss-2018b
#merge with cuffmerge (without reference!) 

cuffmerge --keep-tmp -s $genome --min-isoform-fraction 0 -p 4 -o $working_folder $working_folder/$name_assembly.assemblies_for_merge.txt

 #convert to bed 
file=transcripts
cat  $working_folder/$file.gtf |  awk -v OFS="\t" '!($7=="."){print $0}' > $working_folder/$file.nodots.gtf
/groups/nordborg/projects/cegs/alexandra/software/gtfToGenePred $working_folder/$file.nodots.gtf $working_folder/$file.genepred
/groups/nordborg/projects/cegs/alexandra/software/genePredToBed $working_folder/$file.genepred $working_folder/$file.1.bed
cat $working_folder/$file.1.bed | awk -v OFS="\t"  '{$7=$2;print $0}'| sed -e 's/chloroplast/ChrC/g'| sed -e 's/mitochondria/ChrM/g'  > $working_folder/$file.bed  


#############

echo "transcripts.bed transcript number">  $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.bed >>  $working_folder/linenumbers.txt
echo "transcripts.bed non-overlapping loci number">>  $working_folder/linenumbers.txt
cat $working_folder/transcripts.bed  | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l  >> $working_folder/linenumbers.txt
  
  
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


 ######################
 #1. artefact filtering 
 ######################
 
# 1. remove transcripts with mRNA length <200nt, remove single-exon transcripts <400nt 
cat $working_folder/transcripts.bed | \
awk '{split($11,s,","); total=0; for(x in s){ total+=s[x];} print total,$0}'  | awk -v OFS="\t" '$1>=200 {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13 }' \
  | awk -v OFS="\t" '!(($3-$2)<=400 && $10==1) {print $0 }' > $working_folder/transcripts.f1.bed 

echo "number of transcripts after removing transcripts <200nt and single-exon transcripts <400nt " >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.f1.bed >> $working_folder/linenumbers.txt
echo "number of loci" >> $working_folder/linenumbers.txt
cat $working_folder/transcripts.f1.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l  >> $working_folder/linenumbers.txt

#2. Additional antisense filter 

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


#cat $working_folder/mirror_transcripts_spliced_artefacts_names.txt $working_folder/chimeric_transcripts_artedacts_names.txt $working_folder/potential_AS_artifacts.txt |  uniq > $working_folder/potential_artifacts.txt
cat $working_folder/mirror_transcripts_spliced_artefacts_names.txt $working_folder/potential_AS_artifacts.txt |  uniq > $working_folder/potential_artifacts.txt
 
 
export transcripts_tofiltert_loc=$working_folder/potential_artifacts.txt
export transcripts_loc=$working_folder/transcripts.f1.bed  
export transcripts_filtered=$working_folder/transcripts.f1A.bed  

#cp /net/gmi.oeaw.ac.at/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/Identification_Saturation_Curve/short_Rscripts/filtering_Rscript_1A.r $R_folder

# R script to remove artefacts 
Rscript $R_folder/filtering_Rscript_1A.r  

#write down number of transcripts/loci
echo "number of transcripts after (stringently) removing potential AS and chimeric atrifacts transcripts " >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.f1A.bed >> $working_folder/linenumbers.txt
echo "number of loci after (stringently) removing potential AS and chimeric atrifacts transcripts" >> $working_folder/linenumbers.txt
cat $working_folder/transcripts.f1A.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l  >> $working_folder/linenumbers.txt

wc -l $working_folder/transcripts.f1A.bed

###################
# 3. splitting annotation into PC and non-coding part 
###################
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


#make denovo PC annotation: 
#(there will be many chimeric loci!! so do not take this PC annotation too seriously, use Araport)
f1=$working_folder/transcripts.f1A.bed
intersectBed -u -split -s -a $f1 -b $PC_tair $PC_araport | intersectBed -v -split -s -a stdin -b $tRNA $rRNA $snRNA $snoRNA|sortBed -i stdin > $working_folder/transcripts.f1.denovoPC.with_chimeras.bed

#will need this fuller denovoPC annotation later for excluding PC gene extensions exons 
 
wc -l $working_folder/transcripts.f1.denovoPC.with_chimeras.bed


# Remove chimeric transcript assemblies merging 2 or more genes
# remove chimeric PC-PC PC-pseudogene and pseudo-pseudo transcripts

export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation


 #sortBed -i $annotationfolder/Araport11_pseudogene.201606.bed | mergeBed -s -c 4,5,6 -o first -i stdin > $annotationfolder/Araport11_pseudogene.201606.mergebed.loci.bed

#intersectBed -split -s -v -a $annotationfolder/Araport11_pseudogene.201606.mergebed.loci.bed -b /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.mergeBed.loci.bed |awk -v OFS="\t" '{split($4,a,".");print $1,$2,$3,a[1],$5,$6}'  >  $annotationfolder/Araport11_pseudogene.201606.mergebed.loci.nooverlap_with_Ar11PC.bed

 intersectBed -u -s -split -a $annotationfolder/Araport11_pseudogene.201606.mergebed.loci.nooverlap_with_Ar11PC.bed -b /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes_nooverlapping_genes.bed | wc -l
#0
 cat $annotationfolder/Araport11_pseudogene.201606.mergebed.loci.nooverlap_with_Ar11PC.bed  /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes_nooverlapping_genes.bed | sortBed -i stdin > $working_folder/pseudo_and_PC_for_chimeras_filt.bed

intersectBed -s -f 0.3 -e -F 0.3 -wo -a $working_folder/transcripts.f1.denovoPC.with_chimeras.bed  -b pseudo_and_PC_for_chimeras_filt.bed | awk -v OFS="\t" '{print $4,$16}' | uniq|awk -v OFS="\t" '{print $1}'| uniq -d >  $working_folder/chimeric_transcripts_artefacts_names.txt




export transcripts_tofiltert_loc=$working_folder/chimeric_transcripts_artefacts_names.txt
export transcripts_loc=$working_folder/transcripts.f1.denovoPC.with_chimeras.bed 
export transcripts_filtered=$working_folder/transcripts.f1.denovoPC.without_chimeras.bed 

#cp /net/gmi.oeaw.ac.at/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/Identification_Saturation_Curve/short_Rscripts/filtering_Rscript_1A.r $R_folder

# R script to remove artefacts 
Rscript $R_folder/filtering_Rscript_1A.r  

wc -l $working_folder/transcripts.f1.denovoPC.without_chimeras.bed

#redefine PC loci 

#PC
#Cuffmerge loci names in the $working_folder/transcripts.f1.denovoPC.bed


cat $working_folder/transcripts.f1.denovoPC.without_chimeras.bed | awk -v OFS="\t"  '{print $4,$5}' | awk -v OFS="\t"  '{split($1,s,".");   print s[1]"."s[2]}' |sort -V -k1,1 -k2,2 | uniq  > $working_folder/transcripts.f1.denovoPC.cuff_locinames.bed


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

wc -l $working_folder/loci_table

cat $working_folder/loci_table > $working_folder/transcripts.f1.denovoPC.loci_old-new_names.bed

cat $working_folder/loci_table |awk -v OFS="\t"  '{ print $1,$2,$3,$7,$5,$6}'> $working_folder/denovoPC.loci.bed 
cat $working_folder/denovoPC.loci.bed  > inter 
sortBed -i inter > $working_folder/denovoPC.loci.bed 


cat $working_folder/loci_table  |  awk  -v OFS="\t"  'split($4,s,","){ for(i in s) {print s[i],$7,$7"."i  }}' > $working_folder/denovoPC.loci_transcripts_names_old_new.bed 
wc -l  $working_folder/denovoPC.loci_transcripts_names_old_new.bed

# make annotation with new names 

sort -k 4b,4  $working_folder/transcripts.f1.denovoPC.without_chimeras.bed  > $working_folder/old_names_annotation_sorted
sort -k 1b,1 $working_folder/denovoPC.loci_transcripts_names_old_new.bed  >  $working_folder/new_old_names_sorted
join -1 1 -2 4 $working_folder/new_old_names_sorted $working_folder/old_names_annotation_sorted |  awk -v OFS="\t" '{print $2,$3,$4,$5,$6,$1,$7,$8,$9,$10,$11,$12,$13,$14}' | sort -V -k1,1 -k2,2 > $working_folder/denovoPC.transcripts.old_names_new_names.bed

cat $working_folder/denovoPC.transcripts.old_names_new_names.bed |  awk -v OFS="\t" '{print $3,$4,$5,$2,$7,$8,$9,$10,$11,$12,$13,$14}' | sortBed -i stdin > $working_folder/denovoPC.transcripts.bed

 wc -l $working_folder/denovoPC.transcripts.bed
 
echo "number of denovo protein-coding transcripts" >> $working_folder/linenumbers.txt
wc -l $working_folder/denovoPC.transcripts.bed >> $working_folder/linenumbers.txt

echo "number of denovo protein-coding loci" >> $working_folder/linenumbers.txt
wc -l $working_folder/denovoPC.loci.bed  >> $working_folder/linenumbers.txt

 wc -l $working_folder/denovoPC.loci.bed

 
denovoPC=$working_folder/denovoPC.transcripts.bed


####################
#Define pseudogenes 
####################

#will not allow PC and pseudogenes to share exons in my annotation although there are 45 transcripts in araport pseudogene annotation that have a sense overlap with PC genes 
 #intersectBed -s  -u -a $pseudogenes -b  $PC_araport $PC_tair | wc -l                               
# 46
#do not allow overlap with trna,snRNA snoRNA and rRNA 

intersectBed -u -split -s -a $f1 -b $pseudogenes | intersectBed -v  -s -a stdin -b $denovoPC $PC_araport $PC_tair | intersectBed -v -split -s -a stdin -b $tRNA $rRNA $snRNA $snoRNA| sortBed -i stdin  > $working_folder/transcripts.f1.denovo_pseudogene.bed
wc -l $working_folder/transcripts.f1.denovo_pseudogene.bed

echo "number of denovo pseudogene transcripts" >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.f1.denovo_pseudogene.bed >> $working_folder/linenumbers.txt

echo "number of denovo pseudogene loci" >> $working_folder/linenumbers.txt
cat $working_folder/transcripts.f1.denovo_pseudogene.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l >> $working_folder/linenumbers.txt

#######################################
# define sense overlapping transcripts
#######################################

intersectBed -u -s -a $f1 -b $PC_tair $PC_araport | intersectBed -v -s -split -a stdin -b $PC_tair $PC_araport|intersectBed -v -split -s -a stdin -b $tRNA $rRNA $snRNA $snoRNA | sortBed -i stdin > $working_folder/transcripts.f1.senseoverlapping.bed
wc -l $working_folder/transcripts.f1.senseoverlapping.bed

echo "number of denovo senseoverlapping transcripts" >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.f1.senseoverlapping.bed >> $working_folder/linenumbers.txt
echo "number of denovo senseoverlapping loci" >> $working_folder/linenumbers.txt
cat $working_folder/transcripts.f1.senseoverlapping.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l >> $working_folder/linenumbers.txt



######################
#TEs				##
######################

#TE genes
denovoPC_with_chimeras=$working_folder/transcripts.f1.denovoPC.with_chimeras.bed 

#define TE genes - transcripts that have an exonic overlap with an Araport 11 TE gene
intersectBed -v -s -a $f1 -b $PC_tair $PC_araport $pseudogenes|   intersectBed -v -split -s -a stdin -b  $denovoPC_with_chimeras|intersectBed -v -split -s -a stdin -b $tRNA $rRNA $snRNA $snoRNA |   intersectBed -u -s  -split -a stdin -b  $TE_genes | sortBed -i stdin >  $working_folder/transcripts.TE_genes.bed
wc -l  $working_folder/transcripts.TE_genes.bed

#expressed TE fragments 
#Define TE fragment transcripts (arbitrary separation from lncRNAs - threshold 20% exon sequence)
intersectBed -v -s -a $f1 -b $PC_tair $PC_araport $pseudogenes|   intersectBed -v -split -s -a stdin -b  $denovoPC_with_chimeras| intersectBed -v -split -s -a stdin -b $tRNA $rRNA $snRNA $snoRNA | intersectBed -v -s -split  -a stdin -b  $TE_genes | sortBed -i stdin  >  $working_folder/transcripts.TEs.prelim.bed


#################
# definition of lncRNAs : 
# join -1 4 -2 1 $working_folder/tab3_sorted $working_folder/tab4_sorted |  awk -v OFS="\t" '($13<0.6 &&($13<0.5 || $15<0.8)){print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12}'|sortBed -i stdin > $working_folder/transcripts.denovoNC.f2.bed
##################

bed12ToBed6 -i $working_folder/transcripts.TEs.prelim.bed > $working_folder/transcripts.TEs.prelim.bed6

coverageBed -split -s -a  $working_folder/transcripts.TEs.prelim.bed6 -b $TE_fragments | awk -v OFS="\t" '{arr[$4]+=$8} END {for (i in arr) {print i,arr[i]}}' | sort  -k1,1n > $working_folder/exon_length_covered

coverageBed -split -s -a  $working_folder/transcripts.TEs.prelim.bed6 -b $TE_fragments | awk -v OFS="\t" '{arr[$4]+=$9} END {for (i in arr) {print i,arr[i]}}' | sort -k1,1n > $working_folder/exon_length

sort -k 1b,1  $working_folder/exon_length_covered  > $working_folder/exon_length_covered_sorted
sort -k 1b,1 $working_folder/exon_length  >  $working_folder/exon_length_sorted
join -1 1 -2 1 $working_folder/exon_length_covered_sorted $working_folder/exon_length_sorted |  awk -v OFS="\t" '{print $1, $2/$3}'  > $working_folder/transcripts.TEs.prelim.exon_coverage_by.TE_frags.bed

wc -l  $working_folder/transcripts.TEs.prelim.exon_coverage_by.TE_frags.bed


coverageBed -s -a  $working_folder/transcripts.TEs.prelim.bed -b $TE_fragments | awk -v OFS="\t" '{print $4, $15,$16}'>  $working_folder/transcripts.TEs.prelim.locus_coverage_by.TE_frags.bed


sort -k 4b,4  $working_folder/transcripts.TEs.prelim.bed   > $working_folder/tab1_sorted
sort -k 1b,1 $working_folder/transcripts.TEs.prelim.exon_coverage_by.TE_frags.bed   >  $working_folder/tab2_sorted
join -1 4 -2 1 $working_folder/tab1_sorted $working_folder/tab2_sorted |  awk -v OFS="\t" '{print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13}'  | sort -k 4b,4 > $working_folder/tab3_sorted
sort -k 1b,1 $working_folder/transcripts.TEs.prelim.locus_coverage_by.TE_frags.bed > $working_folder/tab4_sorted

join -1 4 -2 1 $working_folder/tab3_sorted $working_folder/tab4_sorted |  awk -v OFS="\t" '($13>=0.6 ||($13>=0.5 && $15>=0.8)){print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12}'|sortBed -i stdin > $working_folder/transcripts.TEs.bed

wc -l  $working_folder/transcripts.TEs.bed 

echo "number of expressed TEs genes' transcripts expressed" >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.TE_genes.bed >> $working_folder/linenumbers.txt
echo "number of expressed TEs genes' loci" >> $working_folder/linenumbers.txt
cat $working_folder/transcripts.TE_genes.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l >> $working_folder/linenumbers.txt


echo "number of expressed TEs transcripts expressed" >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.TEs.bed >> $working_folder/linenumbers.txt
echo "number of expressed TEs loci" >> $working_folder/linenumbers.txt
cat $working_folder/transcripts.TEs.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l >> $working_folder/linenumbers.txt



# how many TE frag transcripts overlap TE gene transcripts? 
intersectBed -s -u -a $working_folder/transcripts.TEs.bed -b $working_folder/transcripts.TE_genes.bed | wc -l 
#
intersectBed -s -u -a $working_folder/transcripts.TEs.bed -b $working_folder/transcripts.TE_genes.bed |  sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l
# 
intersectBed -split -s -u -a $working_folder/transcripts.TEs.bed -b $working_folder/transcripts.TE_genes.bed | wc -l 
#

#remove those ambiguous TE transcripts to avoid confusions 

intersectBed -s -v -a $working_folder/transcripts.TEs.bed -b $working_folder/transcripts.TE_genes.bed > $working_folder/transcripts.TE_frags.bed
wc -l $working_folder/transcripts.TE_frags.bed



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


# check exon overlap with TE fragments 
 
bed12ToBed6 -i $working_folder/transcripts.denovoNC.f1b.bed > $working_folder/transcripts.denovoNC.f1b.bed6

coverageBed -split -s -a  $working_folder/transcripts.denovoNC.f1b.bed6 -b $TE_fragments | awk -v OFS="\t" '{arr[$4]+=$8} END {for (i in arr) {print i,arr[i]}}' | sort -k1> $working_folder/exon_length_covered

coverageBed -split -s -a  $working_folder/transcripts.denovoNC.f1b.bed6 -b $TE_fragments | awk -v OFS="\t" '{arr[$4]+=$9} END {for (i in arr) {print i,arr[i]}}' | sort -k1> $working_folder/exon_length

sort -k 1b,1  $working_folder/exon_length_covered  > $working_folder/exon_length_covered_sorted
sort -k 1b,1 $working_folder/exon_length  >  $working_folder/exon_length_sorted
join -1 1 -2 1 $working_folder/exon_length_covered_sorted $working_folder/exon_length_sorted |  awk -v OFS="\t" '{print $1, $2/$3}'  > $working_folder/transcripts.denovoNC.f1b.exon_coverage_by.TE_frags.bed
wc -l $working_folder/transcripts.denovoNC.f1b.exon_coverage_by.TE_frags.bed


coverageBed -s -a  $working_folder/transcripts.denovoNC.f1b.bed -b $TE_fragments | awk -v OFS="\t" '{print $4, $15,$16}'>  $working_folder/transcripts.denovoNC.f1b.locus_coverage_by.TE_frags.bed

sort -k 4b,4  $working_folder/transcripts.denovoNC.f1b.bed  > $working_folder/tab1_sorted
sort -k 1b,1 $working_folder/transcripts.denovoNC.f1b.exon_coverage_by.TE_frags.bed   >  $working_folder/tab2_sorted
join -1 4 -2 1 $working_folder/tab1_sorted $working_folder/tab2_sorted |  awk -v OFS="\t" '{print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13}'  | sort -k 4b,4 > $working_folder/tab3_sorted
sort -k 1b,1 $working_folder/transcripts.denovoNC.f1b.locus_coverage_by.TE_frags.bed > $working_folder/tab4_sorted

join -1 4 -2 1 $working_folder/tab3_sorted $working_folder/tab4_sorted |  awk -v OFS="\t" '($13<0.6 &&($13<0.5 || $15<0.8)){print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12}'|sortBed -i stdin > $working_folder/transcripts.denovoNC.f2.bed

wc -l $working_folder/transcripts.denovoNC.f2.bed
cat $working_folder/transcripts.denovoNC.f2.bed | sort -k1,1 -k2,2g |  mergeBed -s -i stdin | wc -l


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
#/groups/nordborg/projects/cegs/alexandra/software/CPC2/CPC2-beta/bin/CPC2.py -i $annotationfolder/Araport11_protein_coding.201606.fa  -o $working_folder/Araport11_protein_coding.201606.CPC_results.txt
#cat $working_folder/Araport11_protein_coding.201606.CPC_results.txt |  awk -v OFS="\t" '($8=="noncoding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/Araport11_protein_coding.201606.CPC_results.noncoding.txt
#cat $working_folder/Araport11_protein_coding.201606.CPC_results.txt |  awk -v OFS="\t" '($8=="coding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/Araport11_protein_coding.201606.CPC_results.coding.txt

#Araport NC
#/groups/nordborg/projects/cegs/alexandra/software/CPC2/CPC2-beta/bin/CPC2.py -i $annotationfolder/Araport11_non_coding.2016016.fa  -o $working_folder/Araport11_non_coding.2016016.CPC_results.txt
#cat $working_folder/Araport11_non_coding.2016016.CPC_results.txt |  awk -v OFS="\t" '($8=="noncoding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/Araport11_non_coding.2016016.CPC_results.noncoding.txt
#cat $working_folder/Araport11_non_coding.2016016.CPC_results.txt |  awk -v OFS="\t" '($8=="coding"){split($1,a,"(");print a[1], $8} ' >  $working_folder/Araport11_non_coding.2016016.CPC_results.coding.txt




wc -l *CPC_results*


# !!!!!!!!!!!!!!!!!filter not just the potentially coding transcript, but the whole locus ( all the transcripts overlapping the transcript with high PC potential)



export CPC_loc=$working_folder/transcripts.denovoNC.f2.CPC_results.noncoding.txt
export CPC_coding_loc=$working_folder/transcripts.denovoNC.f2.CPC_results.coding.txt

export transcripts_loc=$working_folder/transcripts.denovoNC.f2.bed
export transcripts_filt_loc=$working_folder/transcripts.denovoNC.f3.bed
export transcripts_coding_loc=$working_folder/transcripts.denovoNC.f2.potentially_coding.bed

Rscript $R_folder/filtering_Rscript_2.r

#exonic overlap vs any overlap - very tiny difference - just few transcripts (~10)

intersectBed -v -s -a $working_folder/transcripts.denovoNC.f3.bed -b $working_folder/transcripts.denovoNC.f2.potentially_coding.bed  > $working_folder/transcripts.denovoNC.f3a.bed
wc -l $working_folder/transcripts.denovoNC.f3.bed
 
wc -l $working_folder/transcripts.denovoNC.f3a.bed




#remove all ambiguous NC rnas 
#remove transcripts with strand exonic overlap with de novo TE_frags and de novo pseudogenes
intersectBed -v  -split -s -a $working_folder/transcripts.denovoNC.f3a.bed -b $working_folder/transcripts.TE_frags.bed >$working_folder/transcripts.denovoNC.f3b.bed
intersectBed -v  -split -s -a $working_folder/transcripts.denovoNC.f3b.bed -b $working_folder/transcripts.TE_genes.bed >$working_folder/transcripts.denovoNC.f3b1.bed
intersectBed -v -split -s -a $working_folder/transcripts.denovoNC.f3b1.bed -b $working_folder/transcripts.f1.denovo_pseudogene.bed >$working_folder/transcripts.denovoNC.f3c.bed 

wc -l $working_folder/transcripts.denovoNC.f3a.bed 
wc -l $working_folder/transcripts.denovoNC.f3b.bed 
wc -l $working_folder/transcripts.denovoNC.f3b1.bed 
wc -l $working_folder/transcripts.denovoNC.f3c.bed 
sort -k1,1 -k2,2g $working_folder/transcripts.denovoNC.f3c.bed    |  mergeBed -s -i  stdin| wc -l 

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


#tRNAs
intersectBed -s -split -u -a $working_folder/transcripts.bed -b $tRNA | sortBed -i stdin > $working_folder/transcripts.tRNA.bed 
echo "number of de novo transcripts with exonic overlap to a tRNA" >> $working_folder/linenumbers.txt
wc -l  $working_folder/transcripts.tRNA.bed >> $working_folder/linenumbers.txt
sort -k1,1 -k2,2g $working_folder/transcripts.tRNA.bed |  mergeBed -s -i  stdin| wc -l 

# exclude ribosomal RNA and ribosomal RNA extensions 
intersectBed -s -v -a $working_folder/transcripts.denovoNC.f3c.bed -b $working_folder/transcripts.rRNA.bed $working_folder/transcripts.tRNA.bed   > $working_folder/transcripts.denovoNC.f4.bed

echo "number of de novo lncRNA transcripts" >> $working_folder/linenumbers.txt
wc -l $working_folder/transcripts.denovoNC.f4.bed >> $working_folder/linenumbers.txt
echo "number of de novo lncRNA loci " >> $working_folder/linenumbers.txt
sort -k1,1 -k2,2g $working_folder/transcripts.denovoNC.f4.bed  |  mergeBed -s -i  stdin | wc -l >> $working_folder/linenumbers.txt


###############################################################
###############################################################

#rename loci and transcripts for NC 

#rename NC loci 
#
#redefine loci 

#NC
cat $working_folder/transcripts.denovoNC.f4.bed | awk -v OFS="\t"  '{print $4,$5}' | awk -v OFS="\t"  '{split($1,s,".");   print s[1]"."s[2]}' |sort -V -k1,1 -k2,2 | uniq  > $working_folder/transcripts.denovoNC.f4.cuff_locinames.bed
wc -l $working_folder/transcripts.denovoNC.f4.cuff_locinames.bed


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

cat $working_folder/loci_table > $working_folder/transcripts.denovoNC.loci_old-new_names.bed
cat $working_folder/loci_table |awk -v OFS="\t"  '{ print $1,$2,$3,$7,$5,$6}' | sortBed -i > $working_folder/denovoNC.loci.bed 
 
cat $working_folder/loci_table  |  awk  -v OFS="\t"  'split($4,s,","){ for(i in s) {print s[i],$7,$7"."i  }}' > $working_folder/denovoNC.loci_transcripts_names_old_new.bed 


# make annotation with new names 

sort -k 4b,4  $working_folder/transcripts.denovoNC.f4.bed  > $working_folder/old_names_annotation_sorted
sort -k 1b,1 $working_folder/denovoNC.loci_transcripts_names_old_new.bed  >  $working_folder/new_old_names_sorted
join -1 1 -2 4 $working_folder/new_old_names_sorted $working_folder/old_names_annotation_sorted |  awk -v OFS="\t" '{print $2,$3,$4,$5,$6,$1,$7,$8,$9,$10,$11,$12,$13,$14}' | sort -V -k1,1 -k2,2 > $working_folder/denovoNC.transcripts.old_names_new_names.bed

cat $working_folder/denovoNC.transcripts.old_names_new_names.bed |  awk -v OFS="\t" '{print $3,$4,$5,$2,$7,$8,$9,$10,$11,$12,$13,$14}'| sortBed -i stdin > $working_folder/denovoNC.transcripts.bed

wc -l $working_folder/denovoNC.transcripts.bed 
wc -l $working_folder/denovoNC.loci.bed 


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
cat $working_folder/denovoNC.transcripts.bed | grep -w -f $working_folder/lncRNAs.antisense.loci_names.txt > $working_folder/lncRNAs.antisense.transcripts.bed
wc -l $working_folder/lncRNAs.antisense.transcripts.bed

#	4. Antisense to pseudogene  (should not overlap with antisense (to PC) lncRNAs)
intersectBed -s  -v -a $working_folder/denovoNC.loci.bed -b $working_folder/lncRNAs.antisense.loci.bed $miRNA $snRNA $snoRNA | intersectBed -S -u -a stdin -b $pseudogenes   | sortBed -i stdin > $working_folder/lncRNAs.AS_to_pseudogenes.loci.bed
cat $working_folder/lncRNAs.AS_to_pseudogenes.loci.bed | awk -v OFS="\t" '{print $4}' > $working_folder/lncRNAs.AS_to_pseudogenes.loci_names.txt
cat $working_folder/denovoNC.transcripts.bed | grep -w -f $working_folder/lncRNAs.AS_to_pseudogenes.loci_names.txt > $working_folder/lncRNAs.AS_to_pseudogenes.transcripts.bed
wc -l  $working_folder/lncRNAs.AS_to_pseudogenes.transcripts.bed

# 5.Antisense to TE genes 
export TE_genes=$annotationfolder/Araport11_transposable_element_gene.201606.bed
intersectBed -s  -v -a $working_folder/denovoNC.loci.bed -b $working_folder/lncRNAs.antisense.loci.bed $working_folder/lncRNAs.AS_to_pseudogenes.loci.bed $miRNA $snRNA $snoRNA | intersectBed -S -u -a stdin -b $TE_genes | sortBed -i stdin > $working_folder/lncRNAs.AS_to_TE.loci.bed
cat $working_folder/lncRNAs.AS_to_TE.loci.bed | awk -v OFS="\t" '{print $4}' > $working_folder/lncRNAs.AS_to_TE.loci_names.txt
cat $working_folder/denovoNC.transcripts.bed | grep -w -f $working_folder/lncRNAs.AS_to_TE.loci_names.txt > $working_folder/lncRNAs.AS_to_TE.transcripts.bed
wc -l $working_folder/lncRNAs.AS_to_TE.loci.bed


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


#remove transcripts that are too close (and downstream) to the PC genes and pseudogenes (might be trnascription read-through)

bedtools closest -a $working_folder/lncRNAs.intergenic.loci.preliminary.bed -b $PC_tair $PC_araport $pseudogenes -s -D b -iu| awk '($20<100){print $4}' | uniq > $working_folder/lincRNAs.potential_PC_runons.names.txt
wc -l $working_folder/lincRNAs.potential_PC_runons.names.txt


cat $working_folder/lncRNAs.intergenic.loci.preliminary.bed | grep -v -w -f $working_folder/lincRNAs.potential_PC_runons.names.txt >$working_folder/lncRNAs.intergenic.loci.bed 
 
wc -l $working_folder/lncRNAs.intergenic.loci.bed

 
cat $working_folder/lncRNAs.intergenic.loci.bed | awk -v OFS="\t" '{print $4}' > $working_folder/lncRNAs.intergenic.loci_names.txt
cat $working_folder/denovoNC.transcripts.bed | grep -w -f $working_folder/lncRNAs.intergenic.loci_names.txt > $working_folder/lncRNAs.intergenic.transcripts.bed
wc -l $working_folder/lncRNAs.intergenic.transcripts.bed
 
 
cat lncRNAs.intergenic.loci.bed lncRNAs.antisense.loci.bed lncRNAs.AS_to_TE.loci.bed lncRNAs.AS_to_pseudogenes.loci.bed lncRNAs.sn_snoRNA_precursors.loci.bed lncRNAs.miRNA_precursors.loci.bed | sortBed -i > lncRNAs.6types.loci.bed

cat lncRNAs.6types.loci.bed |awk -v OFS="\t" '{print $4}' | uniq | wc -l 
wc -l lncRNAs.6types.loci.bed
wc -l denovoNC.loci.bed

 
#######################
 
cat lncRNAs.antisense.transcripts.bed lncRNAs.AS_to_pseudogenes.transcripts.bed  lncRNAs.AS_to_TE.transcripts.bed lncRNAs.intergenic.transcripts.bed lncRNAs.miRNA_precursors.transcripts.bed  lncRNAs.sn_snoRNA_precursors.transcripts.bed | sortBed -i > lncRNAs.6types.transcripts.bed 
wc -l lncRNAs.6types.transcripts.bed 
wc -l denovoNC.transcripts.bed

 
#known lncRNAs (present in Araport11) 
intersectBed -s -u -a $working_folder/lncRNAs.6types.transcripts.bed  -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed $annotationfolder/Araport11_novel_transcribed_region.201606.bed > $working_folder/lncRNAs.6types.transcripts.known.bed 
intersectBed -s -u -a $working_folder/lncRNAs.6types.loci.bed -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed $annotationfolder/Araport11_novel_transcribed_region.201606.bed > $working_folder/lncRNAs.6types.loci.known.bed 
 
#novel lncRNAs (not present in Araport11)
intersectBed -s -v -a $working_folder/lncRNAs.6types.transcripts.bed  -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed $annotationfolder/Araport11_novel_transcribed_region.201606.bed > $working_folder/lncRNAs.6types.transcripts.novel.bed 
intersectBed -s -v -a $working_folder/lncRNAs.6types.loci.bed -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed $annotationfolder/Araport11_novel_transcribed_region.201606.bed> $working_folder/lncRNAs.6types.loci.novel.bed 

 
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

 known_lnc_loci=`cat  $working_folder/lncRNAs.6types.loci.known.bed | wc -l`
 known_lnc_tr=`cat $working_folder/lncRNAs.6types.transcripts.known.bed | wc -l`
 
 novel_lnc_loci=`cat  $working_folder/lncRNAs.6types.loci.novel.bed | wc -l`
 novel_lnc_tr=`cat $working_folder/lncRNAs.6types.transcripts.novel.bed | wc -l`
 
 echo -e "PC_loci\t"  "PC_tr\t"   "pseudo_loci\t" "pseudo_tr\t" "pot_cod_loci\t"  "pot_cod_tr\t" \
 "lncRNA_6types_loci\t"  "lncRNA_6types_tr\t" "novel_lnc_loci\t" "novel_lnc_tr\t" "known_lnc_loci\t" "known_lnc_tr\t" \
 "AS_loci\t"  "AS_tr\t"  "linc_loci\t"  "linc_tr\t" "AS_to_TE_loci\t"  "AS_to_TE_tr\t" \
 "rRNA_loci\t"  "rRNA_tr\t" "tRNA_loci\t" "tRNA_tr\t" "intronic_loci\t" "intronic_tr\t" \
 "AS_to_pseudo_loci\t"  "AS_to_pseudo_loci_tr\t"  "TEgene_loci\t"  "TEgene_tr\t"  "TEs_loci\t"  "TEs_tr\t" \
 "miRNA_loci\t"  "miRNA_tr\t" "snoRNA_loci\t"  "snoRNA_tr\t"   > $working_folder/$name_assembly.gene_numbers.txt

echo -e  $PC_loci"\t" $PC_tr"\t" $pseudo_loci"\t" $pseudo_tr"\t" $pot_cod_loci"\t" $pot_cod_tr"\t"\
$lncRNA_6types_loci"\t" $lncRNA_6types_tr"\t" $novel_lnc_loci"\t" $novel_lnc_tr"\t" $known_lnc_loci"\t" $known_lnc_tr"\t" \
$AS_loci"\t" $AS_tr"\t"  $linc_loci"\t" $linc_tr"\t"  $AS_to_TE_loci"\t" $AS_to_TE_tr"\t"\
$rRNA_loci"\t" $rRNA_tr"\t" $tRNA_loci"\t" $tRNA_tr"\t" $intronic_loci"\t" $intronic_tr"\t"\
$AS_to_pseudo_loci"\t" $AS_to_pseudo_tr"\t" $TEgene_loci"\t" $TEgene_tr"\t" $TEs_loci"\t" $TEs_tr"\t"\
$miRNA_loci"\t" $miRNA_tr"\t" $snoRNA_loci"\t" $snoRNA_tr"\t">> $working_folder/$name_assembly.gene_numbers.txt


 


