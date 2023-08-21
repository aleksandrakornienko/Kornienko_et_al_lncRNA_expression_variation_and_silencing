#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --qos=medium
#SBATCH --time=15:50:00
#SBATCH --array=1-15
#SBATCH --output=/groups/nordborg/user/aleksandra.kornienko/analyses/logs/array_job_slurm_%A_%a.out

export numberofline=$SLURM_ARRAY_TASK_ID 

samples=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/Public_data/He_et_al_Methylation_free_Arabidopsis/He_et_al_srr_sample_name_RNAseq.txt
wc -l $samples 
#15 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/Public_data/He_et_al_Methylation_free_Arabidopsis/He_et_al_srr_sample_name_RNAseq.txt

working_folder=/groups/nordborg/projects/cegs/alexandra/Public_data/He_et_al

export sample=`cat $samples | sed -n $numberofline,"$numberofline"p  | awk '{print $2}'`
export srr=`cat $samples | sed -n $numberofline,"$numberofline"p| awk '{print $1}'`

fastq1=$srr
fastq1+=_1.fastq.gz

fastq2=$srr
fastq2+=_2.fastq.gz

module load star/2.7.1a-foss-2018b

#align
cd $working_folder
mkdir BAM/$sample
export stargenomedir=/groups/nordborg/user/aleksandra.kornienko/analyses/STAR/STAR2.7.1a/STAR_TAIR10_noSJDB_sizeaccounted
export out=BAM/$sample/$sample.

#95% of Arabidopsis introns are smaller than 3.5kb

STAR --readFilesIn  $fastq1 $fastq2 --readFilesCommand zcat  \
--runThreadN 4 \
--alignIntronMax 6000 \
--alignMatesGapMax 6000 \
--genomeDir $stargenomedir \
--outFilterIntronMotifs RemoveNoncanonical \
--outFilterMismatchNoverReadLmax 0.1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.3 \
--outFilterMultimapNmax 10 \
--outReadsUnmapped None \
--alignSJoverhangMin 8 \
--outSAMattributes NH HI AS nM NM MD jM jI XS \
--outSAMtype BAM SortedByCoordinate  \
--runMode alignReads \
--twopassMode Basic \
--outFileNamePrefix $out


## create indexed bam file
module load samtools/1.9-foss-2018b
samtools index $working_folder/BAM/$sample/$sample.Aligned.sortedByCoord.out.bam  

#check BAM quality 
bam=$working_folder/BAM/$sample/$sample.Aligned.sortedByCoord.out.bam  

#check ribo and chloroplast contamination 
module load bedtools/2.27.1-foss-2018b

bamToBed -i $bam   | coverageBed -a /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Athaliana_full_chromosomes_andRibo.bed.txt -b stdin > $working_folder/BAM/$sample/$sample.Rib_ChrC_contamination.txt
 
module unload  rseqc/2.6.5-foss-2018b-python-2.7.15

#make bigwig 
module load deeptools/3.1.2-foss-2018b-python-2.7.15


bamCoverage -b $bam  --filterRNAstrand forward -o $working_folder/BAM/$sample/$sample.F.bw --binSize=25
bamCoverage -b $bam  --filterRNAstrand reverse -o $working_folder/BAM/$sample/$sample.R.bw --binSize=25
 



ml subread/2.0.0-foss-2018b


#de novo 
 
export featureCountsfolder=/groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/02_expression/featureCounts/20211013_annotation
export saf=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/20211013_annotation.transcripts.saf 
export out=$featureCountsfolder/$sample

featureCounts -p  -T 4  -F SAF -O -s 2  -t exon -g gene_id -a $saf -o $out.txt $bam

# the program does not count multimapping reads! 
#convert to simpler columns
sum=`cat $out.txt| awk -v OFS="\t" '{ sum += $7} END {print sum}'`
cat $out.txt| grep CUFF| awk -v OFS="\t" -v sum="$sum" ' {print $1,$6,$7,$7*1000000*1000/($6*sum)}' > $out.counts_tpm.bed



#Araport+ TE fragments 
export featureCountsfolder=/groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/02_expression/featureCounts/Araport11_with_TE_elements
export saf=/groups/nordborg/user/aleksandra.kornienko/analyses/Araport11.PC_NC_Pseud_TE.longer200nt_and_TEelements.saf
export out=$featureCountsfolder/$sample

featureCounts -p  -T 4  -F SAF -O -s 2  -t exon -g gene_id -a $saf -o $out.txt $bam
# the program does not count multimapping reads! 
#convert to simpler columns
sum=`cat $out.txt| awk -v OFS="\t" '{ sum += $7} END {print sum}'`
cat $out.txt| grep AT| awk -v OFS="\t" -v sum="$sum" ' {print $1,$6,$7,$7*1000000*1000/($6*sum)}' > $out.counts_tpm.bed

 
# combine TPMs 




cd /groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/02_expression/featureCounts/20211013_annotation


samples=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/Public_data/He_et_al_Methylation_free_Arabidopsis/He_et_al_sample_names_RNAseq.txt

echo "gene" > tpm_table
cat RNA_WT_ML.rep2.counts_tpm.bed |  awk -v OFS="\t" '{print $1}' >> tpm_table
while read sample
do 
echo $sample
echo $sample > tpm
cat $sample.counts_tpm.bed | awk -v OFS="\t" '{print $4}' >> tpm
paste tpm_table tpm > inter
cat inter > tpm_table
done <$samples

cat tpm_table > denovo_Oct2021.TPMs.genes.He_et_al.bed


cp /groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/02_expression/featureCounts/20211013_annotation/denovo_Oct2021.TPMs.genes.He_et_al.bed /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation 


cd /groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/02_expression/featureCounts/Araport11_with_TE_elements


samples=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/Public_data/He_et_al_Methylation_free_Arabidopsis/He_et_al_sample_names_RNAseq.txt

echo "gene" > tpm_table
cat RNA_WT_ML.rep2.counts_tpm.bed |  awk -v OFS="\t" '{print $1}' >> tpm_table
while read sample
do 
echo $sample
echo $sample > tpm
cat $sample.counts_tpm.bed | awk -v OFS="\t" '{print $4}' >> tpm
paste tpm_table tpm > inter
cat inter > tpm_table
done <$samples

cat tpm_table > Araport11_with_TE_elements.TPMs.genes.He_et_al.bed


cp /groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/02_expression/featureCounts/Araport11_with_TE_elements/Araport11_with_TE_elements.TPMs.genes.He_et_al.bed /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation 

