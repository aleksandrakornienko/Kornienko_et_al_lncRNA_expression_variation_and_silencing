#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --qos=rapid
#SBATCH --time=00:50:00
#SBATCH --array=1-461
#SBATCH --output=/groups/nordborg/user/aleksandra.kornienko/analyses/logs/array_job_slurm_%A_%a.out

export numberofline=$SLURM_ARRAY_TASK_ID 

ml subread/2.0.0-foss-2018b

export list_samples=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/000_data_preparation/2021/1001G_461accessions_used_for_annotation.txt
export accession=`sed -n $numberofline,"$numberofline"p $list_samples | awk '{print $1}'`
export sample=1001G.$accession
echo $sample
export bamfolder=/groups/nordborg/projects/cegs/alexandra/the1001genome_RNAseq/rerun2020/BAM
export bam=$bamfolder/$accession/$accession.Aligned.sortedByCoord.out.bam 

 
#de novo 
 
export featureCountsfolder=/groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/02_expression/featureCounts/20211013_annotation
export saf=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/20211013_annotation.transcripts.saf 

export out=$featureCountsfolder/$sample

featureCounts  -T 4  -F SAF -O -s 2  -t exon -g gene_id -a $saf -o $out.txt $bam

# the program does not count multimapping reads! 
 
#convert to simpler columns
sum=`cat $out.txt| awk -v OFS="\t" '{ sum += $7} END {print sum}'`
cat $out.txt| grep CUFF| awk -v OFS="\t" -v sum="$sum" ' {print $1,$6,$7,$7*1000000*1000/($6*sum)}' > $out.counts_tpm.bed


#Araport
export featureCountsfolder=/groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/02_expression/featureCounts/Araport11
 
export saf=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11.PC_NC_Pseud_TE.longer200nt.saf

export out=$featureCountsfolder/$sample

featureCounts  -T 4  -F SAF -O -s 2  -t exon -g gene_id -a $saf -o $out.txt $bam

# the program does not count multimapping reads! 
 
#convert to simpler columns
sum=`cat $out.txt| awk -v OFS="\t" '{ sum += $7} END {print sum}'`
cat $out.txt| grep AT| awk -v OFS="\t" -v sum="$sum" ' {print $1,$6,$7,$7*1000000*1000/($6*sum)}' > $out.counts_tpm.bed


 

