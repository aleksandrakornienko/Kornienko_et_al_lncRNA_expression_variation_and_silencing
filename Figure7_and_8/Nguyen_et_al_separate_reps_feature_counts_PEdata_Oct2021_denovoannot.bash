#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30G
#SBATCH --qos=medium
#SBATCH --time=16:00:00
#SBATCH --array=1-60
#SBATCH --output=/groups/nordborg/user/aleksandra.kornienko/analyses/logs/array_job_slurm_%A_%a.out
 
#  wc -l /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/PE_samples_basename.txt
# 8 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/PE_samples_basename.txt
#wc -l /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/SE_sample_basename.txt
# 18 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/SE_sample_basename.txt

# wc -l /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/SE_samples.txt
#71 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/SE_samples.txt

# wc -l /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/PE_samples.txt
#60 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/PE_samples.txt


export numberofline=$SLURM_ARRAY_TASK_ID 



export working_folder=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data
#export list_samples=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/SE_samples.txt
export list_samples=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/PE_samples.txt



export name_accession=`cat $list_samples | sed -n $numberofline,"$numberofline"p  | awk '{print $1}'`
echo $name_accession

 
 
ml subread/2.0.0-foss-2018b

bam=$working_folder/BAM/$name_accession/$name_accession.Aligned.sortedByCoord.out.bam  
export sample=Vu.$name_accession
echo $sample

#de novo 
export featureCountsfolder=/groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/02_expression/featureCounts/20211013_annotation
export saf=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/20211013_annotation.transcripts.saf 
export out=$featureCountsfolder/$sample
featureCounts  -T 4  -F SAF -O -s 2 -p -t exon -g gene_id -a $saf -o $out.txt $bam

# the program does not count multimapping reads! 
 #convert to simpler columns
sum=`cat $out.txt| awk -v OFS="\t" '{ sum += $7} END {print sum}'`
cat $out.txt| grep CUFF| awk -v OFS="\t" -v sum="$sum" ' {print $1,$6,$7,$7*1000000*1000/($6*sum)}' > $out.counts_tpm.bed

#Araport
export featureCountsfolder=/groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/02_expression/featureCounts/Araport11
export saf=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11.PC_NC_Pseud_TE.longer200nt.saf
export out=$featureCountsfolder/$sample
featureCounts  -T 4  -F SAF -O -s 2  -p  -t exon -g gene_id -a $saf -o $out.txt $bam

# the program does not count multimapping reads! 
 
#convert to simpler columns
sum=`cat $out.txt| awk -v OFS="\t" '{ sum += $7} END {print sum}'`
cat $out.txt| grep AT| awk -v OFS="\t" -v sum="$sum" ' {print $1,$6,$7,$7*1000000*1000/($6*sum)}' > $out.counts_tpm.bed


