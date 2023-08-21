#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=20G
#SBATCH --time=6:00:00
#SBATCH --array=1-93
#SBATCH --output=/groups/nordborg/pub/Mirjam/logs/array_job_slurm_%A_%a.out



export numberofline=$SLURM_ARRAY_TASK_ID

module load  rseqc/2.6.5-foss-2018b-python-2.7.15

ml bedtools/2.27.1-foss-2018b

cd  /groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/methylation/additional_samples
#mkdir CG
#mkdir CHG
#mkdir CHH

methcalls_folder=/groups/nordborg/projects/nordborg_common/datasets/alexandra_Ath_bisulfiteseq/methylseq/methylpy
export samples=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/AK_methyl_samplenames.txt



cd /groups/nordborg/projects/nordborg_common/datasets/alexandra_Ath_bisulfiteseq/methylseq/bismark_deduplicated


for numberofline in {1..93}
do 
export sampleN=`sed -n $numberofline,"$numberofline"p $samples | awk '{print $1}'`
export sample_name=`sed -n $numberofline,"$numberofline"p $samples | awk '{print $2}'`
echo $sampleN
echo $sample_name
 bam_stat.py -i $sampleN.R1_val_1_bismark_bt2_pe.deduplicated.bam > /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/bamstat/$sample_name.bamstat
done

export samples=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/AK_methyl_samplenames.txt
cd /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/bamstat/
echo "sample" "read_pair_N"> readN.bisulfiteseq_samples.txt
for numberofline in {1..93}
do 
export sampleN=`sed -n $numberofline,"$numberofline"p $samples | awk '{print $1}'`
export sample_name=`sed -n $numberofline,"$numberofline"p $samples | awk '{print $2}'`
echo $sampleN
echo $sample_name
 
readN=`cat $sample_name.bamstat | grep Read-1 | awk '{print $2}'`
echo $sample_name $readN >> readN.bisulfiteseq_samples.txt
done







