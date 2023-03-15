#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30G
#SBATCH --time=08:00:00
#SBATCH --array=1-{N}
#SBATCH --output=/groups/nordborg/user/aleksandra.kornienko/analyses/logs/array_job_slurm_%A_%a.out

#example script for RNA-seq data alignment

export numberofline=$SLURM_ARRAY_TASK_ID 

module load star/2.7.1a-foss-2018b

export samples=$working_folder/list_of_sample_names.txt 
export name_sample=`sed -n $numberofline,"$numberofline"p $samples | awk '{print $1}'`

#align
export Read1_lane1=$fastq_folder/$name_sample.1.R1.fastq.gz
export Read2_lane1=$fastq_folder/$name_sample.1.R2.fastq.gz
export Read1_lane2=$fastq_folder/$name_sample.2.R1.fastq.gz
export Read2_lane2=$fastq_folder/$name_sample.2.R2.fastq.gz


mkdir $working_folder/BAM/$name_sample
export stargenomedir=/groups/nordborg/user/aleksandra.kornienko/analyses/STAR/STAR2.7.1a/$accession.PacBio.STAR_genome
export out=$working_folder/BAM/$name_sample/$name_sample.

#95% of Arabidopsis introns are smaller than 3.5kb
# align with STAR

STAR --readFilesIn $Read1_lane1,$Read1_lane2 $Read2_lane1,$Read2_lane2  --readFilesCommand zcat \
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
samtools index $working_folder/BAM/$name_sample/$name_sample.Aligned.sortedByCoord.out.bam  

#check BAM quality 
bam=$working_folder/BAM/$name_sample/$name_sample.Aligned.sortedByCoord.out.bam  

#characterize the merged bam file
module load  rseqc/2.6.5-foss-2018b-python-2.7.15
bam_stat.py -i  $bam  >$working_folder/BAM/$name_sample/$name_sample.bam_stat
#check strand specificity
infer_experiment.py -i $bam   -r $annotationfolder/Araport11_protein_coding.201606.bed > $working_folder/BAM/$name_sample/$name_sample.strandspec.txt
#check insert size
inner_distance.py -i $bam -o $working_folder/BAM/$name_sample/$name_sample.innerdist.txt -r $annotationfolder/Araport11_protein_coding.201606.bed

#check ribo and chloroplast contamination 
module load bedtools/2.27.1-foss-2018b
bamToBed -i $bam   | coverageBed -a /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Athaliana_full_chromosomes_andRibo.bed.txt -b stdin > $working_folder/BAM/$name_sample/$name_sample.Rib_ChrC_contamination.txt
module unload  rseqc/2.6.5-foss-2018b-python-2.7.15

#make bigwig 
module load deeptools/3.1.2-foss-2018b-python-2.7.15
bamCoverage -b $bam  --filterRNAstrand forward -o $working_folder/BAM/$name_sample/$name_sample.F.bw --binSize=25
bamCoverage -b $bam  --filterRNAstrand reverse -o $working_folder/BAM/$name_sample/$name_sample.R.bw --binSize=25
 

