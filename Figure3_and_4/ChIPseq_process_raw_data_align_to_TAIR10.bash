#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=20G
#SBATCH --time=6:00:00
#SBATCH --array=1-178
#SBATCH --output=/groups/nordborg/user/aleksandra.kornienko/analyses/logs/array_job_slurm_%A_%a.out


export numberofline=$SLURM_ARRAY_TASK_ID

#load necessary modules
module load sra-toolkit/2.9.6-1-centos_linux64
module load bwa/0.7.17-foss-2018b
module load picard/2.18.27-java-1.8
module load trim_galore/0.6.2-foss-2018b-python-2.7.15
module load samtools/0.1.20-foss-2018b
module load  cutadapt/1.18-foss-2018b-python-3.6.6
module load r/3.5.1-foss-2018b


#set working directory to save space in the script below + make your script more resistant to folder relocations - will need to change fewer things if folders move or file system changes
working_folder=/groups/nordborg/projects/cegs/alexandra/ChIP-seq_2021

#the list of sample names and raw data locations

list=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/All_ChIP_data_samples_name_raw_data_location.txt

#get the sample name from the list
export name=`sed -n $numberofline,"$numberofline"p $list | awk '{print $1}'`
export common_folder=/groups/nordborg/user/aleksandra.kornienko/rawdata/chipseq

export folder1=`sed -n $numberofline,"$numberofline"p $list | awk '{print $4}'`
export folder2=`sed -n $numberofline,"$numberofline"p $list | awk '{print $6}'`
export folder3=`sed -n $numberofline,"$numberofline"p $list | awk '{print $8}'`
export folder4=`sed -n $numberofline,"$numberofline"p $list | awk '{print $10}'`

export unmapped_bam1=`sed -n $numberofline,"$numberofline"p $list | awk '{print $5}'`
export unmapped_bam2=`sed -n $numberofline,"$numberofline"p $list | awk '{print $7}'`
export unmapped_bam3=`sed -n $numberofline,"$numberofline"p $list | awk '{print $9}'`
export unmapped_bam4=`sed -n $numberofline,"$numberofline"p $list | awk '{print $11}'`

echo $name 
echo $folder1/$unmapped_bam1
echo $unmapped_bam2
echo $unmapped_bam3
echo $unmapped_bam4

#java -jar ${EBROOTPICARD}/picard.jar SamToFastq I=$common_folder/$folder1/$unmapped_bam1 FASTQ=$working_folder/fastq/$name.1.fastq 

#java -jar ${EBROOTPICARD}/picard.jar SamToFastq I=$common_folder/$folder2/$unmapped_bam2 FASTQ=$working_folder/fastq/$name.2.fastq  

java -jar ${EBROOTPICARD}/picard.jar SamToFastq I=$common_folder/$folder3/$unmapped_bam3 FASTQ=$working_folder/fastq/$name.3.fastq  

java -jar ${EBROOTPICARD}/picard.jar SamToFastq I=$common_folder/$folder4/$unmapped_bam4 FASTQ=$working_folder/fastq/$name.4.fastq  


# raw reads in (archived .gz) fastq format:
export r1=$working_folder/fastq/$name.1.fastq
export r2=$working_folder/fastq/$name.2.fastq
export r3=$working_folder/fastq/$name.3.fastq
export r4=$working_folder/fastq/$name.4.fastq


#align with star 
ml star/2.7.1a-foss-2018b
 
mkdir $working_folder/BAM/$name
export stargenomedir=/groups/nordborg/user/aleksandra.kornienko/analyses/STAR/STAR2.7.1a/STAR_TAIR10_noSJDB_sizeaccounted
export out=$working_folder/BAM/$name/$name.



STAR --readFilesIn $r1,$r2,$r3,$r4 \
--runThreadN 4 \
--alignIntronMax 5 \
--genomeDir $stargenomedir \
--outFilterMismatchNmax 10 \
--outFilterMultimapNmax 1 \
--outReadsUnmapped None \
 --alignEndsType EndToEnd \
--outSAMtype BAM SortedByCoordinate  \
--runMode alignReads \
--twopassMode Basic \
--outFileNamePrefix $out



## create indexed bam file
module load samtools/1.9-foss-2018b
samtools index $working_folder/BAM/$name/$name.Aligned.sortedByCoord.out.bam  


#check BAM quality 
bam=$working_folder/BAM/$name/$name.Aligned.sortedByCoord.out.bam  

#characterize the bam file
module load  rseqc/2.6.5-foss-2018b-python-2.7.15
bam_stat.py -i  $bam  >$working_folder/BAM/$name/$name.bam_stat
 
#check duplicate number 
read_duplication.py -i $bam  -o $working_folder/BAM/$name/$name
 
module unload  rseqc/2.6.5-foss-2018b-python-2.7.15

#make bigwig 
module load deeptools/3.1.2-foss-2018b-python-2.7.15


bamCoverage -b $bam  -o $working_folder/BAM/$name/$name.bw --binSize=10
