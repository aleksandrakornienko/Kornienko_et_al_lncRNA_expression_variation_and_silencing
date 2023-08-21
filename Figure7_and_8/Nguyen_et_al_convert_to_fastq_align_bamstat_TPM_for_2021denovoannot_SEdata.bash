#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30G
#SBATCH --qos=short
#SBATCH --time=08:00:00
#SBATCH --array=1-71
#SBATCH --output=/groups/nordborg/user/aleksandra.kornienko/analyses/logs/array_job_slurm_%A_%a.out
 
# wc -l InfoDataSingleEnd.txt
# 143 InfoDataSingleEnd.txt
# wc -l /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/InfoDataPairEnd.clean.txt
#60 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/InfoDataPairEnd.clean.txt


export numberofline=$SLURM_ARRAY_TASK_ID 

module load sra-toolkit/2.9.6-1-centos_linux64
module load bwa/0.7.17-foss-2018b
module load picard/2.18.27-java-1.8
module load star/2.7.1a-foss-2018b

#cat /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/InfoDataSingleEnd.txt |grep -v replicate| awk -v OFS="\t" '{print $1".rep"$2,$3}'>/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/InfoDataSingleEnd.clean.txt
#cat /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/InfoDataSingleEnd.clean.txt | awk -v OFS="\t" '{print $1}'|sort|uniq>/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/SE_samples.txt
#cat /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/InfoDataPairEnd.clean.txt |  awk -v OFS="\t" '{print $1".rep"$2,$3}'|awk -v OFS="\t" '{print $1}'|sort|uniq>/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/PE_samples.txt
#wc -l /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/SE_samples.txt
#71 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/SE_samples.txt




export working_folder=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data
export data_folder=/scratch-cbe/users/vu.nguyen/RNA_seq/For_Aleksandra/singleEnd_bam
#export data_folder=/scratch-cbe/users/vu.nguyen/RNA_seq/For_Aleksandra/pairEnd_bam
export sample_table=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/InfoDataSingleEnd.clean.txt
export list_samples=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/SE_samples.txt
#export sample_table=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/Vu_s_data/InfoDataPairEnd.clean.txt


export name_accession=`cat $list_samples | sed -n $numberofline,"$numberofline"p  | awk '{print $1}'`
export unmapped_BAM1=`cat $sample_table | grep -w $name_accession| head -1 | awk '{print $2}'`
export unmapped_BAM2=`cat $sample_table | grep -w $name_accession| head -2 |tail -1| awk '{print $2}'`
echo $name_accession
echo $unmapped_BAM1
echo $unmapped_BAM2

export unmapped_bamfile1=$data_folder/$unmapped_BAM1
export unmapped_bamfile2=$data_folder/$unmapped_BAM2

java -jar ${EBROOTPICARD}/picard.jar SamToFastq I=$unmapped_bamfile1      FASTQ=$working_folder/fastq/$name_accession.1.fastq  

java -jar ${EBROOTPICARD}/picard.jar SamToFastq I=$unmapped_bamfile2      FASTQ=$working_folder/fastq/$name_accession.2.fastq  
#align

export R1_lane1=$working_folder/fastq/$name_accession.1.fastq
export R1_lane2=$working_folder/fastq/$name_accession.2.fastq


mkdir $working_folder/BAM/$name_accession
export stargenomedir=/groups/nordborg/user/aleksandra.kornienko/analyses/STAR/STAR2.7.1a/STAR_TAIR10_noSJDB_sizeaccounted
export out=$working_folder/BAM/$name_accession/$name_accession.

#95% of Arabidopsis introns are smaller than 3.5kb
#!!!!!!!!!!!!!!!!!!!!!!!! changed the multimaps to 10!

STAR --readFilesIn $R1_lane1,$R1_lane2  \
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
samtools index $working_folder/BAM/$name_accession/$name_accession.Aligned.sortedByCoord.out.bam  

#check BAM quality 
bam=$working_folder/BAM/$name_accession/$name_accession.Aligned.sortedByCoord.out.bam  

#characterize the merged bam file
module load  rseqc/2.6.5-foss-2018b-python-2.7.15
bam_stat.py -i  $bam  >$working_folder/BAM/$name_accession/$name_accession.bam_stat
#check strand specificity
infer_experiment.py -i $bam   -r /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.bed > $working_folder/BAM/$name_accession/$name_accession.strandspec.txt
#check insert size
inner_distance.py -i $bam -o $working_folder/BAM/$name_accession/$name_accession.innerdist.txt -r /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.bed

#check ribo and chloroplast contamination 
module load bedtools/2.27.1-foss-2018b
bamToBed -i $bam   | coverageBed -a /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Athaliana_full_chromosomes_andRibo.bed.txt -b stdin > $working_folder/BAM/$name_accession/$name_accession.Rib_ChrC_contamination.txt
 
module unload  rseqc/2.6.5-foss-2018b-python-2.7.15

#make bigwig 
 module load deeptools/3.1.2-foss-2018b-python-2.7.15


bamCoverage -b $bam  --filterRNAstrand forward -o $working_folder/BAM/$name_accession/$name_accession.F.bw --binSize=25
bamCoverage -b $bam  --filterRNAstrand reverse -o $working_folder/BAM/$name_accession/$name_accession.R.bw --binSize=25
 
#calculate TPM by Stringtie in 2021 de novo annotation 

ml stringtie/2.1.5-gcccore-7.3.0 

# de novo annotation  
#Z:\01_POSTDOC\03_Projects\2018_lncRNA_variation_paper\01_lncRNA_identification\01_Identification_pipeline\20210513_cuffmerge_filtering_and_classification_automated_pipeline_cortijo_eracaps_1001g_1001gnew_2stepmerge.bash
export gtf=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20210513_annotation/20210513_annotation.transcripts.gtf
export stringtieExpfolder=$working_folder/stringtie_expression/20210513_annotation

export out=$stringtieExpfolder/$name_accession.STRG.e.txt
export gene_abund=$stringtieExpfolder/$name_accession.gene_abund.txt
stringtie -e -B -G $gtf -A $gene_abund  -p 4 $bam -o  $out
cat $out |  awk -v OFS="\t" '($3=="transcript") {print $1,$4,$5,$7,$12,$14,$16,$18 }' | sed -e 's/\"//g'|sed -e 's/;//g'> $stringtieExpfolder/$name_accession.transcript_abund.txt
cat $stringtieExpfolder/$name_accession.transcript_abund.txt| sort -u  -k5,5  | awk -v OFS="\t"   '{print $5,$8}' >$stringtieExpfolder/$name_accession.transcripts.txt
cat $stringtieExpfolder/$name_accession.gene_abund.txt |grep -v Reference|sort -u  -k1,1  | awk -v OFS="\t"   '{print $1,$9}' >$stringtieExpfolder/$name_accession.genes.txt


#Araport   Z:\01_POSTDOC\03_Projects\2018_lncRNA_variation_paper\02_expression_and_variation\prepare_Araport_and_TAIR10_for_stringtie_expression.bash
export gtf=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11.PC_NC_Pseud_TE.longer200nt.gtf
export stringtieExpfolder=$working_folder/stringtie_expression/Araport11

export out=$stringtieExpfolder/$name_accession.STRG.e.txt
export gene_abund=$stringtieExpfolder/$name_accession.gene_abund.txt
stringtie -e -B -G $gtf -A $gene_abund  -p 4 $bam -o  $out
cat $out |  awk -v OFS="\t" '($3=="transcript") {print $1,$4,$5,$7,$12,$14,$16,$18 }' | sed -e 's/\"//g'|sed -e 's/;//g'> $stringtieExpfolder/$name_accession.transcript_abund.txt
cat $stringtieExpfolder/$name_accession.transcript_abund.txt| sort -u  -k5,5  | awk -v OFS="\t"   '{print $5,$8}' >$stringtieExpfolder/$name_accession.transcripts.txt
cat $stringtieExpfolder/$name_accession.gene_abund.txt |grep -v Reference|sort -u  -k1,1  | awk -v OFS="\t"   '{print $1,$9}' >$stringtieExpfolder/$name_accession.genes.txt


 
 
