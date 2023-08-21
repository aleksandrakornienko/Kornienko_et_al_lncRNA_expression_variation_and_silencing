
cd /groups/nordborg/projects/cegs/alexandra/Public_data/Bagghy

wget -c -P . https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR11780903/SRR11780903.1

wget -c -P . https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR11780904/SRR11780904.1
wget -c -P . https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR11780905/SRR11780905.1
wget -c -P . https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR11780906/SRR11780906.1

wget -c -P . https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR11780907/SRR11780907.1

wget -c -P . https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR11780908/SRR11780908.1




module load sra-toolkit/2.9.6-1-centos_linux64
module load bwa/0.7.17-foss-2018b
module load picard/2.18.27-java-1.8
module load star/2.7.1a-foss-2018b

#convert sra to fastq.gz
fastq-dump.2.9.6 --outdir fastq --gzip --skip-technical  --readids --dumpbase --clip SRR11780903.1
fastq-dump.2.9.6 --outdir fastq --gzip --skip-technical  --readids --dumpbase --clip SRR11780904.1
fastq-dump.2.9.6 --outdir fastq --gzip --skip-technical  --readids --dumpbase --clip SRR11780905.1
fastq-dump.2.9.6 --outdir fastq --gzip --skip-technical  --readids --dumpbase --clip SRR11780906.1
fastq-dump.2.9.6 --outdir fastq --gzip --skip-technical  --readids --dumpbase --clip SRR11780907.1
fastq-dump.2.9.6 --outdir fastq --gzip --skip-technical  --readids --dumpbase --clip SRR11780908.1


#align

fastq=fastq/SRR11780903.1.fastq.gz

sra	sample
SRR11780903	Col_WT.rep1
SRR11780904	Col_WT.rep2
SRR11780905	Col_WT.rep3
SRR11780906	Col_ddm1.rep1
SRR11780907	Col_ddm1.rep2
SRR11780908	Col_ddm1.rep3


name_accession=Col_WT.rep1
mkdir BAM/$name_accession
export stargenomedir=/groups/nordborg/user/aleksandra.kornienko/analyses/STAR/STAR2.7.1a/STAR_TAIR10_noSJDB_sizeaccounted
export out=BAM/$name_accession/$name_accession.

#95% of Arabidopsis introns are smaller than 3.5kb

STAR --readFilesIn  $fastq --readFilesCommand zcat  \
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


working_folder=/groups/nordborg/projects/cegs/alexandra/Public_data/Bagghy

## create indexed bam file
module load samtools/1.9-foss-2018b
samtools index $working_folder/BAM/$name_accession/$name_accession.Aligned.sortedByCoord.out.bam  

#check BAM quality 
bam=$working_folder/BAM/$name_accession/$name_accession.Aligned.sortedByCoord.out.bam  

#check ribo and chloroplast contamination 
module load bedtools/2.27.1-foss-2018b

srun --mem-per-cpu=30G bamToBed -i $bam   | coverageBed -a /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Athaliana_full_chromosomes_andRibo.bed.txt -b stdin > $working_folder/BAM/$name_accession/$name_accession.Rib_ChrC_contamination.txt
 
module unload  rseqc/2.6.5-foss-2018b-python-2.7.15

#make bigwig 
module load deeptools/3.1.2-foss-2018b-python-2.7.15


srun --mem-per-cpu=30G bamCoverage -b $bam  --filterRNAstrand forward -o $working_folder/BAM/$name_accession/$name_accession.F.bw --binSize=25

srun --mem-per-cpu=30G bamCoverage -b $bam  --filterRNAstrand reverse -o $working_folder/BAM/$name_accession/$name_accession.R.bw --binSize=25
 


 
 
module load sra-toolkit/2.9.6-1-centos_linux64
module load bwa/0.7.17-foss-2018b
module load picard/2.18.27-java-1.8
module load star/2.7.1a-foss-2018b

cd /groups/nordborg/projects/cegs/alexandra/Public_data/Bagghy

sra=SRR11780904
name_accession=Col_WT.rep2
#convert sra to fastq.gz
fastq-dump.2.9.6 --outdir fastq --gzip --skip-technical  --readids --dumpbase --clip $sra.1
fastq=fastq/$sra.1.fastq.gz

mkdir BAM/$name_accession
export stargenomedir=/groups/nordborg/user/aleksandra.kornienko/analyses/STAR/STAR2.7.1a/STAR_TAIR10_noSJDB_sizeaccounted
export out=BAM/$name_accession/$name_accession.

#95% of Arabidopsis introns are smaller than 3.5kb

srun --mem-per-cpu=30G STAR --readFilesIn  $fastq --readFilesCommand zcat  \
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


working_folder=/groups/nordborg/projects/cegs/alexandra/Public_data/Bagghy

## create indexed bam file
module load samtools/1.9-foss-2018b

srun --mem-per-cpu=30G samtools index $working_folder/BAM/$name_accession/$name_accession.Aligned.sortedByCoord.out.bam  

#check BAM quality 
bam=$working_folder/BAM/$name_accession/$name_accession.Aligned.sortedByCoord.out.bam  

#check ribo and chloroplast contamination 
module load bedtools/2.27.1-foss-2018b

srun --mem-per-cpu=30G bamToBed -i $bam   | coverageBed -a /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Athaliana_full_chromosomes_andRibo.bed.txt -b stdin > $working_folder/BAM/$name_accession/$name_accession.Rib_ChrC_contamination.txt
 
module unload  rseqc/2.6.5-foss-2018b-python-2.7.15

#make bigwig 
module load deeptools/3.1.2-foss-2018b-python-2.7.15


srun --mem-per-cpu=30G bamCoverage -b $bam  --filterRNAstrand forward -o $working_folder/BAM/$name_accession/$name_accession.F.bw --binSize=25

srun --mem-per-cpu=30G bamCoverage -b $bam  --filterRNAstrand reverse -o $working_folder/BAM/$name_accession/$name_accession.R.bw --binSize=25
 

 
 
 
sra=SRR11780905
name_accession=Col_WT.rep3


sra=SRR11780906
name_accession=Col_ddm1.rep1


sra=SRR11780907
name_accession=Col_ddm1.rep2 

sra=SRR11780908
name_accession=Col_ddm1.rep3



sra	sample
SRR11780903	Col_WT.rep1
SRR11780904	Col_WT.rep2
SRR11780905	Col_WT.rep3
SRR11780906	Col_ddm1.rep1
SRR11780907	Col_ddm1.rep2
SRR11780908	Col_ddm1.rep3



 
module load sra-toolkit/2.9.6-1-centos_linux64
module load bwa/0.7.17-foss-2018b
module load picard/2.18.27-java-1.8
module load star/2.7.1a-foss-2018b

cd /groups/nordborg/projects/cegs/alexandra/Public_data/Bagghy

#convert sra to fastq.gz
srun --mem-per-cpu=30G fastq-dump.2.9.6 --outdir fastq --gzip --skip-technical  --readids --dumpbase --clip $sra.1
fastq=fastq/$sra.1.fastq.gz

#align to TAIR10
mkdir BAM/$name_accession
export stargenomedir=/groups/nordborg/user/aleksandra.kornienko/analyses/STAR/STAR2.7.1a/STAR_TAIR10_noSJDB_sizeaccounted
export out=BAM/$name_accession/$name_accession.

#95% of Arabidopsis introns are smaller than 3.5kb

srun --mem-per-cpu=30G STAR --readFilesIn  $fastq --readFilesCommand zcat  \
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


working_folder=/groups/nordborg/projects/cegs/alexandra/Public_data/Bagghy

## create indexed bam file
module load samtools/1.9-foss-2018b

srun --mem-per-cpu=30G samtools index $working_folder/BAM/$name_accession/$name_accession.Aligned.sortedByCoord.out.bam  

#check BAM quality 
bam=$working_folder/BAM/$name_accession/$name_accession.Aligned.sortedByCoord.out.bam  

#check ribo and chloroplast contamination 
module load bedtools/2.27.1-foss-2018b

srun --mem-per-cpu=30G bamToBed -i $bam   | coverageBed -a /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Athaliana_full_chromosomes_andRibo.bed.txt -b stdin > $working_folder/BAM/$name_accession/$name_accession.Rib_ChrC_contamination.txt
 
module unload  rseqc/2.6.5-foss-2018b-python-2.7.15

#make bigwig 
module load deeptools/3.1.2-foss-2018b-python-2.7.15


srun --mem-per-cpu=30G bamCoverage -b $bam  --filterRNAstrand forward -o $working_folder/BAM/$name_accession/$name_accession.F.bw --binSize=25

srun --mem-per-cpu=30G bamCoverage -b $bam  --filterRNAstrand reverse -o $working_folder/BAM/$name_accession/$name_accession.R.bw --binSize=25
 
 
 
 
 
 
 
 
 
 #@##################### Calculate TPMs
 
 #single end data

working_folder=/groups/nordborg/projects/cegs/alexandra/Public_data/Bagghy
cd $working_folder
ml subread/2.0.0-foss-2018b

export list_samples=samples.txt

while read accession
do

export sample=Bhagy.$accession
echo $sample
export bamfolder=BAM
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

done < samples.txt
 
 
 #combine TPMs
 
 

samples=/groups/nordborg/projects/cegs/alexandra/Public_data/Bagghy/samples.txt

cd /groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/02_expression/featureCounts/Araport11

echo "gene" > tpm_table
cat Bhagy.Col_WT.rep2.counts_tpm.bed |  awk -v OFS="\t" '{print $1}' >> tpm_table
while read sample
do 
echo $sample
echo $sample > tpm
cat Bhagy.$sample.counts_tpm.bed | awk -v OFS="\t" '{print $4}' >> tpm
paste tpm_table tpm > inter
cat inter > tpm_table
done <$samples
cat tpm_table > Araport11.TPMs.genes.Bhagy_WT_ddm1.bed 

cp /groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/02_expression/featureCounts/Araport11/Araport11.TPMs.genes.Bhagy_WT_ddm1.bed  /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation


cd /groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/02_expression/featureCounts/20211013_annotation

echo "gene" > tpm_table
cat Bhagy.Col_WT.rep1.counts_tpm.bed |  awk -v OFS="\t" '{print $1}' >> tpm_table
while read sample
do 
echo $sample
echo $sample > tpm
cat Bhagy.$sample.counts_tpm.bed | awk -v OFS="\t" '{print $4}' >> tpm
paste tpm_table tpm > inter
cat inter > tpm_table
done <$samples
cat tpm_table > denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1.bed


cp /groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/02_expression/featureCounts/20211013_annotation/denovo_Oct2021.TPMs.genes.Bhagy_WT_ddm1.bed /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/02_expression_and_variation 

