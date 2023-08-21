
# Papareddy et al 2020  (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02163-4)

samples_srr=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/Public_data/small_RNA_and_Methyl_Seq_Papareddy_et_al_2020/Ranj_srr_samplename.txt
samples=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/Public_data/small_RNA_and_Methyl_Seq_Papareddy_et_al_2020/Ranj_samplename.txt

module load sra-toolkit/2.9.6-1-centos_linux64
module load bwa/0.7.17-foss-2018b
module load picard/2.18.27-java-1.8
module load star/2.7.1a-foss-2018b
ml cutadapt/1.18-foss-2018b-python-3.6.6
ml bedtools/2.27.1-foss-2018b


export working_folder=/groups/nordborg/projects/cegs/alexandra/Public_data/Ranj


#Z:\01_POSTDOC\03_Projects\Public_data\small_RNA_and_Methyl_Seq_Papareddy_et_al_2020\Ranj_sRNA_cuadapt.bash
#ml cutadapt/1.18-foss-2018b-python-3.6.6

#export srr=`cat $samples_srr | sed -n $numberofline,"$numberofline"p  | awk '{print $1}'`
#export sample=`cat $samples_srr | sed -n $numberofline,"$numberofline"p  | awk '{print $2}'`

#echo $srr
#echo $sample
#cd $working_folder

#convert sra to fastq.gz
#fastq-dump.2.9.6 --outdir fastq --gzip --skip-technical  --readids --dumpbase --clip $srr.1
#cutadapt -a AGATCGGAAGA -m 18 -o $working_folder/fastq/$sample.trimmed.fastq $working_folder/fastq/$srr.1.fastq.gz

# align with STAR

while read sample
do
ml bedtools/2.27.1-foss-2018b
mkdir $working_folder/BAM/$sample
echo $sample
export stargenomedir=/groups/nordborg/user/aleksandra.kornienko/analyses/STAR/STAR2.7.1a/STAR_TAIR10_withSJDB_sizeaccounted
export out=$working_folder/BAM/$sample/$sample.
fastq=$working_folder/fastq/$sample


STAR --runMode alignReads \
--runThreadN 4 \
--runRNGseed 12345 \
--genomeDir $stargenomedir \
--readFilesIn $fastq.trimmed.fastq  \
--outFileNamePrefix $out \
--limitBAMsortRAM 30000000000 \
--alignEndsType Extend5pOfRead1 \
--alignIntronMax 5000 \
--alignSJDBoverhangMin 1 \
--outReadsUnmapped Fastx \
--outSAMtype BAM Unsorted \
--outSAMmultNmax 100 \
--outSAMprimaryFlag AllBestScore \
--outSAMattributes NH HI AS nM NM MD jM jI XS \
--outFilterMultimapNmax 10 \
--outFilterMatchNmin 16 \
--outFilterMatchNminOverLread 0.66 \
--outFilterMismatchNmax 2 \
--outFilterMismatchNoverReadLmax 0.05 \
--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
--twopassMode None \
--quantMode GeneCounts

done < $samples



ml samtools/1.9-foss-2018b
while read sample
do
samtools view $working_folder/BAM/$sample/$sample.Aligned.out.bam |  awk -v OFS="\t" '{print $3,$4,$5,$6,$10}' | sed 's/M//g'|  awk -v OFS="\t" '($4==21||$4==22){print $1,$2,$2+$4,$5,$4}' | grep -v chloroplast | grep -v mitochondria | sortBed -i stdin  > $working_folder/BAM/$sample/$sample.21-22nt.bed

samtools view $working_folder/BAM/$sample/$sample.Aligned.out.bam |  awk -v OFS="\t" '{print $3,$4,$5,$6,$10}' | sed 's/M//g'|  awk -v OFS="\t" '($4==24){print $1,$2,$2+$4,$5,$4}' | grep -v chloroplast | grep -v mitochondria | sortBed -i stdin  > $working_folder/BAM/$sample/$sample.24nt.bed

dir=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation 
export denovo=$dir/20211013_annotation.5genetypes.loci.bed
export araport=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.nochrM_C.and_TEs.bed

genomeCoverageBed -i $working_folder/BAM/$sample/$sample.24nt.bed -g /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/chr_length.txt  -bg > $working_folder/BAM/$sample/$sample.24nt.bedgraph
/groups/nordborg/projects/cegs/alexandra/software/bedGraphToBigWig $working_folder/BAM/$sample/$sample.24nt.bedgraph /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/chr_length.txt $working_folder/BAM/$sample/$sample.24nt.bw

genomeCoverageBed -i $working_folder/BAM/$sample/$sample.21-22nt.bed -g /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/chr_length.txt  -bg > $working_folder/BAM/$sample/$sample.21-22nt.bedgraph
/groups/nordborg/projects/cegs/alexandra/software/bedGraphToBigWig $working_folder/BAM/$sample/$sample.21-22nt.bedgraph /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/chr_length.txt $working_folder/BAM/$sample/$sample.21-22nt.bw

done < $samples




#combine read number for normalization
echo "accession" "total_readN" "uniq_readN" "multimap_readN"> $working_folder/readN_stats_Rnaj_samples.bed
cd $working_folder/BAM
while read sample
do 
echo $sample
readN=`cat $sample/$sample.Log.final.out | grep "Uniquely mapped reads number"| awk '{print $6}'` 
readNtotal=`cat $sample/$sample.Log.final.out | grep "Number of input reads"| awk '{print $6}'` 
readNmulti=`cat $sample/$sample.Log.final.out | grep "Number of reads mapped to multiple loci"| awk '{print $9}'` 
echo $sample $readNtotal $readN $readNmulti 
echo $sample $readNtotal $readN $readNmulti >> $working_folder/readN_stats_Rnaj_samples.bed

done < $samples

cp $working_folder/readN_stats_Rnaj_samples.bed /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/Public_data/small_RNA_and_Methyl_Seq_Papareddy_et_al_2020/


#calculate read coverage for different types of genes
while read sample
do
echo $sample
readN=`cat $working_folder/BAM/$sample/$sample.Log.final.out | grep "Uniquely mapped reads number"| awk '{print $6}'` 

dir=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation 

export denovo=$dir/20211013_annotation.5genetypes.loci.bed
export araport=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.nochrM_C.and_TEs.bed

genomeCoverageBed -i $working_folder/BAM/$sample/$sample.24nt.bed -g /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/chr_length.txt  -d | awk -v OFS="\t" -v readN="$readN" '{print $1,$2-1,$2,$3*1000000/readN}' > $working_folder/BAM/$sample/$sample.24nt.perbase.bed
genomeCoverageBed -i $working_folder/BAM/$sample/$sample.21-22nt.bed -g /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/chr_length.txt  -d | awk -v OFS="\t" -v readN="$readN" '{print $1,$2-1,$2,$3*1000000/readN}' > $working_folder/BAM/$sample/$sample.21-22nt.perbase.bed

bedtools sort -i $denovo | bedtools map -a stdin -b $working_folder/BAM/$sample/$sample.24nt.perbase.bed -c 4 -o mean | sort  -k4,4 > $working_folder/BAM/$sample/$sample.24nt.denovo_Oct2021.coverage.perbase_calc.bed
bedtools sort -i $denovo | bedtools map -a stdin -b $working_folder/BAM/$sample/$sample.21-22nt.perbase.bed -c 4 -o mean | sort  -k4,4 > $working_folder/BAM/$sample/$sample.21-22nt.denovo_Oct2021.coverage.perbase_calc.bed

bedtools sort -i $araport | bedtools map -a stdin -b $working_folder/BAM/$sample/$sample.24nt.perbase.bed -c 4 -o mean | sort  -k4,4 > $working_folder/BAM/$sample/$sample.24nt.Araport11.coverage.perbase_calc.bed
bedtools sort -i $araport | bedtools map -a stdin -b $working_folder/BAM/$sample/$sample.21-22nt.perbase.bed -c 4 -o mean |sort  -k4,4 > $working_folder/BAM/$sample/$sample.21-22nt.Araport11.coverage.perbase_calc.bed
done < $samples



# combine coverage 
cd $working_folder/BAM

cat nrpda3.leaf.r3/nrpda3.leaf.r3.21-22nt.denovo_Oct2021.coverage.perbase_calc.bed| awk -v OFS="\t" '{print "Chr","start","end","gene","score","strand"}'| head -1 > table_21_22
cat nrpda3.leaf.r3/nrpda3.leaf.r3.21-22nt.denovo_Oct2021.coverage.perbase_calc.bed| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' >> table_21_22


cat nrpda3.leaf.r3/nrpda3.leaf.r3.24nt.denovo_Oct2021.coverage.perbase_calc.bed| awk -v OFS="\t" '{print "Chr","start","end","gene","score","strand"}'| head -1 > table_24
cat nrpda3.leaf.r3/nrpda3.leaf.r3.24nt.denovo_Oct2021.coverage.perbase_calc.bed| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' >> table_24


cat nrpda3.leaf.r3/nrpda3.leaf.r3.21-22nt.Araport11.coverage.perbase_calc.bed| awk -v OFS="\t" '{print "Chr","start","end","gene","score","strand"}'| head -1 > table_21_22_ar11
cat nrpda3.leaf.r3/nrpda3.leaf.r3.21-22nt.Araport11.coverage.perbase_calc.bed| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' >> table_21_22_ar11


cat nrpda3.leaf.r3/nrpda3.leaf.r3.24nt.Araport11.coverage.perbase_calc.bed| awk -v OFS="\t" '{print "Chr","start","end","gene","score","strand"}'| head -1 > table_24_ar11
cat nrpda3.leaf.r3/nrpda3.leaf.r3.24nt.Araport11.coverage.perbase_calc.bed| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' >> table_24_ar11



while read sample 
do 
echo $sample

echo $sample > column21_22
echo $sample > column24
echo $sample > column21_22_ar11
echo $sample > column24_ar11

cat $sample/$sample.21-22nt.denovo_Oct2021.coverage.perbase_calc.bed |  awk -v OFS="\t" '{print $7}'  >> column21_22
paste table_21_22 column21_22 > inter 
cat inter > table_21_22

cat $sample/$sample.24nt.denovo_Oct2021.coverage.perbase_calc.bed |  awk -v OFS="\t" '{print $7}'  >> column24
paste table_24 column24 > inter 
cat inter > table_24

cat $sample/$sample.21-22nt.Araport11.coverage.perbase_calc.bed |  awk -v OFS="\t" '{print $7}'  >> column21_22_ar11
paste table_21_22_ar11 column21_22_ar11 > inter 
cat inter > table_21_22_ar11

cat $sample/$sample.24nt.Araport11.coverage.perbase_calc.bed |  awk -v OFS="\t" '{print $7}'  >> column24_ar11
paste table_24_ar11 column24_ar11 > inter 
cat inter > table_24_ar11

done < $samples


cat table_21_22 > /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/Public_data/small_RNA_and_Methyl_Seq_Papareddy_et_al_2020/Ranj.sRNA.21-22nt.denovo_Oct2021.coverage.perbase_calc.normalized.bed 

cat table_24 > /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/Public_data/small_RNA_and_Methyl_Seq_Papareddy_et_al_2020/Ranj.sRNA.24nt.denovo_Oct2021.coverage.perbase_calc.normalized.bed 



cat table_21_22_ar11 > /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/Public_data/small_RNA_and_Methyl_Seq_Papareddy_et_al_2020/Ranj.sRNA.21-22nt.Araport11.coverage.perbase_calc.normalized.bed 

cat table_24_ar11 > /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/Public_data/small_RNA_and_Methyl_Seq_Papareddy_et_al_2020/Ranj.sRNA.24nt.Araport11.coverage.perbase_calc.normalized.bed 























