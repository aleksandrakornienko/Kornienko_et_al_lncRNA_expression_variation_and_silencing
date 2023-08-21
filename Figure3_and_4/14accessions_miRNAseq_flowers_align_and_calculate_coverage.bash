
module load sra-toolkit/2.9.6-1-centos_linux64
module load bwa/0.7.17-foss-2018b
module load picard/2.18.27-java-1.8
module load star/2.7.1a-foss-2018b
ml cutadapt/1.18-foss-2018b-python-3.6.6
ml bedtools/2.27.1-foss-2018b


export working_folder=/groups/nordborg/projects/cegs/alexandra/miRNAseq
export accessions_and_unmappedBam=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/RNA-seq_data/miRNA/14_miRNA_samples_flowers.txt
export unmappedbam_folder=/groups/nordborg/projects/cegs/alexandra/miRNAseq/raw_data





while read name_accession
do 
echo $name_accession
export unmapped_BAM1=`cat $accessions_and_unmappedBam | grep -w  $name_accession | awk '{print $1}'`
echo $unmapped_BAM1

export unmapped_bamfile=$unmappedbam_folder/$unmapped_BAM1
echo $unmapped_bamfile 
# convert unmapped BAM to fastq
java -jar ${EBROOTPICARD}/picard.jar SamToFastq I=$unmapped_bamfile      FASTQ=$working_folder/fastq/$name_accession.fastq  

mkdir $working_folder/BAM/$name_accession
export stargenomedir=/groups/nordborg/user/aleksandra.kornienko/analyses/STAR/STAR2.7.1a/STAR_TAIR10_withSJDB_sizeaccounted
export out=$working_folder/BAM/$name_accession/$name_accession.
fastq=$working_folder/fastq/$name_accession

#trim reads
cutadapt -a AACTGTAGGCACCATCAAT --minimum-length 18 -o $fastq.trimmed.fastq $fastq.fastq

export working_folder=/groups/nordborg/projects/cegs/alexandra/miRNAseq

ml bedtools/2.27.1-foss-2018b

#align reads to TAIR10 using STAR
export stargenomedir=/groups/nordborg/user/aleksandra.kornienko/analyses/STAR/STAR2.7.1a/STAR_TAIR10_withSJDB_sizeaccounted
export out=$working_folder/BAM/$name_accession/$name_accession.
fastq=$working_folder/fastq/$name_accession

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

ml samtools/1.9-foss-2018b

#extract 21-22nt reads
samtools view $working_folder/BAM/$name_accession/$name_accession.Aligned.out.bam |  awk -v OFS="\t" '{print $3,$4,$5,$6,$10}' | sed 's/M//g'|  awk -v OFS="\t" '($4==21||$4==22){print $1,$2,$2+$4,$5,$4}' | grep -v chloroplast | grep -v mitochondria | sortBed -i stdin  > $working_folder/BAM/$name_accession/$name_accession.21-22nt.bed

#extract 24nt reads
samtools view $working_folder/BAM/$name_accession/$name_accession.Aligned.out.bam |  awk -v OFS="\t" '{print $3,$4,$5,$6,$10}' | sed 's/M//g'|  awk -v OFS="\t" '($4==24){print $1,$2,$2+$4,$5,$4}' | grep -v chloroplast | grep -v mitochondria | sortBed -i stdin  > $working_folder/BAM/$name_accession/$name_accession.24nt.bed

#convert bed to bedgraph
genomeCoverageBed -i $working_folder/BAM/$name_accession/$name_accession.24nt.bed -g /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/chr_length.txt  -bg > $working_folder/BAM/$name_accession/$name_accession.24nt.bedgraph
#convert bedgraph to bigwig
/groups/nordborg/projects/cegs/alexandra/software/bedGraphToBigWig $working_folder/BAM/$name_accession/$name_accession.24nt.bedgraph /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/chr_length.txt $working_folder/BAM/$name_accession/$name_accession.24nt.bw

#convert bed to bedgraph
genomeCoverageBed -i $working_folder/BAM/$name_accession/$name_accession.21-22nt.bed -g /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/chr_length.txt  -bg > $working_folder/BAM/$name_accession/$name_accession.21-22nt.bedgraph
#convert bedgraph to bigwig
/groups/nordborg/projects/cegs/alexandra/software/bedGraphToBigWig $working_folder/BAM/$name_accession/$name_accession.21-22nt.bedgraph /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/chr_length.txt $working_folder/BAM/$name_accession/$name_accession.21-22nt.bw

done < /groups/nordborg/projects/cegs/alexandra/miRNAseq/names_accessions_mirnaseq.txt





dir=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation 
export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation
denovo_linc=$dir/lncRNAs.intergenic.loci.bed
denovo_linc_names=$dir/lncRNAs.intergenic.loci.names.txt
denovo_as=$dir/lncRNAs.antisense.loci.bed
denovo_as_names=$dir/lncRNAs.antisense.loci.names.txt
denovo_pc=$dir/denovoPC.loci.bed
denovo_pc_names=$dir/denovoPC.loci.names.txt
denovo_tegenes=$dir/loci.TE_genes.sorted.bed
denovo_tegenes_names=$dir/loci.TE_genes.names.txt
as_to_te=$dir/lncRNAs.AS_to_TE.loci.bed

araport_pc=$annotationfolder/Araport11_protein_coding.201606.genes.nochrM_C.bed

#cat $denovo_linc  $denovo_as $denovo_pc $denovo_tegenes $as_to_te | sortBed -i stdin > $dir/20211013_annotation.5genetypes.loci.bed


export denovo=$dir/20211013_annotation.5genetypes.loci.bed
export araport=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.nochrM_C.and_TEs.bed



# calculate read number for normalisation 
echo "accession" "total_readN" "uniq_readN" "multimap_readN">/groups/nordborg/projects/cegs/alexandra/miRNAseq/readN_stats14samples.bed
cd BAM
while read name_accession
do 
echo $name_accession
cat 10002/10002.Log.final.out | grep "Number of input reads"
readN=`cat $name_accession/$name_accession.Log.final.out | grep "Uniquely mapped reads number"| awk '{print $6}'` 
readNtotal=`cat $name_accession/$name_accession.Log.final.out | grep "Number of input reads"| awk '{print $6}'` 
readNmulti=`cat $name_accession/$name_accession.Log.final.out | grep "Number of reads mapped to multiple loci"| awk '{print $9}'` 
echo $name_accession $readNtotal $readN $readNmulti >> /groups/nordborg/projects/cegs/alexandra/miRNAseq/readN_stats14samples.bed
done < /groups/nordborg/projects/cegs/alexandra/miRNAseq/names_accessions_mirnaseq.txt


cp /groups/nordborg/projects/cegs/alexandra/miRNAseq/readN_stats14samples.bed /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/RNA-seq_data/miRNA/


# calculate sRNA coverage on different kinds of genes

while read name_accession
do 
echo $name_accession

readN=`cat $name_accession/$name_accession.Log.final.out | grep "Uniquely mapped reads number"| awk '{print $6}'` 

#calculate coverage per base for the whole genome (and divide by total read number to normalize)
genomeCoverageBed -i $working_folder/BAM/$name_accession/$name_accession.24nt.bed -g /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/chr_length.txt  -d | awk -v OFS="\t" -v readN="$readN" '{print $1,$2-1,$2,$3*1000000/readN}' > $working_folder/BAM/$name_accession/$name_accession.24nt.perbase.bed
genomeCoverageBed -i $working_folder/BAM/$name_accession/$name_accession.21-22nt.bed -g /groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/chr_length.txt  -d | awk -v OFS="\t" -v readN="$readN" '{print $1,$2-1,$2,$3*1000000/readN}' > $working_folder/BAM/$name_accession/$name_accession.21-22nt.perbase.bed

#calculate average perbase coverage for each gene in denovo annotation
bedtools sort -i $denovo | bedtools map -a stdin -b $working_folder/BAM/$name_accession/$name_accession.24nt.perbase.bed -c 4 -o mean | sort  -k4,4 > $working_folder/BAM/$name_accession/$name_accession.24nt.denovo_Oct2021.coverage.perbase_calc.bed
bedtools sort -i $denovo | bedtools map -a stdin -b $working_folder/BAM/$name_accession/$name_accession.21-22nt.perbase.bed -c 4 -o mean | sort  -k4,4 > $working_folder/BAM/$name_accession/$name_accession.21-22nt.denovo_Oct2021.coverage.perbase_calc.bed

#calculate average perbase coverage for each gene in Araport annotation
bedtools sort -i $araport | bedtools map -a stdin -b $working_folder/BAM/$name_accession/$name_accession.24nt.perbase.bed -c 4 -o mean | sort  -k4,4 > $working_folder/BAM/$name_accession/$name_accession.24nt.Araport11.coverage.perbase_calc.bed
bedtools sort -i $araport | bedtools map -a stdin -b $working_folder/BAM/$name_accession/$name_accession.21-22nt.perbase.bed -c 4 -o mean |sort  -k4,4 > $working_folder/BAM/$name_accession/$name_accession.21-22nt.Araport11.coverage.perbase_calc.bed
done < /groups/nordborg/projects/cegs/alexandra/miRNAseq/names_accessions_mirnaseq.txt



# combine coverage 

cd /groups/nordborg/projects/cegs/alexandra/miRNAseq
cat BAM/10002/10002.21-22nt.denovo_Oct2021.coverage.perbase_calc.bed| awk -v OFS="\t" '{print "Chr","start","end","gene","score","strand"}'| head -1 > table_21_22
cat BAM/10002/10002.21-22nt.denovo_Oct2021.coverage.perbase_calc.bed| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' >> table_21_22


cat BAM/10002/10002.24nt.denovo_Oct2021.coverage.perbase_calc.bed| awk -v OFS="\t" '{print "Chr","start","end","gene","score","strand"}'| head -1 > table_24
cat BAM/10002/10002.24nt.denovo_Oct2021.coverage.perbase_calc.bed| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' >> table_24


cat BAM/10002/10002.21-22nt.Araport11.coverage.perbase_calc.bed| awk -v OFS="\t" '{print "Chr","start","end","gene","score","strand"}'| head -1 > table_21_22_ar11
cat BAM/10002/10002.21-22nt.Araport11.coverage.perbase_calc.bed| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' >> table_21_22_ar11


cat BAM/10002/10002.24nt.Araport11.coverage.perbase_calc.bed| awk -v OFS="\t" '{print "Chr","start","end","gene","score","strand"}'| head -1 > table_24_ar11
cat BAM/10002/10002.24nt.Araport11.coverage.perbase_calc.bed| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' >> table_24_ar11



while read name_accession 
do 
echo $name_accession

echo $name_accession > column21_22
echo $name_accession > column24
echo $name_accession > column21_22_ar11
echo $name_accession > column24_ar11

cat BAM/$name_accession/$name_accession.21-22nt.denovo_Oct2021.coverage.perbase_calc.bed |  awk -v OFS="\t" '{print $7}'  >> column21_22
paste table_21_22 column21_22 > inter 
cat inter > table_21_22

cat BAM/$name_accession/$name_accession.24nt.denovo_Oct2021.coverage.perbase_calc.bed |  awk -v OFS="\t" '{print $7}'  >> column24
paste table_24 column24 > inter 
cat inter > table_24

cat BAM/$name_accession/$name_accession.21-22nt.Araport11.coverage.perbase_calc.bed |  awk -v OFS="\t" '{print $7}'  >> column21_22_ar11
paste table_21_22_ar11 column21_22_ar11 > inter 
cat inter > table_21_22_ar11

cat BAM/$name_accession/$name_accession.24nt.Araport11.coverage.perbase_calc.bed |  awk -v OFS="\t" '{print $7}'  >> column24_ar11
paste table_24_ar11 column24_ar11 > inter 
cat inter > table_24_ar11

done < $working_folder/names_accessions_mirnaseq.txt


cat table_21_22 > /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/RNA-seq_data/miRNA/flowers_14acc.21-22nt.denovo_Oct2021.coverage.perbase_calc.normalized.bed 

cat table_24 > /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/RNA-seq_data/miRNA/flowers_14acc.24nt.denovo_Oct2021.coverage.perbase_calc.normalized.bed 



cat table_21_22_ar11 > /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/RNA-seq_data/miRNA/flowers_14acc.21-22nt.Araport11.coverage.perbase_calc.normalized.bed 

cat table_24_ar11 > /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/RNA-seq_data/miRNA/flowers_14acc.24nt.Araport11.coverage.perbase_calc.normalized.bed 
