#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30G
#SBATCH --time=6:00:00
#SBATCH --array=1-{N}
#SBATCH --output=$outputfolder/array_job_slurm_%A_%a.out


export numberofline=$SLURM_ARRAY_TASK_ID 

export accession=`sed "${numberofline}q;d" $list_of_sample_names`

ml sra-toolkit/2.9.6-1-centos_linux64
ml picard/2.18.27-java-1.8
ml bedtools/2.27.1-foss-2018b
ml samtools/1.9-foss-2018b
#!!! older version of stringtie! better at noticing antisense transcripts 
# increase filtering - remove single exons of <400bp and require TPM>0.5
###
ml stringtie/1.3.5-gcccore-7.3.0
###



##assemble transcriptome using Stringtie

mkdir $working_folder/$stringtie_folder/$name_sample
# Filename of the BAM file to input
bam_file=$bamfolder/$name_sample/$name_sample.Aligned.sortedByCoord.out.bam 

output_file=$working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.gtf          # Filename of the GTF file to output
abundance_file=$working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie_abundances.tsv       # Filename of the output table of estimated TPM values for stringtie transcripts
threads=4                                                  # Number of CPUs (Stringtie default is 10)
minimum_transcript_length=150                                 # Only outputs assembled transcripts that are >= this long
minimum_read_coverage=2                                     # Only outputs assembled transcripts with >= this many reads per nucleotide
minimum_junction_coverage=2.5                                   # Only considers splice junctions with >= this many reads
minimum_junction_overhang=15                                  # At least 1 read should overhang >= this many nt on either side of a splice junction
ref_annotation=/groups/nordborg/projects/cegs/alexandra/Annotation/Arabidopsis_thaliana.TAIR10.40.5chromosomes.gtf


stringtie_command="stringtie \
$bam_file \
-o $output_file \
-A $abundance_file \
-p $threads \
-c $minimum_read_coverage \
-m $minimum_transcript_length \
-j $minimum_junction_coverage \
-a $minimum_junction_overhang \
-x chloroplast,mitochondria --rf \
-G $ref_annotation"

cd $working_folder/$stringtie_folder/$name_sample
echo "$stringtie_command"
eval "$stringtie_command >& $stringtie.log"


# create abundance file for transcripts 
cat $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.gtf | grep -w transcript | awk -v OFS="\t"  -v val1='cov' -v val2='FPKM' -v val3='TPM'  -v val4='reference_id'  '{for (i=1; i<=NF; i++) if ($i==val1) {a=i};for (j=1;j<=NF; j++) if ($j==val2) {b=j};for (k=1; k<=NF; k++) if ($k==val3) {c=k};for (m=1; m<=NF; m++) if ($m==val4) {d=m}; { print $1, $4,$5,$12, 0, $7, $(a+1), $(b+1), $(c+1),$(d+1)} }' | awk -v OFS="\t" '{ if ($10!~/AT/) {$10="noref"}; print $0}'|  sed -e 's/"//g' | sed -e 's/;//g'|sed -e 's/chloroplast/ChrC/g'| sed -e 's/mitochondria/ChrM/g'  > $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie_tr_abundances.tsv


#convert to bed 
cat  $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.gtf |  awk -v OFS="\t" '!($7=="."){print $0}' > $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.nodots.gtf
/groups/nordborg/projects/cegs/alexandra/software/gtfToGenePred $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.nodots.gtf $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.genepred
/groups/nordborg/projects/cegs/alexandra/software/genePredToBed $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.genepred $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.1.bed
cat $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.1.bed |  sed -e 's/chloroplast/ChrC/g'| sed -e 's/mitochondria/ChrM/g'  > $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.bed


#######################################
# initial each-pool artifact filtering 
#######################################

	# main problem - antisense leakage 
#additional points:
	# 1. expression cutoff for single exon transcripts should be higher!	+ length cutoff should be stricter for them! >400nt! 
	# 2. there are highly expressed lncRNAs!!! they also give leakage artifacts!! 
	# 3. Chimeras appear already in stringtie assemblies from each pool. Then when I merge the pool assemblies with cuffmerge the chimeras are kept and the normal gene isoform is trashed by the program


#1. leakage
#2. single exon short (<400nt)! 
#3. transcripts with low expression TPM<0.1  ? (reference protein coding genes that are not expressed)
#4. chimeras? - not a gigantic problem, especially since I give Stringtie guidance.. 


# number 1: transcripts that are too lowly expressed (for example reference genes that are not expressed)
# remove all transcripts whose TPM (see table $name_sample.stringtie_tr_abundances.tsv ) is lower than 0.1 

	 
#find sense-antisense pairs of ANY assembled transcripts and find the antisense transcripts that have a suspiciously low level - <1.2% ( so the lncRNA is expressed at a level lower than expected leakage)
#make a file with "mirror ratios" - TPM ratio between sense and antisense transcript - artifact or not an artifact?
	 #all vs all 
#there needs to be a full exonic overlap with the highly expressed transcript in order to call it a leakage 
		 #use -split option since the transcript has to overlap the exon of the transcript expressed on the other strand to be suspected in leakage. Require that at least 30% of the antisense transcript's exons overlapped the corresponding sense transcript 
intersectBed -split -f 0.3 -S -wo -a $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie_tr_abundances.tsv -b $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie_tr_abundances.tsv |awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$8,$9,$14,$18,$19}' > $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.mirror_ratios.FPKMs.TPMs.bed

# !!! very strict leakage cutoff!!! hypothesis is: 
# 1. there is a basal antisense TPM of 0.06 (FPKM values are taken from stringtie's calculation, 
# 2. leakage should be no more than 1.2% (cutoff I used for choosing accessions) => any lncRNA with TPM < 0.06 + 0.012*TPM(PC) should be leakage artifact => we only keep ~36% of Antisense transcripts!! very strict
# TPM(antisense transcript) must be > 0.06+ 0.012*TPM(sense transcript) 
#TPM(antisense transcript) is column 8, TPM(sense transcript) is column 11 in $name_sample.stringtie.mirror_ratios.FPKMs.TPMs.bed 
#curve(x/(0.06 + 0.012 * x)fits the border the best. Everything below the line is OK. Everything above the line will be discarded 
#checked removed artifacts with this formula (by eye a ~50 loci) - looks very good! The only few residual artifacts can be removed by size selection - single exon transcripts must be >300-400bp

cat $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.mirror_ratios.FPKMs.TPMs.bed \
| awk -v OFS="\t" '($8<($11*0.012+0.06)){print $4}' \
| uniq | awk -v OFS="\t" '{print $0,1}' > $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.potential_leakage_artifacts.names.txt
#cat $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.potential_leakage_artifacts.names.txt |sed  's/^/transcript_id "/;s/$/"/'> $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.potential_leakage_artifacts.names.withquotesforgrep.txt

cat $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.mirror_ratios.FPKMs.TPMs.bed | awk -v OFS="\t" '($8<($11*0.012+0.06)){print $1,$2,$3,$4,$5,$6}' > $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.potential_leakage_artifacts.bed


#find chimeric transcript assemblies 

#make araport 11 pc gene loci annotation (find max and min coordinates for each gene name)
#cat /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.bed |awk -v OFS="\t"  '{split($4,s,"."); print $1,$2,$3,s[1],$5,$6}' | uniq| awk -v OFS="\t"  '{a[$4]+=$3;}END{for (i in a)print i,a[i];}'| sort | head 
#cat /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.bed |awk -v OFS="\t"  '{split($4,s,"."); print $1,$2,$3,s[1],$5,$6}' | uniq> /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes_draft.bed

#intersectBed -s -f 0.3 -e -F 0.3 -wo  -a /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes_draft.bed  -b /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes_draft.bed | awk -v OFS="\t" '($4!=$10){print $4}' >  /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes_senseoverlapping_other_genes.names.txt

#cat /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes_draft.bed | grep -v -w -f /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes_senseoverlapping_other_genes.names.txt > /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes_nooverlapping_genes.bed


#intersectBed -s -f 0.3 -e -F 0.3 -c  -a $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.bed  -b /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes_nooverlapping_genes.bed |awk -v OFS="\t" '{print $1,$4,$10}' | awk -v OFS="\t" '($3>1){print $2,1}'> $working_folder/$stringtie_folder/$name_sample/$name_sample.chimeric_transcripts_artedacts_names.txt

intersectBed -s -f 0.3 -e -F 0.3 -wo -a $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.bed  -b /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes_nooverlapping_genes.bed|awk -v OFS="\t" '{print $4,$16}' | uniq|awk -v OFS="\t" '{print $1}'| uniq -d >  $working_folder/$stringtie_folder/$name_sample/$name_sample.chimeric_transcripts_artedacts_names.txt


#all artefacts 
cat $working_folder/$stringtie_folder/$name_sample/$name_sample.chimeric_transcripts_artedacts_names.txt $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.potential_leakage_artifacts.names.txt > $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.artifacts.names.txt

#filter the stringtie annotation 

#add a column with transcript name for filtering
cat $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.gtf | grep -v "StringTie version" | grep -v "stringtie" | awk -v OFS="\t" '{split($12,a,"\""); print a[2],$0}' > $working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.forfilt.gtf

export bedloc=$working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.bed 
export abundance=$working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie_tr_abundances.tsv 
export gtf_loc=$working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.forfilt.gtf
export artefacts=$working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.artifacts.names.txt
export bedloc_filter=$working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.filtered.bed 
export gtf_loc_filter=$working_folder/$stringtie_folder/$name_sample/$name_sample.stringtie.filtered.gtf



module load r/3.5.1-foss-2018b
Rscript $rscript_folder/stringtie_output_filtering_for_monoaccession_final.r

