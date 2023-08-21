#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=2:00:00
#SBATCH --array=1-125
#SBATCH --output=/groups/nordborg/pub/Mirjam/logs/array_job_slurm_%A_%a.out


export numberofline=$SLURM_ARRAY_TASK_ID

ml deeptools/3.1.2-foss-2018b-python-2.7.15
ml bedtools/2.27.1-foss-2018b


#list of samples
list=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/ChIPsamples_to_process_names_only_antibodies.txt

#extract the sample name
export sample=`awk -v a="$numberofline" 'NR==a' $list`
echo $sample
#determine which input sample to normalize against
export input=`awk -v a="$numberofline" '(NR==a) {split($1,b,"."); print b[1]"."b[2]".INPUT"}' $list`
echo $input
#set working directory
working_folder=/groups/nordborg/projects/cegs/alexandra/ChIP-seq_2021

cd $working_folder/BAM/$sample

bamCompare -b1 $sample.Aligned.sortedByCoord.out.bam -b2 $working_folder/BAM/$input/$input.Aligned.sortedByCoord.out.bam   --outFileName  $sample.log2.input_norm.bedGraph --operation log2 --effectiveGenomeSize 119481543 -p 4 --ignoreDuplicates --outFileFormat bedgraph

bamCompare -b1 $sample.Aligned.sortedByCoord.out.bam -b2 $working_folder/BAM/$input/$input.Aligned.sortedByCoord.out.bam --outFileName  $sample.subtr.input_norm.bedGraph --operation subtract --effectiveGenomeSize 119481543 -p 4 --ignoreDuplicates --outFileFormat bedgraph


bamCompare -b1 $sample.Aligned.sortedByCoord.out.bam -b2  $working_folder/BAM/$input/$input.Aligned.sortedByCoord.out.bam --outFileName   $sample.subtr.input_norm.bw --operation subtract --effectiveGenomeSize 119481543 -p 4 --ignoreDuplicates

bamCompare -b1 $sample.Aligned.sortedByCoord.out.bam -b2  $working_folder/BAM/$input/$input.Aligned.sortedByCoord.out.bam --outFileName   $sample.log2.input_norm.bw --operation log2 --effectiveGenomeSize 119481543 -p 4 --ignoreDuplicates



# calculate chipseq coverage 

dir=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation 



export denovo=$dir/20211013_annotation.5genetypes.loci.bed
# cat /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.nochrM_C.bed /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_TEs.transposable_elements.bed | sortBed -i stdin > /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.nochrM_C.and_TEs.bed

export araport=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.nochrM_C.and_TEs.bed

#TSS +- 200
#cat $denovo | awk -v OFS="\t" '{ if ($6=="+") {print $1,$2-200,$2+200,$4,$5,$6} else {print $1,$3-200,$3+200,$4,$5,$6}}'|awk -v OFS="\t" '{ if ($2<0) {$2=0}; {print $0}}'| sortBed -i stdin > $dir/20211013_annotation.5genetypes.TSS+-200bp.bed 
#cat $araport | awk -v OFS="\t" '{ if ($6=="+") {print $1,$2-200,$2+200,$4,$5,$6} else {print $1,$3-200,$3+200,$4,$5,$6}}'|awk -v OFS="\t" '{ if ($2<0) {$2=0}; {print $0}}'| sortBed -i stdin > /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.nochrM_C.and_TEs.TSS+-200bp.bed 

#TES +-200
#cat $denovo | awk -v OFS="\t" '{ if ($6=="+") {print $1,$3-200,$3+200,$4,$5,$6} else {print $1,$2-200,$2+200,$4,$5,$6}}'|awk -v OFS="\t" '{ if ($2<0) {$2=0}; {print $0}}'| sortBed -i stdin > $dir/20211013_annotation.5genetypes.TES+-200bp.bed 
#cat $araport | awk -v OFS="\t" '{ if ($6=="+") {print $1,$3-200,$3+200,$4,$5,$6} else {print $1,$2-200,$2+200,$4,$5,$6}}'|awk -v OFS="\t" '{ if ($2<0) {$2=0}; {print $0}}'| sortBed -i stdin > /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.nochrM_C.and_TEs.TES+-200bp.bed 

denovo_TSS=$dir/20211013_annotation.5genetypes.TSS+-200bp.bed 
denovo_TES=$dir/20211013_annotation.5genetypes.TES+-200bp.bed 

araport_TSS=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.nochrM_C.and_TEs.TSS+-200bp.bed 

araport_TES=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.nochrM_C.and_TEs.TES+-200bp.bed 


bedtools sort -i $denovo | bedtools map -a stdin -b $working_folder/BAM/$sample/$sample.log2.input_norm.bedGraph -c 4 -o mean  > $working_folder/coverage/$sample.denovo.log2.mean_cov.bed

bedtools sort -i $denovo_TSS | bedtools map -a stdin -b $working_folder/BAM/$sample/$sample.log2.input_norm.bedGraph -c 4 -o mean  > $working_folder/coverage/$sample.denovo_TSS.log2.mean_cov.bed

bedtools sort -i $denovo_TES | bedtools map -a stdin -b $working_folder/BAM/$sample/$sample.log2.input_norm.bedGraph -c 4 -o mean  > $working_folder/coverage/$sample.denovo_TES.log2.mean_cov.bed


bedtools sort -i $araport | bedtools map -a stdin -b $working_folder/BAM/$sample/$sample.log2.input_norm.bedGraph -c 4 -o mean > $working_folder/coverage/$sample.araport.log2.mean_cov.bed

bedtools sort -i $araport_TSS | bedtools map -a stdin -b $working_folder/BAM/$sample/$sample.log2.input_norm.bedGraph -c 4 -o mean > $working_folder/coverage/$sample.araport_TSS.log2.mean_cov.bed

bedtools sort -i $araport_TES | bedtools map -a stdin -b $working_folder/BAM/$sample/$sample.log2.input_norm.bedGraph -c 4 -o mean  > $working_folder/coverage/$sample.araport_TES.log2.mean_cov.bed












