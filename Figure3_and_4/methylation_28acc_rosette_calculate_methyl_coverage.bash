#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=20G
#SBATCH --time=6:00:00
#SBATCH --array=1-93
#SBATCH --output=/groups/nordborg/pub/Mirjam/logs/array_job_slurm_%A_%a.out



export numberofline=$SLURM_ARRAY_TASK_ID


ml bedtools/2.27.1-foss-2018b

cd  /groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/methylation/additional_samples
#mkdir CG
#mkdir CHG
#mkdir CHH

methcalls_folder=/groups/nordborg/projects/nordborg_common/datasets/alexandra_Ath_bisulfiteseq/methylseq/methylpy
export samples=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/2022_methylation_additionalsamples/AK_methyl_samplenames.txt

dir=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation 
export denovo=$dir/20211013_annotation.5genetypes.loci.bed
export denovo_TSS=$dir/20211013_annotation.5genetypes.TSS+-200bp.bed
export denovo_TES=$dir/20211013_annotation.5genetypes.TES+-200bp.bed

 
export araport=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.nochrM_C.and_TEs.bed
export araport_TES=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.nochrM_C.and_TEs.TES+-200bp.bed
export araport_TSS=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.nochrM_C.and_TEs.TSS+-200bp.bed



export sampleN=`sed -n $numberofline,"$numberofline"p $samples | awk '{print $1}'`

export sample_name=`sed -n $numberofline,"$numberofline"p $samples | awk '{print $2}'`


echo $sampleN

echo $sample_name


cd  /groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/methylation/additional_samples
echo $sampleN
index=allc_
index+=$sampleN
##############################
echo $index

gunzip -c $methcalls_folder/$index.tsv.gz > $methcalls_folder/$index.tsv
# take only CG methylation info, reformat the file into bed6-like file 
cat $methcalls_folder/$index.tsv | awk -v OFS="\t"   '{if (($4=="CGG" || $4=="CGA" || $4=="CGT"|| $4=="CGC")) print $1,$2,$2, $4, 0, $3, $5,$6,$7 }'  > CG/$index.CG.bed

#  take only CHG methylation info, reformat the file into bed6-like file 
cat  $methcalls_folder/$index.tsv	| awk -v OFS="\t"   '{if ($4=="CAG" || $4=="CCG" || $4=="CTG") print $1,$2,$2, $4, 0, $3, $5,$6,$7 }'  > CHG/$index.CHG.bed

cat $methcalls_folder/$index.tsv	| awk -v OFS="\t"   '{if ($4=="CAA" || $4=="CAT" || $4=="CAC"|| $4=="CTA"|| $4=="CTT"|| $4=="CTC"|| $4=="CCA"|| $4=="CCT"|| $4=="CCC") print $1,$2,$2, $4, 0, $3, $5,$6,$7 }'  > CHH/$index.CHH.bed

rm   $methcalls_folder/$index.tsv



sortBed -i $denovo |bedtools map -a stdin -b CG/$index.CG.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CG/$sample_name.CG.denovo.bed
sortBed -i $denovo_TSS |bedtools map -a stdin    -b CG/$index.CG.bed -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CG/$sample_name.CG.denovo_TSS.bed
sortBed -i $denovo_TES |bedtools map -a stdin    -b CG/$index.CG.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CG/$sample_name.CG.denovo_TES.bed



sortBed -i $araport |bedtools map -a stdin -b CG/$index.CG.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CG/$sample_name.CG.araport.bed
sortBed -i $araport_TSS |bedtools map -a stdin    -b CG/$index.CG.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CG/$sample_name.CG.araport_TSS.bed
sortBed -i $araport_TES |bedtools map -a stdin    -b CG/$index.CG.bed -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CG/$sample_name.CG.araport_TES.bed


sortBed -i $denovo |bedtools map -a stdin -b CHG/$index.CHG.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CHG/$sample_name.CHG.denovo.bed
sortBed -i $denovo_TSS |bedtools map -a stdin    -b CHG/$index.CHG.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CHG/$sample_name.CHG.denovo_TSS.bed
sortBed -i $denovo_TES |bedtools map -a stdin    -b CHG/$index.CHG.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CHG/$sample_name.CHG.denovo_TES.bed

sortBed -i $araport |bedtools map -a stdin -b CHG/$index.CHG.bed -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CHG/$sample_name.CHG.araport.bed
sortBed -i $araport_TSS |bedtools map -a stdin    -b CHG/$index.CHG.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CHG/$sample_name.CHG.araport_TSS.bed
sortBed -i $araport_TES |bedtools map -a stdin    -b CHG/$index.CHG.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CHG/$sample_name.CHG.araport_TES.bed



sortBed -i $denovo |bedtools map -a stdin -b CHH/$index.CHH.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CHH/$sample_name.CHH.denovo.bed
sortBed -i $denovo_TSS |bedtools map -a stdin    -b CHH/$index.CHH.bed -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CHH/$sample_name.CHH.denovo_TSS.bed
sortBed -i $denovo_TES |bedtools map -a stdin    -b CHH/$index.CHH.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CHH/$sample_name.CHH.denovo_TES.bed

sortBed -i $araport |bedtools map -a stdin -b CHH/$index.CHH.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CHH/$sample_name.CHH.araport.bed
sortBed -i $araport_TSS |bedtools map -a stdin    -b CHH/$index.CHH.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CHH/$sample_name.CHH.araport_TSS.bed
sortBed -i $araport_TES |bedtools map -a stdin    -b CHH/$index.CHH.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CHH/$sample_name.CHH.araport_TES.bed



# calculate average  methylation levels

#CG
cat CG/$index.CG.bed | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }'  > CG/$sample_name.aver  
echo $sample_name > CG/$sample_name.name 
#calculate average CG methylation 
#PC
intersectBed -u -a CG/$index.CG.bed -b $pc | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$sample_name.aver_pc
#linc
intersectBed -u -a CG/$index.CG.bed -b $as | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$sample_name.aver_linc
#AS
intersectBed -u -a CG/$index.CG.bed -b $linc | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$sample_name.aver_as
#te genes
intersectBed -u -a CG/$index.CG.bed -b $te | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$sample_name.aver_tegene

#all tair10 TEs
intersectBed -u -a CG/$index.CG.bed -b $TAIR10_TEs | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$sample_name.aver_te
cat $TAIR10_TEs | grep  LTR_Copia|  intersectBed -u -a CG/$index.CG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$sample_name.aver_LTR_Copia
cat $TAIR10_TEs | grep  RC_Helitron|  intersectBed -u -a CG/$index.CG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$sample_name.aver_RC_Helitron
cat $TAIR10_TEs | grep  DNA_MuDR|  intersectBed -u -a CG/$index.CG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$sample_name.aver_DNA_MuDR
cat $TAIR10_TEs | grep  DNA|  intersectBed -u -a CG/$index.CG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$sample_name.aver_allDNAtes
cat $TAIR10_TEs | grep  LTR_Gypsy|  intersectBed -u -a CG/$index.CG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$sample_name.aver_LTR_Gypsy
cat $TAIR10_TEs | grep  LINE|  intersectBed -u -a CG/$index.CG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$sample_name.aver_LINE
cat $TAIR10_TEs | grep  SINE|  intersectBed -u -a CG/$index.CG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$sample_name.aver_SINE

paste CG/$sample_name.name  CG/$sample_name.aver  CG/$sample_name.aver_pc CG/$sample_name.aver_as CG/$sample_name.aver_linc CG/$sample_name.aver_tegene CG/$sample_name.aver_te CG/$sample_name.aver_LTR_Copia CG/$sample_name.aver_LTR_Gypsy CG/$sample_name.aver_RC_Helitron  CG/$sample_name.aver_DNA_MuDR CG/$sample_name.aver_allDNAtes CG/$sample_name.aver_allDNAtes  CG/$sample_name.aver_LINE  CG/$sample_name.aver_SINE


#CHG
cat CHG/$index.CHG.bed | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$sample_name.aver  
echo $sample_name > CHG/$sample_name.name 
#calculate average CHG methylation 
#PC
intersectBed -u -a CHG/$index.CHG.bed -b $pc | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$sample_name.aver_pc
#linc
intersectBed -u -a CHG/$index.CHG.bed -b $linc | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$sample_name.aver_linc
#AS
intersectBed -u -a CHG/$index.CHG.bed -b $as | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$sample_name.aver_as
#te genes
intersectBed -u -a CHG/$index.CHG.bed -b $te | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$sample_name.aver_tegene
#all tair10 TEs
intersectBed -u -a CHG/$index.CHG.bed -b $TAIR10_TEs | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$sample_name.aver_te
cat $TAIR10_TEs | grep  LTR_Copia|  intersectBed -u -a CHG/$index.CHG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$sample_name.aver_LTR_Copia
cat $TAIR10_TEs | grep  RC_Helitron|  intersectBed -u -a CHG/$index.CHG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$sample_name.aver_RC_Helitron
cat $TAIR10_TEs | grep  DNA_MuDR|  intersectBed -u -a CHG/$index.CHG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$sample_name.aver_DNA_MuDR
cat $TAIR10_TEs | grep  DNA|  intersectBed -u -a CHG/$index.CHG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$sample_name.aver_allDNAtes
cat $TAIR10_TEs | grep  LTR_Gypsy|  intersectBed -u -a CHG/$index.CHG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$sample_name.aver_LTR_Gypsy
cat $TAIR10_TEs | grep  LINE|  intersectBed -u -a CHG/$index.CHG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$sample_name.aver_LINE
cat $TAIR10_TEs | grep  SINE|  intersectBed -u -a CHG/$index.CHG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$sample_name.aver_SINE


paste CHG/$sample_name.name  CHG/$sample_name.aver  CHG/$sample_name.aver_pc CHG/$sample_name.aver_as CHG/$sample_name.aver_linc CHG/$sample_name.aver_tegene CHG/$sample_name.aver_te CHG/$sample_name.aver_LTR_Copia CHG/$sample_name.aver_LTR_Gypsy CHG/$sample_name.aver_RC_Helitron  CHG/$sample_name.aver_DNA_MuDR CHG/$sample_name.aver_allDNAtes CHG/$sample_name.aver_allDNAtes  CHG/$sample_name.aver_LINE  CHG/$sample_name.aver_SINE



#CHH
cat CHH/$index.CHH.bed | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$sample_name.aver  
echo $sample_name > CHH/$sample_name.name 
#calculate average CHH methylation 
#PC
intersectBed -u -a CHH/$index.CHH.bed -b $pc | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$sample_name.aver_pc
#linc
intersectBed -u -a CHH/$index.CHH.bed -b $linc | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$sample_name.aver_linc
#AS
intersectBed -u -a CHH/$index.CHH.bed -b $as | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$sample_name.aver_as
#te genes
intersectBed -u -a CHH/$index.CHH.bed -b $te | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$sample_name.aver_tegene
#all tair10 TEs
intersectBed -u -a CHH/$index.CHH.bed -b $TAIR10_TEs | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$sample_name.aver_te
cat $TAIR10_TEs | grep  LTR_Copia|  intersectBed -u -a CHH/$index.CHH.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$sample_name.aver_LTR_Copia
cat $TAIR10_TEs | grep  RC_Helitron|  intersectBed -u -a CHH/$index.CHH.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$sample_name.aver_RC_Helitron
cat $TAIR10_TEs | grep  DNA_MuDR|  intersectBed -u -a CHH/$index.CHH.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$sample_name.aver_DNA_MuDR
cat $TAIR10_TEs | grep  DNA|  intersectBed -u -a CHH/$index.CHH.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$sample_name.aver_allDNAtes
cat $TAIR10_TEs | grep  LTR_Gypsy|  intersectBed -u -a CHH/$index.CHH.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$sample_name.aver_LTR_Gypsy
cat $TAIR10_TEs | grep  LINE|  intersectBed -u -a CHH/$index.CHH.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$sample_name.aver_LINE
cat $TAIR10_TEs | grep  SINE|  intersectBed -u -a CHH/$index.CHH.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$sample_name.aver_SINE




