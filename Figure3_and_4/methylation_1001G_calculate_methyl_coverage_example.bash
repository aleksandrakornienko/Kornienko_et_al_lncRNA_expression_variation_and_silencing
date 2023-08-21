#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=20G
#SBATCH --qos=rapid
#SBATCH --time=00:50:00
#SBATCH --array=1-444
#SBATCH --output=/groups/nordborg/pub/Mirjam/logs/array_job_slurm_%A_%a.out



export numberofline=$SLURM_ARRAY_TASK_ID


ml bedtools/2.27.1-foss-2018b

cd  /groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/methylation/1001G_data
#mkdir CG
#mkdir CHG
#mkdir CHH


export list_accessions_uniq=/groups/nordborg/projects/cegs/alexandra/2021_lncRNApaper/methylation/accessions_1001G_with_RNAseq_and_methyl_data.txt

export name_accession=`sed -n $numberofline,"$numberofline"p $list_accessions_uniq `


echo $name_accession
index=allc_
index+=$name_accession
##############################
echo $index
# combine all 5 chromosomes, take only CG methylation info, reformat the file into bed6-like file 
#cat  ${index}_"1".tsv   ${index}_"2".tsv    ${index}_"3".tsv ${index}_"4".tsv  ${index}_"5".tsv	| awk -v OFS="\t"   '{if (($4=="CG")) print "Chr"$1,$2,$2, $4, 0, $3, $5,$6,$7 }'  > CG/$index.CG.bed
# combine all 5 chromosomes, take only CHG methylation info, reformat the file into bed6-like file 
#cat  ${index}_"1".tsv   ${index}_"2".tsv    ${index}_"3".tsv ${index}_"4".tsv  ${index}_"5".tsv	| awk -v OFS="\t"   '{if ($4=="CAG" || $4=="CCG" || $4=="CTG") print "Chr"$1,$2,$2, $4, 0, $3, $5,$6,$7 }'  > CHG/$index.CHG.bed
#cat  ${index}_"1".tsv   ${index}_"2".tsv    ${index}_"3".tsv ${index}_"4".tsv  ${index}_"5".tsv	| awk -v OFS="\t"   '{if ($4=="CAA" || $4=="CAT" || $4=="CAC"|| $4=="CTA"|| $4=="CTT"|| $4=="CTC"|| $4=="CCA"|| $4=="CCT"|| $4=="CCC") print "Chr"$1,$2,$2, $4, 0, $3, $5,$6,$7 }'  > CHH/$index.CHH.bed

#rm  ${index}_"1".tsv   ${index}_"2".tsv    ${index}_"3".tsv ${index}_"4".tsv  ${index}_"5".tsv	



sortBed -i $denovo |bedtools map -a stdin -b CG/$index.CG.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CG/$index.CG.denovo.bed

sortBed -i $denovo |bedtools map -a stdin -b CHG/$index.CHG.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CHG/$index.CHG.denovo.bed

sortBed -i $denovo |bedtools map -a stdin -b CHH/$index.CHH.bed  -c 7,8 -o sum -null "22" |awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7/$8}' > CHH/$index.CHH.denovo.bed



# calculate average  methylation levels

#CG
cat CG/$index.CG.bed | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$index.aver  
echo $index > CG/$index.name 
#calculate average CG methylation 
#PC
intersectBed -u -a CG/$index.CG.bed -b $pc | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$index.aver_pc
#linc
intersectBed -u -a CG/$index.CG.bed -b $as | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$index.aver_linc
#AS
intersectBed -u -a CG/$index.CG.bed -b $linc | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$index.aver_as
#te genes
intersectBed -u -a CG/$index.CG.bed -b $te | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$index.aver_tegene

#all tair10 TEs
intersectBed -u -a CG/$index.CG.bed -b $TAIR10_TEs | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$index.aver_te
cat $TAIR10_TEs | grep  LTR_Copia|  intersectBed -u -a CG/$index.CG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$index.aver_LTR_Copia
cat $TAIR10_TEs | grep  RC_Helitron|  intersectBed -u -a CG/$index.CG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$index.aver_RC_Helitron
cat $TAIR10_TEs | grep  DNA_MuDR|  intersectBed -u -a CG/$index.CG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$index.aver_DNA_MuDR
cat $TAIR10_TEs | grep  DNA|  intersectBed -u -a CG/$index.CG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$index.aver_allDNAtes
cat $TAIR10_TEs | grep  LTR_Gypsy|  intersectBed -u -a CG/$index.CG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$index.aver_LTR_Gypsy
cat $TAIR10_TEs | grep  LINE|  intersectBed -u -a CG/$index.CG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$index.aver_LINE
cat $TAIR10_TEs | grep  SINE|  intersectBed -u -a CG/$index.CG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CG/$index.aver_SINE


#CHG
cat CHG/$index.CHG.bed | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$index.aver  
echo $index > CHG/$index.name 
#calculate average CHG methylation 
#PC
intersectBed -u -a CHG/$index.CHG.bed -b $pc | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$index.aver_pc
#linc
intersectBed -u -a CHG/$index.CHG.bed -b $linc | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$index.aver_linc
#AS
intersectBed -u -a CHG/$index.CHG.bed -b $as | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$index.aver_as
#te genes
intersectBed -u -a CHG/$index.CHG.bed -b $te | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$index.aver_tegene
#all tair10 TEs
intersectBed -u -a CHG/$index.CHG.bed -b $TAIR10_TEs | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$index.aver_te
cat $TAIR10_TEs | grep  LTR_Copia|  intersectBed -u -a CHG/$index.CHG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$index.aver_LTR_Copia
cat $TAIR10_TEs | grep  RC_Helitron|  intersectBed -u -a CHG/$index.CHG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$index.aver_RC_Helitron
cat $TAIR10_TEs | grep  DNA_MuDR|  intersectBed -u -a CHG/$index.CHG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$index.aver_DNA_MuDR
cat $TAIR10_TEs | grep  DNA|  intersectBed -u -a CHG/$index.CHG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$index.aver_allDNAtes
cat $TAIR10_TEs | grep  LTR_Gypsy|  intersectBed -u -a CHG/$index.CHG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$index.aver_LTR_Gypsy
cat $TAIR10_TEs | grep  LINE|  intersectBed -u -a CHG/$index.CHG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$index.aver_LINE
cat $TAIR10_TEs | grep  SINE|  intersectBed -u -a CHG/$index.CHG.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHG/$index.aver_SINE



#CHH
cat CHH/$index.CHH.bed | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$index.aver  
echo $index > CHH/$index.name 
#calculate average CHH methylation 
#PC
intersectBed -u -a CHH/$index.CHH.bed -b $pc | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$index.aver_pc
#linc
intersectBed -u -a CHH/$index.CHH.bed -b $linc | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$index.aver_linc
#AS
intersectBed -u -a CHH/$index.CHH.bed -b $as | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$index.aver_as
#te genes
intersectBed -u -a CHH/$index.CHH.bed -b $te | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$index.aver_tegene
#all tair10 TEs
intersectBed -u -a CHH/$index.CHH.bed -b $TAIR10_TEs | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$index.aver_te
cat $TAIR10_TEs | grep  LTR_Copia|  intersectBed -u -a CHH/$index.CHH.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$index.aver_LTR_Copia
cat $TAIR10_TEs | grep  RC_Helitron|  intersectBed -u -a CHH/$index.CHH.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$index.aver_RC_Helitron
cat $TAIR10_TEs | grep  DNA_MuDR|  intersectBed -u -a CHH/$index.CHH.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$index.aver_DNA_MuDR
cat $TAIR10_TEs | grep  DNA|  intersectBed -u -a CHH/$index.CHH.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$index.aver_allDNAtes
cat $TAIR10_TEs | grep  LTR_Gypsy|  intersectBed -u -a CHH/$index.CHH.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$index.aver_LTR_Gypsy
cat $TAIR10_TEs | grep  LINE|  intersectBed -u -a CHH/$index.CHH.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$index.aver_LINE
cat $TAIR10_TEs | grep  SINE|  intersectBed -u -a CHH/$index.CHH.bed -b stdin | awk '{ sum7 += $7;sum8 += $8; n++ } END { if (n > 0) print sum7/sum8; }' > CHH/$index.aver_SINE




