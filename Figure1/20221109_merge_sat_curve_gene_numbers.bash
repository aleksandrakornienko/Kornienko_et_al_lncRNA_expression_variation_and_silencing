# 1001G sat curve 
cd /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/saturation_curve/20211013_annot/


cat   satur_curve.400.accessions.rep.3/satur_curve.400.accessions.rep.3.gene_numbers.txt | head -1>  gene_numbers_1001G_sat_curve.bed


for i in {1..46}
do  
let "accN=$i*10"
echo $accN


for repN in {1..8}
do
replicate=rep.
replicate+=$repN
echo $replicate
export name_assembly=satur_curve.$accN.accessions.$replicate
echo $name_assembly
cat $name_assembly/$name_assembly.gene_numbers.txt| tail -1 >tmp2
echo $accN $replicate > tmp1
paste tmp1 tmp2 >tmp3
cat tmp3 >>  gene_numbers_1001G_sat_curve.bed
done 
done


satur_curve.90.accessions.rep.8.gene_numbers.txt

satur_curve.90.accessions.rep.8

#tissue saturation curve 
cd /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/saturation_curve/eracaps_sat_curve/merge_2021/

cat  accN.16.tissN.4.repN.3/accN.16.tissN.4.repN.3.gene_numbers.txt  | head -1>  gene_numbers_tissuesEC_sat_curve.bed


while read name_assembly 
do  
echo $name_assembly

cat $name_assembly/$name_assembly.gene_numbers.txt| tail -1 >tmp2
echo $name_assembly > tmp1
paste tmp1 tmp2 >tmp3
cat tmp3 >>  gene_numbers_tissuesEC_sat_curve.bed
done  < /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/saturation_curve/eracaps_sat_curve/names_of_points_eracaps_sat_curve.txt

cp gene_numbers_tissuesEC_sat_curve.bed > /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/02_Identification_Saturation_Curve/2021/


#cortijo 

cd /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/saturation_curve/cortijo

samples=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/02_Identification_Saturation_Curve/2021/cortijo.samples.sat.curve.txt

cat  sat_curve.3.pools.cortijo.rep.3/sat_curve.3.pools.cortijo.rep.3.gene_numbers.txt  | head -1>  gene_numbers_Cortijo_sat_curve.bed

while read name_assembly 
do  
echo $name_assembly
cat $name_assembly/$name_assembly.gene_numbers.txt| tail -1 >tmp2
echo $name_assembly > tmp1
paste tmp1 tmp2 >tmp3
cat tmp3 >>  gene_numbers_Cortijo_sat_curve.bed
done  < $samples

cp gene_numbers_Cortijo_sat_curve.bed /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/02_Identification_Saturation_Curve/2021/

