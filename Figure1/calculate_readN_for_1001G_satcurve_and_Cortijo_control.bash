

#cortijo 
export  working_folder=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/saturation_curve/cortijo
#mkdir $working_folder

samples=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/02_Identification_Saturation_Curve/2021/cortijo.samples.sat.curve.txt
pool_readN=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/02_Identification_Saturation_Curve/2021/Cortijo_alignment2021_pool_read_N.txt

echo "sample" "readN" >$working_folder/sat_curve_samples_Cortijocontrol_readN.bed

while read name_assembly
do 
cat $working_folder/$name_assembly/$name_assembly.assemblies_for_merge.txt | awk '{split($1,s,"Cortijo.");  print s[2]}' |  awk '{split($1,s,"/");  print s[1]}' >tmp 
readN=`cat $pool_readN | grep -f tmp |awk -v  OFS="\t" '{sum+=$2;} END{print sum;}' `
echo $name_assembly $readN 
echo $name_assembly $readN>> $working_folder/sat_curve_samples_Cortijocontrol_readN.bed
done <$samples

cp $working_folder/sat_curve_samples_Cortijocontrol_readN.bed /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/02_Identification_Saturation_Curve/2021/



#1001G

cd /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/saturation_curve/shuffled_accessions


#cat rep.1/10.accessions.txt | awk -v  OFS="\t" '{sum+=$2;} END{print sum;}'

echo "sample" "readN" > satur_curve_1001G_points_sum_uniq_readN.bed

for i in {1..46}
do
let "accN=$i*10"
for j in {1..8}
do
readN=`cat rep.$j/$accN.accessions.txt | awk -v  OFS="\t" '{sum+=$2;} END{print sum;}' `
echo $accN.$j $readN 
echo $accN.$j  $readN>> satur_curve_1001G_points_sum_uniq_readN.bed
done 
done 


cp satur_curve_1001G_points_sum_uniq_readN.bed /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/02_Identification_Saturation_Curve/2021/


# eracaps
cd /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/saturation_curve/eracaps_sat_curve/
readNtable=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/000_data_preparation/2021/ERACAPS_uniqreadN.txt 
assemblies_list_folder=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/saturation_curve/eracaps_sat_curve/list_of_assemblies_for_merge
export points_list=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/saturation_curve/eracaps_sat_curve/names_of_points_eracaps_sat_curve.txt

while read point 
do 
echo $point 
cat $assemblies_list_folder/$point.assemblies_for_merge.txt | awk '{split($1,s,"stringtie_2021/");  print s[2]}' |  awk '{split($1,s,"/");  print s[1]}' >tmp 
readN=`cat $readNtable | grep -f tmp |awk -v  OFS="\t" '{sum+=$2;} END{print sum;}' `
echo $point $readN 
echo $point $readN>> sat_curve_samples_ERACAPS_readN.bed
done < $points_list

cp sat_curve_samples_ERACAPS_readN.bed /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/2018_lncRNA_variation_paper/01_lncRNA_identification/02_Identification_Saturation_Curve/2021/


