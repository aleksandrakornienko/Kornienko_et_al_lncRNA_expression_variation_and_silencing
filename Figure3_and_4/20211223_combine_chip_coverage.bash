

cd /groups/nordborg/projects/cegs/alexandra/ChIP-seq_2021/coverage

echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > cov_ar11pc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > cov_pc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand"> cov_linc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > cov_te_table
echo -e  "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > cov_as_table

echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TSS_cov_ar11pc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TSS_cov_pc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand"> TSS_cov_linc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TSS_cov_te_table
echo -e  "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TSS_cov_as_table

echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TES_cov_ar11pc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TES_cov_pc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand"> TES_cov_linc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TES_cov_te_table
echo -e  "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TES_cov_as_table



cat C1.rep1.H1.ar11_pc.log2.mean_cov.bed |  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >> cov_ar11pc_table
cat C1.rep1.H1.ar11_pc_TES.log2.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>TES_cov_ar11pc_table
cat C1.rep1.H1.ar11_pc_TSS.log2.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>TSS_cov_ar11pc_table

cat C1.rep1.H1.AS.log2.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>cov_as_table
cat C1.rep1.H1.as_TES.log2.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >> TES_cov_as_table
cat C1.rep1.H1.as_TSS.log2.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >> TSS_cov_as_table

cat C1.rep1.H1.linc.log2.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>cov_linc_table
cat C1.rep1.H1.linc_TES.log2.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >> TES_cov_linc_table
cat C1.rep1.H1.linc_TSS.log2.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >> TSS_cov_linc_table

cat C1.rep1.H1.PC.log2.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>cov_pc_table
cat C1.rep1.H1.pc_TES.log2.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>TES_cov_pc_table
cat C1.rep1.H1.pc_TSS.log2.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>TSS_cov_pc_table



cat C1.rep1.H1.te_genes.log2.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>cov_te_table
cat C1.rep1.H1.te_genes_TES.log2.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>TES_cov_te_table
cat C1.rep1.H1.te_genes_TSS.log2.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>TSS_cov_te_table


#log2

while read name 
do 
echo $name


echo $name> cov_ar11pc
echo $name> cov_pc
echo $name> cov_linc
echo $name> cov_te
echo $name> cov_as

echo $name> TSS_cov_ar11pc
echo $name> TSS_cov_pc
echo $name> TSS_cov_linc
echo $name> TSS_cov_te
echo $name> TSS_cov_as

echo $name> TES_cov_ar11pc
echo $name > TES_cov_pc
echo $name> TES_cov_linc
echo $name > TES_cov_te
echo $name > TES_cov_as


cat $name.ar11_pc.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>cov_ar11pc
cat $name.ar11_pc_TES.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TES_cov_ar11pc
cat $name.ar11_pc_TSS.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TSS_cov_ar11pc
cat $name.AS.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>cov_as
cat $name.as_TES.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TES_cov_as
cat $name.as_TSS.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TSS_cov_as
cat $name.linc.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>cov_linc
cat $name.linc_TES.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TES_cov_linc
cat $name.linc_TSS.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TSS_cov_linc
cat $name.PC.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>cov_pc
cat $name.pc_TES.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TES_cov_pc
cat $name.pc_TSS.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TSS_cov_pc
cat $name.te_genes.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>cov_te
cat $name.te_genes_TES.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TES_cov_te
cat $name.te_genes_TSS.log2.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TSS_cov_te


paste cov_pc_table cov_pc > inter 
cat inter > cov_pc_table
paste TES_cov_pc_table TES_cov_pc > inter 
cat inter > TES_cov_pc_table
paste TSS_cov_pc_table TSS_cov_pc > inter 
cat inter > TSS_cov_pc_table


paste cov_ar11pc_table cov_ar11pc > inter 
cat inter > cov_ar11pc_table
paste TES_cov_ar11pc_table TES_cov_ar11pc > inter 
cat inter > TES_cov_ar11pc_table
paste TSS_cov_ar11pc_table TSS_cov_ar11pc > inter 
cat inter > TSS_cov_ar11pc_table

paste cov_te_table cov_te > inter 
cat inter > cov_te_table
paste TES_cov_te_table TES_cov_te > inter 
cat inter > TES_cov_te_table
paste TSS_cov_te_table TSS_cov_te > inter 
cat inter > TSS_cov_te_table


paste cov_as_table cov_as > inter 
cat inter > cov_as_table
paste TES_cov_as_table TES_cov_as > inter 
cat inter > TES_cov_as_table
paste TSS_cov_as_table TSS_cov_as > inter 
cat inter > TSS_cov_as_table


paste cov_linc_table cov_linc > inter 
cat inter > cov_linc_table
paste TES_cov_linc_table TES_cov_linc > inter 
cat inter > TES_cov_linc_table
paste TSS_cov_linc_table TSS_cov_linc > inter 
cat inter > TSS_cov_linc_table




done < /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/ChIPsamples_to_process_names_only_antibodies.txt



cat cov_ar11pc_table> 20211223_Chipseq_coverage.ar11pc.log2.bed
cat cov_as_table> 20211223_Chipseq_coverage.as.log2.bed
cat cov_linc_table> 20211223_Chipseq_coverage.linc.log2.bed
cat cov_pc_table> 20211223_Chipseq_coverage.pc.log2.bed
cat cov_te_table> 20211223_Chipseq_coverage.te_genes.log2.bed
cat TES_cov_ar11pc_table> 20211223_Chipseq_coverage.ar11pc_TES.log2.bed
cat TES_cov_as_table> 20211223_Chipseq_coverage.as_TES.log2.bed
cat TES_cov_linc_table> 20211223_Chipseq_coverage.linc_TES.log2.bed
cat TES_cov_pc_table> 20211223_Chipseq_coverage.pc_TES.log2.bed
cat TES_cov_te_table> 20211223_Chipseq_coverage.te_genes_TES.log2.bed
cat TSS_cov_ar11pc_table> 20211223_Chipseq_coverage.ar11pc_TSS.log2.bed
cat TSS_cov_as_table> 20211223_Chipseq_coverage.as_TSS.log2.bed
cat TSS_cov_linc_table> 20211223_Chipseq_coverage.linc_TSS.log2.bed
cat TSS_cov_pc_table> 20211223_Chipseq_coverage.pc_TSS.log2.bed
cat TSS_cov_te_table> 20211223_Chipseq_coverage.te_genes_TSS.log2.bed





#subtract 





cd /groups/nordborg/projects/cegs/alexandra/ChIP-seq_2021/coverage

echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > cov_ar11pc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > cov_pc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand"> cov_linc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > cov_te_table
echo -e  "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > cov_as_table

echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TSS_cov_ar11pc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TSS_cov_pc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand"> TSS_cov_linc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TSS_cov_te_table
echo -e  "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TSS_cov_as_table

echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TES_cov_ar11pc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TES_cov_pc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand"> TES_cov_linc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TES_cov_te_table
echo -e  "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > TES_cov_as_table



cat C1.rep1.H1.ar11_pc.subtr.mean_cov.bed |  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >> cov_ar11pc_table
cat C1.rep1.H1.ar11_pc_TES.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>TES_cov_ar11pc_table
cat C1.rep1.H1.ar11_pc_TSS.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>TSS_cov_ar11pc_table

cat C1.rep1.H1.AS.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>cov_as_table
cat C1.rep1.H1.as_TES.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >> TES_cov_as_table
cat C1.rep1.H1.as_TSS.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >> TSS_cov_as_table

cat C1.rep1.H1.linc.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>cov_linc_table
cat C1.rep1.H1.linc_TES.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >> TES_cov_linc_table
cat C1.rep1.H1.linc_TSS.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >> TSS_cov_linc_table

cat C1.rep1.H1.PC.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>cov_pc_table
cat C1.rep1.H1.pc_TES.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>TES_cov_pc_table
cat C1.rep1.H1.pc_TSS.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>TSS_cov_pc_table



cat C1.rep1.H1.te_genes.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>cov_te_table
cat C1.rep1.H1.te_genes_TES.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>TES_cov_te_table
cat C1.rep1.H1.te_genes_TSS.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >>TSS_cov_te_table


#subtr

while read name 
do 
echo $name


echo $name> cov_ar11pc
echo $name> cov_pc
echo $name> cov_linc
echo $name> cov_te
echo $name> cov_as

echo $name> TSS_cov_ar11pc
echo $name> TSS_cov_pc
echo $name> TSS_cov_linc
echo $name> TSS_cov_te
echo $name> TSS_cov_as

echo $name> TES_cov_ar11pc
echo $name > TES_cov_pc
echo $name> TES_cov_linc
echo $name > TES_cov_te
echo $name > TES_cov_as


cat $name.ar11_pc.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>cov_ar11pc
cat $name.ar11_pc_TES.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TES_cov_ar11pc
cat $name.ar11_pc_TSS.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TSS_cov_ar11pc
cat $name.AS.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>cov_as
cat $name.as_TES.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TES_cov_as
cat $name.as_TSS.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TSS_cov_as
cat $name.linc.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>cov_linc
cat $name.linc_TES.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TES_cov_linc
cat $name.linc_TSS.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TSS_cov_linc
cat $name.PC.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>cov_pc
cat $name.pc_TES.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TES_cov_pc
cat $name.pc_TSS.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TSS_cov_pc
cat $name.te_genes.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>cov_te
cat $name.te_genes_TES.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TES_cov_te
cat $name.te_genes_TSS.subtr.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >>TSS_cov_te


paste cov_pc_table cov_pc > inter 
cat inter > cov_pc_table
paste TES_cov_pc_table TES_cov_pc > inter 
cat inter > TES_cov_pc_table
paste TSS_cov_pc_table TSS_cov_pc > inter 
cat inter > TSS_cov_pc_table


paste cov_ar11pc_table cov_ar11pc > inter 
cat inter > cov_ar11pc_table
paste TES_cov_ar11pc_table TES_cov_ar11pc > inter 
cat inter > TES_cov_ar11pc_table
paste TSS_cov_ar11pc_table TSS_cov_ar11pc > inter 
cat inter > TSS_cov_ar11pc_table

paste cov_te_table cov_te > inter 
cat inter > cov_te_table
paste TES_cov_te_table TES_cov_te > inter 
cat inter > TES_cov_te_table
paste TSS_cov_te_table TSS_cov_te > inter 
cat inter > TSS_cov_te_table


paste cov_as_table cov_as > inter 
cat inter > cov_as_table
paste TES_cov_as_table TES_cov_as > inter 
cat inter > TES_cov_as_table
paste TSS_cov_as_table TSS_cov_as > inter 
cat inter > TSS_cov_as_table


paste cov_linc_table cov_linc > inter 
cat inter > cov_linc_table
paste TES_cov_linc_table TES_cov_linc > inter 
cat inter > TES_cov_linc_table
paste TSS_cov_linc_table TSS_cov_linc > inter 
cat inter > TSS_cov_linc_table




done < /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/ChIPsamples_to_process_names_only_antibodies.txt



cat cov_ar11pc_table> 20211223_Chipseq_coverage.ar11pc.subtr.bed
cat cov_as_table> 20211223_Chipseq_coverage.as.subtr.bed
cat cov_linc_table> 20211223_Chipseq_coverage.linc.subtr.bed
cat cov_pc_table> 20211223_Chipseq_coverage.pc.subtr.bed
cat cov_te_table> 20211223_Chipseq_coverage.te_genes.subtr.bed
cat TES_cov_ar11pc_table> 20211223_Chipseq_coverage.ar11pc_TES.subtr.bed
cat TES_cov_as_table> 20211223_Chipseq_coverage.as_TES.subtr.bed
cat TES_cov_linc_table> 20211223_Chipseq_coverage.linc_TES.subtr.bed
cat TES_cov_pc_table> 20211223_Chipseq_coverage.pc_TES.subtr.bed
cat TES_cov_te_table> 20211223_Chipseq_coverage.te_genes_TES.subtr.bed
cat TSS_cov_ar11pc_table> 20211223_Chipseq_coverage.ar11pc_TSS.subtr.bed
cat TSS_cov_as_table> 20211223_Chipseq_coverage.as_TSS.subtr.bed
cat TSS_cov_linc_table> 20211223_Chipseq_coverage.linc_TSS.subtr.bed
cat TSS_cov_pc_table> 20211223_Chipseq_coverage.pc_TSS.subtr.bed
cat TSS_cov_te_table> 20211223_Chipseq_coverage.te_genes_TSS.subtr.bed





cp 20211223_Chipseq_coverage* /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/2021/coverage





















































cd /groups/nordborg/projects/cegs/alexandra/ChIP-seq/BAM/coverage

echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > cov_pc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand"> cov_linc_table
echo -e "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > cov_te_table
echo -e  "chr\t" "start\t" "end\t" "gene\t" "score\t" "strand" > cov_as_table

cat r8.rep2.H3K9me2.pc.2ndhalf.mean_cov.bed |  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >> cov_pc_table
cat r25.rep1.H3K9me2.as.PROM.log2.mean_cov.bed   |  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  >> cov_as_table

while read name 
do 
echo $name
echo $name > cov_pc
echo $name > cov_as

cat $name.pc.2ndhalf.mean_cov.bed|  awk -v OFS="\t" '{print $7}'  >> cov_pc
cat $name.as.PROM.log2.mean_cov.bed  |  awk -v OFS="\t" '{print $7}'  >> cov_as

paste cov_pc_table cov_pc > inter 
cat inter > cov_pc_table

paste cov_as_table cov_as > inter 
cat inter > cov_as_table

done < /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/chipseq_samples_with_ab.txt


cat cov_as_table > 20201006_cov_asProm_table.log2.bed
cat cov_pc_table > 20201006_cov_pc2ndHalf_table.log2.bed

cp 20201006* /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/



