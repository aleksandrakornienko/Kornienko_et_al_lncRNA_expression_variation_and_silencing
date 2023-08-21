# look for TE pieces (TAIR10 TEs) in loci 
# identify copy number in TAIR10


ml clustalw2/2.1-foss-2018b
ml muscle/3.8.31-foss-2018b
ml blast+/2.8.1-foss-2018b
ml bedtools/2.27.1-foss-2018b


TAIR10=/groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/TAIR10_all.fa
TE_types_sorted=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/TAIR10_TEs_name_7superfam_length.sorted.txt

#prepare genome for blasting 
makeblastdb -in $TAIR10 -dbtype nucl 
    
dir=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation 

export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation
#pc from TAIR10 
export PC_tair=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.protein_coding.gtf
#pc annotation from Araport 
export PC_araport=$annotationfolder/Araport11_protein_coding.201606.bed
export TE_genes=$annotationfolder/Araport11_transposable_element_gene.201606.bed
#TE annotation 
export TE_fragments=$annotationfolder/Araport11_TEs.transposon_fragments.bed 
export TE_elements=$annotationfolder/Araport11_TEs.transposable_elements.bed
#araport nc RNAs 
export AraportNC=$annotationfolder/Araport11_non_coding.2016016.sorted.bed
export AraportNTR=$annotationfolder/Araport11_novel_transcribed_region.201606.bed

denovo_linc=$dir/lncRNAs.intergenic.loci.bed
denovo_linc_names=$dir/lncRNAs.intergenic.loci.names.txt
denovo_as=$dir/lncRNAs.antisense.loci.bed
denovo_as_names=$dir/lncRNAs.antisense.loci.names.txt
denovo_pc=$dir/denovoPC.loci.bed
denovo_pc_names=$dir/denovoPC.loci.names.txt
denovo_tegenes=$dir/loci.TE_genes.sorted.bed
denovo_tegenes_names=$dir/loci.TE_genes.names.txt

#make annotation of genes (loci) 
#cat $annotationfolder/Araport11_protein_coding.201606.bed | sortBed -i stdin| mergeBed -s -split -c 4,6 -o distinct -i stdin  | awk -v  OFS="\t" '{split($4,s,".");print $1,$2,$3,s[1],".",$5}' > $annotationfolder/Araport11_protein_coding.201606.genes.bed
#cat $annotationfolder/Araport11_protein_coding.201606.bed | sortBed -i stdin| mergeBed -s -split -c 4,6 -o distinct -i stdin  | awk -v  OFS="\t" '{split($4,s,".");print s[1]}' > $annotationfolder/Araport11_protein_coding.201606.genes_names.bed

araport_pc=$annotationfolder/Araport11_protein_coding.201606.genes.bed
araport_pc_names=$annotationfolder/Araport11_protein_coding.201606.genes_names.bed


#cat $denovo_linc | awk -v  OFS="\t" '{print $4}' > $denovo_linc_names
#wc -l  $denovo_linc_names
#2246 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/lncRNAs.intergenic.loci.names.txt

#cat $denovo_as | awk -v  OFS="\t" '{print $4}' > $denovo_as_names
#wc -l  $denovo_as_names
#8195 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/lncRNAs.antisense.loci.names.txt

#cat $denovo_pc | awk -v  OFS="\t" '{print $4}' > $denovo_pc_names
#wc -l  $denovo_pc_names
#23676 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/denovoPC.loci.names.txt
  
#cat $denovo_tegenes | awk -v  OFS="\t" '{print $4}' > $denovo_tegenes_names
#wc -l  $denovo_tegenes_names
#2130 /groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/20211013_annotation/loci.TE_genes.names.txt


##################################################################################
# the start of the actual script 
##################################################################################

genes=$denovo_linc
gene_names=$denovo_linc_names
BLASTfolder=/scratch-cbe/users/aleksandra.kornienko/BLAST/linc_strand 

#/scratch-cbe/users/aleksandra.kornienko/BLAST/pc_strand 
#/scratch-cbe/users/aleksandra.kornienko/BLAST/as_strand 
#/scratch-cbe/users/aleksandra.kornienko/BLAST/te_strand 
#/scratch-cbe/users/aleksandra.kornienko/BLAST/ar11_pc_strand 
 

#####
export locus=`sed -n $numberofline,"$numberofline"p $gene_names `
echo $locus
cd $BLASTfolder
   
mkdir $BLASTfolder/$locus

#1. extract the gene sequence in from the TAIR10 genome
#extract mRNA sequence (full gene)
cat $genes | grep -w $locus > $BLASTfolder/$locus/$locus.bed
length=`awk -v  OFS="\t" '{print $3-$2}' $BLASTfolder/$locus/$locus.bed`
bedtools getfasta   -fi $TAIR10 -bed $BLASTfolder/$locus/$locus.bed  -s  -name  > $BLASTfolder/$locus/$locus.fa
#promoter (TSS +/- 200bp)
cat $BLASTfolder/$locus/$locus.bed |  awk -v OFS="\t" '{ print $1,$2-200,$2+200,$4,$5,$6}'|awk -v OFS="\t" '{ if ($2<0) {$2=0}; {print $0}}'> $BLASTfolder/$locus/$locus.TSS_plusminus200bp.TAIR10.bed
bedtools getfasta -s    -name   -fi $TAIR10 -bed $BLASTfolder/$locus/$locus.TSS_plusminus200bp.TAIR10.bed > $BLASTfolder/$locus/$locus.TSS_plusminus200bp.TAIR10.fa
#TES +/- 200bp
cat $BLASTfolder/$locus/$locus.bed |  awk -v OFS="\t" '{ print $1,$3-200,$3+200,$4,$5,$6}'|awk -v OFS="\t" '{ if ($2<0) {$2=0}; {print $0}}'> $BLASTfolder/$locus/$locus.TES_plusminus200bp.TAIR10.bed
bedtools getfasta -s    -name   -fi $TAIR10 -bed $BLASTfolder/$locus/$locus.TES_plusminus200bp.TAIR10.bed > $BLASTfolder/$locus/$locus.TES_plusminus200bp.TAIR10.fa

locusseq_tair10=$BLASTfolder/$locus/$locus.fa


#copy number in TAIR10 
blastn -task blastn -word_size 10 -query  $locusseq_tair10 -db $TAIR10 -strand both -out $BLASTfolder/$locus/$locus.blast.TAIR10.genome.txt -outfmt 7 -evalue 1e-7
# blast fields 
#1 query (lincRNA) 
#2 subject (genome (chromosome))
#3 % identity, 
#4 alignment length,
#5  mismatches, 
#6 gap opens, 
#7 query start, (start within lincRNA locus)
#8 query end, (end within lincRNA locus)
#9 subject start, (chromosome coordinates start)
#10 subject end, (chromosome coordinates end)
#11 evalue, 
#12 bit score

accession=TAIR10
copies=`cat $BLASTfolder/$locus/$locus.blast.TAIR10.genome.txt | grep CUFF | grep -v Query  | awk -v  OFS="\t" '($3>80){print $0,$8-$7}'| awk -v acc="$accession" -v pc="$locus" -v  OFS="\t" '{if ($9<$10) print $2, $9,$10,"besthit."pc".GNM."acc,"0","+",$3*$4,$4,$3; else  print $2, $10,$9,"besthit."pc".GNM."acc,"0","-",$3*$4,$4,$3}'|sortBed -i stdin|  mergeBed  -s -i stdin -d 1500 -c 6,7,8,8,9 -o first,sum,sum,max,max |  awk -v  OFS="\t" '{print $0,$3-$2}'|  awk -v len="$length"  -v  OFS="\t" '($6>0.8*len){print $0}'| sort -k8,8nr| awk -v acc="$accession" -v pc="$locus" -v  OFS="\t" '{if ($6>($3-$2+1)) print $1,$2,$3,pc".GNM."acc".copy"NR,$8,$4; else print $1,$2,$3,pc".GNM."acc".copy"NR,$5/$6,$4 }'  | wc -l`  
echo $copies > $locus.CN_TAIR10
 
cat $BLASTfolder/$locus/$locus.blast.TAIR10.genome.txt | grep CUFF | grep -v Query  | awk -v  OFS="\t" '($3>80){print $0,$8-$7}'| awk -v acc="$accession" -v pc="$locus" -v  OFS="\t" '{if ($9<$10) print $2, $9,$10,"besthit."pc".GNM."acc,"0","+",$3*$4,$4,$3; else  print $2, $10,$9,"besthit."pc".GNM."acc,"0","-",$3*$4,$4,$3}'|sortBed -i stdin|  mergeBed  -s -i stdin -d 1500 -c 6,7,8,8,9 -o first,sum,sum,max,max |  awk -v  OFS="\t" '{print $0,$3-$2}'|  awk -v len="$length"  -v  OFS="\t" '($6>0.8*len){print $0}'| sort -k8,8nr| awk -v acc="$accession" -v pc="$locus" -v  OFS="\t" '{if ($6>($3-$2+1)) print $1,$2,$3,pc".GNM."acc".copy"NR,$8,$4; else print $1,$2,$3,pc".GNM."acc".copy"NR,$5/$6,$4 }'> $BLASTfolder/$locus/$locus.blast.TAIR10.copies.bed 



# find TE pieces 

cat $BLASTfolder/$locus/$locus.fa $BLASTfolder/$locus/$locus.blast.allgenomes.fa > $BLASTfolder/$locus/$locus.blast.allgenomes_and_TAIR.fa

cat $BLASTfolder/$locus/$locus.TSS_plusminus200bp.TAIR10.fa $BLASTfolder/$locus/$locus.TSS_plusminus200bp.besthit.allgenomes.fa >$BLASTfolder/$locus/$locus.TSS_plusminus200bp.besthit.allgenomes_and_TAIR.fa

cat $BLASTfolder/$locus/$locus.TES_plusminus200bp.TAIR10.fa $BLASTfolder/$locus/$locus.TES_plusminus200bp.besthit.allgenomes.fa >$BLASTfolder/$locus/$locus.TES_plusminus200bp.besthit.allgenomes_and_TAIR.fa

makeblastdb -in $BLASTfolder/$locus/$locus.blast.allgenomes_and_TAIR.fa -dbtype nucl  
makeblastdb -in $BLASTfolder/$locus/$locus.TSS_plusminus200bp.besthit.allgenomes_and_TAIR.fa -dbtype nucl  
makeblastdb -in $BLASTfolder/$locus/$locus.TES_plusminus200bp.besthit.allgenomes_and_TAIR.fa -dbtype nucl  

TEs_TAIR10seq=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/TAIR10_TEs.fa

############################
#find TEs in the whole locus
############################ 
blastn -task blastn -word_size 10 -query  $TEs_TAIR10seq -db $BLASTfolder/$locus/$locus.blast.allgenomes_and_TAIR.fa -strand both -out $BLASTfolder/$locus/blast.$locus.TEinserts_each_genome.1.txt -outfmt 7 -evalue 1e-7
# blast fields 
#1 query (TE) 
#2 subject (lincRNA)
#3 % identity, 
#4 alignment length,
#5  mismatches, 
#6 gap opens, 
#7 query start, (start within TE )
#8 query end, (end within TE )
#9 subject start, (end within lincRNA locus)
#10 subject end, (end within lincRNA locus)
#11 evalue, 
#12 bit score


#require >80 sequence identity , but no restriction on length
cat $BLASTfolder/$locus/blast.$locus.TEinserts_each_genome.1.txt | grep AT | grep -v Query |grep -v Database | awk -v  OFS="\t" '($3>80){print $0}'| awk -v acc="$accession" -v pc="$locus" -v  OFS="\t" '{if ( ($10-$9)>0) print $2, $9,$10,$1,$3,"+",$4; else print $2, $10,$9,$1,$3,"-",$4 }' > $BLASTfolder/$locus/blast.$locus.TEinserts_each_genome.2.txt
#add TE type and length
length=`awk -v  OFS="\t" '{print $3-$2}' $BLASTfolder/$locus/$locus.bed` 

cat $BLASTfolder/$locus/blast.$locus.TEinserts_each_genome.2.txt |  sort -k4  > $BLASTfolder/$locus/tmp1
join -1 4 -2 1 $BLASTfolder/$locus/tmp1  $TE_types_sorted|  awk -v len="$length" -v  OFS="\t" '{print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$7/len,$5*$7}' |sortBed -i stdin  > $BLASTfolder/$locus/blast.$locus.TEinserts_each_genome.bed

#columns in $BLASTfolder/$locus/blast.$locus.TEinserts_each_genome.bed:
#1 lincRNA name (besthit name with strand in the corresponding genome)
#2 start of the TE alignment within the lincRNA 
#3 end of the TE alignment within the lincRNA 
#4 TE name (TAIR10 TEs)
#5 % identity (sequnce similarity)
#6 strand of the TE alignment relative to the lincRNA ( +: TE is in the same direction as the lincRNA gene, -: TE sequence is antisence to the lincRNA gene)
#7 alignment length
#8 TE family (7 superfamilies: DNA_other, DNA_MuDR, SINE_LINE, RC_Helitron, LTR_Gypsy, LTR_Copia, Unassigned_NA )
#9 TE length (TAIR10 TEs)
#10 proportion of the lincRNA length covered by this TE (in the sense or antisense direction)
#11 alignment length * % identity 
 
#merge overlapping TE pieces and make a bed file with distinct TE types instead of TE names
cat $BLASTfolder/$locus/blast.$locus.TEinserts_each_genome.bed | mergeBed -s -i stdin -c 6,8 -o first,distinct| awk  -v len="$length" -v  OFS="\t" '{print $1,$2,$3,$5,0,$4}' > $BLASTfolder/$locus/blast.$locus.TEinserts_each_genome.position_on_gene.merged.bed

#columns in $BLASTfolder/$locus/blast.$locus.TEinserts_each_genome.position_on_gene.merged.bed - a proper bed6 file:
#1 lincRNA name (besthit name with strand in the corresponding genome)
#2 start of the merged TE insertion within the lincRNA 
#3 end of the merged TE insertion within the lincRNA 
#4 distinct TE families within the merged TE insertion (1 family - the all the overlapping aligned TEs are of the same type)  (7 superfamilies: DNA_other, DNA_MuDR, SINE_LINE, RC_Helitron, LTR_Gypsy, LTR_Copia, Unassigned_NA )
#5 portion of the lincRNA length covered by this TE piece  
#6 strand (relative to the lincRNA strand )


#remove unneeded files
rm $BLASTfolder/$locus/blast.$locus.TEinserts_each_genome.1.txt $BLASTfolder/$locus/blast.$locus.TEinserts_each_genome.2.txt 


############################
#find TEs in the promoter (TSS +- 200bp)
############################ 
blastn -task blastn -word_size 10 -query  $TEs_TAIR10seq -db $BLASTfolder/$locus/$locus.TSS_plusminus200bp.besthit.allgenomes_and_TAIR.fa -strand both -out $BLASTfolder/$locus/blast.$locus.TSS_plusminus200bp.TEinserts_each_genome.1.txt -outfmt 7 -evalue 1e-7

#require >80 sequence identity , but no restriction on length
cat $BLASTfolder/$locus/blast.$locus.TSS_plusminus200bp.TEinserts_each_genome.1.txt | grep AT | grep -v Query |grep -v Database | awk -v  OFS="\t" '($3>80){print $0}'| awk -v acc="$accession" -v pc="$locus" -v  OFS="\t" '{if ( ($10-$9)>0) print $2, $9,$10,$1,$3,"+",$4; else print $2, $10,$9,$1,$3,"-",$4 }' > $BLASTfolder/$locus/blast.$locus.TSS_plusminus200bp.TEinserts_each_genome.2.txt
#add TE type and length
cat $BLASTfolder/$locus/blast.$locus.TSS_plusminus200bp.TEinserts_each_genome.2.txt |  sort -k4  > $BLASTfolder/$locus/tmp1
join -1 4 -2 1 $BLASTfolder/$locus/tmp1  $TE_types_sorted|  awk -v len=400 -v  OFS="\t" '{print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$7/len,$5*$7}' |sortBed -i stdin  > $BLASTfolder/$locus/blast.$locus.TSS_plusminus200bp.TEinserts_each_genome.bed
 
#merge overlapping TE pieces and make a bed file with distinct TE types instead of TE names
cat $BLASTfolder/$locus/blast.$locus.TSS_plusminus200bp.TEinserts_each_genome.bed | mergeBed -s -i stdin -c 6,8 -o first,distinct| awk  -v  OFS="\t" '{print $1,$2,$3,$5,0,$4}' > $BLASTfolder/$locus/blast.$locus.TSS_plusminus200bp.TEinserts_each_genome.position_on_gene.merged.bed

#remove unneeded files
rm $BLASTfolder/$locus/blast.$locus.TSS_plusminus200bp.TEinserts_each_genome.1.txt $BLASTfolder/$locus/blast.$locus.TSS_plusminus200bp.TEinserts_each_genome.2.txt 



############################
#find TEs in the end of the gene (TES +- 200bp)
############################ 
blastn -task blastn -word_size 10 -query  $TEs_TAIR10seq -db $BLASTfolder/$locus/$locus.TES_plusminus200bp.besthit.allgenomes_and_TAIR.fa -strand both -out $BLASTfolder/$locus/blast.$locus.TES_plusminus200bp.TEinserts_each_genome.1.txt -outfmt 7 -evalue 1e-7

#require >80 sequence identity , but no restriction on length
cat $BLASTfolder/$locus/blast.$locus.TES_plusminus200bp.TEinserts_each_genome.1.txt | grep AT | grep -v Query |grep -v Database | awk -v  OFS="\t" '($3>80){print $0}'| awk -v acc="$accession" -v pc="$locus" -v  OFS="\t" '{if ( ($10-$9)>0) print $2, $9,$10,$1,$3,"+",$4; else print $2, $10,$9,$1,$3,"-",$4 }' > $BLASTfolder/$locus/blast.$locus.TES_plusminus200bp.TEinserts_each_genome.2.txt
#add TE type and length
cat $BLASTfolder/$locus/blast.$locus.TES_plusminus200bp.TEinserts_each_genome.2.txt |  sort -k4  > $BLASTfolder/$locus/tmp1
join -1 4 -2 1 $BLASTfolder/$locus/tmp1  $TE_types_sorted|  awk -v len=400 -v  OFS="\t" '{print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$7/len,$5*$7}' |sortBed -i stdin  > $BLASTfolder/$locus/blast.$locus.TES_plusminus200bp.TEinserts_each_genome.bed
 
#merge overlapping TE pieces and make a bed file with distinct TE types instead of TE names
cat $BLASTfolder/$locus/blast.$locus.TES_plusminus200bp.TEinserts_each_genome.bed | mergeBed -s -i stdin -c 6,8 -o first,distinct| awk  -v  OFS="\t" '{print $1,$2,$3,$5,0,$4}' > $BLASTfolder/$locus/blast.$locus.TES_plusminus200bp.TEinserts_each_genome.position_on_gene.merged.bed

#remove unneeded files
rm $BLASTfolder/$locus/blast.$locus.TES_plusminus200bp.TEinserts_each_genome.1.txt $BLASTfolder/$locus/blast.$locus.TES_plusminus200bp.TEinserts_each_genome.2.txt 


