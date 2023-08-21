# my de novo annotation 

export name_assembly=20210513_annotation
export  working_folder=/groups/nordborg/projects/cegs/alexandra/2020_lncRNA_paper/01_lncRNA_identification/tr_assembly/merge/$name_assembly
 





###################
# 3. splitting annotation into PC and non-coding part 
###################
#use both Araport and TAIR10 for PC annotation 
export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation
#pc from TAIR10 
export PC_tair=$annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.protein_coding.gtf
sortBed -i $annotationfolder/TAIR10/Arabidopsis_thaliana.TAIR10.40.protein_coding.gtf|mergeBed -s -i stdin | wc -l
#27103 PC loci
#pc annotation from Araport 
export PC_araport=$annotationfolder/Araport11_protein_coding.201606.bed
sortBed -i $annotationfolder/Araport11_protein_coding.201606.bed|mergeBed -s -i stdin | wc -l
#27129 loci
#take pseudogene annotation from araport 
export pseudogenes=$annotationfolder/Araport11_pseudogene.201606.bed
sortBed -i  $annotationfolder/Araport11_pseudogene.201606.bed |mergeBed -s -i stdin |  wc -l
#953 pseudogene loci
#take transposon annotation from araport 
export TEs=$annotationfolder/Araport11_transposable_element_gene.201606.bed
mergeBed -s -i  $annotationfolder/Araport11_transposable_element_gene.201606.bed | wc -l
#3897 TE gene loci



#araport nc 
NC_araport=$annotationfolder/Araport11_non_coding.2016016.sorted.bed 

# 1. remove transcripts with mRNA length <200nt
cat $annotationfolder/Araport11_non_coding.2016016.sorted.bed  | awk '{split($11,s,","); total=0; for(x in s){ total+=s[x];} print total,$0}'  | awk -v OFS="\t" '$1>=200 {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13 }'   > $annotationfolder/Araport11.lncRNAs.bed 
 
 wc -l $annotationfolder/Araport11.lncRNAs.bed
4143 /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11.lncRNAs.bed

$annotationfolder/Araport11.lncRNAs.bed 

mergeBed -s -i  $annotationfolder/Araport11.lncRNAs.bed | wc -l
#3685 lncRNA loci

#araport novel transcribed regionbs 
NTR_araport=$annotationfolder/Araport11_novel_transcribed_region.201606.bed
mergeBed -s -i $annotationfolder/Araport11_novel_transcribed_region.201606.bed | wc -l
# 507 NTR loci

############# check novelty of lncRNAs again 

################################

#known lncRNAs (present in Araport11 ncRNA annotation) 
intersectBed -s -u -a $working_folder/denovoNC.transcripts.bed -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed > $working_folder/denovoNC.transcripts.known.bed 
intersectBed -s -u -a $working_folder/denovoNC.loci.bed -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed > $working_folder/denovoNC.loci.known.bed 

#known lncRNAs (present in Araport11 ncRNA and novel transcribed regions annotation)

intersectBed -s -u -a $working_folder/denovoNC.transcripts.bed -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed $annotationfolder/Araport11_novel_transcribed_region.201606.bed > $working_folder/denovoNC.transcripts.known_NTR_and_NC.bed 
intersectBed -s -u -a $working_folder/denovoNC.loci.bed -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed $annotationfolder/Araport11_novel_transcribed_region.201606.bed > $working_folder/denovoNC.loci.known_NTR_and_NC.bed 
 
 wc -l *known*
   797 denovoNC.loci.known.bed
  1000 denovoNC.loci.known_NTR_and_NC.bed
  2211 denovoNC.transcripts.known.bed
  2703 denovoNC.transcripts.known_NTR_and_NC.bed

  
 
 
#novel lncRNAs (not present in Araport11 ncRNA and novel transcribed regions annotation)
intersectBed -s -v -a $working_folder/denovoNC.transcripts.bed -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed > $working_folder/denovoNC.transcripts.novel.bed 
intersectBed -s -v -a $working_folder/denovoNC.loci.bed -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed > $working_folder/denovoNC.loci.novel.bed 

intersectBed -s -v -a $working_folder/denovoNC.transcripts.bed -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed $annotationfolder/Araport11_novel_transcribed_region.201606.bed > $working_folder/denovoNC.transcripts.novel_NTR_and_NC.bed 
intersectBed -s -v -a $working_folder/denovoNC.loci.bed -b $annotationfolder/Araport11_non_coding.2016016.sorted.bed $annotationfolder/Araport11_novel_transcribed_region.201606.bed > $working_folder/denovoNC.loci.novel_NTR_and_NC.bed 


  wc -l *novel*
  11095 denovoNC.loci.novel.bed
  10892 denovoNC.loci.novel_NTR_and_NC.bed
  14416 denovoNC.transcripts.novel.bed
  13924 denovoNC.transcripts.novel_NTR_and_NC.bed

  
  
  