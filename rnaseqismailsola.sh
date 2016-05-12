 #! /usr/bin/env bash

 twoBittoFA mm10.2bit mm10.fa

 bowtie2-build -f mm10.fa mm10

 bowtie2 -x mm10 -1 C17ControlCGATGT_S1_L007_R1_001.fastq.gz \
 -2 C17ControlCGATGT_S1_L007_R2_001.fastq.gz \
 -s C17Control_Aligned.sam > Bowtie2Summary.txt

 bowtie2 -x mm10 -1 C17UpregTGACCA_S4_L007_R1_001.fastq.gz \
 -2 C17UpregTGACCA_S4_L007_R2_001.fastq.gz 
 -s C17Upreg_Aligned.sam > Bowtie2Summary1.txt

 bowtie2 -x mm10 -1 C17DeltaSetCCGTCC_S2_L007_R1_001.fastq.gz \
 -2 C17DeltaSetCCGTCC_S2_L007_R2_001.fastq.gz 
 -s C17DeltaSet_Aligned.sam > Bowtie2Summary2.txt

 samtools view -b -S C17Control_Aligned.sam > C17Control_Aligned.bam
 samtools view -b -S C17Upreg_Aligned.sam > C17Upreg_Aligned.bam
 samtools view -b -S C17DeltaSet_Aligned.sam > C17DeltaSet_Aligned.bam

 samtools sort C17Control_Aligned.bam C17Control_Aligned_sorted.bam
 samtools sort C17Upreg_Aligned.bam C17Upreg_Aligned_sorted.bam
 samtools sort C17DeltaSet_Aligned.C17 C17DeltaSet_Aligned_sorted.bam

 cuffdiff -o C17Upreg_FC/ mm10.refgene.gtf \
 C17Control_Aligned_sorted.bam C17Upreg_Aligned_sorted.bam 

 cuffdiff -o C17DeltaSet_FC/ mm10.refgene.gtf \
 C17Control_Aligned_sorted.bam C17DeltaSet_Aligned_sorted.bam

 cp C17Upreg_FC/gene_exp.diff C17Upreg.diff
 cp C17DeltaSet_FC/gene_exp.diff C17DeltaSet.diff

 cut -f2,6,9 C17Upreg.diff > C17Upreg_FC.diff
 cut -f2,6,9 C17DeltaSet.diff > C17DeltaSet_FC.diff

 grep -w "OK" C17Upreg_FC.diff > C17Upreg_FC_OK.diff
 grep -w "OK" C17DeltaSet_FC.diff > C17DeltaSet_FC_OK.diff

 join C17Upreg_FC_OK.diff C17DeltaSet_FC_OK.diff > Genes_FC.diff

 awk '{print $1 "\t" $3 "\t" $5 "\t" $5-$3}' > Genes_FC_Maxdiff.diff

 sort -K 2n Genes_FC_Maxdiff.diff > Genes_FC_Maxdiff_sorted.diff

 head Genes_FC_Maxdiff_sorted.diff > top_ten_genes.diff

BREAK

# R CODE

 top_ten_genes = read.table("top_ten_genes.diff")

 colnames(top_ten_genes.diff) <- c("gene", "FC_Upreg", "FC_DeltaSet")

 gene_1 <- df[1,]

 ggplot2(gene_1, aes(x=gene, y=FC_Upreg, fill=gene)) \
 + geom_bar(stat="identity") + theme_minimal()



