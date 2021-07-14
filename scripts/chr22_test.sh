#!/bin/bash

# Step 0: Preparation

# Load conda environment
conda activate iso-seq-apa

# Steps 1-45 have already been completed for this example. You will be working with merged Pacbio Iso-Seq data from 8 individuals, restricted to chr22 

# Step 6: Define PolyA Reads 

python scripts/01_filter_for_polyA.py -i data/chr22.bam -o test -f ref/hg19.fa 

# Step 7: Define PAS Sites
bedtools  bamtobed -i test.noMP.bam > test.noMP.bed
bedtools intersect -b ref/hg19_3utr.bed -a test.noMP.bed -s -wa -wb > test.noMP.bed.UTR
bedtools intersect -b ref/hg19_exons.bed -a test.noMP.bed.UTR -s -wa -wb > test.noMP.bed.UTR.EXONS
bedtools intersect -b test.noMP.bed.UTR -a test.noMP.bed -v -wa > test.noMP.bed.INTRON
python scripts/02_noMP_exon_UTR.py -i test.noMP.meta.txt -p test

# Step 8: Annotate

sort -k1,1 -k2,2n test.noMP.restricted.peaks.bed > test.noMP.restricted.peaks.sort.bed
bedtools merge -i test.noMP.restricted.peaks.sort.bed -c 4,5,6 -o collapse,mean,distinct > test.noMP.restricted.peaks.sort.refined.bed
python scripts/03_final_peaks.py  -p test 



/home/ankeetashah/homer/.//bin/annotatePeaks.pl test.noMP.restricted.peaks.sort.refined.score.bed hg19 > test.homer.annotate.BED
python scripts/04_annotate_usage.py -p test -r ref/Refseq2Gene.txt -a annotate 

sort -k1,1 -k2,2n test.FINAL.bed | uniq >  test.FINAL.UNIQ.bed
python scripts/04_annotate_usage.py -p test -r ref/Refseq2Gene.txt -u usage



