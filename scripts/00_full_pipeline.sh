#!/bin/bash

conda activate isoseq-apa

#Step 0: Preparation
python external_scripts/extract_transcript_regions.py  -i ref/hg19.refGene.gtf.gz -o ref/hg19 --gtf

#Step 1: Generate CCS
ccs m54304_190611_190806.subreads.bam m54304_190611_190806.ccs.bam --noPolish --minPasses 1 --numThreads 16
ccs m54304_190613_163530.subreads.bam m54304_190613_163530.ccs.bam --noPolish --minPasses 1 --numThreads 16
ccs m54304_190614_125656.subreads.bam m54304_190614_125656.ccs.bam --noPolish --minPasses 1 --numThreads 16
ccs m54304_190615_091702.subreads.bam m54304_190615_091702.ccs.bam --noPolish --minPasses 1 --numThreads 16

#Step 2: Step 2: Primer Removal & Demux Bams
lima --isoseq --dump-clips --peek-guess -j 24 m54304_190611_190806.ccs.bam data/YG-ISO-SEQ-Barcodes.fasta demux.m54304_190611_190806.fl.bam

#Step 3: Convert to Fastqs
for i in 2 4 6 8 10
do
bamtools convert -format fastq -in demux.5p--YG_GM$i\_3p.bam -out demux.5p--YG_GM$i_3p.bam.fastq
done

#Step 4: Alignment
for i in 2 4 6 8 10
do
minimap2  -ax splice -uf -C5 ref/hg19.fa demux.5p--YG_GM$i_3p.bam.fastq > demux.5p--YG_GM$i_3p.bam
samtools sort demux.5p--YG_GM$i_3p.bam > demux.5p--YG_GM$i_3p.sort.bam 
samtools index demux.5p--YG_GM$i_3p.sort.bam 
done

samtools merge total.merge.bam demux.5p--YG_GM2_3p.sort.bam demux.5p--YG_GM4_3p.sort.bam demux.5p--YG_GM6_3p.sort.bam demux.5p--YG_GM8_3p.sort.bam demux.5p--YG_GM10_3p.sort.bam
samtools sort total.merge.bam total.merge.sort.bam
samtools index total.merge.sort.bam

#Step 6: Define PolyA Reads
python scripts/01_filter_for_polyA.py -i total.merge.sort.bam -o test -f ref/hg19.fa 

#Step 7: Define PAS Sites
cut -f1-6 ref/hg19_3utr.bed > ref/hg19_3utr.1-6.bed
cut -f1-6 ref/hg19_exons.bed > ref/hg19_exons.1-6.bed
bedtools intersect -b ref/hg19_3utr.1-6.bed -a test.noMP.bed -s -wa -wb > test.noMP.bed.UTR
bedtools intersect -b ref/hg19_exons.1-6.bed -a test.noMP.bed.UTR -s -wa -wb > test.noMP.bed.UTR.EXONS
bedtools intersect -b test.noMP.bed.UTR -a test.noMP.bed -v -wa > test.noMP.bed.INTRON
python scripts/02_noMP_exon_UTR.py -i test.noMP.meta.txt -p test

#Step 8: Annotate
sort -k1,1 -k2,2n test.noMP.restricted.peaks.bed > test.noMP.restricted.peaks.sort.bed
bedtools merge -i test.noMP.restricted.peaks.sort.bed -c 4,5,6 -o collapse,mean,distinct > test.noMP.restricted3.peaks.sort.refined.bed
python scripts/03_final_peaks.py  -p test 

annotatePeaks.pl test.noMP.restricted.peaks.sort.refined.score.bed hg19 > test.homer.annotate.BED
python scripts/04_annotate_usage.py -p test -r ref/Refseq2Gene.txt -a 

sort -k1,1 -k2,2n test.FINAL.BED | uniq >  test.FINAL.UNIQ.BED
python scripts/04_annotate_usage.py -p test -r ref/Refseq2Gene.txt -u







