# Polyadenylation site identification using Iso-Seq data

## Step 0: Preparation

- Download [hg19.refGene.gtf.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz) annotations. Move it to the (ref/) (https://github.com/ankeetashah/Benchmarking-APA/tree/main/ref) directory (e.g. ```mv hg19.refGene.gtf.gz ref/.```)
- Download [hg19.fa.gz](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz). Move it to the (ref/) (https://github.com/ankeetashah/Benchmarking-APA/tree/main/ref) directory (e.g. ```mv hg19.fa.gz ref/.```)
- Download the [extract_transcript_regions.py](https://github.com/stephenfloor/extract-transcript-regions) script. We have also included it and its dependencies in the the [external_scripts](https://github.com/ankeetashah/Benchmarking-APA/tree/main/external_scripts) directory.
- Install [Homer](http://homer.ucsd.edu/homer/introduction/install.html). In particular, you will need to be able to run the ```annotatePeaks.pl``` script. 

Create a conda environment, and install the following packages:
```bash
conda create -n isoseq-apa python=2.7
conda activate isoseq-apa
conda install -c bioconda pbccs 
conda install -c bioconda lima 
conda install -c bioconda bamtools 
conda install -c bioconda samtools
conda install -c bioconda bedtools 
conda install -c bioconda pysam
conda install -c bioconda pyfaidx
conda install -c bioconda minimap2
```

The apprpriate output bed files (specifically for 3'UTR and exon annotations) are in the (ref/)(https://github.com/ankeetashah/Benchmarking-APA/tree/main/ref). The following command was used. If one wishes to work with a different reference genome, this step can be re-run accordingly.

```bash
python external_scripts/extract_transcript_regions.py  -i ref/hg19.refGene.gtf.gz -o ref/hg19 --gtf
```

## Step 1: Generate CCS 

For each SMRT cell you will have a ```X.subreads.bam file```, such as ```m54304_190611_190806.subreads.bam```.

Sequencing runs are processed to generate one representative circular consensus sequence (CCS) for every ZMW using ```ccs```. Ensure that you include the ```--noPolish``` flag, otherwise polyA tails from these reads will be clipped.

```
ccs m54304_190611_190806.subreads.bam m54304_190611_190806.ccs.bam --noPolish --minPasses 1 --numThreads 16
ccs m54304_190613_163530.subreads.bam m54304_190613_163530.ccs.bam --noPolish --minPasses 1 --numThreads 16
ccs m54304_190614_125656.subreads.bam m54304_190614_125656.ccs.bam --noPolish --minPasses 1 --numThreads 16
ccs m54304_190615_091702.subreads.bam m54304_190615_091702.ccs.bam --noPolish --minPasses 1 --numThreads 16
```

## Step 2: Primer Removal & Demux Bams

Removal of primers and identification of barcodes is performed using ```lima```. In this case, barcodes file can be found in the [data/](https://github.com/ankeetashah/Benchmarking-APA/tree/main/data) directory. As an example:

```
lima --isoseq --dump-clips --peek-guess -j 24 m54304_190611_190806.ccs.bam data/YG-ISO-SEQ-Barcodes.fasta demux.m54304_190611_190806.fl.bam
```

The barcodes file looks like:
```
cat data/YG-ISO-SEQ-Barcodes.fasta
>5p
AAGCAGTGGTATCAACGCAGAGTACATGGG
>YG-GM1_3p
ctcacagtctgtgtgtGTACTCTGCGTTGATACCACTGCTT
>YG-GM2_3p
ctctcacgagatgtgtGTACTCTGCGTTGATACCACTGCTT
>YG-GM3_3p
cgcgcgtgtgtgcgtgGTACTCTGCGTTGATACCACTGCTT
>YG-GM4_3p
cgcgagagtcgagtgGTACTCTGCGTTGATACCACTGCTT
>YG-GM5_3p
cagctgatatatatgGTACTCTGCGTTGATACCACTGCTT
>YG-GM6_3p
cacatagagatacagaGTACTCTGCGTTGATACCACTGCTT
>YG-GM7_3p
cgcagcgctcgactgtGTACTCTGCGTTGATACCACTGCTT
>YG-GM8_3p
tctgtctcgcgtgtgtGTACTCTGCGTTGATACCACTGCTT
>YG-GM9_3p
ctctgagatagcgcgtGTACTCTGCGTTGATACCACTGCTT
>YG-GM10_3p
tagatatacgtatagGTACTCTGCGTTGATACCACTGCTT
```

Ouptut files will be named in the following way:
```
cat demux.*.bamm
demux.5p--YG_GM10_3p.bam
demux.5p--YG_GM2_3p.bam
demux.5p--YG_GM4_3p.bam
demux.5p--YG_GM6_3p.bam
demux.5p--YG_GM8_3p.bam
```

## Step 3: Convert to Fastqs

As an example:
```
bamtools convert -format fastq -in demux.5p--YG_GM10_3p.bam -out demux.5p--YG_GM8_3p.bam.fastq
```

## Step 4: Alignment 

Using ```minimap2```, map reads to the hg19 reference genome
```
minimap2  -ax splice -uf -C5 ref/hg19.fa demux.5p--YG_GM10_3p.bam.fastq > demux.5p--YG_GM10_3p.bam
samtools sort demux.5p--YG_GM8_3p.bam > demux.5p--YG_GM2_3p.sort.bam 
samtools index demux.5p--YG_GM2_3p.sort.bam 
```

Once you have repeated this step for all samples, merge the bams together:

```
samtools merge total.merge.bam demux.5p--YG_GM2_3p.sort.bam demux.5p--YG_GM4_3p.sort.bam demux.5p--YG_GM6_3p.sort.bam demux.5p--YG_GM8_3p.sort.bam demux.5p--YG_GM10_3p.sort.bam
```

## Step 6: Define polyA reads 

```
python scripts/01_filter_for_poly.py -ib total.merge.bam -o test -f ref/hg19.fa 
```

## Step 7: Define PAS sites

First, truncate the 3'UTR and exon annotation files to be in comprised of 6 columns.
```
cut -f1-6 ref/hg19_3utr.bed > ref/hg19_3utr.1-6.bed
cut -f1-6 ref/hg19_exons.bed > ref/hg19_exons.1-6.bed
```

Remove reads that do not overlap an upstream exon.
```
bedtools intersect -b ref/hg19_3utr.1-6.bed -a test.noMP.bed -s -wa -wb > test.noMP.bed.UTR
bedtools intersect -b ref/hg19_exons.1-6.bed -a test.noMP.bed.UTR -s -wa -wb > test.noMP.bed.UTR.EXONS
bedtools intersect -b test.noMP.bed.UTR -a test.noMP.bed -v -wa > test.noMP.bed.INTRON
python scripts/02_noMP_exon_UTR.py -i test.noMP.meta.txt -p test
```

## Step 8: Annotate

Prepare input files for annotation
```
sort -k1,1 -k2,2n test.noMP.restricted.peaks.bed > test.noMP.restricted.peaks.sort.bed
bedtools merge -i test.noMP.restricted.peaks.sort.bed -c 4,5,6 -o collapse,mean,distinct > test.noMP.restricted3.peaks.sort.refined.bed
python scripts/03_final_peaks.py  -p test 
```

We included RefSeq annotations in the (ref/)(https://github.com/ankeetashah/Benchmarking-APA/tree/main/ref) directory. However, if one is working with a different reference genome, download the approporiate annotation file in the following way (replace ``` -D hg19```):
```
mysql --user=genome -N --host=genome-mysql.cse.ucsc.edu -A -D hg19 -e "select name,name2,strand from refGene" > Refseq2Gene.txt
```

Use the ```annotatePeaks.pl``` script from [Homer](http://homer.ucsd.edu/homer/introduction/install.html) to annotate the putative polyadenylation sites. 
```
annotatePeaks.pl test.noMP.restricted.peaks.sort.refined.score.bed hg19 > test.homer.annotate.BED
python scripts/04_annotate_usage.py -p test -r ref/Refseq2Gene.txt -a 
```

Calculate the usage of polyadenylation sites (i.e. PAUs).
```
sort -k1,1 -k2,2n test.FINAL.BED | uniq >  test.FINAL.UNIQ.BED
python scripts/04_annotate_usage.py -p test -r ref/Refseq2Gene.txt -u 
```

The output file will be formatted as such, with ```chromosome number, start, end, RefSeq annotation_genicfeature, usage, strand, and coverage``` separated by tabs:
```
chr1    2335113 2335213 NM_007033_RER1_3' UTR   0.443037974684  +       70
chr1    2495208 2495308 NM_003820_TNFRSF14_3' UTR       0.978647686833  +       275
chr1    2495429 2495529 NM_003820_TNFRSF14_3' UTR       0.0142348754448 +       4
chr1    3650371 3650471 NM_001204186_TP73_3' UTR        0.466666666667  +       7
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
