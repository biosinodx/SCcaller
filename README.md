# SCcaller
Single Cell Caller (SCcaller) - Identify single nucleotide variations (SNVs) and short insertions and deletions (INDELs) from single cell sequencing data

Version 2.0.0

Updated date: 2019.04.01

Cite us:

Dong X et al. Accurate identification of single-nucleotide variants in whole-genome-amplified single cells. Nat Methods. 2017 May;14(5):491-493. doi: 10.1038/nmeth.4227.

#####
## Author and License

Author: Xiao Dong, Yujue Wang

Email: biosinodx@gmail.com (X.D.), xiao.dong@einstein.yu.edu (X.D.), spsc83@gmail.com (Y.W.)

Licensed under the GNU Affero General Public License version 3 or later

#####
## DEPENDENCIES

python 2.7

python modules os, argparse, sys, subprocess, re, collections, itertools, logging, time, functiontools, random, string, math, numpy, multiprocessing, pysam(0.15.1)

samtools v.1.9+ (Other versions not tested)

#####
## USAGE

###
### I. Basic usage: calling SNV and INDEL from a cell

#### I.a When you have heterozygous SNPs pre-called from bulk DNA of the same subject,

python sccaller_v2.0.0.py \

  --bam cell.bam \ # bam file of a single cell
  
  --fasta ref.fa \ # reference genome in fasta format
  
  --output cell.vcf \ # output vcf file
  
  --snp_type hsnp \ # using heterozygous SNPs pre-called from bulk DNA
  
  --snp_in hsnp.vcf (or bed) \ # vcf or bed file of heterozygous SNPs pre-called from bulk DNA
  
  --cpu_num 8 \ # using 8 cpu threads
  
  --engine samtools # using samtools engine

#### I.b When you do not have heterozygous SNPs pre-called from bulk DNA of the same subject, obtain SNPs from dbSNP or other databases,

python sccaller_v2.0.0.py \

  --bam cell.bam \ # bam file of a single cell
  
  --fasta ref.fa \ # reference genome in fasta format
  
  --output cell.vcf \ # output vcf file
  
  --snp_type dbsnp \ # using SNPs from dbSNP database (or other database)
  
  --snp_in hsnp.vcf (or bed) \ # vcf or bed file containing all SNPs in dbSNP (or other) database
       
  --engine samtools # using samtools engine

### II. Calling somatic mutations not present in bulk DNA

#### II.a Step 1. Calling SNVs and INDELs from a cell together with bulk DNA in input,

python sccaller_v2.0.0.py \

  --bam cell.bam \ # bam file of a single cell
  
  --bulk bulk.bam \ # bam file of bulk DNA
  
  --fasta ref.fa \ # reference genome in fasta format
  
  --output cell.vcf \ # output vcf file
  
  --snp_type hsnp \ # using heterozygous SNPs pre-called from bulk DNA
  
  --snp_in hsnp.vcf (or bed) \ # vcf or bed file of heterozygous SNPs pre-called from bulk DNA
  
  --cpu_num 8 \ # using 8 cpu threads
     
  --engine samtools # using samtools engine

#### II.b Step 2. Filtering

grep '0/1' cell.vcf | grep 'True' | awk '$1<=22 && $7=="." && length($5)==1' | awk -F "[:,]" '$8+$9>=20' > cell.somatic.snv.vcf

grep '0/1' cell.vcf | grep 'True' | awk '$1<=22 && $7=="." && length($5)>1' | awk -F "[:,]" '$8+$9>=20' > cell.somatic.indel.vcf

#####
## RELEASE NOTES

v2.0.0, 2019.04.01, allowing parallele processing, optimizing I/O, optimizing pipeline, output in vcf format, and fixing bugs

v1.21, 2018.08.18, fixing bugs

v1.2, 2017.05.01, allow INDEL calling, release version

v1.1.3, 2017.01.09, users can change the min mapQ, default to 40

v1.1.2, 2016.12.30, fixing bugs

v1.1.1, 2016.12.29, update read_mpileup to consider indels

V1.1, 2016.07.25, fixing bugs, release version

v1.0, 2016.04.26, release version

v0.0.4, 2016.04.26, fixed bugs in readling mpileup file

v0.0.3, 2016.04.22, read_mpilup function returns mindepth fails before returning reference genotype

v0.0.3, 2016.04.22, default mapQ change from 20 to 40

v0.0.2, 2016.04.19, fix bugs - jump mpileup file column not fit problem.

v0.0.1, 2016.03, add likelihood ratio test based on null distribution from the data resampling.
