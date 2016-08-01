# SCcaller
Single Cell Caller (SCcaller) - Identify single nucleotide variations (SNVs) from single cell sequencing data

Version 1.1

Updated date: 2016.07.25

#####
## Author and License

Author: Xiao Dong

Email: biosinodx@gmail.com, xiao.dong@einstein.yu.edu

Licensed under the GNU Affero General Public License version 3 or later

#####
##DEPENDENCIES

python 2.7

python modules os,argparse,sys,string,math,numpy

samtools v.1.3+ (Other versions not tested)

#####
##USAGE

###
###General SNV detection

#####I. Prepare mpileupfile of a single cell bam for analysis

samtools mpileup -C50 -r chr1 -Osf refgenome.fa cell.bam > cell.chr1.mpileup

I highly recommend users generate mpileup files for each chromosome (as above), and analyze them separately.

#####II. Obtain known heterozygous SNP (hsnp) positions

python sccaller.py -a hsnp -i cell.chr1.mpileup --snpin snpin.bed/vcf -o outputheader.chr1

The snpin.bed or snpin.vcf is a file of candidate snp positions in bed (1-based) or vcf format. One can use either public dbSNP SNPs or SNPs identified in bulk cell population sequencing obtained from the same individual.
Of note, this list of known heterozygous SNPs does not have to be precise, because the following analysis is robust to noise at this step and SNVs will be recalled using SCcaller in later steps.

#####III. Calling all potential SNVs

python sccaller.py -a varcall -i cell.chr1.mpileup -s outputheader.chr1.hsnp.bed -o outputheader.chr1

In this step, likelihoods of models in SCcaller are calculated. It will generate a file outputheader.chr1.varcall.bed including the following columns,

Column 1, chromosome

Column 2 & 3, position

Column 4, genotype of reference genome

Column 5, variant

Column 6, # reads of reference allele

Column 7, # reads of variant allele

Column 8, indel filter (other filters may be developed in future updates)

Column 9, allelic amplification bias of single cells

Column 10, likelihoods of artifact variant by sequencing noise, after correcting allelic amplification bias

Column 11, likelihoods of artifact variant by amplification error, after correcting allelic amplification bias

Column 12, likelihoods of heterozygous variant, after correcting allelic amplification bias

Column 13, likelihoods of homozygous variant, after correcting allelic amplification bias

#####IV. Calling SNVs using SCcaller

Generate null distribution of artifacts

python sccaller.py -a cutoff -i outputheader.chr1.varcall.bed

This will add one line to the end of outputheader.varcall.bed file, for example,

"##Cutoffs for L_artifact / L_H1 (if smaller than the following eta, reject H_artifact, accept H1): at alpha=0.01,eta=0.196643099187; alpha=0.05,eta=0.681174408573"

Filter variants at alpha=0.01 (corresponding eta=0.197 in the case above)

awk '$12!=0' outputheader.chr1.varcall.bed | awk '$11/$12<0.197 && $7/($6+$7)>1/8 && $8=="PASS"' > outputheader.chr1.sccaller.bed

The file outputheader.sccaller.bed is the final callset using SCcaller.

###
###Somatic SNV detection

#####I. Prepare mpileupfile of bulk/control bam for as controls

samtools mpileup -C50 -r chr1 -Osf refgenome.fa bulk.bam > bulk.chr1.mpileup

Highly recommend users generate mpileup files for each chromosome, and analyze them separately.

#####II. Detect if variant exist in bulk (by chromosome)

python reasoning.py -m bulk.chr1.mpileup -i cell.chr1.sccaller.bed -o cell-bulk.reason.txt -c chr1

#####III. Keep variants only with 'refgenotype' in cell-bulk.reason.txt file.

#####IV. Users should consider whether to filter out known SNPs (e.g. dbSNP) from the output.

#####
##RELEASE NOTES

v1.0, 2016.04.26, the cleanup and release version of v0.0.4

v0.0.4, 2016.04.26, fixed bugs in readling mpileup file

v0.0.3, 2016.04.22, read_mpilup function returns mindepth fails before returning reference genotype

v0.0.3, 2016.04.22, default mapQ change from 20 to 40

v0.0.2, 2016.04.19, fix bugs - jump mpileup file column not fit problem.

v0.0.1, 2016.03, add likelihood ratio test based on null distribution from the data resampling.
