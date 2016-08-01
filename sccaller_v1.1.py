### Single Cell Caller (SCcaller) - Identify single nucleotide variations (SNV) from single cell sequencing data
# Copyright (C) 2016  Dong, et. al.
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero General Public License for more details.

### updates
# V1.1, 2016.07.25, fixing bugs
# v1.0, 2016.04.26, cleanup and release version of v0.0.4
# v0.0.4, 2016.04.26, fixed bugs in readling mpileup file
# v0.0.3, 2016.04.22, read_mpilup function returns mindepth fails before returning reference genotype
# v0.0.3, 2016.04.22, default mapQ change from 20 to 40
# v0.0.2, 2016.04.19, fix bugs - jump mpileup file column not fit problem.
# v0.0.1, 2016.03, add likelihood ratio test based on null distribution from the data resampling.

import os
import argparse
import sys
import string
import math
import numpy as np

#Parse commandline arguments
parser=argparse.ArgumentParser(description="sccaller.py, v1.0, Xiao Dong, biosinodx@gmail.com, xiao.dong@einstein.yu.edu")
parser.add_argument("-a", "--analysis", type=str,required=True, help="{hsnp, varcall, cutoff}")

parser.add_argument("-i", "--input", type=str,required=True, help="The input file, required. {hsnp, varcall} mpileup file generated using 'samtools mpileup -C50 -Osf refgenome.fa input.bam > input.mpileup', samtools version 1.3 was tested good; {cutoff}: bedfile generated from {varcall}")
parser.add_argument("-o","--output",type=str,help="The output header, required for {hsnp, varcall}")
parser.add_argument("-s","--gss",type=str,help="Input known heterozygous snps (hSNP) - required for {varcall} output from {hsnp} for bias estimation")
parser.add_argument("-w","--wkdir",type=str,default='./',help="Work dir. Default: ./")

parser.add_argument("--indel",type=str,default='off',help="indel call {off, on}. Default: off, can't turn on at the moment")
parser.add_argument("--min",type=int,default=10,help="min. # reads. Default: 10")
parser.add_argument("--minvar",type=int,default=4,help="min. # variant supporting reads. Default: 4")

parser.add_argument("--RD",type=int,default=20,help="Known heterogous SNP call required read depths. Default: 20. Recommand: 10,15,20, depending on average read depth")
parser.add_argument("--snp",type=str,default='control',help="Known snp type input for known heterogous SNP call {dbsnp, control}. Default: control.")
parser.add_argument("--snpin",type=str,default='',help="Candidate snp input file for known heterogous call. file type: bed (1-based) or vcf. NOTE: should be sorted in the same order as the mpileup file.")

parser.add_argument("--lamb",type=int,default=10000,help="lambda for bias estimation. Default=10000.")
parser.add_argument("--bias",type=float,default=0.75,help="Default theta (bias) for SNVs whose theta cannot be estimated. Default=0.75.")
parser.add_argument("--null",type=float,default=0.03,help="Allelic fraction of artifacts. Default=0.03.")

args=parser.parse_args()

### set work dir
os.chdir(args.wkdir)

### indel call
indel_call=False
#if args.indel=='on':
#	indel_call=True

### known standard heterogeneity call read depth
gd_rd=args.RD
# gd_rd=20;#for testing

### minvar. # variant supporting reads
min_depth=args.min

var_min=args.minvar
#var_min=4;#for testing

### snp input for known heterogous call {dbsnp, control}
snptype=args.snp
#snptype='control'

### lambda
lamb=args.lamb
#lamb=10000

### mpile up # colunms expected
mpu_ncol=8

### read one mpileup line
def read_mpileup(rm_line, rm_indel=indel_call, rm_minvar=var_min, illumina_score='1.3+', rm_minmapQ=40, rm_mindepth=min_depth):
	rm_line=rm_line.split('\n')[0]
	rm_line=rm_line.split('\t')
	if not rm_indel:
		if ('-' in rm_line[4]) or ('+' in rm_line[4]) or ('*' in rm_line[4]):
			return(['indel'])
	rm_out=rm_line[0:3]
	rm_tmp1=rm_line[4].split('^');# remove read starting code '^]' in mpileup genotype column
	rm_tmp2=rm_tmp1[0]
	for i in rm_tmp1[1:]:
		rm_tmp2=rm_tmp2+i[1:]
	rm_line[4]=rm_tmp2
	rm_line[4]=rm_line[4].replace('$','');# remove read ending code '$' in mpileup genotype column
	pick=[]
	for i in range(0, int(rm_line[3])):
		if ord(rm_line[6][i])-33 >= rm_minmapQ :
			pick.append(i)
	tmp=rm_line[7].split(',')
	t1='';t2=[];t3=[];t4=[];
	phredscore=33;
	for i in pick:
		t1=t1+rm_line[4][i]
		t2.append(ord(rm_line[5][i])-phredscore)
		t3.append(ord(rm_line[6][i])-33)
		t4.append(int(tmp[i]))
	t1=t1.upper()
	tmp=t1.replace(',','')
	tmp=tmp.replace('.','')
	tmp=tmp.upper()
	maxvarp=['A','C','G','T']
	maxvar=[tmp.count('A'),tmp.count('C'),tmp.count('G'),tmp.count('T')]
	tmp3=max(maxvar)
	for i in range(0,len(maxvar)):
		if tmp3==maxvar[i]:
			break
	rm_out.append(maxvarp[i])
	rm_out.append(len(t1)-maxvar[i])
	rm_out.append(maxvar[i])
	if len(t1) < rm_mindepth:
		return(['lessmindepth'])
	if maxvar[i] < rm_minvar:
		return(['refgenotype']);#return reference genome type and for skipping this line
	for rm_i in range(0, len(t1)):
		rm_tmp1='-'
		if t1[rm_i].isupper():
			rm_tmp1='+'
		if  ord(rm_line[6][rm_i])-33 >= rm_minmapQ:
			rm_out.append([t1[rm_i].upper(), rm_tmp1, t2[rm_i], t3[rm_i], t4[rm_i]])
	return(rm_out)

### known heterogenous site call for input bias estimates, using mpileup file and candiddate known (dbsnp) or real known (control) input
def golden_hetero(gh_mpinfile, gh_snpinfile, gh_outfile, gh_snptype=snptype, gh_readdepth=gd_rd):
	gh_mpin = open(gh_mpinfile)
	gh_snpin = open(gh_snpinfile)
	gh_out = open(gh_outfile,'w')
	while True:
		gh_mpline=gh_mpin.readline()
		gh_snpline=gh_snpin.readline()
		if len(gh_mpline)==0 or len(gh_snpline)==0:
			break
		while gh_mpline[0]=='#' and len(gh_mpline.split('\t'))!=mpu_ncol:
			gh_mpline=gh_mpin.readline()
			if len(gh_mpline)==0:
				break
		if len(gh_mpline)==0:
			break
		while gh_snpline[0]=='#':
			gh_snpline=gh_snpin.readline()
			if len(gh_snpline)==0:
				break
		if len(gh_snpline)==0:
			break
		while gh_mpline.split('\t')[0] != gh_snpline.split('\t')[0] or len(gh_mpline.split('\t'))<mpu_ncol:
			if gh_mpline.split('\t')[0] < gh_snpline.split('\t')[0]:
				gh_mpline=gh_mpin.readline()
			elif gh_mpline.split('\t')[0] > gh_snpline.split('\t')[0]:
				gh_snpline=gh_snpin.readline()
			elif len(gh_mpline.split('\t'))<mpu_ncol:
				gh_mpline=gh_mpin.readline()
			elif len(gh_snpline.split('\t'))<mpu_ncol:
				gh_snpline=gh_snpin.readline()
			if len(gh_mpline)==0 or len(gh_snpline)==0:
				break
		if len(gh_mpline)==0 or len(gh_snpline)==0:
			break
		while gh_mpline.split('\t')[1] != gh_snpline.split('\t')[1] or len(gh_mpline.split('\t'))<mpu_ncol:
			if int(gh_mpline.split('\t')[1]) < int(gh_snpline.split('\t')[1]):
				gh_mpline=gh_mpin.readline()
			elif int(gh_mpline.split('\t')[1]) > int(gh_snpline.split('\t')[1]):
				gh_snpline=gh_snpin.readline()
			elif len(gh_mpline.split('\t'))<mpu_ncol:
				gh_mpline=gh_mpin.readline()
			elif len(gh_snpline.split('\t'))<mpu_ncol:
				gh_snpline=gh_snpin.readline()			
			if len(gh_mpline)==0 or len(gh_snpline)==0:
				break
		if len(gh_mpline)==0 or len(gh_snpline)==0:
			break
		gh_x = read_mpileup(rm_line=gh_mpline, rm_indel=False, rm_minvar=0, rm_minmapQ=20)
		if len(gh_x)==1:
			continue
		gh_xref=gh_x[2]
		gh_xalt=[]
		gh_xaltfinal=''
		gh_xrefcount=0
		gh_xaltcount=0
		# current known standard snp allele read counts are based on all mapped reads without considering base quality
		for i in gh_x[6:]:
			if i[0]==gh_xref or i[0]=='.' or i[0]==',':
				gh_xrefcount=gh_xrefcount+1
			else:
				gh_xaltcount=gh_xaltcount+1
				if not (i[0] in gh_xalt):
					gh_xalt.append(i[0])
					gh_xaltfinal=gh_xaltfinal+i[0]
		if gh_snptype=='dbsnp' and (gh_xaltcount<2 or gh_xrefcount<2 or (gh_xrefcount + gh_xaltcount < gh_readdepth)) :
			continue
		gh_l=gh_x[0] + '\t' + gh_x[1] + '\t' + gh_x[1] + '\t' + gh_x[2] + '\t' + gh_xaltfinal + '\t' + str(gh_xrefcount) + '\t' + str(gh_xaltcount) + '\n'
		gh_out.write(gh_l)
	gh_mpin.close()
	gh_snpin.close()
	gh_out.close()


### returning snp for testing and also update gd_pos for tracking, gt_snp=gd_snp
def gd_tracking(gt_chr, gt_pos, gt_i, gt_snp, gt_lambda=lamb):
	gt_chr = str(gt_chr)
	gt_lambda = int(gt_lambda)
	gt_pos = int(gt_pos)
	gt_out=[]
	gt_chrpass='No'
	gt_pospass='No'
	gt_start=-1
	while (gt_chrpass!='Passed' or gt_pospass!='Passed') and gt_i < len(gt_snp):
		#print(gt_snp[gt_i])
		if gt_snp[gt_i][0]==gt_chr:
			gt_chrpass='Passing'
		elif gt_chrpass=='Passing' and gt_snp[gt_i][0]!=gt_chr:
			gt_chrpass='Passed'
			break
		if gt_snp[gt_i][0]==gt_chr and (int(gt_snp[gt_i][1]) <= gt_pos+gt_lambda and int(gt_snp[gt_i][1]) >= gt_pos-gt_lambda):
			gt_out.append(gt_snp[gt_i])
			if gt_start == -1:
				gt_start = gt_i
		elif gt_snp[gt_i][0]==gt_chr and int(gt_snp[gt_i][1]) > gt_pos+gt_lambda :
			gt_pospass='Passed'
			break
		gt_i = gt_i + 1
	# gd_pos=gt_start
	if gt_start!=-1:
		gt_out.append(gt_start)
	else:
		gt_out.append(gt_i)
	return(gt_out)

### kernal for bias esimator, EQ for Epanechnikov quadratic kernel
def bias_kernal(bk_x0,bk_xi,bk_lambda,bk_method='EQ'):
	bk_t=float(abs(bk_x0 - bk_xi))/bk_lambda
	#print(bk_t)
	if abs(bk_t)<1:
		bk_d = 3.0 / 4.0 * (1-bk_t*bk_t)
	else:
		bk_d = 0.0
	return(bk_d)

### bias estimator, be_goldensnp is the output goldensnp file in memory, be_goldensnp from gd_tracking output
def bias_estimator(be_pos, be_goldensnp, be_lambda=lamb,be_default=args.bias):
	be_pos=int(be_pos)
	be_kwy=[]
	be_kw=[]
	for i in be_goldensnp:
		be_tmp=bias_kernal(be_pos, int(i[1]), bk_lambda=be_lambda, bk_method='EQ')
		#print(be_tmp)
		be_tmp1=float(i[5])+float(i[6])
		if be_tmp1<=0:
			continue
		be_tmp2=max([float(i[5]),float(i[6])])/(float(i[5])+float(i[6]))
		be_kwy.append(be_tmp*be_tmp1*be_tmp2)
		be_kw.append(be_tmp*be_tmp1)
	#Nadaraya-Watson kernel-weighted average
	if len(be_kwy)>0 and sum(be_kw)>0 :
		be_a = sum(be_kwy) / sum(be_kw)
	# return args.bias when no neighboring heterozygous base
	else:
		be_a = be_default
	return(be_a)

### P_b_GG: log Probability value
### r - reference base (not in use), m - alternative base, f = {0, 0.125, bias (or 1-bias depending on mut>ref or ref>mut), 1} for {ref/ref, ref/artifacts, ref/mut; mut/mut}
def P_b_GG(bp,r,m,f):
	e=10**(-bp[2]/10.0)
	if bp[0]==m:
		a=f * (1-e) + (1-f) * e/3
	elif bp[0]==',' or bp[0]=='.' or bp[0]==r:
		a=f*e/3+(1-f)*(1-e)
	else:
		a=e/3
	return(math.log(a,10))

### sc caller, wait to write, sc_candidate is the read_mpileup output, sc_bias is the bias estimates from bias_estimator
def sc_caller(sc_candidate, sc_bias, sc_artifact=args.null):
	ref=sc_candidate[2]
	mut=sc_candidate[3]
	f=sc_bias
	RR=0.0;RA=0.0;RM=0.0;MM=0.0
	if sc_candidate[4] > sc_candidate[5]:
		f=1-f
	for i in sc_candidate[6:]:
		RR=RR+P_b_GG(i, ref, mut, sc_artifact);# ref/ empiracal artifact; if sc_artifact set as 0.0, than assumes that all errors came from sequencing noise
		RA=RA+P_b_GG(i, ref, mut, 0.125);# ref / artifact
		RM=RM+P_b_GG(i, ref, mut, f);# ref/mut
		MM=MM+P_b_GG(i, ref, mut, 1.0);# mut/mut
	return([RR, RA, RM, MM])

### var_calling total, vc_mpileupfile - mpileupfile; vc_snpfile from formated known hetero snp;
def var_calling(vc_mpileupfile, vc_snpfile, vc_outfile, vc_lambda=lamb,vc_indelneighbor=11):
	# read known snps
	gd_snpfile=open(vc_snpfile)
	gd_snpin=gd_snpfile.read()
	gd_snpfile.close()
	gd_snpin=gd_snpin.split('\n')
	gd_snp=[]
	for i in range(0, len(gd_snpin)-1):
		tmp=gd_snpin[i].split('\t')
		gd_snp.append(tmp)
	gd_snpin=None
	# perform variant calling
	f=open(vc_mpileupfile)
	outfile=open(vc_outfile, 'w')
	pos=0
	indel=[False for i in range(0, vc_indelneighbor)];#track indel
	while True:
		l=f.readline()
		if len(l)==0:
			break
		if len(l.split('\t'))!=mpu_ncol:
			continue
		x=read_mpileup(l)
		indelupdate=False;#track indel
		if x[0]=='indel':#track indel
			indelupdate=True;#track indel
		del indel[0];#track indel
		indel.append(indelupdate);#track indel
		if len(x)<4 :
			continue
		#print(l)
		tmp=gd_tracking(x[0], x[1], pos, gd_snp)
		#print(tmp)
		pos=tmp[len(tmp)-1]
		tmp=tmp[0:(len(tmp)-1)]
		tmp1=bias_estimator(x[1], tmp)
		#print(tmp1)
		out=x[0:6]
		#####out.append(sc_filter()) add filter results here
		filters=[]
		if any(indel):
			filters.append('indel_neighbor')
		if len(filters)==0:
			filterline='PASS'
		else:
			filterline=filters[0]
			if len(filters)>1:
				filterline=[filterline+','+i for i in filters[1:]]
		out.append(filterline)
		#####finish filter
		out.append(tmp1)
		#print(x)
		out=out + [10**i for i in sc_caller(x, tmp1)]
		outline=out[0]
		outline=outline + '\t' + str(out[1])
		for i in out[1:]:
			outline=outline+'\t' + str(i)
		outline=outline+'\n'
		outfile.write(outline)

def cutoff(co_infile, co_repeattime=100000, co_artifact=0.125):
	### read theta & depth in the real data
	xfile=open(co_infile)
	depth=[]
	theta=[]
	while len(theta)<=co_repeattime:
		xline=xfile.readline()
		if len(xline)==0:
			break
		if xline[0]=='#':
			continue
		y=xline.split('\t')
		if float(y[8])==0:
			continue
		depth.append(int(y[5])+int(y[6]))
		theta.append(float(y[8]))
	xfile.close()
	### sampling according to H_filter, llr=-log10(L0/L1)
	LLR=[]
	for i in range(0, len(depth)):
		f_artifact=0.125*theta[i]/0.5
		alt=np.random.binomial(depth[i], f_artifact)
		ref=depth[i]-alt
		L_filter=(1-f_artifact)**ref * ((f_artifact)**alt)
		## random select major allele
		major=np.random.binomial(1,.5)
		if major==0:
			f_true=0.5*0.5/theta[i]
		if major==1:
			f_true=0.5*theta[i]/0.5
		L_true=(1-f_true)**ref * ((f_true)**alt)
		## if L_filter/true is 0, assign a very small value
		if L_filter==0:
			L_filter=10**-100
		if L_true==0:
			L_true=10**-100
		LLR.append(-math.log10(L_filter / L_true))
	LLR=np.array(LLR)
	co_001=np.percentile(LLR,99)
	co_005=np.percentile(LLR,95)
	return([co_001, co_005])

def main(analysis):
	if analysis=='hsnp':
		golden_hetero(gh_mpinfile=args.input, gh_snpinfile=args.snpin, gh_outfile=args.output+'.hsnp.bed', gh_snptype=args.snp, gh_readdepth=gd_rd)
	if analysis=='varcall':
		var_calling(vc_mpileupfile=args.input, vc_snpfile= args.gss, vc_outfile=args.output+'.varcall.bed', vc_lambda=lamb)
	if analysis=='cutoff':
		x1=cutoff(co_infile=args.input)
		x=open(args.input,'a')
		#print(x1)
		line='##Cutoffs for L_artifact / L_H1 (if smaller than the following eta, reject H_artifact, accept H1): at alpha=0.01,eta=' + str(10**-x1[0]) + '; alpha=0.05,eta=' + str(10**-x1[1]) + '\n'
		x.write(line); x.close()

main(analysis=args.analysis)

