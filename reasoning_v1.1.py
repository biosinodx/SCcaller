### check bulk cell population sequencing file for germline variants, indel reads, and low sequencing depth for calling somatic mutations
# Copyright (C) 2016  Dong, et. al.
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero General Public License for more details.

import os
import argparse
import sys
import string
import math

parser=argparse.ArgumentParser(description="reasoning.py, v1.0, Xiao Dong, biosinodx@gmail.com")

parser.add_argument("-i", "--input", type=str,required=True, help="The input bedfile")
parser.add_argument("-o","--output",type=str,required=True, help="The output file")
parser.add_argument("-m","--mpu",type=str,required=True, help="Input mpileup file")
parser.add_argument("-c","--chr",type=str,required=True, help="Input chr")


args=parser.parse_args()

filein=open(args.input)
filempl=open(args.mpu)
fileoout=open(args.output,'w')

indel_call=False
var_min=1
min_depth=20

### read one mpileup line
def read_mpileup(rm_line, rm_indel=indel_call, rm_minvar=var_min, illumina_score='1.3+', rm_minmapQ=20, rm_mindepth=min_depth):
	rm_line=rm_line.split('\n')[0]
	rm_line=rm_line.split('\t')
	if not rm_indel:
		if ('-' in rm_line[4]) or ('+' in rm_line[4]) or ('*' in rm_line[4]):
			return(['indel'])
	#rm_line=rm_line.split('\n')[0]
	#rm_line=rm_line.split('\t')
	rm_out=rm_line[0:3]
	#rm_line[4]=rm_line[4].replace('^]','');# remove read starting code '^]' in mpileup genotype column
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
	elif maxvar[i] < rm_minvar:
		return(['refgenotype']);#return reference genome type and for skipping this line
	else:
		return(['varreads'])

d1=filein.read()
d1=d1.split('\n')
d1p=[]
chrom=args.chr
for i in d1:
	if len(i)>0:
		if i.split('\t')[0]==chrom:
			a=int(i.split('\t')[1])
			d1p.append(a)

pointer=0; a=0
while True:
	if pointer>=len(d1p):
		break
	if a < d1p[pointer]:
		x=filempl.readline()
		if len(x)==0:
			break
		a=int(x.split('\t')[1])
	if a < d1p[pointer]:
		continue
	elif a == d1p[pointer]:
		abc=read_mpileup(x)
		out=chrom + '\t'+ str(a)
		for i in abc:
			out=out + '\t' + str(i)
		out = out + '\n'
		fileoout.write(out)
		pointer = pointer + 1
	else:
		out = 'nocoverageincontrol\n'
		fileooute.write(out)
		pointer = pointer + 1

filein.close()
filempl.close()
fileoout.close()

