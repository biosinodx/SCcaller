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
# v2.0.0, 2019.04.01, allowing parallel processing, optimizing I/O, optimizing pipeline, output in vcf format, and fixing bugs
# v1.21, 2018.08.18, fixing bugs
# v1.2, 2017.05.01, allowing INDEL calling
# v1.1.3, 2017.01.09, users can change the min mapQ, default to 40
# v1.1.2, 2016.12.30, fixing bugs
# v1.1.1, 2016.12.29, updating read_mpileup to consider indels
# v1.1, 2016.07.25, fixing bugs
# v1.0, 2016.04.26, cleanup and release version of v0.0.4
# v0.0.4, 2016.04.26, fixing bugs in readling mpileup file
# v0.0.3, 2016.04.22, read_mpilup function returns mindepth fails before returning reference genotype
# v0.0.3, 2016.04.22, default mapQ change from 20 to 40
# v0.0.2, 2016.04.19, fixing bugs - jump mpileup file column not fit problem.
# v0.0.1, 2016.03, adding likelihood ratio test based on null distribution from the data resampling.

import sys
import argparse
import os
from subprocess import Popen, PIPE
import multiprocessing
import math
import re
from collections import Counter
from itertools import compress
import logging
import time
from functools import wraps
import pysam  # 0.15.1
import copy
import random
import numpy as np

description_info = "Description: SCcaller, v2.0.0; Xiao Dong, biosinodx@gmail.com, xiao.dong@einstein.yu.edu; Yujue Wang, spsc83@gmail.com\n"

para_dict_m = {"-b, --bam": "bamfile of a single cell",
               "-f, --fasta": "fasta file of reference genome",
               "-s, --snp_in": "Candidate snp input file, either from dbsnp data or heterozygous snp (hsnp) data of the bulk, for known heterogous call. file type: bed (1-based) or vcf.",
               "-t, --snp_type": "SNP type for --snp_in. It could be either \"dbsnp\" or \"hsnp\". When choosing dbsnp, --bulk bulk_bamfile is required.",
               "-o, --output": "output file name"}

para_dict_o = {
    "-d, --wkdir": "work dir. Default: ./",
    "-l, --lamb": "lambda for bias estimation. Default=10000",
    "--bias": "default theta (bias) for SNVs whose theta cannot be estimated. Default=0.75",
    "--minvar": "min. num. variant supporting reads. Default: 4",
    "--mapq": "min. mapQ. Default: 40",
    "--min_depth": "min. reads. Default: 10",
    "--RD": "min. read depth of known heterogous SNP called from bulk when choosing -t dbsnp. Default: 20. Recommand: 10,15,20, depending "
            "on average read depth",
    "--null": "min. allelic fraction considered. Default=0.03",
    "-e, --engine": "pileup engine. samtools or pysam. Default: pysam",
    "-w, --work_num": "num. splits per chromosome for multi-process computing. Default: 100",
    "-n, --cpu_num": "num. processes. Default: 1",
    "--head": "first chromosome as sorted as in fasta file to analyze (1-based). Default: the first chr. in the fasta",
    "--tail": "last chromosome as sorted as in fasta file to analyze (1-based). Default: the last chr. in the fasta",
    "--format": "output file format. bed or vcf. Default: vcf",
    # "--coverage": "use \"--coverage\" to generate the coverage file at the same time",
    "--bulk": "bamfile of bulk DNA sequencing",
    "--bulk_min_var": "min. num. variant supporting reads for bulk. Default: 1",
    "--bulk_min_mapq": "min. mapQ for bulk. Default: 20",
    "--bulk_min_depth": "min. reads for bulk. Default: 20",
    "-h, --help": "Help"}

# expected mpileup minimum columns
MPU_COL_NUM = 8

# regular expression of indels
INSERTION = "\+[0-9]+[ACGTNacgtn]+"
DELETION = "\-[0-9]+[ACGTNacgtn]+"
INDEL = "\+[0-9]+[ACGTNacgtn]+|\-[0-9]+[ACGTNacgtn]+"
PHREDSCORE = 33  #
INVALIDNUM = -1

CFGNAME = "sys.cfg"
SECTIONPARA = "PARA"
SECTIONSYS = "SYS"
KEYINTERFACEON = "interface_on"

INVALIDCHROMOSOMENAME = "XXX"
DATAENDNAME = "End"
DONEMSG = "I'm done"
INVALIDCOORDINATE1 = -1
DATAENDNUM = -2
INDEXHEAD = 0

VARCALLFORMAT = "bed"
VCFFORMAT = "vcf"
ENGINESAMTOOLS = "samtools"
ENGINEPYSAM = "pysam"
MULTIPLEGENOTYPE = "multiple-genotype"  # type: str
NOTENOUGHVARIANTS = "No.variants<"
WORKVAR = 1
WORKCOVERAGE = 2
WORKVCF = 3
STOPSIGN = -1

DEFAULTBASEQUALITY = 30  # type: int
CUTOFFNUM = 100000

VERSION = "2.0.0"


# should_analyze = False
# q7 = multiprocessing.Manager().Queue()


class DataInQ7:
    def __init__(self, name, time_float):
        # type: (str, float) -> None
        self.name = name  # type: str
        self.time_float = time_float  # type: float


def send_name_time(should_send, q7):
    def send(f):
        @wraps(f)
        def decorated(*args, **kwargs):
            if should_send:
                start_ = time.time()
            ret_ = f(*args, **kwargs)
            if should_send:
                q7.put(DataInQ7(f.__name__, time.time() - start_), block=False)
            return ret_

        return decorated

    return send


# @send_name_time(should_send=should_analyze, q7=q7)


class VcfInfo:
    def __init__(self, gt, ad, bi, gq, pl, qual, genotype_num, genotype_class):
        """

        :type qual: str
        :type pl: str
        :type gq: str
        :type bi: float
        :type genotype_num: int
        :type ad: str
        :type so: str
        :type gt: str
        """
        self.gt = gt

        self.ad = ad
        self.bi = bi
        self.gq = gq
        self.pl = pl
        self.qual = qual
        self.genotype_num = genotype_num
        self.genotype_class = genotype_class


class OutLineStruct:
    def __init__(self, name, pos,
                 ref, var,
                 ref_num, var_num,
                 bias, sn,
                 ae, he,
                 ho, vcf_info, so, var_all):
        # type: (str, int, str, str, int, int, float, float, float, float, float, VcfInfo) -> OutLineStruct
        self.name = name
        self.pos = pos
        self.ref = ref
        self.var = var
        self.ref_num = ref_num  # 6
        self.var_num = var_num  # 7
        self.bias = bias
        self.sn = sn  # 10
        self.ae = ae  # 11
        self.he = he  # 12
        self.ho = ho  # 13
        self.so = so
        self.vcf_info = vcf_info
        self.var_all = var_all


class DataInQ5:
    def __init__(self, outline_struct, worker_id, work_type, coverage):
        self.outlineStruct = outline_struct  # type: OutLineStruct
        self.coverage = coverage  # type: list  # [name str, start int, stop int]
        self.worker_id = worker_id  # type: int
        self.work_type = work_type


class DataInQ2:
    chromosome_name = None  # type: str
    coordinate_1 = None  # type: int
    refcount = None  # type: int
    altcount = None  # type: int

    def __init__(self, chromosome_name, coordinate_1, refcount, altcount):
        # type: (str, int, int, int) -> DataInQ2
        self.chromosome_name = chromosome_name  # type: str
        self.coordinate_1 = coordinate_1  # type: int
        self.refcount = refcount  # type: int
        self.altcount = altcount  # type: int


class ReadInfo:
    def __init__(self, read_base, base_quality):
        # type: (str, str) -> None
        self.read_base = read_base  # type: str
        self.base_quality = base_quality  # type: int


class PileupDataStruct:
    def __init__(self, chromosome_name, coordinate_1, reference_base, variant, reference_allele_num, variant_allele_num,
                 genotype_num, read_info_list, so=None, genotype_class=-1, variant_all=""):
        # type: (str, int, str, str, int, int, int, list[ReadInfo], str, int) -> object
        self.chromosome_name = chromosome_name  # type: str
        self.coordinate_1 = coordinate_1  # type: int
        self.reference_base = reference_base  # type: str
        self.variant = variant  # type: str
        self.reference_allele_num = reference_allele_num  # type: int
        self.variant_allele_num = variant_allele_num  # type: int
        self.genotype_num = genotype_num  # type: int
        self.read_info_list = read_info_list  # type: list[ReadInfo]
        self.so = so  # type: str  # reasoning  bulk ref -- "True"   bulk var -- "False"  else -- "NA"
        self.genotype_class = genotype_class  # type: int  # 0 unknown 1:0,0 0,1 1,1      2:1,1 1,2 2,2
        self.variant_all = variant_all


def my_print(msg, name):
    pass
    # if name in ["main", "control"]:
    #     print "{}".format(msg)


# @send_name_time(should_send=should_analyze, q7=q7)
def handle_paras():
    """
    handle the parameters.
    :return: argparse.Namespace
    """

    def get_usage():
        """
        usage description
        :return: usage description string
        """
        return "Usage: " + os.path.basename(sys.argv[0]) + \
               "[-h] [-d WKDIR] [-l LAMB] [--bias NUM] [--minvar NUM]" + \
               " [--mapq NUM] [--min_depth NUM] [--RD NUM] [--null NUM]" + \
               " [--bulk BAM] [--bulk_min_depth NUM] [--bulk_min_mapq NUM] [--bulk_min_var NUM]" + \
               " [--format {bed,vcf}] [--head NUM] [--tail NUM] [-e {pysam, samtools}] [--cpu_num NUM] [-w NUM] [-n NUM]" + \
               " -t {dbsnp,hsnp} -b BAM -f FASTA -s " + \
               "SNP_IN -o OUTPUT\n\n"

    def is_arg_file_exists(args):
        """
        check if the bam, fasta, vcf file in arguments exits
        :param args: arguments
        :return:
        true files exit. It is OK to continue.
        false file doesn't exit. Should not continue.
        """
        should_continue = True

        if not os.path.exists(args.bam):
            print os.path.basename(sys.argv[0]) + ": error: file [" + args.bam + "] doesn't exist."
            should_continue = False

        if not os.path.exists(args.fasta):
            print os.path.basename(sys.argv[0]) + ": error: file [" + args.fasta + "] doesn't exist."
            should_continue = False

        if not os.path.exists(args.snp_in):
            print os.path.basename(sys.argv[0]) + ": error: file [" + args.snp_in + "] doesn't exist."
            should_continue = False
        return should_continue

    def get_argument_info(m, o):
        # type: (dict, dict) -> str
        """
        get the description of arguments
        :param m: dict of mandatory arguments
        :param o: dict of optional arguments
        :return: description string of arguments
        """
        argument_info = "\nArguments:\n"
        data_list = list(m.items())
        data_list.sort()
        # print "aaaa=" + str(data_list)
        for key, value in data_list:
            argument_info += "{0:16}{1}\n".format(key + ":", value)
        argument_info += "\nOptional arguments:\n"

        data_list = list(o.items())
        data_list.sort()

        for key, value in data_list:
            argument_info += "{0:16}{1}\n".format(key + ":", value)
        return argument_info

    if len(sys.argv) == 1:
        print get_usage() + description_info + get_argument_info(para_dict_m, para_dict_o)
        exit(0)

    parser = argparse.ArgumentParser(description="sccaller v1.3", add_help=False)

    parser.add_argument("-b", "--bam", type=str, help="bam file name")

    parser.add_argument("-f", "--fasta", type=str, help="fasta file name")

    parser.add_argument("-s", "--snp_in", type=str, help="Candidate snp input file for known heterogous call. "
                                                         "file type: bed (1-based) or vcf. NOTE: should be sorted "
                                                         "in the same order as the mpileup file.")
    parser.add_argument("-o", "--output", type=str, help='output file name')
    parser.add_argument("-t", "--snp_type", type=str, choices=["dbsnp", "hsnp"],
                        help="Known snp type input for known heterogous SNP call. ")
    parser.add_argument("-d", "--wkdir", type=str, default="./", help="Work dir. Default: ./")
    parser.add_argument("-h", "--help", action="store_true")
    parser.add_argument("-l", "--lamb", type=int, default=10000, help="lambda for bias estimation. Default=10000.")
    parser.add_argument("--bias", type=float, default=0.75, help="Default theta (bias) for SNVs whose theta cannot be "
                                                                 "estimated. Default=0.75.")
    parser.add_argument("--minvar", type=int, default=4, help="min. # variant supporting reads. Default: 4")
    parser.add_argument("--mapq", type=int, default=40, help="min. mapQ. Default: 40")
    parser.add_argument("--min_depth", type=int, default=10, help="min. # reads. Default: 10")
    parser.add_argument("--RD", type=int, default=20,
                        help="Known heterogous SNP call required read depths. Default: 20. Recommand: 10,15,20, "
                             "depending on average read depth")
    parser.add_argument("--null", type=float, default=0.03, help="Allelic fraction of artifacts. Default=0.03.")
    parser.add_argument("-e", "--engine", type=str, choices=[ENGINEPYSAM, ENGINESAMTOOLS],
                        default=ENGINEPYSAM)
    parser.add_argument("-w", "--work_num", type=int, default=100)
    parser.add_argument("-n", "--cpu_num", type=int, default=1)

    parser.add_argument("--head", type=int, default=1)
    parser.add_argument("--tail", type=int, default=-1)
    parser.add_argument("--format", type=str, choices=[VARCALLFORMAT, VCFFORMAT],
                        default=VCFFORMAT)
    parser.add_argument("--coverage", action="store_true", default=False)
    parser.add_argument("--bulk", type=str, default="")
    parser.add_argument("--bulk_min_var", type=int, default=1)
    parser.add_argument("--bulk_min_mapq", type=int, default=20)
    parser.add_argument("--bulk_min_depth", type=int, default=20)

    args = parser.parse_args()

    os.chdir(args.wkdir)
    if args.help:
        print "{0}{1}{2}".format(get_usage(), description_info, get_argument_info(para_dict_m, para_dict_o))
        exit(0)

    # check mandatory arguments
    if args.bam:
        para_dict_m.pop("-b, --bam")
    if args.fasta:
        para_dict_m.pop("-f, --fasta")
    if args.snp_in:
        para_dict_m.pop("-s, --snp_in")
    if args.output:
        para_dict_m.pop("-o, --output")
    if args.snp_type:
        para_dict_m.pop("-t, --snp_type")

    if len(para_dict_m) > 0:
        for key, value in para_dict_m.items():
            print os.path.basename(sys.argv[0]) + ": error: argument " + key + " is required"
        exit(0)
    if args.snp_type == "dbsnp" and args.bulk == "":
        print os.path.basename(sys.argv[0]) + ": When choosing dbsnp, --bulk bulk_bamfile is required."
        exit(0)

    # check mandatory files
    if not is_arg_file_exists(args):
        exit(0)
    # check result file
    if os.path.exists(args.output):
        print os.path.basename(sys.argv[
                                   0]) + ": error: output file: " + args.output + " already exists. Please delete it manually or try a different file name."
        exit(0)

    if os.path.exists(get_my_filename(args.output, ".coverage", os.path.dirname(args.output) + "/")):
        os.remove(get_my_filename(args.output, ".coverage", os.path.dirname(args.output) + "/"))
    if os.path.exists(get_my_filename(args.output, ".reasoning", os.path.dirname(args.output) + "/")):
        os.remove(get_my_filename(args.output, ".reasoning", os.path.dirname(args.output) + "/"))
    return args


# @send_name_time(should_send=should_analyze, q7=q7)
def parse_indel(indel_str, indel_list_out):
    # type: (str, list) -> int
    """
    parse indel string, return the indel string and length of indel string
    :param indel_str:indel string (the '+' and '-' in indel string doesn't affect the result)
    :param indel_list_out:indel string (should be clear before)
    :return:
    >0 return indel string length
    =0 in case of indel and svn
    <0 in case of invalid format of indel
    eg:
        rr = []
        ret = parse_indel("+2AAC",rr)   #it will print:
        print ret                       #0
        print rr                        #['+', '2', 'A', 'A']

        rr = []
        ret = parse_indel("+2AA",rr)    #it will print:
        print ret                       #4
        print rr                        #['+', '2', 'A', 'A']
    """
    i = 0  # type: int
    j = 0  # type: int
    len_indel_str = len(indel_str)
    for i in xrange(len_indel_str):
        if indel_str[i] == "+" or indel_str[i] == "-":
            continue
        else:
            break
    # more than 1 '+' or '-'
    if i > 1:
        return -1
    buf = indel_str[i:]
    len_buf = len(buf)
    for j in xrange(len_buf):
        if buf[j].isalpha():
            break

    if len_indel_str - i - j > int(buf[:j]):
        indel_list_out.extend(indel_str[0:i + j + int(buf[:j])])
        return 0
    elif len_indel_str - i - j < int(buf[:j]):
        return -2
    indel_list_out.extend(indel_str[0:i + j + int(buf[:j])])
    return int(buf[:j]) + j + i


# @send_name_time(should_send=should_analyze, q7=q7)
def remove_head_end(str_in, readbase_len, str_map_q, name, pos):
    # type: (str, int, str, str, int) -> str

    def get_head_MQ(str_read_base, head_index, str_map_q):
        """
        calculate the mapQ of read head (^) from readbase string without I
        :param str_read_base:
        :param head_index:
        :param str_map_q:
        :return:
        """
        # logging.debug("str_read_base={0} head_index={1} str_map_q={2}".format(str_read_base,head_index,str_map_q))
        counter = 0
        len_read_base = len(str_read_base)
        for i in xrange(len_read_base):
            if i == head_index:
                break
            if str_read_base[i] in [".", ",", "A", "T", "G", "C"]:
                if i == 0:
                    counter += 1
                else:
                    if str_read_base[i - 1] != "^":
                        counter += 1
        # logging.debug("result={}".format(str_map_q[counter]))
        return str_map_q[counter]

    str1 = str_in
    # remove '$' which is not head mapQ
    str1 = re.sub("(?<!\^)\$", "", str1)

    str2 = re.sub("\^\S{1}", "", str1)
    missing_num = readbase_len - len(str2.replace("I", ""))
    if missing_num == 0:
        return str2
    # look for pattern like '^*^'
    total_found_num = 0
    while 1:
        tmp = re.findall("\^\S{1}\^", str1)
        if len(tmp) == 0:
            break
        total_found_num += 1
        # fill in the lost data
        head_index = str1.find(tmp[0])
        str1 = str1[:head_index + 1] + get_head_MQ(str1, head_index, str_map_q) + str1[head_index + 1:]

    if total_found_num == missing_num:
        return re.sub("\^\S{1}", "", str1)

    num_found = 0
    # look for pattern like '^'
    for i in xrange(len(str1) - 1):
        if str1[i] == "^" and str1[i + 1] != get_head_MQ(str1, i, str_map_q):
            num_found += 1
            str1 = str1[:i + 1] + get_head_MQ(str1, i, str_map_q) + str1[i + 1:]

    total_found_num += num_found
    if total_found_num == missing_num:
        return re.sub("\^\S{1}", "", str1)
    logging.critical("Can not handle this read base right now: "
                     "name={0} pos={1} readbase_len={2} readbase={3} str_map_q={4}".format(name, pos,
                                                                                           readbase_len,
                                                                                           str_in, str_map_q))
    return "???"


def compress_read_base(read_base, my_filter):
    """
    compress read base but leave the I (indel)
    :type read_base: str
    :type my_filter: list
    """
    counter = 0
    for i in read_base:
        if i == "I":
            my_filter.insert(counter, True)
        counter += 1
    return "".join(list(compress(read_base, my_filter)))


def get_reference_variant_allele_num(read_base):
    v_num = 0
    r_num = 0
    for i in xrange(len(read_base)):
        if read_base[i] in [".", ","]:
            r_num += 1
        if read_base[i] in ["A", "T", "C", "G"]:
            v_num += 1
        if i > 0 and read_base[i] in ["I", "i"] and read_base[i - 1] in [".", ","]:
            r_num -= 1
            v_num += 1
    return [r_num, v_num]


def rebuild_read_base_list(read_base_without_i, indel_name_list, indel_count_list):
    """
    add indel to the end of the readbase without I, return as list
    gather
    :param read_base_without_i:
    :param indel_name_list:
    :param indel_count_list:
    :return:
    """
    read_base_list = []
    read_base_list.extend(list(read_base_without_i))
    tmp = [[indel_name_list[list_index] for i in xrange(indel_count_list[list_index])]
           for list_index in xrange(len(indel_count_list))]
    for i in tmp:
        read_base_list.extend(i)
    # [read_base_list.extend([indel_name_list[list_index] for i in xrange(indel_count_list[list_index])])
    #  for list_index in xrange(len(indel_count_list))]

    return read_base_list


def read_mpileup(mpileup_list, rm_minvar, rm_minmap_q, rm_mindepth, is_gh, worker_id):
    """
    screenfor mutations, filter out alignments that do not meet the requirements of mapQ, and convert the data format
    :param is_gh:
    :type is_gh: bool
    :type rm_mindepth: int
    :type rm_minmap_q: int
    :type rm_minvar: int
    :type mpileup_list: list
    :param mpileup_list: 1 line data of mpileup
    :param rm_minvar: variant supporting reads.
    :param rm_minmap_q: minimun mapQ
    :param rm_mindepth: minimun read depth
    :return:
    not empty PileupDataStruct: mutation data
    -1: return reference genome type and for skipping this line, but should be count in coverage
    []: not mutation, should not be count in coverage
    """
    rm_indel = ""  # type: list
    if "*" in mpileup_list[4]:
        return []  # indel read overlaps
    indel_name_list = []
    indel_count_list = []
    if ("-" in mpileup_list[4]) or ("+" in mpileup_list[4]):
        indel_list = re.findall(INDEL, mpileup_list[4])
        if indel_list:
            tmp = Counter(indel_list).most_common()
            if len(tmp) > 1:
                return []  # multiple different indel calls
            indel_name_list = map(lambda x: x[0], tmp)
            indel_count_list = map(lambda x: x[1], tmp)
            rm_indel = tmp[0][0]
            result = []
            if parse_indel(rm_indel, result) == 0:  # indel reads followed by a SNV read, e.g. ....-2AAAA.....
                return []
            mpileup_list[4] = mpileup_list[4].replace(rm_indel, "I")
    # remove head (^) and tail ($)
    mpileup_list[4] = remove_head_end(mpileup_list[4], int(mpileup_list[3]), mpileup_list[6], mpileup_list[0],
                                      int(mpileup_list[1]))

    # filter out alignments that do not meet the requirements of mapQ
    data_filter = [ord(i) - PHREDSCORE >= rm_minmap_q for i in mpileup_list[6]]
    if not is_gh:
        base_quality = map(lambda x: ord(x) - PHREDSCORE,
                           list(compress(mpileup_list[5], data_filter)))  # type: list

    read_base_with_i = compress_read_base(mpileup_list[4], data_filter)  # type: str
    read_base_with_d = "".join(
        [mpileup_list[4][i] if data_filter[i] else "D" for i in xrange(len(mpileup_list[4]))])  # type: str
    read_base_without_i = re.sub("I", "", read_base_with_i)  # type: str
    if len(read_base_without_i) < rm_mindepth:
        return []

    if is_gh:
        return read_base_without_i

    maxvarp = ["A", "C", "G", "T", "I"]
    maxvar = [read_base_with_i.count(element) for element in maxvarp]
    i_index_max = maxvar.index(max(maxvar))
    if 4 == i_index_max:
        if len(rm_indel) == 0:
            logging.critical("has I but no rm_indel name = {0} pos = {1}".format(mpileup_list[0], mpileup_list[1]))
        maxvarp[4] = rm_indel

    if maxvar[i_index_max] < rm_minvar:
        return -1  # return reference genome type and for skipping this line
    r_num, v_num = get_reference_variant_allele_num(read_base_with_d)
    num_of_i = len(read_base_with_i) - len(read_base_without_i)
    read_base_list_final = rebuild_read_base_list(read_base_without_i, indel_name_list, indel_count_list)
    base_quality.extend([DEFAULTBASEQUALITY for i in xrange(num_of_i)])  # type: list

    rm_out = PileupDataStruct(mpileup_list[0],  # name str
                              mpileup_list[1],  # pos  int
                              mpileup_list[2],  # reference_base  str
                              maxvarp[i_index_max],  # variant  str
                              r_num,  # reference_allele_num  int
                              v_num,  # variant_allele_num  int
                              1,  # genotype_num  int
                              map(lambda x: ReadInfo(read_base_list_final[x], base_quality[x]),
                                  xrange(len(read_base_list_final))),
                              1)  # genotype class
    return rm_out


def choose(candidate_index_list, current_basis, choice_out, basis_filter):
    """

    choose indexs from the candidate_index_list according to current_basis
    :param candidate_index_list: the length equal to the number of 'True's in basis_filter
    :type values: list[int]
    """
    values = list(compress(current_basis, basis_filter))  # compressed basis
    tmp = Counter(values).most_common()
    value_list = map(lambda x: x[0], tmp)  # value
    count_list = map(lambda x: x[1], tmp)  # count of value
    sorted_value_list = copy.copy(value_list)
    sorted_value_list.sort(reverse=True)  # descending sorted values
    candidate = []
    if len(choice_out) >= 2:
        return candidate
    if count_list[value_list.index(sorted_value_list[0])] == 1:  # Only one of the maximum
        choice_out.append(candidate_index_list[values.index(sorted_value_list[0])])
        basis_filter[candidate_index_list[values.index(sorted_value_list[0])]] = False

        if len(choice_out) < 2:
            if count_list[value_list.index(sorted_value_list[1])] == 1:  # only one second largest value
                choice_out.append(candidate_index_list[values.index(sorted_value_list[1])])
                basis_filter[candidate_index_list[values.index(sorted_value_list[1])]] = False
            else:
                candidate.extend(map(lambda y: candidate_index_list[y],
                                     filter(lambda x: values[x] == sorted_value_list[1], xrange(len(values)))))

    else:
        candidate.extend(map(lambda y: candidate_index_list[y],
                             filter(lambda x: values[x] == sorted_value_list[0], xrange(len(values)))))
    for i in xrange(len(basis_filter)):
        if i not in candidate:
            basis_filter[i] = False
    return candidate


def choose_random(candidate, num):
    result = []
    for i in xrange(num):
        # index_random = random.randint(0, len(candidate) - 1)
        # result.append(candidate[index_random])
        # del candidate[index_random]
        result.append(candidate[0])
        del candidate[0]
    return result


def select_variants(candidate_in, basis, should_log):
    """
    Select 2 indexes from candidate indexes according to basis.
    If basis doesn't work, choose 2 randomly.
    :param candidate_in: the candidate original list
    :param basis: [reads list, mapq list, baseq list] reads list, mapq list, baseq list has the same length as original list
    :return: [index1, index2]
    """

    choice = []
    candidate_index = [i for i in xrange(len(candidate_in))]
    basis_filter = [True for i in candidate_in]
    for current_basis in basis:
        candidate_index = choose(candidate_index, current_basis, choice, basis_filter)
        if len(choice) >= 2:
            break
    if len(choice) < 2:
        choice.extend(choose_random(candidate_index, 2 - len(choice)))
    return choice[:2]


def get_variant_info(v_name_list, basis, read_base_with_i, should_log, geno_class):
    # type: (list[str], list[list[int]], str, str) -> list
    """
    calculate the variant's name_str, ref_num and alt_num
    if variant is more than 2, choose 2 according to basis. If basis doesn't work, choose randomly.
    :param basis: [reads, map_q, base_q]
    :param v_name_list:
    :param read_base_with_i:
    :return: [variant_str, ref_num, alt_num]
    """

    def get_reference_variant_allele_num_vcf(read_base):
        v_num = 0
        r_num = 0
        for i in xrange(len(read_base)):
            if read_base[i] in [".", ","]:
                r_num += 1
            if read_base[i] in ["A", "T", "C", "G"] or "+" in read_base[i] or "-" in read_base[i]:
                v_num += 1
        return [r_num, v_num]

    if not v_name_list:  # is empty
        return ["", 0, 0]

    if len(v_name_list) == 1:
        r_num, v_num = get_reference_variant_allele_num_vcf(read_base_with_i)
        return [v_name_list[0], r_num, v_num]

    if len(v_name_list) >= 2:

        index1, index2 = select_variants(v_name_list, basis, should_log)
        if geno_class == 2:
            return ["{0},{1}".format(v_name_list[index1], v_name_list[index2]),
                    basis[0][index1],
                    basis[0][index2]]
        else:
            ref_num = read_base_with_i.count(".") + read_base_with_i.count(",")

            return [v_name_list[index1], ref_num, basis[0][index1]]


def split_readbase(readbase_with_i, indel_list):
    """
    ignore the ',' and '.' before indel
    ..+1A ---> ['.', '+1A']
    :param readbase_with_i:
    :param indel_list:
    :return:
    """
    tmp_list = []
    index = 0
    for element in readbase_with_i:
        if element != "I":
            tmp_list.append(element)
        else:
            tmp_list[-1] = indel_list[index]
            index += 1
    # logging.debug("22")
    return tmp_list


def sort_name(v_name_list, v_count_list, v_map_q_list, v_base_q_list):
    def bigger_than_second(count1, mapq1, baseq1, count2, mapq2, baseq2):
        if count1 > count2:
            return True
        if count1 < count2:
            return False
        if mapq1 > mapq2:
            return True
        if mapq1 < mapq2:
            return False
        if baseq1 >= baseq2:
            return True
        return False

    def exchange(v_name_list, v_count_list, v_map_q_list, v_base_q_list, index1, index2):
        tmp_name = v_name_list[index1]
        tmp_count = v_count_list[index1]
        tmp_mapq = v_map_q_list[index1]
        tmp_baseq = v_base_q_list[index1]

        v_name_list[index1] = v_name_list[index2]
        v_count_list[index1] = v_count_list[index2]
        v_map_q_list[index1] = v_map_q_list[index2]
        v_base_q_list[index1] = v_base_q_list[index2]

        v_name_list[index2] = tmp_name
        v_count_list[index2] = tmp_count
        v_map_q_list[index2] = tmp_mapq
        v_base_q_list[index2] = tmp_baseq

    for i in xrange(len(v_name_list)):
        for j in xrange(len(v_name_list)):
            index = len(v_name_list) - 1 - j - i
            if index == 0:
                break
            if bigger_than_second(v_count_list[index], v_map_q_list[index], v_base_q_list[index],
                                  v_count_list[index - 1], v_map_q_list[index - 1], v_base_q_list[index - 1]):
                exchange(v_name_list, v_count_list, v_map_q_list, v_base_q_list, index, index - 1)
    return v_name_list


def read_mpileup_vcf(mpileup_list, rm_minvar, rm_minmap_q, rm_mindepth):
    """
    screenfor mutations for vcf format output, filter out alignments that do not meet the requirements of mapQ, and convert the data format
    :type rm_mindepth: int
    :type rm_minmap_q: int
    :type rm_minvar: int
    :type mpileup_list: list
    :param mpileup_list: one line of mpileup file
    :param rm_minvar: variant supporting reads.
    :param rm_minmap_q: mapQ
    :param rm_mindepth: minimum reads
    :return:
    not empty PileupDataStruct: data
    -1: return reference genome type and for skipping this line, but should be count in coverage
    []: should be skip.
    """

    def get_q_list(v_name_list, v_count_list, read_base_list_with_i, map_q, base_quality):
        """
        calculate the mapq sum and baseq sum of all kinds of variants.
        indel is calculated according to default
        :param v_name_list: variant's names
        :param v_count_list: variant's reads number
        :param read_base_without_i: read_base without indel
        :param map_q: mapq list
        :param base_quality: baseq list
        :param default_map_q: default mapq for indel
        :param default_base_q: default baseq for indel
        :return: [[mapq sum of all kinds of variants], [baseq sum of all kinds of variants]]
        """
        map_q_list = []
        base_q_list = []
        for i in xrange(len(v_name_list)):
            q_filter = [j == v_name_list[i] for j in read_base_list_with_i]
            map_q_list.append(sum(list(compress(map_q, q_filter))))
            base_q_list.append(sum(list(compress(base_quality, q_filter))))

        return [map_q_list, base_q_list]

    counter = 0
    if "*" in mpileup_list[4]:
        return []  # indel read overlaps
    indel_name_list = []
    indel_count_list = []
    indel_list = []
    if ("-" in mpileup_list[4]) or ("+" in mpileup_list[4]):
        indel_list = re.findall(INDEL, mpileup_list[4])
        if indel_list:
            while 1:  # check data: eliminate the invalid indel like "+9AA", split the indel and the mutation together, for example, convert +2AACT to +2AA
                if counter == len(indel_list):
                    break
                result = []
                ret = parse_indel(indel_list[counter], result)
                if ret < 0:
                    indel_list.remove(indel_list[counter])
                elif ret == 0:
                    indel_list[counter] = "".join(result)
                    counter += 1
                else:
                    counter += 1
            indel_count_list = Counter(indel_list).most_common()
            indel_name_list = map(lambda x: x[0], indel_count_list)
            indel_count_list = map(lambda x: x[1], indel_count_list)
            for i in indel_name_list:
                mpileup_list[4] = mpileup_list[4].replace(i, "I")

    # remove head(^) and tail($)
    mpileup_list[4] = remove_head_end(mpileup_list[4], int(mpileup_list[3]), mpileup_list[6], mpileup_list[0],
                                      int(mpileup_list[1]))
    # logging.debug("111 readbase = {0} indel_list={1}".format(mpileup_list[4], indel_list))
    readbase_list_with_i = split_readbase(mpileup_list[4], indel_list)
    # logging.debug("222")
    # screen according to mapQ
    map_q_filter = [ord(i) - PHREDSCORE >= rm_minmap_q for i in mpileup_list[6]]
    base_quality = map(lambda x: ord(x) - PHREDSCORE, list(compress(mpileup_list[5], map_q_filter)))
    map_q = map(lambda x: ord(x) - PHREDSCORE, list(compress(mpileup_list[6], map_q_filter)))

    read_base_list_with_i = list(compress(readbase_list_with_i, map_q_filter))
    # read_base_with_i = compress_read_base(mpileup_list[4], map_q_filter)
    # read_base_with_d = "".join([mpileup_list[4][i] if map_q_filter[i] else "D" for i in
    #                             xrange(len(
    #                                 mpileup_list[4]))])  # Replace the compressed part with "D" to calculate allele num
    read_base_list_without_i = list(
        compress(read_base_list_with_i, ["+" not in i and "-" not in i for i in read_base_list_with_i]))
    # read_base_without_i = re.sub("I", "", read_base_with_i)
    if len(read_base_list_with_i) < rm_mindepth:
        return []  # lessmindepth
    v_name_list = ["A", "C", "G", "T"]
    v_count_list = [read_base_list_without_i.count(element) for element in v_name_list]

    name_filter = [i > 0 for i in v_count_list]
    v_name_list = list(compress(v_name_list, name_filter))
    v_count_list = list(compress(v_count_list, name_filter))

    v_name_list.extend(indel_name_list)
    v_count_list.extend(indel_count_list)

    if not v_count_list:
        return -1
    # if max(v_count_list) < rm_minvar:
    #     return -1  # reference genome type, skip this line

    [v_map_q_list, v_base_q_list] = get_q_list(v_name_list, v_count_list,
                                               read_base_list_with_i, map_q,
                                               base_quality)

    geno_class = get_geno_class(read_base_list_with_i, map_q, base_quality, v_count_list, v_map_q_list,
                                v_base_q_list)

    [variant_str, r_num, v_num] = get_variant_info(v_name_list, [v_count_list, v_map_q_list,
                                                                 v_base_q_list],
                                                   read_base_list_with_i,
                                                   False, geno_class)
    sort_name(v_name_list, v_count_list, v_map_q_list, v_base_q_list)

    var_all_str = ",".join(v_name_list)
    if mpileup_list[1] == 16647426:  # todo
        logging.debug("geno_class={0} r_num={1} v_num={2}".format(geno_class, r_num, v_num))

    return PileupDataStruct(mpileup_list[0],  # name str
                            mpileup_list[1],  # pos  int
                            mpileup_list[2],  # reference_base  str
                            variant_str,  # variant  str
                            r_num,  # reference_allele_num  int
                            v_num,  # variant_allele_num  int
                            len(v_name_list),  # genotype_num  int
                            map(lambda x: ReadInfo(list(readbase_list_with_i)[x], base_quality[x]),
                                xrange(len(readbase_list_with_i))),
                            None,
                            geno_class,  # genotype class
                            var_all_str)


def second_large_num(num_list):
    # type: (list) -> int
    tmp_list = sorted(num_list)
    return tmp_list[-2]


def get_geno_class(read_base_list_with_i, map_q, base_quality, v_count_list, v_map_q_list, v_base_q_list):
    # type: (list, list, list, list, list, list, list) -> int

    if len(v_count_list) < 2:
        return 1
    second_num = second_large_num(v_count_list)
    ref_num = read_base_list_with_i.count(".") + read_base_list_with_i.count(",")
    if ref_num > second_num:
        return 1
    if ref_num < second_num:
        return 2

    # ref_num == second_large_num
    read_base_filter = [i in [".", ","] for i in read_base_list_with_i]
    sum_ref_map_q = sum(list(compress(map_q, read_base_filter)))
    v_map_q = v_map_q_list[v_count_list.index(second_num)]
    if sum_ref_map_q > v_map_q:
        return 1
    if sum_ref_map_q < v_map_q:
        return 2

    # sum_ref_map_q == v_map_q
    sum_ref_base_q = sum(list(compress(base_quality, read_base_filter)))
    v_base_q = v_base_q_list[v_count_list.index(second_num)]
    if sum_ref_base_q >= v_base_q:
        return 1
    if sum_ref_base_q < v_base_q:
        return 2


# @send_name_time(should_send=should_analyze, q7=q7)
def golden_hetero(mpileup_list, min_map_q, min_depth, snp_type, read_depth, worker_id):
    """
    handle the data after vcf screen, store the result in q2
    :param mpileup_list:
    :param read_depth: Known heterogous SNP call required read depths. RD
    :param snp_type: Known snp type input for known heterogous SNP call
    :param min_depth: minimum reads
    :param min_map_q: mapQ
    :return: None  should skip
    """
    read_base = read_mpileup(mpileup_list, 0, min_map_q, min_depth, True, worker_id)  # type: str
    if not read_base or read_base == -1:
        return None
    gh_xrefcount = len(filter(lambda x: x in [mpileup_list[2], ".", ","], read_base))
    gh_xaltcount = len(read_base) - gh_xrefcount

    if snp_type == 'dbsnp' and (gh_xaltcount < 2 or gh_xrefcount < 2 or
                                (gh_xrefcount + gh_xaltcount < read_depth)):
        return None
    return DataInQ2(mpileup_list[0],  # type: str
                    mpileup_list[1],  # type: int
                    gh_xrefcount,
                    gh_xaltcount)


def window_data_one_chromosome(center_pos, data_list, lamb, is_list_ended, current_pos, worker_id):
    """
    window_data from one chromosome data
    :type current_pos: int
    :param current_pos: current handling position.
    :param is_list_ended:
    :rtype: BigForewordList
    :type lamb: int
    :type data_list: BigForewordList
    :type center_pos: int
    :param center_pos: pos of data in q3
    :param data_list: data after GH in q2(DataInQ2)
    :param lamb: half of window width
    :return:
        None can not window yet.
        others data windowed
    """
    if data_list.is_empty():  # is empty
        if is_list_ended:
            return []
        else:
            if current_pos - center_pos > lamb:
                return []
            return None

    # now data_list is not empty, looking for left edge
    if data_list.get_last_element().coordinate_1 < center_pos - lamb:
        if is_list_ended:
            return []
        if current_pos - center_pos > lamb:
            return []
        return None
    counter = 0  # type: int
    len_data_list = data_list.len()
    while 1:
        if counter >= len_data_list:
            if is_list_ended:
                return []
            else:
                if current_pos - center_pos > lamb:
                    return []
                return None
        if data_list.get_current_element().coordinate_1 >= center_pos - lamb:
            left_edge = data_list.get_current_element()  # type: DataInQ2
            if left_edge.coordinate_1 > center_pos + lamb:
                return []
            break
        else:
            data_list.move_foreword()
        counter += 1

    # right edge
    if data_list.get_current_element().coordinate_1 > center_pos + lamb:
        return []
    counter = 0  # type: int
    while 1:
        if counter == len_data_list:
            return []
        if data_list.get_element(-1 - counter).coordinate_1 <= center_pos + lamb:
            if is_list_ended or counter > 0:
                right_edge = data_list.get_element(-1 - counter)
                break
            else:
                if current_pos - center_pos > lamb:
                    right_edge = data_list.get_last_element()  # type: DataInQ2
                    break
                else:
                    return None
        else:
            counter += 1

    # double check edges
    if int(right_edge.coordinate_1) < int(left_edge.coordinate_1):
        return []
    return data_list.filter(lambda x: left_edge.coordinate_1 <= x.coordinate_1 <= right_edge.coordinate_1)


def get_gq(ra_u, rm_u, mm_u):
    # type: (float, float, float) -> str
    while 1:
        value_list = [10 ** ra_u, 10 ** rm_u, 10 ** mm_u]
        if (sum(value_list) - min(value_list)) != 0:
            break
        else:
            ra_u += 1
            rm_u += 1
            mm_u += 1

    value = 1 - max(value_list) / (sum(value_list) - min(value_list))
    if value == 0:
        return "150"
    else:
        return str(int(round(-10 * math.log10(value))))


def round_3(value):
    """
    prevent loss of precision, retaining to 3 decimal places
    :param value:
    :return:
    """
    ret = round(value * 1000) / 1000.0
    return ret


def get_pl(rr_u, ra_u, rm_u, mm_u):
    result = "{3},{0},{1},{2}".format(int(round(-10 * ra_u)),
                                      int(round(-10 * rm_u)),
                                      int(round(-10 * mm_u)),
                                      int(round(-10 * rr_u)))
    return result


def differential(q2_source, data_in_q3_curr, q5, lamb,
                 default_bias, artifact, worker_id,
                 q2_list_is_end, current_pos, my_format):
    """
    window the data in q2, send result to q5.
    :type my_format: str
    :param my_format:
    :type q2_list_is_end: bool
    :param q2_list_is_end:
    :type worker_id: int
    :param worker_id:
    :param artifact:
    :param default_bias:
    :param q5:
    :type current_pos: int
    :param current_pos: current data pos from source
    :type data_in_q3_curr: PileupDataStruct
    :param data_in_q3_curr:
    :param q2_source: data queue to be windowed
    :type q2_source: BigForewordList
    :param lamb: half of window width
    :type lamb: int
    :return:
    0 data is windowed, send result to W
    2 can not window yet
    """

    tracked_data = window_data_one_chromosome(data_in_q3_curr.coordinate_1, q2_source, lamb, q2_list_is_end,
                                              current_pos, worker_id)  # type: List[DataInQ2]
    if tracked_data is None:
        return 2
    bias = bias_estimator(data_in_q3_curr.coordinate_1, tracked_data, lamb, default_bias)
    [rr_u, ra_u, rm_u, mm_u] = sc_caller(data_in_q3_curr, bias, artifact)  # 10 11 12 13

    rr = 10 ** rr_u  # type: float
    ra = 10 ** ra_u  # type: float
    rm = 10 ** rm_u  # type: float
    mm = 10 ** mm_u  # type: float

    if data_in_q3_curr.coordinate_1 == 16647426:  # todo
        logging.debug(
            "pos={6} rr_u={7} ra_u={0} rm_u={1} mm_u={2} rr={8} ra={3} rm={4} mm={5}".format(ra_u, rm_u, mm_u, ra,
                                                                                             rm, mm,
                                                                                             data_in_q3_curr.coordinate_1,
                                                                                             rr_u,
                                                                                             rr))

    if my_format == VARCALLFORMAT:
        vcf_info = None
    else:

        gq = get_gq(ra_u, rm_u, mm_u)  # type: str

        vcf_info = VcfInfo("",
                           "{0},{1}".format(data_in_q3_curr.reference_allele_num,
                                            data_in_q3_curr.variant_allele_num),  # AD
                           round_3(bias),  # BI
                           gq,  # GQ
                           get_pl(rr_u, ra_u, rm_u, mm_u),  # type:str  # PL
                           gq,  # QUAL
                           data_in_q3_curr.genotype_num,
                           data_in_q3_curr.genotype_class)

    outline = OutLineStruct(data_in_q3_curr.chromosome_name,  # type: str
                            data_in_q3_curr.coordinate_1,  # type: int
                            data_in_q3_curr.reference_base,  # type: str
                            data_in_q3_curr.variant,  # type: str
                            data_in_q3_curr.reference_allele_num,  # type: int
                            data_in_q3_curr.variant_allele_num,  # type: int
                            bias,  # type: float
                            rr,  # sn type:float
                            ra,  # ae type:float
                            rm,  # he type:float
                            mm,  # ho type:float
                            vcf_info,
                            data_in_q3_curr.so,
                            data_in_q3_curr.variant_all)
    q5.put(DataInQ5(outline, worker_id, WORKVAR if my_format == VARCALLFORMAT else WORKVCF, []),
           block=False)
    return 0


def read_bulk_mpileup(mpileup_list, rm_minvar, rm_minmapQ, rm_mindepth):
    if "*" in mpileup_list[4]:
        return "indel"
    if ("-" in mpileup_list[4]) or ("+" in mpileup_list[4]):
        if re.findall(INDEL, mpileup_list[4]):
            return "indel"

    # remove head(^) and tail($)
    mpileup_list[4] = remove_head_end(mpileup_list[4], int(mpileup_list[3]), mpileup_list[6], mpileup_list[0],
                                      int(mpileup_list[1]))
    # mapQ
    map_q_filter = [ord(i) - PHREDSCORE >= rm_minmapQ for i in mpileup_list[6]]
    read_base = compress_read_base(mpileup_list[4], map_q_filter)
    maxvar = [read_base.count('A'), read_base.count('C'), read_base.count('G'), read_base.count('T')]

    if len(read_base) < rm_mindepth:
        return "lessmindepth"
    elif max(maxvar) < rm_minvar:
        return "refgenotype"  # return reference genome type and for skipping this line
    else:
        return "varreads"


def get_so_source(bulk_pileup_source, bulk_minvar, bulk_min_mapq, bulk_mindepth):
    # type: (GeneratorExit, int, int, int) -> GeneratorExit
    """
    When pos in bulk greater than data_pos, should not continue to read bulk.
    should use the current bulkpos to compare the next data_pos.
    Build the generator, retain the bulkpos data, and achieve the above functions
    :param bulk_pileup_source: bulk generator
    :return:
        Datapos is not in bulk or bulk has no data: noCoverageInControl
        Datapos in bulk does not have enough depth: lessmindepth
        Datapos is ref in bulk: refgenotype
        Datapos is var in bulk: varreads
        -1 Pysam crushed and should be recalculated
    """
    should_read = True
    pos = -1
    while 1:
        if should_read:
            bulk_pileup_list = bulk_pileup_source.next()
            if bulk_pileup_list == -1:  # pysam crashed
                logging.info("Use samtools engine instead!")
                yield -1
            if bulk_pileup_list is None:
                while 1:
                    pos = yield "noCoverageInControl2"

        if int(bulk_pileup_list[1]) < pos:
            should_read = True
            continue
        elif int(bulk_pileup_list[1]) > pos:
            should_read = False
            pos = yield "noCoverageInControl"

        else:
            bulk_pileup_list[4] = bulk_pileup_list[4].upper()
            should_read = True
            pos = yield read_bulk_mpileup(bulk_pileup_list, bulk_minvar, bulk_min_mapq, bulk_mindepth)


# @send_name_time(should_send=should_analyze, q7=q7)
def bias_estimator(pos, tracked_data, lamb, default):
    """
    calculate bias
    :type lamb: int
    :rtype: float
    :type tracked_data: list[DataInQ2]
    :type default: float
    :type pos: int
    :param pos: center of window
    :param tracked_data: windowed data
    :param lamb: half of window width
    :param default: default value
    :return:
    """

    def bias_kernel(bk_x0, bk_xi, lamb):
        """
        function D in paper
        :rtype: float
        :type bk_x0: int
        :param bk_x0: target position
        :param bk_xi: one of windowed positions
        :param lamb: half of window width
        :return:
        """
        # bk_t = float(abs(bk_x0 - bk_xi)) / lamb
        if -lamb < bk_x0 - bk_xi < lamb:
            return 0.75 * (1 - (float((bk_x0 - bk_xi)) / lamb) ** 2)
        else:
            return 0.0

    be_kwy = []
    be_kw = []
    for i in tracked_data:
        be_tmp1 = float(i.refcount) + float(i.altcount)  # type: float
        if be_tmp1 <= 0:
            continue
        be_tmp = bias_kernel(int(pos), int(i.coordinate_1), lamb)  # K in the formula
        be_tmp2 = float(max(i.refcount, i.altcount)) / be_tmp1
        be_kwy.append(be_tmp * be_tmp1 * be_tmp2)
        be_kw.append(be_tmp * be_tmp1)
    # Nadaraya-Watson kernel-weighted average
    if len(be_kwy) > 0 and sum(be_kw) > 0:
        be_a = sum(be_kwy) / sum(be_kw)
    # return args.bias when no neighboring heterozygous base
    else:
        be_a = default
    return be_a


def sc_caller(sc_candidate, sc_bias, sc_artifact):
    """
    Calculate RR, RA, RM, MM from data in q3 and bias.
    :type sc_bias: float
    :type sc_candidate: PileupDataStruct
    :param sc_candidate: data in q3
    :param sc_bias: Estimated bias in this position
    :param sc_artifact:
    :return:
    [RR, RA, RM, MM]
    """

    def P_b_GG(bp, ref, mut, f, pos):
        """
        Formula (7),calculate lg Probability value
        :type bp: ReadInfo
        :type ref: str
        :type mut: str
        :type f: float
        :rtype: float
        :param bp: data of 1 read
        :param ref: reference base
        :param mut: mismatch
        :param f: {0, 0.125, bias (or 1-bias depending on mut>ref or ref>mut), 1} for {ref/ref, ref/artifacts, ref/mut; mut/mut}
        :return: lg Probability value
        """
        e = 10 ** (-bp.base_quality / 10.0)

        if bp.read_base == mut:
            a = f * (1 - e) + (1 - f) * e / 3
        elif bp.read_base in [',', '.', ref]:
            a = f * e / 3 + (1 - f) * (1 - e)
        else:
            a = e / 3
        ret = math.log(a, 10)
        # if pos == 16211035:
        #     logging.debug("readbase={0} baseq={1} mut={2} ref={3} f={4} ret={5} e={6}".format(bp.read_base,
        #                                                                                       bp.base_quality,
        #                                                                                       mut,
        #                                                                                       ref,
        #                                                                                       f,
        #                                                                                       ret,
        #                                                                                       e)) # todo
        return ret

    if sc_candidate.coordinate_1 == 16647426:
        logging.debug("genotype_num = {0} reference_base={1} "
                      "variant={2} ref_allele_num={3} var_allele_num={4}".format(sc_candidate.genotype_num,
                                                                                 sc_candidate.reference_base,
                                                                                 sc_candidate.variant,
                                                                                 sc_candidate.reference_allele_num,
                                                                                 sc_candidate.variant_allele_num))  # todo
    # logging.debug("11 genotype_num={0} variant={1}".format(sc_candidate.genotype_num, sc_candidate.variant))
    # if sc_candidate.genotype_num == 1:
    #     ref = sc_candidate.reference_base  # type: str
    #     mut = sc_candidate.variant  # type: str
    # else:
    #     ref, mut = sc_candidate.variant.split(",")

    if "," in sc_candidate.variant:
        ref, mut = sc_candidate.variant.split(",")
    else:
        ref = sc_candidate.reference_base  # type: str
        mut = sc_candidate.variant  # type: str

    f = sc_bias  # type: float

    # if sc_candidate.coordinate_1 == 16211035:
    #     logging.debug("ref_allele_num={0} var_allele_num={1} bias={2}".format(sc_candidate.reference_allele_num,
    #                                                                           sc_candidate.variant_allele_num,
    #                                                                           sc_bias))  # todo
    #     logging.debug("len read_info_list={}".format(len(sc_candidate.read_info_list)))
    # for i in sc_candidate.read_info_list:
    #     logging.debug("base_quality={0} read_base={1}".format(i.base_quality, i.read_base))
    if sc_candidate.reference_allele_num > sc_candidate.variant_allele_num:
        f = 1 - f
    # if sc_candidate.coordinate_1 == 16053791:  # todo
    #     logging.debug("f={}".format(f))

    rr_u = sum(map(lambda x: P_b_GG(x, ref, mut, sc_artifact, sc_candidate.coordinate_1), sc_candidate.read_info_list))
    # if sc_candidate.coordinate_1 == 16211035:  # todo
    #     logging.debug("rr_u={}".format(rr_u))
    ra_u = sum(map(lambda x: P_b_GG(x, ref, mut, 0.125, sc_candidate.coordinate_1), sc_candidate.read_info_list))
    # if sc_candidate.coordinate_1 == 16211035:  # todo
    #     logging.debug("ra_u={}".format(ra_u))
    rm_u = sum(map(lambda x: P_b_GG(x, ref, mut, f, sc_candidate.coordinate_1), sc_candidate.read_info_list))
    # if sc_candidate.coordinate_1 == 16211035:  # todo
    #     logging.debug("rm_u={}".format(rm_u))
    mm_u = sum(map(lambda x: P_b_GG(x, ref, mut, 1.0, sc_candidate.coordinate_1), sc_candidate.read_info_list))
    # if sc_candidate.coordinate_1 == 16211035:  # todo
    #     logging.debug("mm_u={}".format(rr_u))
    result = [rr_u, ra_u, rm_u, mm_u]

    return result


def packup_send_vcf(current_data, q5, worker_id, vcf_info):
    outline = OutLineStruct(name=current_data.chromosome_name,
                            pos=current_data.coordinate_1,
                            ref=current_data.reference_base,
                            var=current_data.variant,
                            ref_num=current_data.reference_allele_num,
                            var_num=current_data.variant_allele_num,
                            vcf_info=vcf_info)
    q5.put(DataInQ5(outline, worker_id, WORKVCF, []), block=False)
    return


def get_my_filename(output, suffix_str, prefix_str):
    base = os.path.basename(output)
    if prefix_str == "/":
        prefix_str = "./"
    return prefix_str + os.path.splitext(base)[0] + suffix_str


def calculate_eta(list_var_buf, list_var_tag):
    """
    calculate eta
    :param list_var_buf: list[list[OutLineStruct]]
    :param list_var_tag: list[bool]
    :return:
    """
    allele_num_list = []
    bias_list = []
    for i in xrange(len(list_var_tag)):
        if not list_var_tag[i]:
            break
        else:
            allele_num_list.extend([j.var_num + j.ref_num for j in list_var_buf[i]])
            bias_list.extend([j.bias for j in list_var_buf[i]])
    allele_num_list = allele_num_list[:CUTOFFNUM]
    bias_list = bias_list[:CUTOFFNUM]
    LLR = []
    np.random.seed(42)
    for i in range(0, len(allele_num_list)):
        f_artifact = 0.125 * bias_list[i] / 0.5
        alt = np.random.binomial(allele_num_list[i], f_artifact)
        ref = allele_num_list[i] - alt
        L_filter = (1 - f_artifact) ** ref * ((f_artifact) ** alt)
        ## random select major allele
        major = np.random.binomial(1, .5)
        if major == 0:
            f_true = 0.5 * 0.5 / bias_list[i]
        if major == 1:
            f_true = 0.5 * bias_list[i] / 0.5
        L_true = (1 - f_true) ** ref * ((f_true) ** alt)
        ## if L_filter/true is 0, assign a very small value
        if L_filter == 0:
            L_filter = 10 ** -100
        if L_true == 0:
            L_true = 10 ** -100

        LLR.append(-math.log10(L_filter / L_true) if L_filter / L_true != 0 else 10000)
    LLR = np.array(LLR)
    co_001 = np.percentile(LLR, 99)
    result = 10 ** -co_001
    return result


def write_result(q5, var_file_output, worker_num, bulk, name, head, tail, min_var):
    try:
        write_result_inner(q5, var_file_output, worker_num, bulk, name, head, tail, min_var)
    except:
        logging.critical("Exception!")


def write_result_inner(q5, var_file_output, worker_num, bulk, name, head, tail, min_var):
    """
    Receive the data in q5, organize it, and write the result file.
    :return:
    """

    def get_current_cutoff_num(list_var_buf, list_var_tag):
        cutoff_counter = 0
        for i in xrange(len(list_var_tag)):
            if not list_var_tag[i]:
                break
            else:
                cutoff_counter += len(list_var_buf[i])
        return cutoff_counter

    def merge_coverage(list_coverage_buf):
        # list_coverage_buf: [[coverage1, coverage2...]]
        result_list = []  # [coverage1, coverage2...]  coverage: [name str ,start int ,stop int]
        for i in xrange(len(list_coverage_buf)):
            if not result_list:  # is empty
                result_list.extend(
                    list_coverage_buf[i])  # [coverage1, coverage2...]  coverage: [name str ,start int ,stop int]
                continue
            if not list_coverage_buf[i]:
                continue
            if list_coverage_buf[i][0][1] == result_list[-1][2] + 1:
                result_list[-1][2] = list_coverage_buf[i][0][2]
                result_list.extend(list_coverage_buf[i][1:])
            else:
                result_list.extend(list_coverage_buf[i])
        return result_list

    eta = -1
    logging.info("Waiting for data coming")
    current_handling_worker = 1
    list_var_buf = []  # type: list[list[DataInQ5]]
    list_var_tag = []  # worker is done
    list_coverage_buf = []  # type: list  # [[worker1 coverage1, coverage2...],[worker2 coverage1, coverage2...]...]  coverage: [name str ,start int ,stop int]
    file_reasoning = get_my_filename(var_file_output, "_{0:0>2d}to{1:0>2d}.reasoning".format(head, tail),
                                     os.path.dirname(var_file_output) + "/")

    for i in xrange(worker_num):
        list_var_buf.append([])  # Initialize data for each worker
        list_var_tag.append(False)
        list_coverage_buf.append([])
    with open(var_file_output + ".body", "a") as fp_out, open(file_reasoning, "a") as fp_reasoning, open(
            var_file_output + ".eta", "a") as fp_eta:
        while 1:
            msg_q5 = q5.get(block=True)  # type: DataInQ5

            if msg_q5.worker_id == 0:  # from main
                # try to write coverage file
                merged_coverage_list = merge_coverage(
                    list_coverage_buf)  # type: list # [coverage1, coverage2...]  coverage: [name str ,start int ,stop int]
                if merged_coverage_list:
                    with open(get_my_filename(var_file_output, "_{0:0>2d}to{1:0>2d}.coverage".format(head, tail),
                                              os.path.dirname(var_file_output) + "/"),
                              "a") as fp:
                        fp.write("\n".join(
                            map(lambda z: "\t".join(z), map(lambda x: map(lambda y: str(y), x), merged_coverage_list))))
                if eta == -1:
                    if get_current_cutoff_num(list_var_buf, list_var_tag) > 0:
                        eta = calculate_eta(
                            [[i.outlineStruct for i in one_worker_data] for one_worker_data in list_var_buf],
                            list_var_tag)
                        for i in xrange(worker_num):
                            if list_var_buf[i]:
                                write_do([element.outlineStruct for element in list_var_buf[i]], msg_q5.work_type,
                                         i + 1, fp_out, eta, bulk, fp_reasoning, min_var)
                            else:
                                logging.info("\x1b[1;33mworker{0} write len = 0\x1b[0m".format(i + 1))
                    else:
                        logging.info("\x1b[1;33mNothing need to be written.\x1b[0m".format(i + 1))
                fp_eta.write("##contig=<ID={0},eta={1}>\n".format(name, eta))
                break
            if msg_q5.work_type == WORKVAR or msg_q5.work_type == WORKVCF:

                # record tag and data
                if msg_q5.outlineStruct == DONEMSG:
                    list_var_tag[msg_q5.worker_id - 1] = True
                    logging.info("\x1b[1;33mchr {0} worker{1} done\x1b[0m".format(name, msg_q5.worker_id))
                else:
                    list_var_buf[msg_q5.worker_id - 1].append(msg_q5)  # record data

                if msg_q5.outlineStruct == "pysam crashed":
                    logging.info("\x1b[1;33mworker{} pysam crashed. Data cleaned\x1b[0m".format(msg_q5.worker_id))
                    list_var_tag[msg_q5.worker_id - 1] = False
                    list_var_buf[msg_q5.worker_id - 1] = []
                    list_coverage_buf[msg_q5.worker_id - 1] = []
                # Try to calculate eta, write data
                if msg_q5.outlineStruct == DONEMSG:
                    if get_current_cutoff_num(list_var_buf, list_var_tag) >= CUTOFFNUM and eta == -1:
                        eta = calculate_eta(
                            [[i.outlineStruct for i in one_worker_data] for one_worker_data in list_var_buf],
                            list_var_tag)
                    if eta != -1:
                        for i in xrange(worker_num):
                            if i < current_handling_worker - 1:
                                continue
                            if list_var_tag[i]:
                                if list_var_buf[i]:
                                    write_do([element.outlineStruct for element in list_var_buf[i]], msg_q5.work_type,
                                             i + 1, fp_out, eta, bulk, fp_reasoning, min_var)
                                else:
                                    logging.info("\x1b[1;33mchr {1} worker{0} write len = 0\x1b[0m".format(i + 1, name))
                                current_handling_worker += 1
                            else:
                                break
            else:  # coverage
                list_coverage_buf[msg_q5.worker_id - 1].append(msg_q5.coverage)
    if bulk == "" or msg_q5.work_type == WORKVCF:
        os.remove(file_reasoning)
    logging.info("\x1b[1;33mQuit!\x1b[0m")
    return 0


def write_do(data_list, work_type, worker_id, fp_out, eta, bulk, fp_reasoning, min_var):
    """

    :type data_list: list[OutLineStruct]
    """

    def pack_vcf_str(outline, eta, min_var):
        """
        Assemble data of vcf according to the outline
        :type outline: OutLineStruct
        """

        def transform_vcf_so(so_in):
            if so_in == "refgenotype":
                return "True"
            elif so_in == "varreads":
                return "False"
            else:
                return "NA"

        def get_gt(outline, eta):
            AE = outline.ae / eta
            SN = 2 * outline.sn
            HE = outline.he
            HO = outline.ho
            likelihood_list = [AE, SN, HE, HO]

            if outline.pos == 16647426:
                logging.debug("eta={4} [AE, SN, HE, HO]=[{0}, {1}, {2}, {3}]".format(AE, SN, HE, HO, eta))  # todo
            max_likelihood = max(likelihood_list)
            if max_likelihood in [AE, SN]:
                result_str = "0/0"
            elif HE == max_likelihood:
                if 7 * outline.var_num > outline.ref_num:
                    result_str = "0/1"
                else:
                    result_str = "0/0"
            else:
                result_str = "1/1"

            if outline.vcf_info.genotype_class == 2:
                if result_str == "0/0":
                    result_str = "1/1"
                elif result_str == "0/1":
                    result_str = "1/2"
                else:
                    result_str = "1/1"
            return result_str

        # def get_gt(outline, eta):
        #     """
        #     Calculate gt according to eta
        #     :param outline:
        #     :param eta:
        #     :return:
        #     """
        #
        #     if outline.pos == 16053791:
        #         logging.debug("pos={0} eta={1} ae={2} he={3} ho={4} "
        #                       "sn={7} var_num={5} ref_num={6}".format(outline.pos,
        #                                                               eta,
        #                                                               outline.ae,
        #                                                               outline.he,
        #                                                               outline.ho,
        #                                                               outline.var_num,
        #                                                               outline.ref_num,
        #                                                               outline.sn))
        #     max_value = max(outline.ho, outline.ae, outline.he)
        #     if (outline.he != 0 and outline.ae / outline.he < eta) \
        #             and 7 * outline.var_num > outline.ref_num \
        #             and 2 * outline.sn < outline.he == max_value:
        #         if outline.vcf_info.genotype_class == 1:
        #             genotype_str = "0/1"
        #             # if outline.vcf_info.genotype_num > 1:
        #             #     genotype_str = "1/2"
        #         else:
        #             genotype_str = "1/2"
        #         if outline.pos == 16053791:
        #             logging.debug("16114636 AAA")
        #     elif outline.ho > outline.ae:
        #         if outline.vcf_info.genotype_class == 1:
        #             genotype_str = "1/1"
        #         else:
        #             genotype_str = "1/1"
        #         if outline.pos == 16053791:
        #             logging.debug("16114636 BBB")
        #     else:
        #         if outline.he != max_value:
        #
        #             if outline.ho == max_value:
        #                 genotype_str = "1/1"
        #             if outline.ae == max_value:
        #                 genotype_str = "0/0"
        #         else:
        #             if outline.vcf_info.genotype_class == 1:
        #                 genotype_str = "0/1"
        #             else:
        #                 genotype_str = "1/2"
        #         if outline.pos == 16053791:
        #             logging.debug("16114636 CCC")
        #
        #     return genotype_str

        outline.vcf_info.gt = get_gt(outline, eta)
        format_str = "GT:SO:AD:BI:GQ:PL" if outline.so != "" else "GT:AD:BI:GQ:PL"  # format
        if outline.so != "":
            cell1_str = "{0}:{1}:{2}:{3}:{4}:{5}".format(outline.vcf_info.gt, transform_vcf_so(outline.so),
                                                         outline.vcf_info.ad, outline.vcf_info.bi,
                                                         outline.vcf_info.gq, outline.vcf_info.pl)
        else:
            cell1_str = "{0}:{1}:{2}:{3}:{4}".format(outline.vcf_info.gt, outline.vcf_info.ad,
                                                     outline.vcf_info.bi, outline.vcf_info.gq,
                                                     outline.vcf_info.pl)
        if outline.pos == 16647426:
            logging.debug("pos={0} genotype_num={1} "
                          "ref={2} ref_num={3} "
                          "var={4} var_num={5} "
                          "genotype_class={6}".format(outline.pos, outline.vcf_info.genotype_num,
                                                      outline.ref, outline.ref_num,
                                                      outline.var, outline.var_num,
                                                      outline.vcf_info.genotype_class))  # todo

        def get_comment(var_str, genotype_num, min_var, var_num):
            comment = ""
            if len(var_str.split(",")) < genotype_num:
                comment = comment + MULTIPLEGENOTYPE
            if var_num < min_var:
                if len(comment) == 0:
                    comment = "{0}{1}".format(NOTENOUGHVARIANTS, min_var)
                else:
                    comment = comment + ",{0}{1}".format(NOTENOUGHVARIANTS, min_var)
            if len(comment) == 0:
                comment = "."
            return comment

        result = "\t".join([outline.name,  # name
                            str(outline.pos),  # position
                            ".",  # ID
                            outline.ref,  # ref
                            # outline.var,
                            outline.var_all,  # alt
                            outline.vcf_info.qual,  # qual
                            # MULTIPLEGENOTYPE if len(outline.var.split(",")) < outline.vcf_info.genotype_num else ".",
                            get_comment(outline.var, outline.vcf_info.genotype_num, min_var,
                                        outline.var_num if outline.vcf_info.genotype_class == 1 else outline.ref_num),
                            # filter
                            "NS=1",  # info
                            format_str,
                            cell1_str])
        return result

    def pack_varcall_str(outline):
        # type: (OutLineStruct) -> str
        return "{0}\t{1}\t{1}\t{2}\t{3}\t{4}\t" \
               "{5}\tPASS\t{6}\t{7}\t{8}\t{9}\t{10}".format(outline.name, outline.pos, outline.ref,
                                                            outline.var, outline.ref_num, outline.var_num,
                                                            outline.bias, outline.sn, outline.ae,
                                                            outline.he, outline.ho)

    if work_type == WORKVAR:
        my_result_filter = [i.he != 0 and i.ae / i.he < eta and 7 * i.var_num > i.ref_num and 2 * i.sn < i.he for i in
                            data_list]
        data_list = list(compress(data_list, my_result_filter))
        logging.debug("\x1b[1;33mworker{0} write len = {1}\x1b[0m".format(worker_id, len(data_list)))
        fp_out.write("\n".join(map(lambda x: pack_varcall_str(x), data_list)))
        if len(data_list) > 0:
            fp_out.write("\n")
        if bulk != "" and data_list:
            fp_reasoning.write("\n".join(
                map(lambda x: "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(x.name, x.pos, x.ref, x.var, x.ref_num,
                                                                         x.var_num, x.so),
                    data_list)))
            fp_reasoning.write("\n")
    elif work_type == WORKVCF:
        logging.info("\x1b[1;33mworker{0} write len = {1}\x1b[0m".format(worker_id, len(data_list)))
        fp_out.write("\n".join(map(lambda x: pack_vcf_str(x, eta, min_var), data_list)))
        fp_out.write("\n")


def data_generator(bam, fasta, name, start, stop, engine, is_bulk):
    # type: (str, str, str, int, int, str) -> GeneratorExit
    def data_generator_pysam(bam_file, fasta_file, name, start_pos, stop, is_bulk):
        """

        :param stop:
        :param name:
        :param start_pos:
        :type fasta_file: str
        :type bam_file: str
        """
        read_bases_list = []

        bam_file = pysam.AlignmentFile(bam_file, "rb")
        fasta_file = pysam.FastaFile(fasta_file)

        str_ref = fasta_file.fetch(name, start_pos - 1, stop + 1)

        my_arg = {"fastafile": fasta_file, "stepper": "samtools", "adjust_capq_threshold": 50, "contig": name,
                  "start": start_pos, "stop": stop, "min_mapping_quality": 0 if is_bulk else 40}
        for pileup_column in bam_file.pileup(**my_arg):
            pos = pileup_column.reference_pos + 1
            if pos > stop:
                break
            if pos < start_pos:
                continue
            try:
                read_bases_list = pileup_column.get_query_sequences(mark_matches=True, mark_ends=True, add_indels=True)
            except Exception as e:
                logging.debug("Pysam crashed! Unexpected Error: {}".format(e))
                yield -1
            read_bases = ''.join(read_bases_list)
            if len(read_bases) == 0:
                read_bases = "*"

            base_qualities = "".join(map(lambda x: chr(int(x) + 33), pileup_column.get_query_qualities()))
            if len(base_qualities) == 0:
                base_qualities = "*"

            mapq = "".join(map(lambda x: chr(int(x) + 33), pileup_column.get_mapping_qualities()))
            if len(mapq) == 0:
                mapq = "*"
            result = [name,
                      pos,
                      str_ref[pos - start_pos],
                      str(pileup_column.get_num_aligned()),
                      read_bases,
                      base_qualities,
                      mapq]
            yield result
        yield None

    def data_generator_samtools(bam, fasta, name, start, stop, is_bulk):

        def read_line_from_process(process_handle, ch, spliter, columns, data_out):
            # type: (Popen, str, str, int, list) -> int
            """
            Read data from the process's stdout, filter out the comments, split by spliter, and check the format
            :param process_handle: handle of the process
            :param ch: comment character
            :param spliter:
            :param columns: input. Greater than 0: column number of the data. Others, Unlimited
            :param data_out: output data (should clear buffer before using)
            :return: 0, Success. Others, processes have been terminated, and stdout can't read the data.
            """

            def read_line_from_file(from_file, ch, spliter, columns, data_out):
                # type: (file, str, str, int, list) -> int
                """
                Read data from file, filter out comments, and split by spliter, check format
                :param from_file: file pointer
                :param ch: comment character
                :param spliter:
                :param columns: input. Greater than 0: column number of the data. Others, Unlimited
                :param data_out: output data (should clear buffer before using)
                :return: 0, Success. Others, end of file.
                eg:
                data = []
                with open("test") as file:
                ret = read_data_from_file(file,"#","\t",8,data)
                """
                while 1:
                    buf = from_file.readline()
                    if len(buf) == 0:
                        return -1
                    if buf[0] == ch:
                        continue
                    buf = buf.strip("\n")
                    buf2 = buf.split(spliter)
                    if columns > 0:
                        if len(buf2) != columns:
                            continue
                        else:
                            break
                    else:
                        break
                buf2[1] = int(buf2[1])
                data_out.extend(buf2)
                return 0

            while 1:
                buf = []
                if read_line_from_file(process_handle.stdout, ch, spliter, columns, buf) != 0:
                    if not process_handle.poll() is None:  # samtools has been terminated
                        return -1
                    else:
                        continue
                else:
                    data_out.extend(buf)
                    return 0

        if is_bulk:
            cmd_str = "samtools mpileup -C50  -f {0} -s {1} -r {2}:{3}-{4}".format(fasta, bam, name, start, stop)
        else:
            cmd_str = "samtools mpileup -C50  -f {0} -q 40 -s {1} -r {2}:{3}-{4}".format(fasta, bam, name, start, stop)
        process_samtools = Popen([cmd_str], shell=True, stdout=PIPE)
        while 1:
            mpileup_list = []
            if read_line_from_process(process_samtools, "#", "\t", 7, mpileup_list) != 0:
                break
            else:
                yield mpileup_list
        yield None

    if engine == ENGINESAMTOOLS:
        return data_generator_samtools(bam, fasta, name, start, stop, is_bulk)
    else:
        return data_generator_pysam(bam, fasta, name, start, stop, is_bulk)


def calculate_coverage(name, pos, coverage_list):
    """
    format of coverage_list is [[name str, start int, stop int]]
    :return:
    0 Pos has joined the coverage_list
    1 Pos and front discontinuity
    :rtype: int
    :type pos: int
    :type name: str
    :param pos:
    :param name:
    :type coverage_list: list
    """
    if not coverage_list:  # is empty
        coverage_list.append([name, pos, pos])
        return 0
    # not empty
    if pos == coverage_list[-1][2] + 1:
        coverage_list[-1][2] = pos
        return 0
    else:
        coverage_list.append([name, pos, pos])
        return 1


def control(my_args, list_vcf, q5, name, start, stop, worker_id, v_rows):
    try:
        control_inner(my_args, list_vcf, q5, name, start, stop, worker_id, v_rows)
    except:
        logging.critical("Exception! name={0} start={1} stop={2} workerid={3}".format(name, start, stop, worker_id))


def control_inner(my_args, list_vcf, q5, name, start, stop, worker_id, v_rows):
    """
    :param stop: Unexpanded region
    :type my_args: object
    :type list_vcf: list
    :type name: str
    :type q5: object
    :type worker_id: int
    :type stop: int
    :type start: int
    """
    random.seed(10)
    q3_list = BigForewordList([])  # type: BigForewordList
    q2_list = BigForewordList([])  # type: BigForewordList
    q2_list_is_end = False
    coverage_list = []  # list of [name str, start int, stop int]
    counter = 0

    main_index = 0
    logging.info("\x1b[1;32mworker_id={0}\x1b[0m begin! name={1} "
                 "head={2} tail={3} len={4} engine={5}".format(worker_id, name, start,
                                                               stop, stop - start, my_args.engine))
    my_vcf_list = BigForewordList(list_vcf)
    copy_vcf_list = copy.copy(list_vcf)
    # Expanding the edge
    if worker_id != 1:
        my_start = start - my_args.lamb
    else:
        my_start = start

    if worker_id != my_args.work_num:
        my_stop = stop + my_args.lamb
    else:
        my_stop = stop
    total_work = my_stop - my_start
    # data source
    pileup_source = data_generator(my_args.bam, my_args.fasta, name, my_start, my_stop, my_args.engine, False)
    # bulk source
    if my_args.bulk != "":
        bulk_pileup_source = data_generator(my_args.bulk, my_args.fasta, name, start, stop, my_args.engine, True)
        so_source = get_so_source(bulk_pileup_source, my_args.bulk_min_var, my_args.bulk_min_mapq,
                                  my_args.bulk_min_depth)
        so_source.next()
    row = (worker_id - 1) % (v_rows - 2) + 2
    col = (worker_id - 1) / (v_rows - 2) * 21 + 1
    my_print("\x1b[{2};{3}Hw{0:<3d}      0/{1:<6}   ".format(worker_id, total_work, row, col),
             "control")

    while 1:
        mpileup_list = pileup_source.next()

        if mpileup_list == -1:  # pysam crashed
            logging.info("Use samtools engine instead!")
            q5.put(DataInQ5("pysam crashed", worker_id, WORKVAR, []), block=False)
            my_args.engine = ENGINESAMTOOLS
            control(my_args, copy_vcf_list, q5, name, start, stop, worker_id, v_rows)
            return
        if mpileup_list is None:
            if my_args.coverage and len(coverage_list) > 0:
                q5.put(DataInQ5(DONEMSG, worker_id, WORKCOVERAGE, coverage_list[0]), block=False)
                del coverage_list[0]
            break

        mpileup_list[4] = mpileup_list[4].upper()
        # append the data into q3_list
        if start <= mpileup_list[1] < stop:  # [start, stop)
            if mpileup_list[1] == 16432988:  # todo
                logging.debug("mpileup_list = {}".format(mpileup_list))

            if my_args.format == VARCALLFORMAT:
                data_in_q3_curr = read_mpileup(copy.copy(mpileup_list), my_args.minvar,
                                               my_args.mapq, my_args.min_depth, False,
                                               worker_id)  # type: PileupDataStruct

            else:
                data_in_q3_curr = read_mpileup_vcf(copy.copy(mpileup_list), my_args.minvar,
                                                   my_args.mapq, my_args.min_depth)  # type: PileupDataStruct

            # if mpileup_list[1] == 16114636:  # todo
            #     logging.debug("genotype_class={}".format(data_in_q3_curr.genotype_class))
            #     logging.debug("len read info = {}".format(len(data_in_q3_curr.read_info_list)))
            if data_in_q3_curr and data_in_q3_curr != -1:

                if my_args.bulk != "":
                    #  calculate so

                    result = so_source.send(data_in_q3_curr.coordinate_1)
                    if result == -1:
                        logging.info("Use samtools engine instead!")
                        q5.put(DataInQ5("pysam crashed", worker_id, WORKVAR, []), block=False)
                        my_args.engine = ENGINESAMTOOLS
                        control(my_args, copy_vcf_list, q5, name, start, stop, worker_id, v_rows)
                        return
                    else:
                        data_in_q3_curr.so = result
                else:
                    data_in_q3_curr.so = ""

                q3_list.append(data_in_q3_curr)

            # calculate coverage
            if my_args.bulk == "":
                if my_args.coverage and data_in_q3_curr:
                    if calculate_coverage(name, mpileup_list[1], coverage_list) == 1:
                        q5.put(DataInQ5("", worker_id, WORKCOVERAGE, coverage_list[0]), block=False)
                        del coverage_list[0]
            else:
                if my_args.coverage and \
                        data_in_q3_curr and \
                        ((data_in_q3_curr != -1 and data_in_q3_curr.so in ["", "varreads", "refgenotype"]) or
                         (data_in_q3_curr == -1)):
                    if calculate_coverage(name, mpileup_list[1], coverage_list) == 1:
                        q5.put(DataInQ5("", worker_id, WORKCOVERAGE, coverage_list[0]), block=False)
                        del coverage_list[0]
        # handle the head data in q3_list
        if not q3_list.is_empty():
            ret = differential(q2_list, q3_list.get_current_element(), q5, my_args.lamb, my_args.bias, my_args.null,
                               worker_id, q2_list_is_end, mpileup_list[1], my_args.format)
            if ret == 0:
                q3_list.move_foreword()
        counter += 1
        if counter == 10000:
            counter = 0
            main_index += 1
            my_print("\x1b[{3};{4}Hw{0:<3d} "
                     "{1:>6}/{2:<6}   ".format(worker_id, main_index * 10000, total_work,
                                               row, col),
                     "control")
        # vcf screen
        if not q2_list_is_end:
            while not my_vcf_list.is_empty():
                if my_vcf_list.get_current_element() < mpileup_list[1]:
                    my_vcf_list.move_foreword()
                else:
                    break
            if mpileup_list[1] == my_vcf_list.get_current_element():
                ret = golden_hetero(mpileup_list, 20, my_args.min_depth, my_args.snp_type,
                                    my_args.RD, worker_id)  # type: DataInQ2

                if ret is not None:
                    q2_list.append(ret)
            if my_vcf_list.is_empty() or mpileup_list[1] >= my_vcf_list.get_last_element():
                q2_list_is_end = True

        else:
            continue

    while not q3_list.is_empty():
        ret = differential(q2_list, q3_list.get_current_element(), q5, my_args.lamb, my_args.bias, my_args.null,
                           worker_id, True, my_stop, my_args.format)
        if ret == 0:
            q3_list.move_foreword()
    q5.put(DataInQ5(DONEMSG, worker_id,
                    WORKVAR if my_args.format == VARCALLFORMAT else WORKVCF, []), block=False)

    my_print("\x1b[{1};{2}H\x1b[1;32mw{0:<3d} I'm done!\x1b[0m       ".format(worker_id, row, col),
             "control")
    logging.info("\x1b[1;35mworker_id={0}\x1b[0m Quit".format(worker_id))
    return


def load_vcf(name, vcf_info, vcf_file_name):
    """
    Read the position information of the specified name chromosome from the vcf file, and store it as a list
    :param name: chromosome name
    :param vcf_info:
    :param vcf_file_name: vcf file name
    :return: POS information of vcf (list[int])
    """
    vcf_list = []
    curr_info = filter(lambda x: x[0] == name, vcf_info)
    if not curr_info:  # is empty
        return []
    else:
        with open(vcf_file_name) as fp:
            for i in curr_info:
                fp.seek(i[1], 0)
                buf = fp.read(i[2])
                buf_list = buf.splitlines(False)
                del buf
                vcf_list.extend(map(lambda y: int(y[1]), map(lambda x: x.split("\t"), buf_list)))
                del buf_list
        return vcf_list


def parse_vcf(vcf_file_name, need_check):
    """
    Parse the vcf file into [[name,head,length],[name,head,length]...] format
    Here head, length is the file pointer
    :param vcf_file_name:
    :return: the vcf infomation [[name,head,length],[name,head,length]...]
    """
    has_shown_info = False
    info_name = "{}.catalog".format(os.path.splitext(vcf_file_name)[0])
    if os.path.exists(info_name) and os.path.getmtime(info_name) > os.path.getmtime(vcf_file_name):
        # read vcf catalog
        with open(info_name, "r") as fp:
            return map(lambda x: map(lambda y: y if x.split(",").index(y) == 0 else int(y), x.split(",")),
                       fp.read().split(";"))
    else:
        # parse vcf
        result = []
        name = ""
        counter = 0
        head = 0
        with open(vcf_file_name, "r") as fp:
            while 1:
                text = fp.readline()
                if text == "":
                    result.append([name, head, fp.tell() - head])
                    break
                if text[0] == "#" or text == "\n":
                    counter += len(text)
                    continue
                tmp_list = text.split("\t")
                if name != tmp_list[0]:
                    if name == "":
                        head = fp.tell() - len(text)
                    else:
                        tail = fp.tell() - len(text)
                        result.append([name, head, tail - head])
                        head = tail
                    name = tmp_list[0]
                if need_check and not has_shown_info and len(tmp_list) > 8:
                    tmp_str = "\t".join(tmp_list[9:])
                    tmp_list = re.findall("0\\|0|0/0|1\\|1|1/1|2/2|2\\|2", tmp_str)
                    if len(tmp_list) > 0:
                        logging.info(">>>Please confirm the input VCF only contains heterozygote loci in bulk!!<<<")
                        has_shown_info = True

        # with open(info_name, "w") as fp:
        #     fp.write(";".join(map(lambda z: ",".join(z), map(lambda x: map(lambda y: str(y), x), result))))
        return result


def main():
    def parse_fasta(fasta_file_name):
        """
        Parse the fasta file into the format [[name,head,tail],[name,head,tail]...]
        :param fasta_file_name:
        :return:
        """
        # check catalog
        info_name = "{}.catalog".format(os.path.splitext(fasta_file_name)[0])
        if os.path.exists(info_name) and os.path.getmtime(info_name) > os.path.getmtime(fasta_file_name):
            # if False:
            with open(info_name, "r") as fp:
                return map(lambda x: map(lambda y: y if x.split(",").index(y) == 0 else int(y), x.split(",")),
                           fp.read().split(";"))
        else:  # parse fasta
            with open(fasta_file_name) as fp:
                fasta = fp.read()
            #  split into different chromosomes
            list_fasta = fasta.split(">")
            del fasta
            my_result = []

            for target in list_fasta:
                if target == "":
                    continue
                first_line = target[:target.find("\n")]

                if not (" " in first_line):
                    name = first_line
                else:
                    first_line_len = len(first_line)
                    for j in xrange(first_line_len):
                        if first_line[j] == " ":
                            break
                    name = first_line[:j]

                target = target[target.find("\n") + 1:]
                target = target.replace("\n", "")
                target_len = len(target)
                for j in xrange(target_len):
                    if target[j] != "N":
                        break
                head = j + 1

                for j in xrange(target_len):
                    if target[-1 - j] != "N":
                        break
                tail = target_len - j
                my_result.append([name, head, tail])
            # with open(info_name, "w") as fp:
            #     fp.write(";".join(map(lambda z: ",".join(z), map(lambda x: map(lambda y: str(y), x), my_result))))
            return my_result

    def write_vcf_head(output_file, my_format, fasta_file, bulk, minvar):
        """
        write the vcf file head
        :param output_file:
        :param my_format:
        :param fasta_file:
        :return:
        """
        if "./" in fasta_file:
            fasta_file = "{0}/{1}".format(os.getcwd(), fasta_file[2:])
        elif "/" not in fasta_file:
            fasta_file = "{0}/{1}".format(os.getcwd(), fasta_file)

        head_str = "##fileformat=VCFv4.3\n" \
                   "{0}\n" \
                   "##source=SCcallerV{3}\n" \
                   "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n" \
                   "##FILTER=<ID={2},Description=\"Multiple genotype\">\n" \
                   "##FILTER=<ID={4}{5},Description=\"Number of variant supporting reads < {5}\">\n" \
                   "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" \
                   "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"" \
                   "Allelic depths for the ref and alt alleles in the order listed\">\n" \
                   "##FORMAT=<ID=BI,Number=1,Type=Float," \
                   "Description=\"Amplification Bias\">\n" \
                   "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" \
                   "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"" \
                   "sequencing noise, amplification artifact, heterozygous SNV and " \
                   "homozygous SNV respectively\">\n".format(time.strftime("##fileDate=%Y%m%d", time.localtime()),
                                                             fasta_file,
                                                             MULTIPLEGENOTYPE,
                                                             VERSION,
                                                             NOTENOUGHVARIANTS,
                                                             minvar)
        if bulk != "":
            head_str = head_str + "##FORMAT=<ID=SO,Number=1,Type=String,Description=\"Whether it is a somatic mutation.\">\n"
        if my_format == VCFFORMAT:
            with open(output_file, "a") as fp_out:
                fp_out.write("{0}".format(head_str))
        return

    if float(sys.version[:3]) != 2.7:
        print "CRITICAL: Python version must be 2.7!\n"
        sys.exit(1)

    # handle parameter
    my_args = handle_paras()
    log_format = "%(asctime)s %(levelname)-8s %(process)d [%(funcName)s]%(message)s(%(filename)s line%(lineno)d)"
    logging.basicConfig(
        filename=get_my_filename(my_args.output, "_{0:0>2d}to{1:0>2d}.log".format(my_args.head, my_args.tail),
                                 "sc_"), level=logging.DEBUG, format=log_format,
        filemode="w")
    my_print("parsing vcf...", "main")
    logging.info("Welcome to use SCcaller 2.0.0")
    logging.info("now parsing vcf...")
    vcf_info = parse_vcf(my_args.snp_in, my_args.snp_type == "hsnp")

    my_print("handling fasta...", "main")
    logging.info("now parsing fasta...")
    fasta_info = parse_fasta(my_args.fasta)

    q5 = multiprocessing.Manager().Queue()

    if my_args.tail == -1 or my_args.tail > len(fasta_info):
        my_args.tail = len(fasta_info)
    # if my_args.cpu_num == -1:
    #    my_args.cpu_num = multiprocessing.cpu_count()
    #    logging.info("num of CPUs in this computer = {}".format(my_args.cpu_num))
    logging.info("args:{}".format(my_args))

    write_vcf_head(my_args.output, my_args.format, my_args.fasta, my_args.bulk, my_args.minvar)

    my_print("\x1b[?25l", "main")
    # Start the operation process
    for j in xrange(len(fasta_info)):
        if j + 1 < my_args.head:
            continue
        if j + 1 > my_args.tail:
            break
        process_pool = multiprocessing.Pool(processes=my_args.cpu_num)
        # v_rows, v_columns = map(lambda x: int(x), os.popen('stty size', 'r').read().split())
        v_rows = 26

        # fasta_info[j][2] = 16647426 + 20000  # tail # todo
        # fasta_info[j][1] = 16647426 - 20000  # head
        logging.info("loading vcf...")
        list_vcf = load_vcf(fasta_info[j][0], vcf_info, my_args.snp_in)
        my_print("loading vcf...", "main")
        logging.debug("name = {0} list vcf len = {1}".format(fasta_info[j][0], len(list_vcf)))
        my_print("\x1b[2J\x1b[0;0HSCcaller v2.0.0 is handling chromosome {}...\x1b[0J".format(fasta_info[j][0]), "main")

        # Calculate the amount of tasks for each process
        step = int(math.ceil((fasta_info[j][2] - fasta_info[j][1]) / float(my_args.work_num)))
        if step < my_args.lamb:
            logging.info(
                "work_num = {0} is too large for {1}. Use 1 instead. lamb = {2}".format(my_args.work_num,
                                                                                        fasta_info[j][0],
                                                                                        my_args.lamb))
            my_args.work_num = 1
            step = int(fasta_info[j][2] - fasta_info[j][1])

        # Start the write file process
        proc_write_result = multiprocessing.Process(target=write_result,
                                                    args=(
                                                        q5, my_args.output, my_args.work_num, my_args.bulk,
                                                        fasta_info[j][0],
                                                        my_args.head,
                                                        my_args.tail,
                                                        my_args.minvar))
        proc_write_result.daemon = True
        proc_write_result.start()
        for i in xrange(my_args.work_num):
            head = fasta_info[j][1] + i * step
            if i != 0:
                if (head - my_args.lamb) < 0:
                    print "lamb = {0} is too large. fasta info: name = {1} from {2} to {3}.".format(my_args.lamb,
                                                                                                    fasta_info[j][0],
                                                                                                    fasta_info[j][1],
                                                                                                    fasta_info[j][2])
                    logging.critical(
                        "lamb = {0} is too large. fasta info: name = {1} from {2} to {3}.".format(my_args.lamb,
                                                                                                  fasta_info[j][0],
                                                                                                  fasta_info[j][1],
                                                                                                  fasta_info[j][2]))
                    exit()
            tail = fasta_info[j][1] + i * step + step
            process_pool.apply_async(control, (my_args,
                                               filter(lambda x: head - my_args.lamb <= x <= tail + my_args.lamb,
                                                      list_vcf),
                                               q5,
                                               fasta_info[j][0],
                                               head, tail, i + 1,
                                               v_rows))

        del list_vcf
        process_pool.close()
        process_pool.join()

        # Exit the W process
        q5.put(DataInQ5(None, 0, WORKVAR if my_args.format == VARCALLFORMAT else WORKVCF, None),
               block=True)
        proc_write_result.join()
    extend_file(extend_file(my_args.output, my_args.output + ".eta",
                            "##reference=file://{0}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCELL001\n".format(my_args.fasta)),
                my_args.output + ".body", "")

    my_print("\x1b[{};0H\x1b[?25h".format(my_args.work_num + 2), "main")
    logging.info("W quit. All done.")
    return


def extend_file(file_a, file_b, str_in):
    with open(file_a, "a") as fp_out, open(file_b, "r") as fp_in:
        data = fp_in.read()
        fp_out.write(data)
        fp_out.write(str_in)
    os.remove(file_b)
    return file_a


class BigForewordList:
    """
    the list only consider the location after the current
    """
    my_list = None
    pos = 0

    def __init__(self, the_list):
        self.my_list = the_list
        self.pos = 0

    # def condition_jump(self, condition=lambda x: True):
    #     if self.my_list:
    #         ge = (x for x in xrange(len(self.my_list[self.pos:])) if condition(self.my_list[x]))
    #         try:
    #             tmp = ge.next()  # type: int
    #             if tmp < len(self.my_list):
    #                 self.pos = tmp
    #         except StopIteration:
    #             self.pos = len(self.my_list)-1
    def contain(self, element):
        """
        Determine if the element is in the List
        :param element:
        :return:
        """
        if self.my_list == [] or self.pos == len(self.my_list) - 1:
            return False
        if element not in self.my_list:
            return False
        return element in self.my_list[self.pos:]

    def filter(self, func):
        if self.my_list:
            return filter(func, self.my_list[self.pos:])
        else:
            return []

    def get_element(self, index):
        # type: (int) -> object
        if not self.my_list:
            return None
        if index >= 0:
            if self.pos + index < len(self.my_list):
                return self.my_list[self.pos + index]
            else:
                return None
        else:
            if len(self.my_list) + index >= self.pos:
                return self.my_list[index]
            else:
                return None

    def len(self):
        if not self.my_list:
            return 0
        else:
            return len(self.my_list) - self.pos

    def append(self, element):
        self.my_list.append(element)

    def display_list(self):
        print self.my_list

    def get_foreword_list(self):
        if self.my_list and self.pos < len(self.my_list):
            return self.my_list[self.pos:]
        else:
            return []

    def get_current_element(self):
        if self.my_list and self.pos < len(self.my_list):
            return self.my_list[self.pos]
        else:
            return None

    def get_last_element(self):
        if self.my_list:
            return self.my_list[-1]
        else:
            return None

    def move_foreword(self):
        if self.pos <= len(self.my_list) - 1:
            self.pos += 1
        if self.pos == 10000:
            del self.my_list[:10000]
            self.pos = 0

    def is_empty(self):
        if self.my_list:
            return self.pos == len(self.my_list)
        else:
            return True


if __name__ == "__main__":
    start = time.time()
    main()
    print "duration = {}s".format(str(time.time() - start))
