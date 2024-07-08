#!/usr/bin/env python

import os
import pyfastx
import operator
import random
import argparse
from collections import defaultdict

import utils


# Template Swithcing Oligos Sequence. [16-bp cell barcode][10/12-bp UMI][TSO].
TSO = "TTTCTTATATGGG"
UMI_10X_LEN = 10


def gen_sgr_tenX_dict(fq2, whitelist_10X):
    sgr_tenX = {}
    count_dict = defaultdict(int)
    whitelist_10X_fh = open(whitelist_10X, 'r')
    
    fq2 = pyfastx.Fastx(fq2)
    for (name, seq, qual) in fq2:
        sgr_barcode = name.split(':')[0]
        count_dict[sgr_barcode] += 1

    count_dict = dict(sorted(count_dict.items(), key=operator.itemgetter(1), reverse=True))
        
    for sgr_barcode in count_dict:
        sgr_tenX[sgr_barcode] = whitelist_10X_fh.readline().strip()

    # Add invalid barcode
    for sgr_barcode, barcode_10X in sgr_tenX.items():
        if barcode_10X == '':
            sgr_tenX[sgr_barcode] = "AAAA" + ''.join(random.choice("ATCG") for _ in range(12))

    return sgr_tenX


def convert_seq(sgr_tenX, barcode_sgr, umi_sgr):

    barcode_10X = sgr_tenX[barcode_sgr]

    umi_len_sgr = len(umi_sgr)
    if umi_len_sgr > UMI_10X_LEN:
        umi_10X = umi_sgr[:UMI_10X_LEN]
    elif umi_len_sgr < UMI_10X_LEN:
        umi_10X = umi_sgr + "C" * (UMI_10X_LEN - umi_len_sgr)
    else:
        umi_10X = umi_sgr

    new_seq1 = barcode_10X + umi_10X + TSO
    new_qual1 = 'F' * len(new_seq1)

    return new_seq1, new_qual1
    
    
def write_fq1(fq2, sample):
    out_fq1 = utils.openfile(f"convert_fq/{sample}_S1_L001_R1_001.fastq.gz")

    fq2 = pyfastx.Fastx(fq2)
    for (name, seq, qual) in fq2:
        attrs = name.split(':')
        sgr_barcode, sgr_umi = attrs[0], attrs[1]
        new_seq1, new_qual1 = convert_seq(sgr_barcode, sgr_umi)
        out_fq1.write(f'@{name}\n{new_seq1}\n+\n{new_qual1}\n')

    out_fq1.close() 


def gzip_fq2(fq2, sample):
    out_fq2_file = f"convert_fq/{sample}_S1_L001_R2_001.fastq.gz"
    cmd = f"gzip -c {fq2} > {out_fq2_file}"
    os.system(cmd)
    

def dump_tenX_sgr_barcode_json(sgr_tenX):
    tenX_sgr = {}
    barcode_convert_json = f"barcode_convert.json"
    for sgr, tenX in sgr_tenX.items():
        tenX_sgr[tenX] = sgr

    utils.dump_dict_to_json(tenX_sgr, barcode_convert_json) 


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--fq2", required=True)
    parser.add_argument("--assets_dir", required=True)
    args = parser.parse_args()
    
    whitelist_10X = f"{args.assets_dir}/whitlist/737K-august-2016.txt"
    sgr_tenX = gen_sgr_tenX_dict(args.fq2, whitelist_10X)
    write_fq1(args.fq2, args.sample)
    gzip_fq2(args.fq2, args.sample)
    dump_tenX_sgr_barcode_json(sgr_tenX)