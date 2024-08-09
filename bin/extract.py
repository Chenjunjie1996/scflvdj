#!/usr/bin/env python

import argparse

import parse_protocol
import pyfastx
import itertools
import random
import utils
from __init__ import ASSAY

SPLIT_N_CHUNKS = 4


class Auto(parse_protocol.Auto):
    def __init__(
        self,
        fq1_list,
        sample,
    ):
        super().__init__(fq1_list, sample)

    def seq_protocol(self, seq):
        """
        Returns: protocol or None

        >>> runner = Auto([], "fake_sample")
        """
        for protocol in self.protocol_dict:
            if self.is_protocol(seq, protocol):
                return protocol


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--fq1", required=True)
    parser.add_argument("--fq2", required=True)
    parser.add_argument("--assets_dir", required=True)
    parser.add_argument("--protocol", required=True)
    args = parser.parse_args()

    fq1_list = args.fq1.split(",")
    fq2_list = args.fq2.split(",")
    # protocol
    protocol_dict = parse_protocol.get_protocol_dict(args.assets_dir)
    p = args.protocol

    fn = f"{args.sample}.{ASSAY}.protocol.stats.json"
    utils.write_json({"Protocol": p}, fn)

    pmeta = protocol_dict[p]
    pattern_dict = pmeta["pattern_dict"]
    raw_list, mismatch_list = parse_protocol.get_raw_mismatch(pmeta["bc"], 1)

    # out_fq
    out_fq_fn = {x: f"{args.sample}_R{x}.fq.gz" for x in [1, 2]}
    outdict = {k: utils.openfile(v, "wt") for k, v in out_fq_fn.items()}

    raw_reads = 0
    valid_reads = 0
    corrected_reads = 0
    total_umi_set = list(itertools.product('ATCGN', repeat=9))
    total_umi_set = set("".join(x) for x in total_umi_set)
    
    # 遍历第一个文库，查找umi
    fq1 = pyfastx.Fastx(fq1_list[0])
    fq2 = pyfastx.Fastx(fq2_list[0])
    used_umi_set = set()
    for (name1, seq1, qual1), (name2, seq2, qual2) in zip(fq1, fq2):
        raw_reads += 1
        bc_list = [utils.rev_compl(seq1[x]) for x in pattern_dict["C"]][::-1]
        valid, corrected, corrected_seq = parse_protocol.check_seq_mismatch(bc_list, raw_list, mismatch_list)
        if valid:
            valid_reads += 1
            if corrected:
                corrected_reads += 1
            umi = parse_protocol.get_seq_str(seq1, pattern_dict["U"])
            bc = corrected_seq
            read_name = f"{bc}:{umi}:{raw_reads}"
            qual1 = "F" * len(bc + umi)
            outdict[1].write(utils.fastq_str(read_name, bc + umi, qual1))
            outdict[2].write(utils.fastq_str(read_name, seq2, qual2))
            used_umi_set.add(umi)
    
    # 用未出现在第一个文库的umi替换第二个文库的umi，得到替换umi的字典
    diff_umi_list = list(total_umi_set - used_umi_set)
    random.shuffle(diff_umi_list)
    umi_dict = {}
    fq1 = pyfastx.Fastx(fq1_list[1])
    for (name1, seq1, qual1) in fq1:
        umi = parse_protocol.get_seq_str(seq1, pattern_dict["U"])
        if umi not in umi_dict:
            umi_dict[umi] = diff_umi_list.pop()
    
    # 遍历第二个文库
    fq1 = pyfastx.Fastx(fq1_list[1])
    fq2 = pyfastx.Fastx(fq2_list[1])

    for (name1, seq1, qual1), (name2, seq2, qual2) in zip(fq1, fq2):
        raw_reads += 1
        bc_list = [utils.rev_compl(seq1[x]) for x in pattern_dict["C"]][::-1]
        valid, corrected, corrected_seq = parse_protocol.check_seq_mismatch(bc_list, raw_list, mismatch_list)
        if valid:
            valid_reads += 1
            if corrected:
                corrected_reads += 1
            umi = parse_protocol.get_seq_str(seq1, pattern_dict["U"])
            new_umi = umi_dict[umi]
            bc = corrected_seq
            read_name = f"{bc}:{new_umi}:{raw_reads}"
            qual1 = "F" * len(bc + umi)
            outdict[1].write(utils.fastq_str(read_name, bc + new_umi, qual1))
            outdict[2].write(utils.fastq_str(read_name, seq2, qual2))

    outdict[1].close()
    outdict[2].close()

    fn = f"{args.sample}.{ASSAY}.extract.stats.json"
    metrics = {"Raw Reads": raw_reads}
    metrics["Valid Reads"] = utils.get_frac(valid_reads / raw_reads)
    metrics["Corrected Barcodes"] = utils.get_frac(corrected_reads / valid_reads)
    utils.write_json(metrics, fn)
