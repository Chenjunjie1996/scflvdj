#!/usr/bin/env python

import argparse
import random
from collections import defaultdict

import parse_protocol
import pyfastx
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
    barcode_read_counter = defaultdict(int)
    for fq1, fq2 in zip(fq1_list, fq2_list):
        fq1 = pyfastx.Fastx(fq1)
        fq2 = pyfastx.Fastx(fq2)

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
                barcode_read_counter[bc] += 1
                if barcode_read_counter[bc] <= 40000:
                    outdict[1].write(utils.fastq_str(read_name, bc + umi, qual1))
                    outdict[2].write(utils.fastq_str(read_name, seq2, qual2))

    outdict[1].close()
    outdict[2].close()

    fn = f"{args.sample}.{ASSAY}.extract.stats.json"
    metrics = {"Raw Reads": raw_reads}
    metrics["Valid Reads"] = utils.get_frac(valid_reads / raw_reads)
    metrics["Corrected Barcodes"] = utils.get_frac(corrected_reads / valid_reads)
    utils.write_json(metrics, fn)

    # split fastq
    barcode_list = list(barcode_read_counter.keys())
    random.shuffle(barcode_list)
    barcode_num = len(barcode_list)
    step = barcode_num // SPLIT_N_CHUNKS
    split_barcodes = [barcode_list[i:i+step] for i in range(0, barcode_num, step)]
    split_barcodes = [set(i) for i in split_barcodes]

    fq1_list = [utils.openfile(f"temp/temp{i}_R1.fq.gz", "wt") for i in range(SPLIT_N_CHUNKS)]
    fq2_list = [utils.openfile(f"temp/temp{i}_R2.fq.gz", "wt") for i in range(SPLIT_N_CHUNKS)]
    fq1 = pyfastx.Fastx(f"{args.sample}_R1.fq.gz")
    fq2 = pyfastx.Fastx(f"{args.sample}_R2.fq.gz")
    for (name1, seq1, qual1), (name2, seq2, qual2) in zip(fq1, fq2):
        bc = name1.split(':')[0]
        for i in range(SPLIT_N_CHUNKS):
            if bc in split_barcodes[i]:
                fq1_list[i].write(utils.fastq_str(name1, seq1, qual1))
                fq2_list[i].write(utils.fastq_str(name2, seq2, qual2))
                break

    for i in range(SPLIT_N_CHUNKS):
        fq1_list[i].close()
        fq2_list[i].close()