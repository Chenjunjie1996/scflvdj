#!/usr/bin/env python

import argparse

import pandas as pd
import pysam
import utils
from __init__ import ASSAY
from summarize import get_vdj_metric


def match_contig_csv(sample, seqtype, contig_csv, match_barcode):
    """
    Output filter contig annotation matched with sc-RNA barcodes.
    """
    if seqtype == "BCR":
        chains = ["IGH", "IGL", "IGK"]
        paired_groups = ["IGK_IGH", "IGL_IGH"]
    else:
        chains = ["TRA", "TRB"]
        paired_groups = ["TRA_TRB"]

    df_match = pd.read_csv(contig_csv)
    df_match = df_match[df_match["barcode"].isin(match_barcode)]
    df_match.to_csv(f"{sample}_matched_contig.csv", sep=",", index=False)
    matched_cell_num = len(set(df_match["barcode"]))
    data_dict = get_vdj_metric(df_match, chains, paired_groups)
    data_dict.update({"Cells Match with ScRNA-seq": matched_cell_num})
    fn = f"{args.sample}.{ASSAY}.match.stats.json"
    utils.write_json(data_dict, fn)


def match_contig_fasta(sample, contig_fasta, match_barcode):
    """
    Output filter contig fasta matched with sc-RNA barcodes.
    """
    fasta_file = pysam.FastxFile(contig_fasta)
    match_fasta_file = open(f"{sample}_matched_contig.fasta", "w")
    for entry in fasta_file:
        name = entry.name
        attrs = name.split("_")
        cb = "_".join(attrs[:3])
        if cb in match_barcode:
            new_name = cb + "_" + attrs[-2] + "_" + attrs[-1]
            seq = entry.sequence
            match_fasta_file.write(f">{new_name}\n{seq}\n")
    match_fasta_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--contig_csv", required=True)
    parser.add_argument("--contig_fasta", required=True)
    parser.add_argument("--seqtype", required=True)
    parser.add_argument("--match_barcode_file", help="File containing matched barcodes", required=True)
    args = parser.parse_args()

    match_barcode = set(utils.read_one_col(args.match_barcode_file))
    match_contig_csv(args.sample, args.seqtype, args.contig_csv, match_barcode)
    match_contig_fasta(args.sample, args.contig_fasta, match_barcode)
