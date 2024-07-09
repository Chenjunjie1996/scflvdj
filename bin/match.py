#!/usr/bin/env python

import pandas as pd
import pysam
import argparse

import utils
from __init__ import ASSAY

def gen_vdj_metric(seqtype):
    """
    Add vdj metrics.
    """
    def get_vj_spanning_pair():
        """
        Get Productive V-J Spanning_Pair metric from annotation file
        Return productive chain pair number. eg: TRA/TRB or IGH/IGL, IGH/IGK.
        """
        df_productive = df[df["productive"]]
        
        if seqtype == "BCR":
            df_chain_heavy = df_productive[(df_productive["chain"] == "IGH")] 
            df_chain_light = df_productive[(df_productive["chain"] == "IGL") | (df_productive["chain"] =="IGK")]
        else:
            df_chain_heavy = df_productive[df_productive["chain"] == "TRA"]
            df_chain_light = df_productive[df_productive["chain"] == "TRB"]

        for _df in [df_chain_heavy, df_chain_light]:
            _df.drop_duplicates(["barcode"], inplace=True)

        vj_spanning_pair_cells = pd.merge(df_chain_heavy, df_chain_light, on="barcode", how="inner")
        
        return vj_spanning_pair_cells.shape[0]

    df = pd.read_csv("matched_contig_annotations.csv")
    data_dict = {}
    
    if seqtype == "BCR":
        chains = ["IGH", "IGL", "IGK"]
        chain_pairs = ["IGK_IGH", "IGL_IGH"]
    else:
        chains = ["TRA", "TRB"]
        chain_pairs = ["TRA_TRB"]
    
    cell_nums = len(set(df["barcode"]))
    data_dict.update({"Cells Match with ScRNA-seq": cell_nums})
    data_dict = {"Cells With Productive V-J Spanning Pair": utils.get_frac(get_vj_spanning_pair() / cell_nums)}

    for pair in chain_pairs:
        chain1, chain2 = pair.split("_")[0], pair.split("_")[1]
        cbs1 = set(df[(df["productive"]) & (df["chain"]==chain1)].barcode)
        cbs2 = set(df[(df["productive"]) & (df["chain"]==chain2)].barcode)
        paired_cbs = len(cbs1.intersection(cbs2))
        data_dict.update(
            {f"Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair": utils.get_frac(paired_cbs / cell_nums)}
        )

    for chain in chains:
        value = len(set(df[df["chain"] == chain].barcode))
        data_dict.update({f"Cells With {chain} Contig": utils.get_frac(value / cell_nums)})

        value = len(set(df[(df["chain"] == chain) & (df["cdr3"] != "None")].barcode))
        data_dict.update({f"Cells With CDR3-annotated {chain} Contig": utils.get_frac(value / cell_nums)})

        value = len(set(df[(df["full_length"]) & (df["chain"] == chain)].barcode))
        data_dict.update({f"Cells With V-J Spanning {chain} Contig": utils.get_frac(value / cell_nums)})

        value = len(set(df[(df["full_length"]) & (df["productive"]) & (df["chain"] == chain)].barcode))
        data_dict.update({f"Cells With Productive {chain} Contig": utils.get_frac(value / cell_nums)})

    return data_dict


def gen_matched_result(annot_csv, contig_fasta, match_cell_barcodes):

    df_annotation = pd.read_csv(annot_csv)
    df_match = df_annotation[df_annotation.barcode.isin(match_cell_barcodes)]
    df_match.to_csv("matched_contig_annotations.csv", sep=",", index=False)

    fasta_fh = pysam.FastxFile(contig_fasta)
    match_fasta = open("filtered_contig.fasta", "w")
    for entry in fasta_fh:
        name = entry.name
        attrs = name.split("_")
        cb = "_".join(attrs[:3])
        if cb in match_cell_barcodes:
            new_name = cb + "_" + attrs[-2] + "_" + attrs[-1]
            seq = entry.sequence
            match_fasta.write(f">{new_name}\n{seq}\n")
    match_fasta.close()


def gen_matched_clonotypes(clonotype_csv):

    raw_clonotypes= pd.read_csv(clonotype_csv, sep=",", index_col=None)
    raw_clonotypes.drop(["frequency", "proportion"], axis=1, inplace=True)
    df_match = pd.read_csv("matched_contig_annotations.csv")
    df_match = df_match[df_match["productive"]]
        
    # Count frequency and proportion
    df_match = df_match.rename(columns={"raw_clonotype_id":"clonotype_id"})\
        .dropna(subset=["clonotype_id"]).groupby("clonotype_id")["barcode"].nunique().to_frame()\
            .reset_index().rename(columns={"barcode": "frequency"})\
                .sort_values("clonotype_id", key=lambda x: x.str.lstrip("clonotype").astype(int))
    df_match["proportion"] = df_match["frequency"] / df_match["frequency"].sum()
        
    df_match = pd.merge(df_match, raw_clonotypes, on="clonotype_id")
    df_match.to_csv("matched_clonotypes.csv", sep=",", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seqtype", required=True)
    parser.add_argument("--clonotype_csv", required=True)
    parser.add_argument("--annot_csv", required=True)
    parser.add_argument("--contig_fasta", required=True)
    args = parser.parse_args()

    if args.seqtype == "BCR":
        chains = ["IGH", "IGL", "IGK"]
        paired_groups = ["IGK_IGH", "IGL_IGH"]
    else:
        chains = ["TRA", "TRB"]
        paired_groups = ["TRA_TRB"]
    
    match_barcode = set(utils.read_one_col(args.match_barcode_file))
    gen_matched_result(args.annot_csv, args.contig_fasta, match_barcode)
    gen_matched_clonotypes(args.clonotype_csv)
    data_dict = gen_vdj_metric(args.seqtype)
    fn = f"{args.sample}.{ASSAY}.match.stats.json"
    utils.write_json(data_dict, fn)
    