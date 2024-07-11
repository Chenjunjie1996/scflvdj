#!/usr/bin/env python

import json
import pandas as pd
import pysam
import os
import argparse
from collections import defaultdict

import utils
from __init__ import ASSAY


MAX_CELL = 2 * 10**5


def get_umi_count(bcs, umis, cell_barcodes, sample):
    """
    Args:
        bcs: raw barcodes
        umis: umi count
        cell_barcodes: cell barcodes
    """
    a = [(umi, bc) for umi, bc in zip(umis, bcs) if umi > 0]
    a.sort(reverse=True)
    cell_barcodes = set(cell_barcodes)
    plot_data = {}
    n = len(a)
    first_noncell = n - 1
    for i, (umi, bc) in enumerate(a):
        if bc not in cell_barcodes:
            first_noncell = i
            break
    print(f"first non-cell barcode rank: {first_noncell}")
    last_cell = 0
    for i in range(min(n - 1, MAX_CELL), -1, -1):
        bc = a[i][1]
        if bc in cell_barcodes:
            last_cell = i
            break
    pure = sample + ".cells.pure" + f"({first_noncell}/{first_noncell}, 100%)"
    bg_cells = n - first_noncell
    bg = sample + ".cells.background" + f"(0/{bg_cells}, 0%)"
    plot_data[pure] = {}
    plot_data[bg] = {}
    for i in range(first_noncell + 1):  # for plot
        plot_data[pure][i + 1] = int(a[i][0])

    n_mix = last_cell - first_noncell + 1
    if n_mix != 0:
        n_total = len(cell_barcodes)
        n_mix_cell = n_total - first_noncell
        mix_rate = round(n_mix_cell / n_mix * 100, 2)
        mix = sample + ".cells.mix" + f"({n_mix_cell}/{n_mix}, {mix_rate}%)"
        plot_data[mix] = {}
        for i in range(first_noncell, last_cell + 1):
            plot_data[mix][i + 1] = int(a[i][0])

    for i in range(last_cell + 1, min(MAX_CELL, n), 10):
        plot_data[bg][i + 1] = int(a[i][0])
    # do not record every umi count
    for i in range(MAX_CELL, n, 1000):
        plot_data[bg][i + 1] = int(a[i][0])
    return plot_data


def barcode_rank_plot(sample, df_annotation, bam, tenx_sgr):
    dic_umi = defaultdict(set)
        
    with pysam.AlignmentFile(bam) as fh:
        for read in fh:
            cb = read.get_tag('CB')
            umi = read.get_tag('UB')
        dic_umi[cb].add(umi)

    df_umi = pd.DataFrame()
    df_umi['barcode'] = list(dic_umi.keys())
    df_umi['UMI'] = [len(dic_umi[i]) for i in dic_umi]
    df_umi = df_umi.sort_values(by='UMI', ascending=False)
    cbs = set(df_annotation['barcode'])
    df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if x in cbs else 'UB')
    df_umi['barcode'] = df_umi['barcode'].apply(lambda x : tenx_sgr[x.split('-')[0]])
    df_umi.to_csv(f"{sample}.count.txt", sep='\t', index=False)
    
    plot_data = get_umi_count(df_umi["barcode"], df_umi["UMI"], cbs, sample)

    return plot_data
 
 
def gen_vdj_metric(metrics_csv, df_annotation, seqtype):
    
    df = pd.read_csv(metrics_csv, index_col=None)
    metrics_dict = df.T.to_dict()[0]

    if seqtype == "BCR":
        chains = ["IGH", "IGL", "IGK"]
        chain_pairs = ["IGK_IGH", "IGL_IGH"]
    else:
        chains = ["TRA", "TRB"]
        chain_pairs = ["TRA_TRB"]

    data_dict = {}
    for k, v in metrics_dict.items():
        metrics_dict[k] = str(v).replace(',', '').replace('%', '')
    # mapping
    mapping_metrics_list = ["Reads Mapped to Any V(D)J Gene"]
    for chain in chains:
        name = f"Reads Mapped to {chain}"
        mapping_metrics_list.append(name)
        
    for name in mapping_metrics_list:
        data_dict.update({name: float(metrics_dict[name])})

    # cells
    cell_metrics_list = [
        "Estimated Number of Cells",
        "Fraction Reads in Cells",
        "Mean Read Pairs per Cell",
        "Mean Used Read Pairs per Cell"
    ]
    for name in cell_metrics_list:
        data_dict.update({name: float(metrics_dict[name])})
        
    for chain in chains:
        name = f"Median used {chain} UMIs per Cell"
        median_value = df_annotation[df_annotation['chain']==chain]["umis"].median()
        if median_value == median_value:
            value = int(median_value)
        else:
            value = 0
        data_dict.update({name: value})

    fn = f"{args.sample}.{ASSAY}.summarize.stats.json"
    utils.write_json(data_dict, fn)

    # annotation
    data_dict = {}
    annotation_metrics_list = ["Cells With Productive V-J Spanning Pair"]
    for pair in chain_pairs:
        chain1, chain2 = pair.split("_")[0], pair.split("_")[1]
        name = f"Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair"
        annotation_metrics_list.append(name)
           
    for chain in chains:
        name_list = [
            f"Cells With {chain} Contig",
            f"Cells With CDR3-annotated {chain} Contig",
            f"Cells With V-J Spanning {chain} Contig",
            f"Cells With Productive {chain} Contig"
        ]
        annotation_metrics_list.extend(name_list)

    for name in annotation_metrics_list:
        data_dict.update({name: float(metrics_dict[name])})
    
    fn = f"{args.sample}.{ASSAY}.annotation.stats.json"
    utils.write_json(data_dict, fn)
    
    
def convert_barcode_to_sgr(sample, tenx_sgr, df_annotation, contig_fasta, clonotype_csv):
        
    df_annotation["barcode"] = df_annotation["barcode"].apply(lambda x: tenx_sgr[x.split("-")[0]])
    df_annotation["contig_id"] = df_annotation["contig_id"].apply(lambda x: 
        tenx_sgr[x.split("-")[0]]+ "_" + x.split("_")[1] +"_" + x.split("_")[2])

    df_annotation.to_csv(f"{sample}_filtered_contig.csv", sep=",", index=False)
    
    tenx_fasta_fh = pysam.FastxFile(contig_fasta)
    with open(f"{sample}_filtered_contig.fasta", "w") as f:
        for entry in tenx_fasta_fh:
            name = entry.name
            seq = entry.sequence
            attrs = name.split("_")
            new_name = tenx_sgr[attrs[0].split("-")[0]] + "_" + attrs[1] + "_" + attrs[2]
            f.write(f">{new_name}\n{seq}\n")
        
    cmd = f"cp {clonotype_csv} {sample}_clonotypes.csv"
    os.system(cmd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--seqtype", required=True)
    parser.add_argument("--barcode_convert_json", required=True)
    parser.add_argument("--clonotype_csv", required=True)
    parser.add_argument("--annot_csv", required=True)
    parser.add_argument("--contig_fasta", required=True)
    parser.add_argument("--metrics_csv", required=True)
    args = parser.parse_args()

    df_annotation = pd.read_csv(args.annot_csv, sep=",", index_col=None)
    with open(args.barcode_convert_json, "r") as f:
        tenx_sgr = json.load(f)
   
    convert_barcode_to_sgr(args.sample, tenx_sgr, df_annotation, args.contig_fasta, args.clonotype_csv)
    gen_vdj_metric(args.metrics_csv, df_annotation, args.seqtype)
    
    plot_data = barcode_rank_plot(args.sample, df_annotation, args.bam, tenx_sgr)
    fn = f"{args.sample}.{ASSAY}.umi_count.json"
    utils.write_json(plot_data, fn) 
