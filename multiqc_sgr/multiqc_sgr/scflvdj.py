import json
import logging
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph

ASSAY = "scflvdj"
# Initialise the logger
log = logging.getLogger("multiqc")


def get_int(x):
    return str(int(x))


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name=ASSAY,
            anchor=ASSAY,
        )
        log.info(f"Running module: {self.name}")

        stat_data = self.parse_json(self.name, "stats")
        umi_count_data = self.parse_json(self.name, "umi_count")
        if all(
            len(x) == 0
            for x in [
                stat_data,
                umi_count_data,
            ]
        ):
            raise ModuleNoSamplesFound

        # Basic Stats Table
        self.general_stats_table(stat_data)

        # barcode rank plot
        self.add_section(
            name="Barcode Rank", anchor=f"{ASSAY}_barcode_rank", plot=self.barcode_rank_plot(umi_count_data)
        )

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

    def parse_json(self, assay, seg):
        data_dict = defaultdict(dict)
        n = 0
        for f in self.find_log_files(f"{assay}/{seg}"):
            log.debug(f"Found file: {f['fn']}")
            n += 1
            parsed_data = json.loads(f["f"])
            if parsed_data is not None:
                x = f["s_name"]
                s_name = x[: x.find(f".{assay}")]
                if s_name in data_dict:
                    log.info(f"Duplicate sample name found! Update: {s_name}")
                self.add_data_source(f, s_name=s_name, section=seg)
                data_dict[s_name].update(parsed_data)

        data_dict = self.ignore_samples(data_dict)

        log.info(f"Found {n} {assay} {seg} reports")
        # Write parsed report data to a file
        self.write_data_file(data_dict, f"multiqc_{assay}_{seg}")
        return data_dict

    def general_stats_table(self, summary_data):
        headers = {
            "Protocol": {
                "title": "Protocol",
                "description": "Predefined pattern of barcode and UMI",
                "scale": "purple",
                "hidden": True,
            },
            "Raw Reads": {
                "title": "Raw Reads",
                "description": "Number of reads in the input file",
                "scale": "blue",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "Valid Reads": {
                "title": "Valid Reads",
                "description": "Percent of reads with valid barcodes",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Corrected Barcodes": {
                "title": "Corrected Barcodes",
                "description": "Percent of corrected barcodes",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Reads Mapped To Any V(D)J Genes": {
                "title": "Reads Mapped To Any V(D)J Genes",
                "description": "Percent of reads mapped to any V(D)J genes",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Estimated Number of Cells": {
                "title": "Number of Cells",
                "description": "Estimated number of cells",
                "scale": "blue",
                "format": "{:,.0f}",
            },
            "Fraction of Reads in Cells": {
                "title": "Reads in Cells",
                "description": "Percent of unique reads in cells",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Mean Read Pairs per Cell": {
                "title": "Mean Reads",
                "description": "Mean number of reads per cell",
                "scale": "blue",
                "format": "{:,.0f}",
            },
            "Mean Used Read Pairs per Cell": {
                "title": "Mean Used Read",
                "description": "Mean number of reads used per cell",
                "scale": "blue",
                "format": "{:,.0f}",
            },
            "Median TRA UMIs per Cell": {
                "title": "Median TRA UMIs per Cell",
                "description": "Median number of TRA umis per cell",
                "scale": "green",
                "format": "{:,.0f}",
            },
            "Median TRB UMIs per Cell": {
                "title": "Median TRB UMIs per Cell",
                "description": "Median number of TRB umis per cell",
                "scale": "green",
                "format": "{:,.0f}",
            },
            "Median IGH UMIs per Cell": {
                "title": "Median IGH UMIs per Cell",
                "description": "Median number of IGH umis per cell",
                "scale": "green",
                "format": "{:,.0f}",
            },
            "Median IGK UMIs per Cell": {
                "title": "Median IGK UMIs per Cell",
                "description": "Median number of IGK umis per cell",
                "scale": "green",
                "format": "{:,.0f}",
            },
            "Median IGL UMIs per Cell": {
                "title": "Median IGL UMIs per Cell",
                "description": "Median number of IGL umis per cell",
                "scale": "green",
                "format": "{:,.0f}",
            },
            "Cells With Productive V-J Spanning Pair": {
                "title": "Cells With Productive V-J Spanning Pair",
                "description": "Cells With Productive V-J Spanning Pair",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Cells With Productive V-J Spanning (TRA, TRB) Pair": {
                "title": "Cells With Productive V-J Spanning (TRA, TRB) Pair",
                "description": "Cells With Productive V-J Spanning (TRA, TRB) Pair",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With Productive V-J Spanning (IGK, IGH) Pair": {
                "title": "Cells With Productive V-J Spanning (IGK, IGH) Pair",
                "description": "Cells With Productive V-J Spanning (IGK, IGH) Pair",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With Productive V-J Spanning (IGL, IGH) Pair": {
                "title": "Cells With Productive V-J Spanning (IGL, IGH) Pair",
                "description": "Cells With Productive V-J Spanning (IGL, IGH) Pair",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With TRA Contig": {
                "title": "Cells With TRA Contig",
                "description": "Cells With TRA Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With CDR3-annotated TRA Contig": {
                "title": "Cells With CDR3-annotated TRA Contig",
                "description": "Cells With CDR3-annotated TRA Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With V-J Spanning TRA Contig": {
                "title": "Cells With V-J Spanning TRA Contig",
                "description": "Cells With V-J Spanning TRA Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With Productive TRA Contig": {
                "title": "Cells With Productive TRA Contig",
                "description": "Cells With Productive TRA Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Cells With TRB Contig": {
                "title": "Cells With TRB Contig",
                "description": "Cells With TRB Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With CDR3-annotated TRB Contig": {
                "title": "Cells With CDR3-annotated TRB Contig",
                "description": "Cells With CDR3-annotated TRB Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With V-J Spanning TRB Contig": {
                "title": "Cells With V-J Spanning TRB Contig",
                "description": "Cells With V-J Spanning TRB Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With Productive TRB Contig": {
                "title": "Cells With Productive TRB Contig",
                "description": "Cells With Productive TRB Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Cells With IGH Contig": {
                "title": "Cells With IGH Contig",
                "description": "Cells With IGH Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With CDR3-annotated IGH Contig": {
                "title": "Cells With CDR3-annotated IGH Contig",
                "description": "Cells With CDR3-annotated IGH Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With V-J Spanning IGH Contig": {
                "title": "Cells With V-J Spanning IGH Contig",
                "description": "Cells With V-J Spanning IGH Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With Productive IGH Contig": {
                "title": "Cells With Productive IGH Contig",
                "description": "Cells With Productive IGH Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Cells With IGL Contig": {
                "title": "Cells With IGL Contig",
                "description": "Cells With IGL Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With CDR3-annotated IGL Contig": {
                "title": "Cells With CDR3-annotated IGL Contig",
                "description": "Cells With CDR3-annotated IGL Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With V-J Spanning IGL Contig": {
                "title": "Cells With V-J Spanning IGL Contig",
                "description": "Cells With V-J Spanning IGL Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With Productive IGL Contig": {
                "title": "Cells With Productive IGL Contig",
                "description": "Cells With Productive IGL Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Cells With IGK Contig": {
                "title": "Cells With IGK Contig",
                "description": "Cells With IGK Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With CDR3-annotated IGK Contig": {
                "title": "Cells With CDR3-annotated IGK Contig",
                "description": "Cells With CDR3-annotated IGK Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With V-J Spanning IGK Contig": {
                "title": "Cells With V-J Spanning IGK Contig",
                "description": "Cells With V-J Spanning IGK Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Cells With Productive IGK Contig": {
                "title": "Cells With Productive IGK Contig",
                "description": "Cells With Productive IGK Contig",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
        }
        self.general_stats_addcols(summary_data, headers=headers)

    def barcode_rank_plot(self, umi_count_data):
        plot_data = {}
        colors = {}
        for sample in umi_count_data:
            for sub in umi_count_data[sample]:
                cur = umi_count_data[sample][sub]
                if not cur:
                    continue
                new = {}
                for k, v in cur.items():
                    new[int(k)] = v
                plot_data[sub] = new
                if "pure" in sub:
                    colors[sub] = "darkblue"
                elif "mix" in sub:
                    colors[sub] = "lightblue"
                elif "background" in sub:
                    colors[sub] = "lightgray"

        # Config for the plot
        pconfig = {
            "id": f"{ASSAY}_barcode_rank_plot",
            "title": f"{ASSAY}: Barcode Rank",
            "ylab": "UMI counts",
            "xlab": "Barcode Rank",
            "yLog": True,
            "xLog": True,
            "colors": colors,
            "ymin": 0,
            "height": 750,
        }

        return linegraph.plot(plot_data, pconfig)
