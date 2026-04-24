"""Main functions for the OmniBenchmark module."""

from pathlib import Path
import gzip

import pandas as pd
import scanpy as sc


def process_data(args):
    """Process data using parsed command-line arguments.

    Args:
        args: Parsed arguments from argparse containing:
            - output_dir: Output directory path
            - name: Module name
            - rawdata_h5ad: Input files for rawdata.h5ad (CLI: --rawdata.h5ad)

    Note: Input IDs with dots (e.g., 'data.raw') are converted to underscores
          in Python variable names (e.g., 'data_raw') but preserve dots in CLI args.
    """
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # print(f"Processing module: {args.name}")

    # Access stage inputs
    rawdata_h5ad_files = args.rawdata_h5ad
    # print(f"  rawdata.h5ad: {rawdata_h5ad_files}")

    # Read inputs
    input_h5ad = rawdata_h5ad_files[0]
    adata = sc.read_h5ad(input_h5ad)

    # Calculate percentage of mitochondrial genes per cell
    adata.var["mt"] = adata.var_names.str.startswith("MT-")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], layer="counts", 
        inplace=True, log1p=False, percent_top=None,
    )

    # Filtering
    if args.filter_type == "manual":
        qc = pd.DataFrame(adata.uns["qc_thresholds"])

        def get_threshold(metric, column):
            return float(qc.loc[qc["metric"] == metric, column].iloc[0])

        keep = (
            (adata.obs["n_genes_by_counts"] >= get_threshold("nFeature", "min"))
            & (adata.obs["n_genes_by_counts"] <= get_threshold("nFeature", "max"))
            & (adata.obs["pct_counts_mt"] < get_threshold("percent.mt", "max"))
            & (adata.obs["total_counts"] <= get_threshold("nCount", "max"))
        )
    else:
        raise ValueError(f"Unknown filter_type: {args.filter_type}")
    
    # Write a simple output file
    output_file = output_dir / f"{args.name}_cellids.txt.gz"
    
    with gzip.open(output_file, "wt") as f:
        for cell_id in adata.obs_names[keep]:
            f.write(f"{cell_id}\n")   