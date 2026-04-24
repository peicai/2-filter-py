"""Main functions for the OmniBenchmark module."""

from pathlib import Path


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

    print(f"Processing module: {args.name}")

    # Access stage inputs
    rawdata_h5ad_files = args.rawdata_h5ad
    print(f"  rawdata.h5ad: {rawdata_h5ad_files}")

    # TODO: Implement your processing logic here
    # Example: Read inputs, process, write outputs

    # Write a simple output file
    output_file = output_dir / f"{args.name}_result.txt"
    with open(output_file, 'w') as f:
        f.write(f"Processed module: {args.name}\n")
        f.write(f"rawdata.h5ad: {len(rawdata_h5ad_files)} file(s)\n")

    print(f"Results written to: {output_file}")
