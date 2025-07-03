#!/usr/bin/env python
"""
edit_sites2loose.py
===================
Convert Sheriff's edit_sites.bed to the 5-column 'loose.tsv' format
required by pairwise_homology.py.

Sheriff outputs edit sites as a 4-column BED file with intervals that
represent windows around potential cut sites. This script:
1. Collapses each interval to its center point (1bp)
2. Removes genome build prefixes (e.g., 'hg38_chr1' → 'chr1')
3. Adds a strand column (set to '.' since Sheriff doesn't track strand)
4. Reorders columns to match the expected format

Input format (Sheriff edit_sites.bed):
--------------------------------------
hg38_chr1    1000    1100    site_001
hg38_chr1    2000    2100    site_002
...

Output format (loose.tsv):
--------------------------
chr    start    end    strand    site_id
chr1   1050     1051   .         site_001
chr1   2050     2051   .         site_002
...

Usage:
------
python edit_sites2loose.py \\
    --bed results/edit_sites.bed \\
    --out results/loose.tsv
"""
import pandas as pd
import argparse
import re
import pathlib as pl

def parse_args():
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(
        description="Convert Sheriff edit_sites.bed to loose.tsv format"
    )
    p.add_argument("--bed", required=True, 
                   help="Path to Sheriff edit_sites.bed file")
    p.add_argument("--out", default="loose.tsv",
                   help="Output TSV file for pairwise_homology.py (default: loose.tsv)")
    return p.parse_args()

def main():
    """Main conversion workflow."""
    args = parse_args()
    
    # Read the 4-column BED file from Sheriff
    print(f"Reading Sheriff edit sites from: {args.bed}")
    df = pd.read_csv(args.bed, sep="\t", header=None,
                     names=["chr", "start", "end", "site_id"])
    
    print(f"Loaded {len(df)} edit sites")
    
    # 1. Collapse intervals to their center point
    # Sheriff outputs windows (e.g., ±50bp around cut site)
    # We need just the central coordinate
    centre = ((df["start"] + df["end"]) // 2).astype(int)
    df["start"] = centre
    df["end"] = centre + 1  # BED format is half-open [start, end)
    
    # 2. Clean chromosome names
    # Sheriff prefixes chromosomes with genome build (e.g., 'hg38_chr1')
    # Remove these prefixes to match standard reference FASTA headers
    df["chr"] = df["chr"].str.replace(r"^hg\d+_", "", regex=True)
    df["chr"] = df["chr"].str.replace(r"^mm\d+_", "", regex=True)  # Also handle mouse
    
    # 3. Add strand column
    # Sheriff collapses forward and reverse strand information
    # We set to '.' (unknown) and let pairwise_homology.py try both strands
    df["strand"] = "."
    
    # 4. Reorder columns to match expected format
    df = df[["chr", "start", "end", "strand", "site_id"]]
    
    # Write output
    df.to_csv(args.out, sep="\t", index=False)
    
    # Print summary statistics
    print(f"\nConversion complete!")
    print(f"Output file: {pl.Path(args.out).resolve()}")
    print(f"Total sites: {len(df):,}")
    print(f"Unique chromosomes: {df['chr'].nunique()}")
    print(f"\nFirst few rows:")
    print(df.head())
    
    # Check for potential issues
    if df["chr"].str.contains("_").any():
        print("\nWARNING: Some chromosome names still contain underscores.")
        print("This might indicate non-standard chromosome names or scaffolds.")
        print("Affected chromosomes:")
        print(df[df["chr"].str.contains("_")]["chr"].unique())

if __name__ == "__main__":
    main()