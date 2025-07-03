#!/usr/bin/env python
"""
sheriff2homology.py
===================
Convert Sheriff's relaxed TSV output to the format expected by pairwise_homology.py.

IMPORTANT: This script expects Sheriff output from RELAXED mode:
    sheriff call -i dedup.bam --min_cells 1 --no_blacklist -o loose.tsv

The relaxed mode keeps all potential sites including:
- Single-cell events (--min_cells 1)
- Repeat-masked regions (--no_blacklist)

This gives us maximum sensitivity for homology scoring. We apply our own
filtering AFTER attaching homology scores, PAM info, etc.

Input format (Sheriff loose.tsv):
---------------------------------
chr    start    end    site_id    cell_count    ...other columns...

Output format (candidates.tsv):
-------------------------------
chr    start    end    strand    site_id

Usage:
------
python sheriff2homology.py \\
    --sheriff loose.tsv \\
    --out candidates.tsv
"""
import pandas as pd
import argparse
import sys
from pathlib import Path

def parse_args():
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(
        description="Convert Sheriff relaxed output to homology scorer format"
    )
    p.add_argument("--sheriff", required=True, 
                   help="Sheriff TSV output from relaxed mode (loose.tsv)")
    p.add_argument("--out", default="candidates.tsv",
                   help="Output TSV for pairwise_homology.py (default: candidates.tsv)")
    p.add_argument("--min-cells", type=int, default=0,
                   help="Optional: filter sites with fewer than N cells (default: 0, keep all)")
    return p.parse_args()

def validate_sheriff_output(df):
    """Check if the input looks like Sheriff output."""
    required_cols = {"chr", "start", "end"}
    
    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        print(f"ERROR: Sheriff output missing required columns: {missing}", file=sys.stderr)
        print(f"Found columns: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)
    
    # Check for common mistake - using strict mode output
    if "cell_count" in df.columns and df["cell_count"].min() > 1:
        print("WARNING: All sites have cell_count > 1.", file=sys.stderr)
        print("This might be strict mode output. For homology scoring, use:", file=sys.stderr)
        print("  sheriff call -i dedup.bam --min_cells 1 --no_blacklist -o loose.tsv", file=sys.stderr)
        print("", file=sys.stderr)

def main():
    """Main conversion workflow."""
    args = parse_args()
    
    # Read Sheriff output
    print(f"Reading Sheriff output from: {args.sheriff}")
    df = pd.read_csv(args.sheriff, sep="\t")
    
    # Validate input
    validate_sheriff_output(df)
    
    print(f"Loaded {len(df)} sites from Sheriff")
    
    # Apply optional cell count filter
    if args.min_cells > 0 and "cell_count" in df.columns:
        before = len(df)
        df = df[df["cell_count"] >= args.min_cells]
        print(f"Filtered to {len(df)} sites with >= {args.min_cells} cells (removed {before - len(df)})")
    
    # Create site_id if not present
    if "site_id" not in df.columns:
        print("Creating site_id column...")
        df["site_id"] = [f"site_{i:06d}" for i in range(len(df))]
    
    # Clean chromosome names (remove any prefixes)
    df["chr"] = df["chr"].str.replace(r"^(hg|mm)\d+_", "", regex=True)
    
    # Add strand column (Sheriff doesn't track strand after merging)
    df["strand"] = "."
    
    # Select and reorder columns
    output_cols = ["chr", "start", "end", "strand", "site_id"]
    df_out = df[output_cols].copy()
    
    # Write output
    df_out.to_csv(args.out, sep="\t", index=False)
    
    # Print summary
    print(f"\nConversion complete!")
    print(f"Output file: {Path(args.out).resolve()}")
    print(f"Total sites: {len(df_out):,}")
    print(f"Chromosomes: {', '.join(sorted(df_out['chr'].unique()[:5]))}", end="")
    if df_out['chr'].nunique() > 5:
        print(f" ... and {df_out['chr'].nunique() - 5} more")
    else:
        print()
    
    # Show example
    print(f"\nFirst few rows:")
    print(df_out.head())
    
    # Reminder about workflow
    print("\n" + "="*60)
    print("NEXT STEPS:")
    print("1. Run pairwise_homology.py on this output:")
    print(f"   python pairwise_homology.py --candidates {args.out} \\")
    print("          --guides guides.fa --fasta hg38.fa --out homology.tsv")
    print("\n2. Merge with Sheriff metadata and filter as needed")
    print("="*60)

if __name__ == "__main__":
    main()