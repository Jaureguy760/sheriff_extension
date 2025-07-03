#!/usr/bin/env python
"""
pairwise_homology.py
====================
Guide–locus local alignment with seed-mismatch counting.

This script performs seed-aware alignment of CRISPR guides to candidate sites
identified by Sheriff. It uses local Smith-Waterman alignment to find the best
match for each guide at each candidate site, considering both strands.

Key features:
- Searches ±100bp around each candidate cut site
- Counts total mismatches and seed region (positions 1-8) mismatches
- Identifies PAM presence (NGG/NAG) immediately 3' of the 20-mer
- Returns the best-scoring alignment for each site-guide pair

Usage:
------
python pairwise_homology.py \\
    --guides guides.fa \\
    --candidates sheriff/loose.tsv \\
    --fasta hg38.fa \\
    --out homology.tsv
"""
import argparse
from collections import namedtuple
from pathlib import Path

from Bio import pairwise2, SeqIO
from Bio.Seq import Seq
from pyfaidx import Fasta
import pandas as pd

# --------------------------------------------------------------------------- #
# Helper dataclass to store alignment hits
Hit = namedtuple(
    "Hit",
    "guide_id strand offset score mm_total mm_seed aln_query aln_sub"
)

# --------------------------------------------------------------------------- #
def revcomp(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    return str(Seq(seq).reverse_complement())

def load_guides(fasta_file: str):
    """Load guide sequences from a FASTA file.
    
    Parameters:
    -----------
    fasta_file : str
        Path to FASTA file containing 20nt guide sequences
        
    Returns:
    --------
    list of Bio.SeqRecord objects
    """
    return [rec for rec in SeqIO.parse(fasta_file, "fasta")]

def load_candidates(tsv: str):
    """Load candidate sites from Sheriff's loose.tsv format.
    
    Parameters:
    -----------
    tsv : str
        Path to TSV file with columns: chr, start, end, strand, site_id
        
    Returns:
    --------
    pandas.DataFrame with candidate sites
    """
    df = pd.read_csv(tsv, sep="\t")
    required = {"chr", "start", "end", "site_id"}
    assert required.issubset(df.columns), f"TSV must contain {required}"
    return df

def pam_in_window(window: str) -> bool:
    """Check if NGG or NAG PAM is present in the last 3 nucleotides.
    
    Parameters:
    -----------
    window : str
        DNA sequence to check (PAM should be in positions -3 to end)
        
    Returns:
    --------
    bool indicating PAM presence
    """
    # Check NGG / NAG immediately 3' of the 20-mer
    return ("GG" in window[-3:]) or ("AG" in window[-3:])

# --------------------------------------------------------------------------- #
def main(args):
    """Main workflow: align guides to candidate sites and score homology."""
    
    # Load reference genome with uppercase sequences
    fa = Fasta(args.fasta, as_raw=True, sequence_always_upper=True)
    
    # Load guide sequences and candidate sites
    guides = load_guides(args.guides)
    cand = load_candidates(args.candidates)
    
    print(f"Loaded {len(guides)} guides and {len(cand)} candidate sites")

    # Open output file and write header
    with open(args.out, "w") as out:
        header = ["site_id", "guide_id", "strand", "genome_offset",
                  "mm_total", "mm_seed", "score", "pam_present",
                  "aln_guide", "aln_locus"]
        out.write("\t".join(header) + "\n")

        # Process each candidate site
        for idx, (_, site) in enumerate(cand.iterrows()):
            if idx % 100 == 0:
                print(f"Processing site {idx+1}/{len(cand)}")
                
            chrom, start, end = site["chr"], site["start"], site["end"]
            center = (start + end) // 2
            
            # Extract ±100bp window around the cut site
            window = fa[chrom][center-100:center+100]
            best_hit = None

            # Try aligning each guide
            for g in guides:
                guide_seq = str(g.seq).upper()
                
                # Check both strands
                for strand, seq in [("+", window), ("-", revcomp(window))]:
                    # Slide 20nt window across the sequence
                    for offset in range(len(seq) - 20):
                        sub = seq[offset:offset+20]
                        
                        # Perform local alignment
                        # Scoring: match=+1, mismatch=-1, gap_open=-0.8, gap_extend=-0.5
                        aln = pairwise2.align.localms(
                            guide_seq, sub,
                            1, -1, -0.8, -0.5,
                            one_alignment_only=True
                        )[0]

                        # Count mismatches
                        mm_tot = sum(b1 != b2 for b1, b2 in zip(aln.seqA, aln.seqB))
                        mm_seed = sum(b1 != b2 for b1, b2 in zip(aln.seqA[:8], aln.seqB[:8]))

                        # Store hit information
                        # offset-100 converts to genome-relative position
                        hit = Hit(g.id, strand, offset - 100, aln.score,
                                  mm_tot, mm_seed, aln.seqA, aln.seqB)

                        # Keep best-scoring hit
                        if (best_hit is None) or (hit.score > best_hit.score):
                            best_hit = hit

            # Write best hit for this site
            if best_hit:
                # Check for PAM presence at the expected position
                # PAM should be at positions 17-19 (0-indexed) relative to window start
                pam_flag = pam_in_window(window[best_hit.offset+100+17:best_hit.offset+100+20])
                
                out.write("\t".join(map(str, [
                    site["site_id"], best_hit.guide_id, best_hit.strand,
                    best_hit.offset, best_hit.mm_total, best_hit.mm_seed,
                    round(best_hit.score, 2), int(pam_flag),
                    best_hit.aln_query, best_hit.aln_sub
                ])) + "\n")
    
    print(f"Done! Results written to {args.out}")

# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Score guide-locus homology with seed mismatches")
    p.add_argument("--guides", required=True, 
                   help="FASTA file of 20-nt protospacer sequences")
    p.add_argument("--candidates", required=True, 
                   help="Sheriff loose.tsv with candidate cut sites")
    p.add_argument("--fasta", required=True, 
                   help="Reference genome FASTA (must be indexed with samtools faidx)")
    p.add_argument("--out", default="homology.tsv",
                   help="Output TSV with homology scores (default: homology.tsv)")
    main(p.parse_args())