# CRISPR-Sheriff → Homology Extension Tutorial  
*In-silico off-target triage with seed-mismatch scoring*

---

## Table of Contents
1. [Overview](#overview)
2. [Project Layout](#project-layout)
3. [Installation](#installation)
4. [Quick Start](#quick-start)
5. [Step-by-Step Tutorial](#step-by-step-tutorial)
   - [Step 0: Run Sheriff](#step-0-run-the-upstream-sheriff-workflow)
   - [Step 1: Convert edit sites](#step-1-convert-sheriff-edit-sites-to-loosetsv)
   - [Step 2: Score homology](#step-2-run-pairwise_homologypy)
   - [Step 3: Merge results](#step-3-merge-homology-scores-with-other-evidence)
6. [Understanding the Scripts](#understanding-the-scripts)
7. [Output Format](#output-format)
8. [Troubleshooting](#troubleshooting)
9. [FAQ](#faq)

---

## Overview

This repository provides tools to extend [CRISPR-Sheriff](https://github.com/pinellolab/sheriff) with local sequence homology scoring. It helps identify and prioritize potential off-target sites by:

- **Seed-aware alignment**: Counts mismatches in the critical seed region (positions 1-8)
- **PAM detection**: Identifies canonical NGG/NAG PAMs near candidate sites
- **Strand-specific scoring**: Tests both DNA strands for best alignment
- **Integration-ready**: Outputs merge with Sheriff's results and other off-target predictors

⚠️ **IMPORTANT**: This pipeline requires Sheriff to be run in **relaxed mode** to capture all potential sites:
```bash
sheriff call -i dedup.bam --min_cells 1 --no_blacklist -o loose.tsv
```

The workflow consists of three main scripts:
1. `sheriff2homology.py` - Converts Sheriff's relaxed TSV output to candidate format
2. `edit_sites2loose.py` - Alternative converter for Sheriff's BED output (legacy)
3. `pairwise_homology.py` - Performs guide-to-genome alignment and scoring

---

## Project Layout

```
homology-extension/
├── environment.yml          # Conda environment specification
├── README.md               # This file
├── scripts/
│   ├── sheriff2homology.py # Sheriff TSV → candidates.tsv (PRIMARY)
│   ├── edit_sites2loose.py # Sheriff BED → loose.tsv (legacy)
│   └── pairwise_homology.py # Main homology scoring script
├── data/                   # Example data files
│   ├── guides.example.fa
│   ├── sheriff_edit_sites.example.bed
│   └── test_genome.fa
├── tests/
│   └── smoke_test.sh       # Basic functionality test
└── .gitignore
```

---

## Installation

### Prerequisites
- Conda or Mamba package manager
- Access to a reference genome (e.g., hg38.fa)
- Sheriff output files (edit_sites.bed)

### Setup

```bash
# Clone the repository
git clone https://github.com/YOUR-LAB/homology-extension.git
cd homology-extension

# Create and activate the conda environment
conda env create -f environment.yml
conda activate sheriff-homology

# Run the smoke test to verify installation
bash tests/smoke_test.sh
```

---

## Quick Start

```bash
# 1. Run Sheriff in RELAXED mode (critical!)
sheriff call -i dedup.bam --min_cells 1 --no_blacklist -o loose.tsv

# 2. Convert Sheriff output to candidates format
python scripts/sheriff2homology.py \
    --sheriff loose.tsv \
    --out candidates.tsv

# 3. Score homology for your guides
python scripts/pairwise_homology.py \
    --guides guides.fa \
    --candidates candidates.tsv \
    --fasta /path/to/hg38.fa \
    --out homology_scores.tsv
```

---

## Step-by-Step Tutorial

### Step 0: Run Sheriff in RELAXED mode

Before using these tools, you need Sheriff's output files. **CRITICAL**: Run Sheriff in relaxed mode to capture all potential sites for homology scoring.

#### Sheriff Modes Explained

| Mode | Command | Use Case |
|------|---------|----------|
| **Strict** (default) | `sheriff call -i dedup.bam -o strict.tsv` | Publication figures; drops low-cell-count events |
| **Relaxed** ✅ | `sheriff call -i dedup.bam --min_cells 1 --no_blacklist -o loose.tsv` | Homology scoring; keeps ALL sites for re-ranking |

```bash
# For homology scoring, ALWAYS use relaxed mode:
sheriff call -i dedup.bam --min_cells 1 --no_blacklist -o loose.tsv
```

Why relaxed mode?
- `--min_cells 1`: Keeps single-cell events that might be real off-targets
- `--no_blacklist`: Retains repeat-masked regions for comprehensive scoring
- We apply our own filtering AFTER attaching homology scores

Sheriff produces:
- `loose.tsv` - All candidate sites (relaxed mode output)
- `edit_site_info.txt` - Detailed metadata per site
- Various per-cell edit count files

### Step 1: Convert Sheriff output to candidates format

The `pairwise_homology.py` script expects a specific 5-column format. Use the appropriate converter:

#### Option A: From Sheriff TSV (RECOMMENDED)
```bash
python scripts/sheriff2homology.py \
    --sheriff loose.tsv \
    --out candidates.tsv
```

#### Option B: From Sheriff BED (legacy workflow)
```bash
python scripts/edit_sites2loose.py \
    --bed results/edit_sites.bed \
    --out candidates.tsv
```

**What this does:**
- Centers each BED interval to a single base pair (the predicted cut site)
- Removes genome build prefixes (e.g., `hg38_chr1` → `chr1`)
- Adds a strand column (set to `.` for unknown)
- Reorders columns to: `chr start end strand site_id`

**Example transformation:**
```
# Input (Sheriff BED):
hg38_chr1    1000    1100    site_001

# Output (loose.tsv):
chr1    1050    1051    .    site_001
```

### Step 2: Run `pairwise_homology.py`

Now score each guide against each candidate site:

```bash
python scripts/pairwise_homology.py \
    --guides data/guides.fa \
    --candidates results/loose.tsv \
    --fasta /path/to/hg38.fa \
    --out results/homology_scores.tsv
```

**What this does for each site × guide combination:**

1. **Extracts genomic context**: Pulls ±100bp around the cut site
2. **Aligns on both strands**: Uses Smith-Waterman local alignment
3. **Counts mismatches**: Total and seed region (positions 1-8)
4. **Checks PAM presence**: Looks for NGG/NAG at the expected position
5. **Reports best match**: Keeps highest-scoring alignment

**Runtime:** ~0.2 seconds per 100 sites per guide (pure Python implementation)

### Step 3: Merge homology scores with other evidence

Combine the homology scores with Sheriff's cell-level data and other off-target predictors:

```python
import pandas as pd

# Load all data sources
homology = pd.read_csv("results/homology_scores.tsv", sep="\t")
sheriff  = pd.read_csv("results/edit_site_info.txt", sep="\t")
offinder = pd.read_csv("results/cas_offinder.txt", sep="\t")  # If available

# Merge on site_id
combined = homology.merge(sheriff, on="site_id", how="left")
combined = combined.merge(offinder, on="site_id", how="left")

# Filter by homology criteria
filtered = combined[
    (combined["mm_total"] <= 4) &      # Max 4 total mismatches
    (combined["mm_seed"] <= 1) &       # Max 1 seed mismatch
    (combined["pam_present"] == 1)     # Canonical PAM required
]

# Sort by evidence strength
filtered = filtered.sort_values(
    ["score", "cell_count"],           # Homology score, then cell support
    ascending=[False, False]
)

filtered.to_csv("results/prioritized_offtargets.tsv", sep="\t", index=False)
```

---

## Understanding the Scripts

### `edit_sites2loose.py`

This adapter script handles format conversion between Sheriff and the homology scorer.

**Key operations:**
```python
# Collapse window to center point
centre = ((df["start"] + df["end"]) // 2).astype(int)
df["start"] = centre
df["end"] = centre + 1  # BED is half-open

# Clean chromosome names
df["chr"] = df["chr"].str.replace(r"^hg\d+_", "", regex=True)

# Add unknown strand
df["strand"] = "."
```

### `pairwise_homology.py`

The main scoring engine that performs local alignment.

**Core algorithm:**
```python
# For each candidate site
for site in candidates:
    # Extract ±100bp window
    window = genome[chr][center-100:center+100]
    
    # Try each guide
    for guide in guides:
        # Check both strands
        for strand in ["+", "-"]:
            # Slide 20nt window
            for offset in range(len(seq) - 20):
                # Local alignment
                alignment = pairwise2.align.localms(
                    guide_seq, genomic_seq,
                    match=1, mismatch=-1, 
                    gap_open=-0.8, gap_extend=-0.5
                )
                
                # Count mismatches
                mm_total = count_mismatches(alignment)
                mm_seed = count_mismatches(alignment[:8])
                
                # Keep best scoring
                if score > best_score:
                    best_hit = current_hit
```

**Scoring parameters:**
- Match: +1
- Mismatch: -1
- Gap open: -0.8
- Gap extend: -0.5

---

## Output Format

The `homology_scores.tsv` file contains:

| Column | Description |
|--------|-------------|
| `site_id` | Sheriff's site identifier |
| `guide_id` | Guide name from FASTA header |
| `strand` | Best-scoring strand (+/-) |
| `genome_offset` | Position relative to cut site (-100 to +100) |
| `mm_total` | Total mismatches across 20nt |
| `mm_seed` | Mismatches in seed region (pos 1-8) |
| `score` | Smith-Waterman alignment score |
| `pam_present` | 1 if NGG/NAG found, 0 otherwise |
| `aln_guide` | Aligned guide sequence |
| `aln_locus` | Aligned genomic sequence |

**Example row:**
```
site_001  gRNA_1  +  -5  2  0  17.2  1  GAGTCCGAGCAGAAGAAGAA  GAGTCCGAGCAGAAGAAGAA
```

---

## Troubleshooting

### Common Issues

**Issue: "Chromosome not found in reference"**
- Check that chromosome names match between your files
- The converter removes `hg38_` prefixes - ensure your reference uses standard names (chr1, not hg38_chr1)

**Issue: "No PAM detected for on-target site"**
- Verify your guides include only the 20nt protospacer (not the PAM)
- Check that the reference genome matches your experimental organism

**Issue: Script runs slowly**
- Consider parallelizing by splitting the candidate file
- For >10k sites, consider using a cluster with GNU parallel

### Validation Checks

| Check | Expected Result |
|-------|-----------------|
| On-target alignment | mm_total ≤ 1, mm_seed = 0, pam_present = 1 |
| Known off-target | Appropriate mm_total/mm_seed for validated site |
| Runtime | <1 second per 1000 site-guide pairs |

---

## FAQ

**Q: My guides are 19nt or 21nt, not 20nt. What should I do?**

A: Modify line containing `range(len(seq) - 20)` in `pairwise_homology.py`:
```python
# Change 20 to your guide length
for offset in range(len(seq) - 19):  # For 19nt guides
```

**Q: Can I use a different scoring matrix?**

A: Yes, modify the `pairwise2.align.localms()` parameters:
```python
# Current: match=1, mismatch=-1, gap_open=-0.8, gap_extend=-0.5
aln = pairwise2.align.localms(
    guide_seq, sub,
    2, -3, -5, -2,  # Your custom scores
    one_alignment_only=True
)
```

**Q: How do I filter results?**

A: Common filtering criteria:
```python
# Stringent (likely off-targets)
df[(df["mm_total"] <= 3) & (df["mm_seed"] == 0) & (df["pam_present"] == 1)]

# Permissive (possible off-targets)
df[(df["mm_total"] <= 5) & (df["mm_seed"] <= 2)]

# Exclude non-canonical PAMs
df[df["pam_present"] == 1]
```

**Q: Can I run this on mouse (mm10) or other genomes?**

A: Yes! Just provide the appropriate reference FASTA. The edit_sites2loose.py script already handles mm10 prefixes.

**Q: How does this compare to Cas-OFFinder or CRISPOR?**

A: This tool is complementary:
- **Cas-OFFinder**: Fast genome-wide search with mismatches/bulges
- **CRISPOR**: Aggregates multiple off-target predictors
- **This tool**: Focuses on Sheriff-identified sites with detailed alignment

Use all three for comprehensive off-target assessment.

---

## Citation

If you use this tool, please cite:
- The Sheriff paper (when published)
- BioPython: Cock et al., 2009
- This repository

---

## License

MIT License - See LICENSE file for details

---

## Support

For issues or questions:
1. Check the [FAQ](#faq) and [Troubleshooting](#troubleshooting) sections
2. Open an issue on GitHub
3. Contact: your-email@institution.edu

Happy off-target hunting! 🧬🔬