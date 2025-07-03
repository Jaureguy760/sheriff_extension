# 🚀 Quick Start Guide for CRISPR Off-Target Analysis

## ⚡ 5-Minute Setup

```bash
# 1. Clone and enter
git clone https://github.com/Jaureguy760/sheriff_extension.git
cd sheriff_extension

# 2. Create environment
conda env create -f environment.yml
conda activate sheriff-homology

# 3. Test installation
bash tests/smoke_test.sh
```

## 🎯 The ONE Critical Rule

**ALWAYS run Sheriff in relaxed mode for this pipeline:**

```bash
# ✅ CORRECT - Relaxed mode
sheriff call -i dedup.bam --min_cells 1 --no_blacklist -o loose.tsv

# ❌ WRONG - Default strict mode (misses off-targets!)
sheriff call -i dedup.bam -o strict.tsv
```

## 📋 Complete Workflow

### Step 1: Generate Sheriff Results (Relaxed Mode)
```bash
sheriff call -i your_dedup.bam \
    --min_cells 1 \
    --no_blacklist \
    -o sheriff_loose.tsv
```

### Step 2: Convert to Candidate Format
```bash
python scripts/sheriff2homology.py \
    --sheriff sheriff_loose.tsv \
    --out candidates.tsv
```

### Step 3: Score Homology
```bash
python scripts/pairwise_homology.py \
    --guides your_guides.fa \
    --candidates candidates.tsv \
    --fasta /path/to/hg38.fa \
    --out homology_scores.tsv
```

### Step 4: Filter Results (Python)
```python
import pandas as pd

# Load results
hom = pd.read_csv("homology_scores.tsv", sep="\t")

# Filter high-confidence off-targets
filtered = hom[
    (hom["mm_total"] <= 3) &     # Max 3 total mismatches
    (hom["mm_seed"] <= 1) &      # Max 1 seed mismatch  
    (hom["pam_present"] == 1)    # Has NGG/NAG PAM
]

# Save filtered results
filtered.to_csv("high_priority_offtargets.tsv", sep="\t", index=False)
print(f"Found {len(filtered)} high-priority off-targets to validate")
```

## 🔍 Understanding Your Results

### Output Columns Explained
- `mm_total`: Total mismatches (0 = perfect match)
- `mm_seed`: Mismatches in critical seed region (positions 1-8)
- `score`: Alignment score (higher = better match)
- `pam_present`: 1 if NGG/NAG PAM detected
- `strand`: Which DNA strand matches best

### Interpretation Guide
| mm_total | mm_seed | PAM | Risk Level | Action |
|----------|---------|-----|------------|---------|
| 0-1 | 0 | ✓ | 🔴 HIGH | Validate immediately |
| 2-3 | 0-1 | ✓ | 🟡 MEDIUM | Validate if possible |
| 4+ | 2+ | ✓ | 🟢 LOW | Monitor only |
| Any | Any | ✗ | 🟢 LOW | Non-canonical, unlikely |

## 🚨 Common Mistakes to Avoid

1. **Using Sheriff strict mode** - You'll miss real off-targets!
2. **Wrong reference genome** - Make sure versions match (hg38, mm10, etc.)
3. **Including PAM in guide sequences** - Use 20nt protospacer only
4. **Forgetting to index reference** - Run `samtools faidx reference.fa` first

## 💡 Pro Tips

### Parallel Processing
```bash
# Split candidates for parallel processing
split -l 1000 candidates.tsv chunk_

# Run multiple jobs
for chunk in chunk_*; do
    python scripts/pairwise_homology.py \
        --guides guides.fa \
        --candidates $chunk \
        --fasta hg38.fa \
        --out ${chunk}.homology.tsv &
done
wait

# Merge results
cat chunk_*.homology.tsv | grep -v "^site_id" > all_homology_body.tsv
head -1 chunk_aa.homology.tsv > all_homology.tsv
cat all_homology_body.tsv >> all_homology.tsv
```

### Quick Validation Check
```bash
# Your on-target site should show up as perfect match
grep "your_target_site_id" homology_scores.tsv
# Expected: mm_total=0, mm_seed=0, pam_present=1
```

## 📊 Typical Results

For a typical genome-editing experiment:
- Sheriff relaxed mode: 5,000-50,000 sites
- After homology scoring: Same number (all sites scored)
- High-priority off-targets: 10-100 sites
- Validated true off-targets: 2-20 sites

## 🆘 Getting Help

1. **Check the test worked**: `bash tests/smoke_test.sh`
2. **Verify Sheriff mode**: Look for sites with cell_count=1 in your TSV
3. **Check guide format**: Should be exactly 20nt, no PAM
4. **Reference genome**: Chromosome names should match (chr1, not 1 or hg38_chr1)

## 📚 Next Steps

Once comfortable with the basics:
1. Read `docs/TECHNICAL_DETAILS.md` for algorithm details
2. Customize filtering criteria for your needs
3. Integrate with Cas-OFFinder results
4. Add custom blacklists/whitelists
5. Generate publication-ready figures

Remember: **Start with relaxed Sheriff mode** - you can always filter later, but you can't recover missed sites!