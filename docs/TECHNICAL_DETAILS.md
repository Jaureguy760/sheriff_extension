# Technical Problem Analysis

## The Core Challenge

CRISPR off-target detection requires balancing two competing needs:
1. **Empirical evidence**: What actually got edited in cells (Sheriff)
2. **Sequence prediction**: What could theoretically be edited (homology/OFFinder)

Neither approach alone is sufficient:
- Empirical-only: Misses low-frequency events, provides no mechanistic insight
- Prediction-only: Many false positives, no cellular validation

## Why Sheriff's Default Mode Falls Short

### The Filtering Problem

```python
# Sheriff default behavior (pseudocode)
def sheriff_strict_mode(sites):
    sites = filter(lambda s: s.cell_count >= 2, sites)  # Loses single-cell events
    sites = filter(lambda s: s not in blacklist, sites)  # Loses repeat regions
    return sites  # May have discarded real off-targets!
```

### Real Example
```
Guide: GAGTCCGAGCAGAAGAAGAA

Scenario 1: Perfect off-target in heterochromatin
- Location: Repeat-masked region
- Cell count: 1 (hard to access)
- Sheriff strict: ❌ Filtered out
- Homology: 0 mismatches, perfect PAM
- Reality: Dangerous off-target missed!

Scenario 2: Degenerate match in active promoter  
- Location: Open chromatin
- Cell count: 50
- Sheriff strict: ✓ Reported
- Homology: 6 mismatches, no PAM
- Reality: Likely false positive from DNA damage
```

## Our Solution Architecture

### 1. Data Flow Design

```
                    ┌─────────────────┐
                    │   Sheriff Call  │
                    │  (Relaxed Mode) │
                    └────────┬────────┘
                             │
                    ┌────────▼────────┐
                    │  sheriff2homology│
                    │   Standardize   │
                    └────────┬────────┘
                             │
                ┌────────────┴────────────┐
                │                         │
       ┌────────▼────────┐      ┌────────▼────────┐
       │ pairwise_homology│      │  Cas-OFFinder   │
       │  Seed scoring   │      │ Bulge detection │
       └────────┬────────┘      └────────┬────────┘
                │                         │
                └────────────┬────────────┘
                             │
                    ┌────────▼────────┐
                    │  Merge & Filter │
                    │  (User-defined) │
                    └─────────────────┘
```

### 2. Algorithm Details

#### Homology Scoring Algorithm
```python
def score_site_guide_pair(site_seq, guide_seq):
    best_score = -inf
    
    # Try both strands
    for strand in ['+', '-']:
        seq = site_seq if strand == '+' else reverse_complement(site_seq)
        
        # Slide guide along sequence
        for offset in range(len(seq) - 20):
            target = seq[offset:offset+20]
            
            # Local alignment with specific penalties
            alignment = pairwise2.align.localms(
                guide_seq, target,
                match=1,       # Reward matches
                mismatch=-1,   # Penalize mismatches
                gap_open=-0.8, # Slightly discourage gaps
                gap_extend=-0.5
            )
            
            # Count mismatches with seed emphasis
            mm_total = count_mismatches(alignment)
            mm_seed = count_mismatches(alignment[:8])  # PAM-proximal region
            
            # Check PAM (NGG/NAG)
            pam = seq[offset+20:offset+23]
            pam_valid = pam[1:3] in ['GG', 'AG']
            
            # Update best hit
            if alignment.score > best_score:
                best_score = alignment.score
                best_hit = Hit(strand, offset, mm_total, mm_seed, pam_valid)
    
    return best_hit
```

#### Key Design Decisions

1. **Local vs Global Alignment**: Local allows for small indels
2. **Scoring Matrix**: Tuned for CRISPR (gaps less common than mismatches)
3. **Seed Region**: Positions 1-8 are PAM-proximal, most critical
4. **PAM Flexibility**: Checks NGG and NAG, expandable for other Cas variants

### 3. Performance Considerations

```python
# Naive approach: O(sites × guides × genome)
# Our approach: O(sites × guides × window)

window_size = 200  # ±100bp around cut site
genome_size = 3e9
speedup = genome_size / window_size  # ~15,000,000x faster!
```

### 4. Memory Efficiency

- Uses pyfaidx for lazy FASTA loading (no genome in memory)
- Streaming processing (one site at a time)
- Suitable for laptop-scale computation

## Integration Points

### Input Formats

1. **Sheriff TSV** (preferred)
```
chr  start  end  site_id  cell_count  ...
chr1 1000   1001 site_001 1           ...
```

2. **Sheriff BED** (legacy)
```
chr1  950   1050  site_001
```

### Output Format
```
site_id  guide_id  strand  offset  mm_total  mm_seed  score  pam_present  aln_guide  aln_locus
site_001 gRNA_1    +       -5      2         0        17.2   1            GAGTCC...  GAGTCC...
```

### Merge-Ready Design
- Common `site_id` key across all files
- Tab-separated for easy pandas/R integration
- Preserves all information for custom filtering

## Validation Strategy

### 1. Positive Controls
- Known on-target should have mm_total=0, mm_seed=0
- Published off-targets should match expected mismatch counts

### 2. Negative Controls  
- Random genomic sites should have high mismatch counts
- Scrambled guides should show poor scores genome-wide

### 3. Concordance Checks
- High-cell-count Sheriff sites should have good homology
- Perfect homology sites should eventually show editing (if accessible)

## Common Pitfalls Addressed

1. **Chromosome naming**: Automatically strips prefixes (hg38_chr1 → chr1)
2. **Strand bias**: Always checks both strands
3. **PAM assumptions**: Configurable PAM checking
4. **Edge cases**: Handles sites near chromosome ends
5. **Version conflicts**: Minimal dependencies, pinned versions

## Extensibility

The modular design allows easy additions:
- New scoring matrices for Cas12/Cas13
- Machine learning rescoring
- Epigenetic accessibility integration
- Batch cloud processing

This architecture ensures no off-target is missed while providing the tools to intelligently prioritize validation efforts.