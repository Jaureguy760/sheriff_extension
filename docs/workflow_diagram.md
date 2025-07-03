# Sheriff → Homology Extension Workflow

## Complete Pipeline Flow

```mermaid
graph TD
    A[Sequencing Data<br/>dedup.bam] --> B{Sheriff Mode?}
    
    B -->|Strict Mode<br/>DEFAULT| C[sheriff call -i dedup.bam<br/>-o strict.tsv]
    B -->|Relaxed Mode<br/>REQUIRED ✓| D[sheriff call -i dedup.bam<br/>--min_cells 1 --no_blacklist<br/>-o loose.tsv]
    
    C --> E[❌ Too Few Sites<br/>Missing real off-targets]
    D --> F[✓ All Potential Sites<br/>Ready for scoring]
    
    F --> G[sheriff2homology.py]
    G --> H[candidates.tsv]
    
    H --> I[pairwise_homology.py<br/>+ guides.fa<br/>+ reference.fa]
    I --> J[homology.tsv]
    
    J --> K[Merge & Filter]
    F --> K
    L[Cas-OFFinder results] --> K
    M[PAM scan results] --> K
    
    K --> N[Final Filtered<br/>Off-Target List]
    
    style D fill:#90EE90
    style F fill:#90EE90
    style E fill:#FFB6C1
```

## Why Relaxed Mode is Critical

### Sheriff Mode Comparison

| Aspect | Strict Mode (Default) | Relaxed Mode (Required) |
|--------|----------------------|-------------------------|
| Command | `sheriff call -i dedup.bam` | `sheriff call -i dedup.bam --min_cells 1 --no_blacklist` |
| Cell threshold | Drops sites with <2 cells | Keeps ALL sites (≥1 cell) |
| Blacklist | Applies repeat mask | No blacklist filtering |
| Output sites | ~100s-1000s | ~1000s-10000s |
| Use case | Publication figures | Homology scoring input |

### Example Impact

```bash
# Strict mode might show:
site_001    chr1:1000    cell_count=5    ✓ Kept
site_002    chr1:2000    cell_count=1    ❌ Filtered out
site_003    chr1:3000    cell_count=1    ❌ Filtered out

# Relaxed mode shows ALL:
site_001    chr1:1000    cell_count=5    ✓
site_002    chr1:2000    cell_count=1    ✓  <- Could be real off-target!
site_003    chr1:3000    cell_count=1    ✓  <- Perfect sequence match?
```

The single-cell sites (site_002, site_003) might have:
- Perfect guide homology (0 mismatches)
- Strong PAM sequences
- Co-occurrence with known off-targets

**We can only discover this by starting with relaxed mode!**

## File Format Evolution

```
Sheriff BAM → Sheriff TSV (relaxed) → candidates.tsv → homology.tsv → filtered results
           ↓
    [Alternative: Sheriff BED → edit_sites2loose.py → candidates.tsv]
```