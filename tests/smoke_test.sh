#!/usr/bin/env bash
set -euo pipefail

echo "[i] Running smoke test for homology-extension..." >&2
echo "[i] Testing the RELAXED Sheriff workflow" >&2
echo "" >&2

# Check if we're in the right directory
if [ ! -f "scripts/pairwise_homology.py" ]; then
    echo "Error: Must run from homology-extension root directory" >&2
    exit 1
fi

# Create temp directory for test outputs
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

echo "[i] Creating test genome index..." >&2
# Create a simple faidx index for the test genome
cd data
samtools faidx test_genome.fa 2>/dev/null || python -c "
import sys
with open('test_genome.fa.fai', 'w') as f:
    f.write('chr1\t3060\t6\t60\t61\n')
    f.write('chr2\t2074\t3078\t60\t61\n')
"
cd ..

echo "[i] Step 1a: Converting Sheriff TSV (relaxed mode) to candidates..." >&2
python scripts/sheriff2homology.py \
    --sheriff data/sheriff_loose.example.tsv \
    --out "$TMPDIR/candidates_from_tsv.tsv"

# Check output exists
if [ ! -f "$TMPDIR/candidates_from_tsv.tsv" ]; then
    echo "FAIL: sheriff2homology.py did not create output file" >&2
    exit 1
fi

echo "[i] Step 1b: Testing legacy BED conversion..." >&2
python scripts/edit_sites2loose.py \
    --bed data/sheriff_edit_sites.example.bed \
    --out "$TMPDIR/candidates_from_bed.tsv"

# Check output exists
if [ ! -f "$TMPDIR/candidates_from_bed.tsv" ]; then
    echo "FAIL: edit_sites2loose.py did not create output file" >&2
    exit 1
fi

# Check output format
LINES=$(wc -l < "$TMPDIR/candidates_from_tsv.tsv")
if [ "$LINES" -ne "6" ]; then  # 5 data rows + 1 header
    echo "FAIL: Expected 6 lines in candidates.tsv, got $LINES" >&2
    exit 1
fi

echo "[i] Step 2: Running homology scoring (using TSV-derived candidates)..." >&2
python scripts/pairwise_homology.py \
    --guides data/guides.example.fa \
    --candidates "$TMPDIR/candidates_from_tsv.tsv" \
    --fasta data/test_genome.fa \
    --out "$TMPDIR/homology.tsv"

# Check output exists
if [ ! -f "$TMPDIR/homology.tsv" ]; then
    echo "FAIL: pairwise_homology.py did not create output file" >&2
    exit 1
fi

# Check that we have results
if ! grep -q "site_" "$TMPDIR/homology.tsv"; then
    echo "FAIL: No results found in homology.tsv" >&2
    exit 1
fi

# Basic validation - check header
EXPECTED_HEADER="site_id	guide_id	strand	genome_offset	mm_total	mm_seed	score	pam_present	aln_guide	aln_locus"
ACTUAL_HEADER=$(head -n1 "$TMPDIR/homology.tsv")
if [ "$ACTUAL_HEADER" != "$EXPECTED_HEADER" ]; then
    echo "FAIL: Unexpected header format" >&2
    echo "Expected: $EXPECTED_HEADER" >&2
    echo "Got:      $ACTUAL_HEADER" >&2
    exit 1
fi

echo "" >&2
echo "[i] Demonstrating Sheriff mode differences:" >&2
echo "  Relaxed mode sites: $(grep -c site_ data/sheriff_loose.example.tsv)" >&2
echo "  Note: 3 of 5 sites have cell_count=1 (would be filtered in strict mode)" >&2
echo "" >&2
echo "[✓] Smoke test passed!" >&2
echo "" >&2
echo "Test outputs:" >&2
echo "  - $TMPDIR/candidates_from_tsv.tsv (from Sheriff TSV)" >&2
echo "  - $TMPDIR/candidates_from_bed.tsv (from Sheriff BED)" >&2
echo "  - $TMPDIR/homology.tsv" >&2
echo "" >&2
echo "To inspect results:" >&2
echo "  head $TMPDIR/homology.tsv" >&2
echo "" >&2
echo "Remember: Always use 'sheriff call --min_cells 1 --no_blacklist' for this pipeline!" >&2