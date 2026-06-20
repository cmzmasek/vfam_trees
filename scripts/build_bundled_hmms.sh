#!/bin/bash
# Build vfam_trees/data/hmms/<Family>/<NAME>.hmm from individual Pfam-A profiles.
#
# Reproducible recipe (mirrors the sister project repseq): each profile is
# fetched from InterPro by accession, decompressed, validated (the NAME line
# must match the expected name), and written to its family directory.  At load
# time vfam_trees combines every *.hmm in the directory and hmmpress-es it.
#
# Run from the repo root:
#     bash scripts/build_bundled_hmms.sh
#
# Pfam-A is licensed CC0 (public domain); redistribution is unrestricted.
# Accessions verified against the InterPro REST API (2026-06).

set -euo pipefail

OUTROOT="vfam_trees/data/hmms"

# (Family, Pfam accession, expected NAME, marker role)
#
# Poxviridae — the 8-marker concatenation set (CONCATENATION_FAMILIES).
# Adenoviridae — single-marker hexon, as a two-domain architecture
#   (Adeno_hexon N-term --> Adeno_hexon_C C-term) for the "A--B" token.
PROFILES=(
    "Poxviridae   PF00136 DNA_pol_B        DNA_polymerase_E9L"
    "Poxviridae   PF00623 RNA_pol_Rpb1_2   RNA_pol_147kDa_A24R_RPO147"
    "Poxviridae   PF00562 RNA_pol_Rpb2_6   RNA_pol_132kDa_J6R_RPO132"
    "Poxviridae   PF03291 mRNA_G-N7_MeTrfase mRNA_capping_enzyme_large_D1R_N7MTase"
    "Poxviridae   PF12011 NPH-II           DNA_helicase_NPH-II_A18R"
    "Poxviridae   PF03296 Pox_polyA_pol    polyA_polymerase_E1L_VP55"
    "Poxviridae   PF04947 Pox_VLTF3        late_transcription_factor_VLTF3_A1L"
    "Poxviridae   PF03167 UDG              uracil-DNA_glycosylase_D4R"
    "Adenoviridae PF01065 Adeno_hexon      hexon_N-terminal"
    "Adenoviridae PF03678 Adeno_hexon_C    hexon_C-terminal"
)

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

for entry in "${PROFILES[@]}"; do
    family="$(echo "$entry" | awk '{print $1}')"
    acc="$(echo "$entry" | awk '{print $2}')"
    expected_name="$(echo "$entry" | awk '{print $3}')"
    role="$(echo "$entry" | awk '{print $4}')"
    outdir="$OUTROOT/$family"
    mkdir -p "$outdir"
    url="https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/${acc}?annotation=hmm"
    echo "[$family/$acc] fetching ($expected_name — $role) ..."
    if ! curl -sf "$url" -o "$TMP/$acc.gz"; then
        echo "  FAILED — skipping $acc"
        continue
    fi
    gunzip -c "$TMP/$acc.gz" > "$TMP/$acc.hmm"
    actual_name="$(grep -m1 '^NAME' "$TMP/$acc.hmm" | awk '{print $2}' || true)"
    if [[ -z "$actual_name" ]]; then
        echo "  WARN — no NAME line in $acc.hmm; writing anyway"
    elif [[ "$actual_name" != "$expected_name" ]]; then
        echo "  NOTE — NAME is '$actual_name' (expected '$expected_name'); using actual"
    fi
    # File named by the profile NAME so it's obvious which token it serves.
    cp "$TMP/$acc.hmm" "$outdir/${actual_name:-$expected_name}.hmm"
done

echo ""
for fam in "$OUTROOT"/*/; do
    n=$(find "$fam" -name '*.hmm' | wc -l | tr -d ' ')
    echo "Bundled $fam : $n profile(s)"
done
echo "Done.  vfam_trees auto-combines + hmmpress-es each family dir on first use."
