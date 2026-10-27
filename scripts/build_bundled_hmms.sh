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
# Accessions verified against the InterPro REST API (2026-06; RdRp set 2026-10-20;
# RT and Rep sets 2026-10-21; Mamonoviridae 2026-10-26; NCLDV 2026-10-27).

set -euo pipefail

OUTROOT="vfam_trees/data/hmms"

# (Family, Pfam accession, expected NAME, marker role)
#
# Poxviridae — the 8-marker concatenation set (CONCATENATION_FAMILIES).
# Adenoviridae — single-marker hexon, as a two-domain architecture
#   (Adeno_hexon N-term --> Adeno_hexon_C C-term) for the "A--B" token.
# Mamonoviridae — single-marker DNA polymerase B (same profile as Poxviridae's).
# NCLDV — shared set for the Asfarviridae / Pithoviridae HMM concatenations.
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
    # Mamonoviridae — single-marker family-B DNA polymerase (2026-10-26).
    "Mamonoviridae PF00136 DNA_pol_B       DNA_polymerase_B"
    # NCLDV — HMM marker set shared by the Asfarviridae and Pithoviridae
    #   concatenations (2026-10-27); accessions checked against InterPro.
    "NCLDV PF00136 DNA_pol_B               DNA_polymerase_B"
    "NCLDV PF04451 Capsid_NCLDV            major_capsid_protein_(ASFV_p72)"
    "NCLDV PF00623 RNA_pol_Rpb1_2          RNA_polymerase_largest_subunit"
    "NCLDV PF00562 RNA_pol_Rpb2_6          RNA_polymerase_second_subunit"
    "NCLDV PF08706 D5_N                    D5-like_primase-helicase"
    "NCLDV PF04947 Pox_VLTF3               late_transcription_factor_VLTF-3"
    "NCLDV PF01331 mRNA_cap_enzyme         capping_enzyme_guanylyltransferase"
    "NCLDV PF03291 mRNA_G-N7_MeTrfase      capping_enzyme_N7-methyltransferase"
    "NCLDV PF00521 DNA_topoisoIV           topoisomerase_II"
    # RdRp — RNA-dependent RNA polymerase profiles for the RNA-virus families
    #   (RDRP_FAMILIES); all tokens OR'd, accessions verified 2026-10-20.
    "RdRp PF00680 RdRP_1                    picorna-/calici-/potyvirus-like_RdRp"
    "RdRp PF00978 RdRP_2                    alphavirus-like_supergroup_RdRp"
    "RdRp PF00998 RdRP_3                    flavi-/tombus-like_and_HCV_RdRp"
    "RdRp PF02123 RdRP_4                    toti-/luteo-/rotavirus_RdRp"
    "RdRp PF07925 RdRP_5                    reovirus_lambda3_RdRp"
    "RdRp PF00946 Mononeg_RNA_pol           Mononegavirales_L"
    "RdRp PF04196 Bunya_RdRp                bunyavirus_L"
    "RdRp PF04197 Birna_RdRp_palm           birnavirus_VP1_palm"
    "RdRp PF06317 Arena_RNA_pol             arenavirus_L"
    "RdRp PF00602 Flu_PB1                   influenza_PB1"
    "RdRp PF00972 Flavi_NS5                 flavivirus_NS5_fingers_palm"
    "RdRp PF20483 Flavi_NS5_thumb           flavivirus_NS5_thumb"
    "RdRp PF03431 RNA_replicase_B           levivirus_replicase"
    "RdRp PF05919 Mitovir_RNA_pol           mitovirus_RdRp"
    "RdRp PF05788 Orbi_VP1                  orbivirus_VP1"
    "RdRp PF06478 CoV_RPol_N                coronavirus_nsp12_N-terminal"
    "RdRp PF22152 Permu_RdRp_palm           permutotetravirus_palm"
    "RdRp PF27669 RNA-dir_pol_phytoreovirus phytoreovirus_RdRp"
    "RdRp PF27845 Seadorna_RNAP             seadornavirus_VP1"
    "RdRp PF22212 CPV_RdRP_pol_dom          cypovirus_polymerase_domain"
    "RdRp PF08467 Luteo_P1-P2               luteovirus_P1-P2_replicase"
    "RdRp PF21561 L_thumb_ring_vir          bunyavirus_L_thumb_ring"
    "RdRp PF12603 L_PA-C-like               bunyavirus_L_PA-C-like"
    # RT — reverse transcriptase of the reverse-transcribing families
    #   (RT_FAMILY_PROFILES); accessions verified 2026-10-21.
    "RT PF00078 RVT_1                       retro-/caulimo-/hepadnavirus_and_Ty3_RT"
    "RT PF07727 RVT_2                       Ty1/copia_RT"
    "RT PF06817 RVT_thumb                   RT_thumb"
    "RT PF06815 RVT_connect                 RT_connection"
    # Rep — CRESS-DNA replication-associated protein (REP_FAMILIES).
    "Rep PF02407 Viral_Rep                  circo-/nanovirus-like_Rep_endonuclease"
    "Rep PF00910 RNA_helicase               SF3_helicase_of_Rep"
    "Rep PF00799 Gemini_AL1                 geminivirus-like_Rep_catalytic"
    "Rep PF08283 Gemini_AL1_M               geminivirus_Rep_central"
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
