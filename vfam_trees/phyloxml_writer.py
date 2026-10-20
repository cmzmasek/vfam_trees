"""Generate PhyloXML output with vipr: metadata properties."""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path
from xml.dom import minidom

from datetime import datetime, timezone

from Bio import Phylo

from . import __version__
from .leaf_fields import (
    clean_value, decimal_year, leaf_species, location_fields, own_taxon, ranked_name,
)
from .logger import get_logger
from .taxonomy import _infer_rank

log = get_logger(__name__)

PHYLOXML_NS = "http://www.phyloxml.org"
VIPR_PREFIX = "vipr"

# Branch lengths of every vfam_trees tree are substitutions per site.  The
# token is the one Archaeopteryx and Archaeopteryx.js use for a divergence
# tree: it labels their scale bar and clock-plot rate, and keeps a tree with
# tip dates from being drawn on a calendar axis.
BRANCH_LENGTH_UNIT = "subs/site"

# Accession sources other than NCBI, keyed by metadata "source": URL template
# for the record page (None: no public page).
_ACCESSION_URLS = {
    "pathoplexus": "https://pathoplexus.org/seq/{acc}",
    "user": None,
}

# Valid rank values per the PhyloXML schema enumeration.
_PHYLOXML_VALID_RANKS = frozenset([
    "domain", "superkingdom", "kingdom", "subkingdom", "branch",
    "infrakingdom", "superphylum", "phylum", "subphylum", "infraphylum",
    "microphylum", "superdivision", "division", "subdivision", "infradivision",
    "superclass", "class", "subclass", "infraclass", "superlegion", "legion",
    "sublegion", "infralegion", "supercohort", "cohort", "subcohort",
    "infracohort", "magnorder", "superorder", "order", "suborder",
    "infraorder", "superfamily", "family", "subfamily", "supertribe", "tribe",
    "subtribe", "infratribe", "genus", "subgenus", "superspecies", "species",
    "subspecies", "variety", "varietas", "subvariety", "form", "subform",
    "cultivar", "strain", "section", "subsection", "unknown", "other",
])

# Fields written as vipr: properties (species and taxon_id are handled via
# <taxonomy>; location via vipr:Country / vipr:Location below)
METADATA_FIELDS = [
    ("host",            "xsd:string", "Host"),
    ("collection_date", "xsd:string", "Collection_Date"),
    ("strain",          "xsd:string", "Strain"),
    # outbreak tag for externally-sourced tips (additional_data_sources); absent
    # on ordinary records, so the empty-value guard below skips it for them.
    ("outbreak",        "xsd:string", "outbreak"),
]


_YEAR_RE = re.compile(r'\b(1[89]\d{2}|20\d{2})\b')


def _extract_year(date_str: str) -> str:
    """Return the 4-digit year from a collection_date string, or ''."""
    if not date_str:
        return ""
    m = _YEAR_RE.search(date_str)
    return m.group(1) if m else ""


def write_phyloxml(
    newick_path: Path,
    output_xml: Path,
    id_map: dict[str, str],
    leaf_metadata: dict[str, dict],
    family: str,
    tree=None,
    phylogeny_name: str | None = None,
    phylogeny_detail: str | None = None,
    confidence_type: str = "SH_aLRT",
    leaf_colors: dict[str, str] | None = None,
    aligned_seqs: dict[str, str] | None = None,
    seq_type: str = "",
) -> None:
    """Generate a PhyloXML file from a Newick tree with vipr: metadata properties.

    Args:
        newick_path: annotated Newick file (short IDs, internal nodes labelled)
        output_xml: output PhyloXML path
        id_map: short_id → display_name
        leaf_metadata: short_id → metadata dict
        family: viral family name (used in phylogeny description)
        tree: pre-parsed BioPython tree (avoids lossy Newick roundtrip if supplied)
        aligned_seqs: short_id → aligned sequence string; when provided, each leaf
            gets a <mol_seq is_aligned="true"> child inside its <sequence> element
        seq_type: "protein" or "nucleotide"; sets the type= attribute on <sequence>
    """
    if tree is None:
        # Fall back to parsing from file — bootstrap on internal nodes may be
        # merged into the name label by BioPython's Newick writer
        tree = list(Phylo.parse(str(newick_path), "newick"))[0]
        # Normalize bootstrap values to 0–100 if needed (e.g. FastTree outputs 0–1)
        confidences = [c.confidence for c in tree.find_clades() if c.confidence is not None]
        if confidences and max(confidences) <= 1.0:
            for clade in tree.find_clades():
                if clade.confidence is not None:
                    clade.confidence = round(clade.confidence * 100)
    tree.ladderize(reverse=True)
    tree.name = family

    root = ET.Element("phyloxml", xmlns=PHYLOXML_NS)
    # Attribute order follows the Archaeopteryx writer: rooted,
    # branch_length_unit, type, rerootable.
    phylogeny_el = ET.SubElement(root, "phylogeny", rooted="true",
                                 branch_length_unit=BRANCH_LENGTH_UNIT,
                                 rerootable="false")
    name_el = ET.SubElement(phylogeny_el, "name")
    name_el.text = phylogeny_name or family
    desc_el = ET.SubElement(phylogeny_el, "description")
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    desc_el.text = (
        f"vfam_trees v{__version__} phylogenetic tree for {family} "
        f"— generated {timestamp}"
        + (f" | {phylogeny_detail}" if phylogeny_detail else "")
    )

    _write_clade(phylogeny_el, tree.root, id_map, leaf_metadata, confidence_type,
                 leaf_colors or {}, aligned_seqs, seq_type)

    output_xml.parent.mkdir(parents=True, exist_ok=True)
    _write_pretty_xml(root, output_xml)
    log.debug("PhyloXML written to %s", output_xml)


def _write_clade(
    parent_el: ET.Element,
    clade,
    id_map: dict[str, str],
    leaf_metadata: dict[str, dict],
    confidence_type: str = "SH_aLRT",
    leaf_colors: dict[str, str] | None = None,
    aligned_seqs: dict[str, str] | None = None,
    seq_type: str = "",
) -> ET.Element:
    clade_el = ET.SubElement(parent_el, "clade")

    # PhyloXML schema requires this element order:
    # name, branch_length, confidence, width, color, ..., taxonomy, property, clade

    is_leaf = clade.is_terminal()

    if is_leaf:
        short_id = clade.name or ""
        display_name = id_map.get(short_id, short_id)
        name_el = ET.SubElement(clade_el, "name")
        name_el.text = display_name
    else:
        # Internal node — use label as taxon annotation if present,
        # but suppress it entirely for same-species pairs (redundant label).
        label = clade.name or ""
        if label and not _is_same_species_pair(clade, leaf_metadata):
            name_el = ET.SubElement(clade_el, "name")
            name_el.text = label

    # branch_length and confidence must follow name
    if clade.branch_length is not None:
        bl_el = ET.SubElement(clade_el, "branch_length")
        bl_el.text = str(clade.branch_length)

    if clade.confidence is not None:
        conf_el = ET.SubElement(clade_el, "confidence", type=confidence_type)
        conf_el.text = str(clade.confidence)

    # taxonomy + properties follow confidence
    if is_leaf:
        short_id = clade.name or ""
        meta = leaf_metadata.get(short_id, {})
        # <taxonomy> is the leaf's own NCBI taxon: its id, NCBI's current name
        # for that id (GenBank organism strings can be stale, splitting one
        # taxon into several names), and its rank when NCBI gives one the
        # phyloXML schema allows.  Without a ranked lineage, the organism.
        species = clean_value(meta.get("species", ""))
        taxon_id = meta.get("taxon_id", "")
        own = own_taxon(meta)
        sci_name = own.get("name") or species
        own_rank = (own.get("rank") or "").lower()
        if sci_name or taxon_id:
            taxonomy_el = ET.SubElement(clade_el, "taxonomy")
            if taxon_id:
                id_el = ET.SubElement(taxonomy_el, "id", provider="ncbi")
                id_el.text = str(taxon_id)
            if sci_name:
                sci_el = ET.SubElement(taxonomy_el, "scientific_name")
                sci_el.text = sci_name
            if own_rank in _PHYLOXML_VALID_RANKS and own_rank not in ("other", "unknown"):
                rank_el = ET.SubElement(taxonomy_el, "rank")
                rank_el.text = own_rank
        # <sequence> element with accession, title, and (optionally) mol_seq
        accession = meta.get("accession", "")
        seq_name = meta.get("seq_name", "")
        aligned_seq = (aligned_seqs or {}).get(short_id, "")
        if accession or seq_name or aligned_seq:
            mol_type = "protein" if seq_type == "protein" else ("dna" if seq_type else "")
            seq_attrs = {"type": mol_type} if mol_type else {}
            seq_el = ET.SubElement(clade_el, "sequence", **seq_attrs)
            source = (meta.get("source") or "ncbi").lower()
            if accession:
                acc_el = ET.SubElement(seq_el, "accession", source=source)
                acc_el.text = accession
            if seq_name:
                seqname_el = ET.SubElement(seq_el, "name")
                seqname_el.text = seq_name
            if aligned_seq:
                molseq_el = ET.SubElement(seq_el, "mol_seq", is_aligned="true")
                molseq_el.text = aligned_seq
            url_template = _ACCESSION_URLS.get(source)
            if accession and url_template:
                uri_el = ET.SubElement(seq_el, "uri", desc=f"{source} record")
                uri_el.text = url_template.format(acc=accession)

        # <date>: the collection date as a decimal year, which is where both
        # Archaeopteryx viewers read tip dates (root-to-tip clock plot).  The
        # raw GenBank string goes in <desc>; no minimum/maximum, which the
        # viewers would draw as age bars.
        raw_date = meta.get("collection_date", "")
        year_value = decimal_year(raw_date)
        if year_value is not None:
            date_el = ET.SubElement(clade_el, "date", unit="year")
            ET.SubElement(date_el, "desc").text = str(raw_date).strip()
            ET.SubElement(date_el, "value").text = f"{year_value:.5f}"

        # Remaining fields as vipr: properties
        for field, datatype, ref_name in METADATA_FIELDS:
            # Only real values: not "unknown" (extract_metadata's placeholder),
            # "n/a", or INSDC's "missing: ..." terms.
            value = clean_value(meta.get(field, ""))
            if value:
                _add_property(clade_el, f"{VIPR_PREFIX}:{ref_name}", datatype, value)

        # Country (the part of an INSDC geo_loc_name before ':') is what the
        # viewers group and colour by; the full Location only when it says more.
        country, location_detail = location_fields(meta)
        if country:
            _add_property(clade_el, "vipr:Country", "xsd:string", country)
        if location_detail:
            _add_property(clade_el, "vipr:Location", "xsd:string", location_detail)

        year = _extract_year(meta.get("collection_date", ""))
        if year:
            _add_property(clade_el, "vipr:Year", "xsd:string", year)

        # Taxonomic rank properties: the viewers' clades-by-rank read a rank
        # from a property's local name.  Species is the ICTV species from the
        # ranked lineage, matching the internal LCA labels.
        species_val = leaf_species(meta)
        if species_val:
            _add_property(clade_el, "vipr:Species", "xsd:string", species_val)
        for rank, ref in (("genus", "vipr:Genus"), ("subgenus", "vipr:Subgenus"),
                          ("subfamily", "vipr:Subfamily")):
            name = ranked_name(meta, rank)
            if name:
                _add_property(clade_el, ref, "xsd:string", name)

        # Genus-based font color for node label
        if leaf_colors:
            hex_color = leaf_colors.get(short_id, "")
            if hex_color:
                _add_property(clade_el, "style:font_color", "xsd:token", hex_color)
    else:
        label = clade.name or ""
        if label and not _is_same_species_pair(clade, leaf_metadata):
            # _taxonomy_rank is set by taxonomy.annotate_tree ('' = NCBI says
            # "no rank"); infer from the name only when it was never set.
            rank = getattr(clade, "_taxonomy_rank", None)
            if rank is None:
                rank = _infer_rank(label)
            taxonomy_el = ET.SubElement(clade_el, "taxonomy")
            taxid = getattr(clade, "_taxonomy_id", "")
            if taxid:
                id_el = ET.SubElement(taxonomy_el, "id", provider="ncbi")
                id_el.text = str(taxid)
            sci_el = ET.SubElement(taxonomy_el, "scientific_name")
            sci_el.text = label
            if rank:
                rank_el = ET.SubElement(taxonomy_el, "rank")
                rank_el.text = _valid_rank(rank)

    # child clades must come last
    for child in clade.clades:
        _write_clade(clade_el, child, id_map, leaf_metadata, confidence_type,
                     leaf_colors, aligned_seqs, seq_type)

    return clade_el


def _valid_rank(rank: str) -> str:
    """Return rank if valid per PhyloXML schema, else 'other'."""
    return rank if rank in _PHYLOXML_VALID_RANKS else "other"


def _is_same_species_pair(clade, leaf_metadata: dict[str, dict]) -> bool:
    """Return True if this node has exactly two leaf children of the same species.

    Such nodes get a redundant species-level LCA label that adds no information
    beyond what the two leaves already carry: the same vipr:Species (ICTV
    species), or the same organism, which their <taxonomy> carries.
    """
    children = clade.clades
    if len(children) != 2:
        return False
    if not all(c.is_terminal() for c in children):
        return False
    metas = [leaf_metadata.get(c.name, {}) for c in children]
    for field in (leaf_species, lambda m: clean_value(m.get("species", ""))):
        a, b = (field(m) for m in metas)
        if a and a == b:
            return True
    return False


def _add_property(
    parent: ET.Element,
    ref: str,
    datatype: str,
    value: str,
    applies_to: str = "node",
) -> ET.Element:
    prop = ET.SubElement(
        parent,
        "property",
        ref=ref,
        datatype=datatype,
        applies_to=applies_to,
    )
    prop.text = value
    return prop


def _write_pretty_xml(root: ET.Element, path: Path) -> None:
    raw = ET.tostring(root, encoding="unicode", xml_declaration=False)
    reparsed = minidom.parseString(raw)
    pretty = reparsed.toprettyxml(indent="  ", encoding="UTF-8")
    # minidom adds <?xml ...?> header
    with open(path, "wb") as f:
        f.write(pretty)
