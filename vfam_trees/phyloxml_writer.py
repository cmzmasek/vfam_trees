"""Generate PhyloXML output with vipr: metadata properties."""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path
from xml.dom import minidom

from datetime import datetime, timezone

from Bio import Phylo

from . import __version__
from .logger import get_logger
from .taxonomy import _infer_rank

log = get_logger(__name__)

PHYLOXML_NS = "http://www.phyloxml.org"
VIPR_PREFIX = "vipr"

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

# Fields written as vipr: properties (species and taxon_id are handled via <taxonomy>)
METADATA_FIELDS = [
    ("host",            "xsd:string", "Host"),
    ("collection_date", "xsd:string", "Collection_Date"),
    ("location",        "xsd:string", "Location"),
    ("strain",          "xsd:string", "Strain"),
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
) -> None:
    """Generate a PhyloXML file from a Newick tree with vipr: metadata properties.

    Args:
        newick_path: annotated Newick file (short IDs, internal nodes labelled)
        output_xml: output PhyloXML path
        id_map: short_id → display_name
        leaf_metadata: short_id → metadata dict
        family: viral family name (used in phylogeny description)
        tree: pre-parsed BioPython tree (avoids lossy Newick roundtrip if supplied)
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
    phylogeny_el = ET.SubElement(root, "phylogeny", rooted="true", rerootable="false")
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
                 leaf_colors or {})

    output_xml.parent.mkdir(parents=True, exist_ok=True)
    _write_pretty_xml(root, output_xml)
    log.info("PhyloXML written to %s", output_xml)


def _write_clade(
    parent_el: ET.Element,
    clade,
    id_map: dict[str, str],
    leaf_metadata: dict[str, dict],
    confidence_type: str = "SH_aLRT",
    leaf_colors: dict[str, str] | None = None,
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
        # Species and taxon_id go into the proper <taxonomy> element.
        # Treat the "unknown" placeholder (emitted by extract_metadata) as absent.
        species = meta.get("species", "")
        if species == "unknown":
            species = ""
        taxon_id = meta.get("taxon_id", "")
        if species or taxon_id:
            taxonomy_el = ET.SubElement(clade_el, "taxonomy")
            if taxon_id:
                id_el = ET.SubElement(taxonomy_el, "id", provider="ncbi")
                id_el.text = str(taxon_id)
            if species:
                sci_el = ET.SubElement(taxonomy_el, "scientific_name")
                sci_el.text = species
        # <sequence> element with accession and title
        accession = meta.get("accession", "")
        seq_name = meta.get("seq_name", "")
        if accession or seq_name:
            seq_el = ET.SubElement(clade_el, "sequence")
            if accession:
                acc_el = ET.SubElement(seq_el, "accession", source="ncbi")
                acc_el.text = accession
            if seq_name:
                seqname_el = ET.SubElement(seq_el, "name")
                seqname_el.text = seq_name

        # Remaining fields as vipr: properties
        for field, datatype, ref_name in METADATA_FIELDS:
            value = meta.get(field, "")
            # Skip the placeholder "unknown" emitted by extract_metadata — emit the
            # property only when the underlying GenBank field carried a real value.
            if value and value != "unknown":
                _add_property(clade_el, f"{VIPR_PREFIX}:{ref_name}", datatype, str(value))

        year = _extract_year(meta.get("collection_date", ""))
        if year:
            _add_property(clade_el, "vipr:Year", "xsd:string", year)

        # Taxonomic rank properties for downstream colorization
        species_val = meta.get("species", "")
        if species_val and species_val != "unknown":
            _add_property(clade_el, "vipr:Species", "xsd:string", species_val)
        lineage = meta.get("lineage_ranked", [])
        _taxa: dict[str, str] = {"genus": "", "subgenus": "", "subfamily": ""}
        for _entry in lineage:
            if isinstance(_entry, dict):
                _rank = (_entry.get("rank") or "").lower()
                if _rank in _taxa:
                    _taxa[_rank] = _entry.get("name", "")
        if _taxa["genus"]:
            _add_property(clade_el, "vipr:Genus", "xsd:string", _taxa["genus"])
        if _taxa["subgenus"]:
            _add_property(clade_el, "vipr:Subgenus", "xsd:string", _taxa["subgenus"])
        if _taxa["subfamily"]:
            _add_property(clade_el, "vipr:Subfamily", "xsd:string", _taxa["subfamily"])

        # Genus-based font color for node label
        if leaf_colors:
            hex_color = leaf_colors.get(short_id, "")
            if hex_color:
                _add_property(clade_el, "style:font_color", "xsd:token", hex_color)
    else:
        label = clade.name or ""
        if label and not _is_same_species_pair(clade, leaf_metadata):
            rank = getattr(clade, "_taxonomy_rank", None) or _infer_rank(label)
            taxonomy_el = ET.SubElement(clade_el, "taxonomy")
            sci_el = ET.SubElement(taxonomy_el, "scientific_name")
            sci_el.text = label
            if rank:
                rank_el = ET.SubElement(taxonomy_el, "rank")
                rank_el.text = _valid_rank(rank)

    # child clades must come last
    for child in clade.clades:
        _write_clade(clade_el, child, id_map, leaf_metadata, confidence_type, leaf_colors)

    return clade_el


def _valid_rank(rank: str) -> str:
    """Return rank if valid per PhyloXML schema, else 'other'."""
    return rank if rank in _PHYLOXML_VALID_RANKS else "other"


def _is_same_species_pair(clade, leaf_metadata: dict[str, dict]) -> bool:
    """Return True if this node has exactly two leaf children of the same species.

    Such nodes get a redundant species-level LCA label that adds no information
    beyond what the two leaves already carry.
    """
    children = clade.clades
    if len(children) != 2:
        return False
    if not all(c.is_terminal() for c in children):
        return False
    species = [
        leaf_metadata.get(c.name, {}).get("species", "")
        for c in children
    ]
    return bool(species[0]) and species[0] == species[1]


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
