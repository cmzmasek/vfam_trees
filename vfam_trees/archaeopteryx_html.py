"""Write a self-contained HTML page that opens a PhyloXML tree in Archaeopteryx.js.

The tree is inlined in the page (so it opens from a local file, where a
``fetch`` of a sibling file is refused); the viewer and its libraries load
from jsDelivr at pinned versions.  Archaeopteryx.js throws on an unknown
config key, so the version is pinned and the config kept minimal.

The page opens coloured by Genus (or by Species on a one-genus tree), as the
PDF/PNG figures are, when that field is one the viewer would accept without a
warning; otherwise the viewer's own choice stands.  Left to itself, the viewer
opens a small tree on its highest-entropy field, which there can be a raw date
string (any field with at most 20 values qualifies).
"""
from __future__ import annotations

import html
import json
import xml.etree.ElementTree as ET
from pathlib import Path

from . import __version__
from .logger import get_logger

log = get_logger(__name__)

ARCHAEOPTERYX_JS_VERSION = "3.28.0"

_CDN = "https://cdn.jsdelivr.net/npm"
_SCRIPTS = [
    f"{_CDN}/d3@7.9.0/dist/d3.min.js",
    f"{_CDN}/sax@1.6.1/lib/sax.js",
    # Optional: the PDF entry in the viewer's Download menu.
    f"{_CDN}/jspdf@4.2.1/dist/jspdf.umd.min.js",
    f"{_CDN}/svg2pdf.js@2.8.1/dist/svg2pdf.umd.min.js",
    f"{_CDN}/phyloxml@1.1.7/phyloxml.js",
    f"{_CDN}/archaeopteryx@{ARCHAEOPTERYX_JS_VERSION}/forester.js",
    f"{_CDN}/archaeopteryx@{ARCHAEOPTERYX_JS_VERSION}/archaeopteryx.js",
]

_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="generator" content="{generator}">
<title>{title}</title>
{scripts}
<style>
  html, body {{ margin: 0; height: 100%; background: #fff;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }}
  #phylogram1 {{ position: absolute; top: 0; bottom: 0; left: 0; right: 0; }}
  #aptx-error {{ margin: 16px; color: #b00020; white-space: pre-wrap; }}
</style>
</head>
<body>
<div id="phylogram1"></div>
<div id="aptx-error"></div>
<script>
var TREE_FILE = {file_name_js};
var TREE_DATA = {tree_data_js};
window.addEventListener('load', function () {{
  var err = document.getElementById('aptx-error');
  if (typeof archaeopteryx === 'undefined') {{
    err.textContent = 'Could not load Archaeopteryx.js from cdn.jsdelivr.net. '
      + 'This page needs an internet connection to show the tree; the tree itself '
      + 'is in ' + TREE_FILE + ' next to this page.';
    return;
  }}
  // The view rides in the URL hash, so a copied link reopens the tree as left.
  var config = {{ enableDynamicSizing: true }};
  var OPEN_WITH = {open_with_js};
  if (OPEN_WITH) {{ config.initialVisualization = OPEN_WITH; }}
  var view = archaeopteryx.decodeViewState(window.location.hash);
  if (view) {{ config.view = view; }}
  config.onViewChange = function (state, encoded) {{
    history.replaceState(null, '', window.location.pathname + window.location.search
      + (encoded ? '#' + encoded : ''));
  }};
  try {{
    archaeopteryx.launchArchaeopteryx('#phylogram1', TREE_FILE, TREE_DATA, config);
  }} catch (e) {{
    err.textContent = 'Error: ' + (e && e.message ? e.message : e);
  }}
}});
</script>
</body>
</html>
"""


# Archaeopteryx.js 3.28.0 accepts a forced initialVisualization silently only
# when the field has at most this many values (LEGEND_MAX_ROWS) and is on at
# least two thirds of the tips.
_MAX_OPENING_VALUES = 20
_OPENING_FIELDS = ("vipr:Genus", "vipr:Species")


def _opening_visualization(tree_text: str) -> str | None:
    """The Color-by id to open with ("prop:vipr:Genus"), or None for the viewer's own."""
    ns = "{http://www.phyloxml.org}"
    try:
        root = ET.fromstring(tree_text.encode("utf-8"))
    except ET.ParseError:
        return None
    tips = [c for c in root.iter(f"{ns}clade") if c.find(f"{ns}clade") is None]
    for ref in _OPENING_FIELDS:
        values = [p.text.strip() for c in tips for p in c.findall(f"{ns}property")
                  if p.get("ref") == ref and p.text and p.text.strip()]
        distinct = {v.lower() for v in values}
        if 2 <= len(distinct) <= _MAX_OPENING_VALUES and len(values) * 3 >= len(tips) * 2:
            return f"prop:{ref}"
    return None


def _js_string(text: str) -> str:
    """A JavaScript string literal safe to place inside a <script> element."""
    # json.dumps gives a valid JS literal (ensure_ascii also escapes U+2028/9);
    # "</" is escaped so the XML can never close the script element.
    return json.dumps(text, ensure_ascii=True).replace("</", "<\\/")


def write_archaeopteryx_html(xml_path: Path, html_path: Path, title: str) -> None:
    """Write ``html_path``: a page showing the PhyloXML at ``xml_path``."""
    tree_text = Path(xml_path).read_text(encoding="utf-8")
    scripts = "\n".join(f'<script src="{src}"></script>' for src in _SCRIPTS)
    page = _TEMPLATE.format(
        title=html.escape(title),
        scripts=scripts,
        generator=f"vfam_trees v{__version__}",
        file_name_js=_js_string(Path(xml_path).name),
        tree_data_js=_js_string(tree_text),
        open_with_js=json.dumps(_opening_visualization(tree_text)),
    )
    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(page, encoding="utf-8")
    log.debug("Archaeopteryx.js page written to %s", html_path)
