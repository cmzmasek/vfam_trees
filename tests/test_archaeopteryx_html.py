"""Tests for vfam_trees.archaeopteryx_html."""
import json
import re

from vfam_trees.archaeopteryx_html import (
    ARCHAEOPTERYX_JS_VERSION, _SCRIPTS, write_archaeopteryx_html,
)

_XML = ('<?xml version="1.0" encoding="UTF-8"?>\n<phyloxml xmlns="http://www.phyloxml.org">'
        '<phylogeny rooted="true"><clade><name>Côte d\'Ivoire </script> "x"</name>'
        '</clade></phylogeny></phyloxml>\n')


def _write(tmp_path, title="Ebola <tree> & co"):
    xml = tmp_path / "Ebola_tree_100.xml"
    xml.write_text(_XML, encoding="utf-8")
    out = tmp_path / "Ebola_tree_100.html"
    write_archaeopteryx_html(xml, out, title=title)
    return out.read_text(encoding="utf-8")


def test_tree_inlined_and_recoverable(tmp_path):
    page = _write(tmp_path)
    m = re.search(r"var TREE_DATA = (.*);\n", page)
    assert json.loads(m.group(1).replace("<\\/", "</")) == _XML


def test_xml_cannot_close_the_script(tmp_path):
    page = _write(tmp_path)
    assert page.count("</script>") == len(_SCRIPTS) + 1


def test_title_escaped(tmp_path):
    assert "<title>Ebola &lt;tree&gt; &amp; co</title>" in _write(tmp_path)


def test_pinned_versions_and_load_order(tmp_path):
    page = _write(tmp_path)
    srcs = re.findall(r'<script src="([^"]+)"', page)
    assert all("@" in s for s in srcs)
    names = [s.rsplit("/", 1)[1] for s in srcs]
    assert names.index("d3.min.js") < names.index("phyloxml.js") < names.index("forester.js") \
        < names.index("archaeopteryx.js")
    assert f"archaeopteryx@{ARCHAEOPTERYX_JS_VERSION}/archaeopteryx.js" in page


def test_launch_uses_file_name(tmp_path):
    page = _write(tmp_path)
    assert 'var TREE_FILE = "Ebola_tree_100.xml";' in page
    assert "archaeopteryx.launchArchaeopteryx('#phylogram1', TREE_FILE, TREE_DATA, config)" in page


def _tree_xml(genera, species=None, extra_tips=0):
    tips = []
    for i, g in enumerate(genera):
        props = f'<property ref="vipr:Genus" datatype="xsd:string" applies_to="node">{g}</property>'
        if species:
            props += (f'<property ref="vipr:Species" datatype="xsd:string" '
                      f'applies_to="node">{species[i]}</property>')
        tips.append(f"<clade><name>t{i}</name>{props}</clade>")
    tips += [f"<clade><name>u{i}</name></clade>" for i in range(extra_tips)]
    return ('<phyloxml xmlns="http://www.phyloxml.org"><phylogeny rooted="true"><clade>'
            + "".join(tips) + "</clade></phylogeny></phyloxml>")


def _open_with(tmp_path, text):
    xml = tmp_path / "t.xml"
    xml.write_text(text, encoding="utf-8")
    out = tmp_path / "t.html"
    write_archaeopteryx_html(xml, out, title="t")
    return json.loads(re.search(r"var OPEN_WITH = (.*);\n", out.read_text()).group(1))


def test_opens_with_genus(tmp_path):
    assert _open_with(tmp_path, _tree_xml(["A", "B", "A"])) == "prop:vipr:Genus"


def test_one_genus_opens_with_species(tmp_path):
    xml = _tree_xml(["A", "A", "A"], species=["s1", "s2", "s1"])
    assert _open_with(tmp_path, xml) == "prop:vipr:Species"


def test_too_many_genera_left_to_viewer(tmp_path):
    assert _open_with(tmp_path, _tree_xml([f"G{i}" for i in range(21)])) is None


def test_sparse_genus_left_to_viewer(tmp_path):
    # genus on 2 of 4 tips: under two thirds
    assert _open_with(tmp_path, _tree_xml(["A", "B"], extra_tips=2)) is None
