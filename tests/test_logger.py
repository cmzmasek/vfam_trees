"""Tests for vfam_trees.logger.

The pipeline relies on Python's hierarchical-logger propagation so that
module-level loggers (``vfam_trees.fetch``, ``vfam_trees.concat``,
``vfam_trees.pipeline_concat``, …) write to the per-family log file.  A
silent regression here would cause the entire concat-mode log stream to
disappear from ``<family>.log`` — see v1.2.4 for the original incident.
"""
from __future__ import annotations

import logging

from vfam_trees.logger import get_logger, setup_logger


def test_module_level_logger_propagates_to_parent_file_handler(tmp_path):
    """A log from `vfam_trees.<module>` lands in the file attached to
    the `vfam_trees` parent logger."""
    log_file = tmp_path / "family.log"
    setup_logger("vfam_trees", log_file=log_file, level="INFO")
    try:
        get_logger("vfam_trees.fetch").info("FETCH-MARKER-MSG-1")
        get_logger("vfam_trees.concat").info("CONCAT-MARKER-MSG-2")
        get_logger("vfam_trees.pipeline_concat").info("PIPELINE-MARKER-MSG-3")
        # Family-named child also propagates
        get_logger("vfam_trees.Adenoviridae").info("FAMILY-MARKER-MSG-4")
        # Force the file handler to flush
        for h in logging.getLogger("vfam_trees").handlers:
            h.flush()
        log_text = log_file.read_text()
        assert "FETCH-MARKER-MSG-1" in log_text
        assert "CONCAT-MARKER-MSG-2" in log_text
        assert "PIPELINE-MARKER-MSG-3" in log_text
        assert "FAMILY-MARKER-MSG-4" in log_text
    finally:
        # Detach handlers so they don't leak to other tests.
        parent = logging.getLogger("vfam_trees")
        for h in list(parent.handlers):
            h.close()
            parent.removeHandler(h)


def test_debug_level_propagates_to_file_when_level_is_debug(tmp_path):
    """File handler is configured at DEBUG, but the parent logger's level
    gates DEBUG records — so DEBUG only reaches the file when the level
    is set to DEBUG."""
    log_file = tmp_path / "family.log"
    setup_logger("vfam_trees", log_file=log_file, level="DEBUG")
    try:
        get_logger("vfam_trees.fetch").debug("DEBUG-FETCH-MARKER")
        for h in logging.getLogger("vfam_trees").handlers:
            h.flush()
        assert "DEBUG-FETCH-MARKER" in log_file.read_text()
    finally:
        parent = logging.getLogger("vfam_trees")
        for h in list(parent.handlers):
            h.close()
            parent.removeHandler(h)


def test_setup_clears_old_handlers_on_reconfigure(tmp_path):
    """Re-running setup_logger swaps the file handler — useful when one
    process runs the pipeline for several families sequentially."""
    log_file_a = tmp_path / "a.log"
    log_file_b = tmp_path / "b.log"
    setup_logger("vfam_trees", log_file=log_file_a, level="INFO")
    get_logger("vfam_trees.fetch").info("MSG-FOR-A")
    setup_logger("vfam_trees", log_file=log_file_b, level="INFO")
    get_logger("vfam_trees.fetch").info("MSG-FOR-B")
    try:
        for h in logging.getLogger("vfam_trees").handlers:
            h.flush()
        assert "MSG-FOR-A" in log_file_a.read_text()
        # After reconfigure, MSG-FOR-B goes only to b.log
        assert "MSG-FOR-B" not in log_file_a.read_text()
        assert "MSG-FOR-B" in log_file_b.read_text()
    finally:
        parent = logging.getLogger("vfam_trees")
        for h in list(parent.handlers):
            h.close()
            parent.removeHandler(h)
