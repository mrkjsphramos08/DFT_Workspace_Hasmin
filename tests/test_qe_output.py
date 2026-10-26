"""Tests for the ``lco`` package.

Run with::

    python -m pytest tests/
"""
import os
from pathlib import Path

import pytest

from lco.qe_output import parse_scf_history, final_accuracy

REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_OUT = REPO_ROOT / "05_metals" / "Pt1" / "Adsorbed" / "3_layers" / "lili" / "OUTPUT_adsorbed_Pt1_latest_lili_new_template..out"


@pytest.mark.skipif(not EXAMPLE_OUT.exists(), reason="example QE output not present in repo")
def test_parse_scf_history_returns_dict():
    history = parse_scf_history(EXAMPLE_OUT)
    assert isinstance(history, dict)
    assert "converged" in history
    assert "iterations" in history
    assert "job_done" in history
    assert "final_energy" in history


@pytest.mark.skipif(not EXAMPLE_OUT.exists(), reason="example QE output not present in repo")
def test_final_accuracy_returns_number_or_none():
    history = parse_scf_history(EXAMPLE_OUT)
    acc = final_accuracy(history)
    assert acc is None or isinstance(acc, float)