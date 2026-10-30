"""API parity: each ported function has the R function's arguments, in the
same order, with the same defaults. R signatures are dumped by make_golden.R."""

import inspect
import json
from pathlib import Path

import pytest

import heatwave3 as hw3

SIGS = Path(__file__).resolve().parent / "data" / "r_signatures.json"

# Functions not yet ported (later phases).
NOT_YET = set()

# ggplot2 layers in R, matplotlib functions in Python: Python-native signatures.
PYTHON_NATIVE = {"geom_flame3", "geom_lolli3"}

# Python-only trailing arguments (the matplotlib Axes to draw on; Zarr output).
PYTHON_ONLY = {"event_line3": {"ax"}, "plot_metric3": {"ax"}, "ts2clm3": {"zarr_out"},
               "detect_event3": {"zarr_out"}, "detect3": {"zarr_out"}}

# R argument names that are Python keywords: R name -> Python name.
RENAMED = {("detect_blob3", "return"): "return_"}

# R defaults that cannot be literal in Python, and their documented Python form.
EXCEPTIONS = {
    ("detect_event3", "minDuration2"): ("minDuration", None),
    ("detect_event3", "maxGap2"): ("maxGap", None),
    ("detect3", "minDuration2"): ("minDuration", None),
    ("detect3", "maxGap2"): ("maxGap", None),
    ("detect_blob3", "return_"): ('c("event", "daily")', ("event", "daily")),
    ("detect_event_array3", "outputs"): ('c("event", "intensity", "category")',
                                         ("event", "intensity", "category")),
}


def r_default_to_python(text):
    """Translate a deparsed R literal default to its Python equivalent."""
    if text == "<required>":
        return inspect.Parameter.empty
    simple = {"NULL": None, "TRUE": True, "FALSE": False}
    if text in simple:
        return simple[text]
    if text.startswith("c(") and text.endswith(")"):  # match.arg choices: first is default
        return json.loads(text[2:-1].split(",")[0].strip())
    if text.startswith('"'):
        return json.loads(text)
    num = text.rstrip("L")
    return int(num) if num.lstrip("-").isdigit() else float(num)


pytestmark = pytest.mark.skipif(not SIGS.exists(), reason="R signatures not available")
R_SIGS = json.loads(SIGS.read_text()) if SIGS.exists() else {}


def test_every_r_export_is_ported_or_pending():
    assert set(R_SIGS) - set(hw3.__all__) == NOT_YET


@pytest.mark.parametrize("name", sorted(set(R_SIGS) - NOT_YET - PYTHON_NATIVE))
def test_signature_matches_r(name):
    r_args = {RENAMED.get((name, k), k): v for k, v in (R_SIGS[name] or {}).items()}
    # R's ... corresponds to Python's **kwargs
    params = {("..." if p.kind is p.VAR_KEYWORD else k): p
              for k, p in inspect.signature(getattr(hw3, name)).parameters.items()
              if k not in PYTHON_ONLY.get(name, ())}
    assert list(params) == list(r_args), f"{name}: argument names or order differ"
    for arg, r_text in r_args.items():
        py_default = params[arg].default
        if (name, arg) in EXCEPTIONS:
            assert (r_text, py_default) == EXCEPTIONS[(name, arg)], f"{name}({arg})"
        else:
            assert py_default == r_default_to_python(r_text), \
                f"{name}({arg}): R default {r_text!r}, Python {py_default!r}"
