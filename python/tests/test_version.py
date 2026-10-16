"""The Python and R packages are pinned to the same version, because both have
the same features: version in pyproject.toml equals Version in DESCRIPTION."""

import re
from importlib.metadata import version
from pathlib import Path

import pytest

import heatwave3 as hw3

DESCRIPTION = Path(__file__).resolve().parents[2] / "DESCRIPTION"


@pytest.mark.skipif(not DESCRIPTION.exists(), reason="DESCRIPTION not available")
def test_python_version_is_pinned_to_the_r_package():
    r_version = re.search(r"^Version:\s*(\S+)", DESCRIPTION.read_text(), re.MULTILINE).group(1)
    assert version("heatwave3") == r_version
    assert hw3.__version__ == r_version
    assert hw3._core.hw3_version()["version"] == r_version  # compiled into the C++ core
