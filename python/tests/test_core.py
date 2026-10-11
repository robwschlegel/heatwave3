"""Smoke tests for the compiled C++ core (heatwave3._core)."""

from pathlib import Path

import numpy as np
import pytest

import heatwave3
from heatwave3 import _core

SST = Path(__file__).resolve().parents[2] / "inst" / "extdata" / "sst_test.nc"


def test_version_matches_package_metadata():
    v = _core.hw3_version()
    assert v["package"] == "heatwave3"
    assert v["version"] == heatwave3.__version__
    assert v["cpp_standard"] == "C++17"


def test_threads_round_trip():
    old = _core.hw3_get_threads()
    try:
        assert _core.hw3_set_threads(2) == 2
        assert _core.hw3_get_threads() == 2
    finally:
        _core.hw3_set_threads(old)


@pytest.mark.skipif(not SST.exists(), reason="sst_test.nc not available")
def test_climatology_and_events_round_trip(tmp_path):
    clim = str(tmp_path / "x_clim.nc")
    events = str(tmp_path / "x_events.nc")
    _core.hw3_compute_clim(str(SST), clim, "", ["1982-01-01", "2011-12-31"])
    cd = _core.hw3_read_clim_nc(clim)
    assert cd["seas"].shape == (cd["nlon"] * cd["nlat"] * 366,)
    assert np.isfinite(cd["thresh"]).any()

    _core.hw3_detect_events(str(SST), clim, events, "", category=True)
    ed = _core.hw3_read_event_nc(events)
    assert ed["nevents"] > 0
    assert (ed["duration"] >= 5).all()
    assert ((ed["category"] >= 1) & (ed["category"] <= 4)).all()


def test_core_errors_become_python_exceptions(tmp_path):
    with pytest.raises(RuntimeError, match="NetCDF error"):
        _core.hw3_read_clim_nc(str(tmp_path / "missing.nc"))
