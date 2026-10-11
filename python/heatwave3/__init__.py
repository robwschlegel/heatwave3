"""Detect marine heatwaves and cold-spells directly in gridded NetCDF data.

Python interface to the heatwave3 C++ core, shared with the R package of the
same name. Function names, arguments and outputs match the R package.
"""

from importlib.metadata import version as _version

from . import _core
from .blob import detect_blob3
from .category import block_average3, category3, category_daily3
from .climatology import ts2clm3
from .events import detect3, detect_event3, exceedance3
from .export import hw3_export
from .plot import event_line3, geom_flame3, geom_lolli3, plot_metric3
from .threads import getHW3threads, setHW3threads

__version__ = _version("heatwave3")

__all__ = [
    "__version__",
    "block_average3",
    "category3",
    "category_daily3",
    "detect3",
    "detect_blob3",
    "detect_event3",
    "event_line3",
    "exceedance3",
    "geom_flame3",
    "geom_lolli3",
    "getHW3threads",
    "hw3_export",
    "plot_metric3",
    "setHW3threads",
    "ts2clm3",
]
