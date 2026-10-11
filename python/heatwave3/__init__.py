"""Detect marine heatwaves and cold-spells directly in gridded NetCDF data.

Python interface to the heatwave3 C++ core, shared with the R package of the
same name.
"""

from importlib.metadata import version as _version

from . import _core

__version__ = _version("heatwave3")

__all__ = ["__version__"]
