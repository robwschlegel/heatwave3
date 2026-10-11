# heatwave3

Detect marine heatwaves (MHWs) and marine cold-spells (MCSs) directly in
gridded NetCDF data, following the Hobday et al. (2016, 2018) definitions.

`heatwave3` is the Python interface to the same C++ core as the
[heatwave3 R package](https://robwschlegel.github.io/heatwave3/): climatology,
event detection, categorisation and blob detection run in multi-threaded C++
with direct NetCDF I/O.

This package is in early development (alpha). The public API is not yet
available.
