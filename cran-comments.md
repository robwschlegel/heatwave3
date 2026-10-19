## Submission

This is a new submission.

heatwave3 detects marine heatwaves and cold-spells on gridded 'NetCDF' data
(Hobday et al. 2016, 2018). The algorithms are implemented in C++ with direct
calls to the netCDF-C library, which is a system requirement.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

The words flagged as possibly misspelled in DESCRIPTION are correct:
"Hobday" (author name), "MHWs" (marine heatwaves), "MCS" (marine cold-spells),
"gridded", "NetCDF", and the British spellings "categorisation" and
"visualisation".

## Test environments

* local macOS (arm64), R 4.6.1
* GitHub Actions: macOS, Ubuntu and Windows, R release
* win-builder: R-devel and R-release
* mac-builder: R-release

## Threads

The C++ core uses std::thread. The default is half of the available cores,
capped at 2 when _R_CHECK_LIMIT_CORES_ is set and at OMP_THREAD_LIMIT when that
is set; all examples, tests and vignettes use at most 2 threads.
