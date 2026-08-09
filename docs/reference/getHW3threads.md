# Get or set the number of threads used by heatwave3

`getHW3threads()` returns the current thread count. `setHW3threads()`
sets it and returns the new value.

## Usage

``` r
getHW3threads()

setHW3threads(threads = 0L)
```

## Arguments

- threads:

  Integer. Number of threads to use. `0` resets to the default (50
  percent of cores).

## Value

Integer: the current (or newly set) thread count.

## Details

The default is 50 percent of available cores, which can be overridden by
the environment variable `R_HEATWAVE3_NUM_THREADS`. The per-function
`n_threads` parameter, when greater than zero, takes precedence over the
package-level setting for that call.

## Examples

``` r
getHW3threads()
#> [1] 4
old <- setHW3threads(2)
getHW3threads()
#> [1] 2
setHW3threads(0)  # reset to default
#> [1] 4
```
