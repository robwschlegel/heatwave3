<div id="main" class="col-md-9" role="main">

# Get or set the number of threads used by heatwave3

<div class="ref-description section level2">

`getHW3threads()` returns the current thread count. `setHW3threads()`
sets it and returns the new value.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
getHW3threads()

setHW3threads(threads = 0L)
```

</div>

</div>

<div class="section level2">

## Arguments

-   threads:

    Integer. Number of threads to use. `0` resets to the default (50
    percent of cores).

</div>

<div class="section level2">

## Value

Integer: the current (or newly set) thread count.

</div>

<div class="section level2">

## Details

The default is 50 percent of available cores, which can be overridden by
the environment variable `R_HEATWAVE3_NUM_THREADS`. The per-function
`n_threads` parameter, when greater than zero, takes precedence over the
package-level setting for that call.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
getHW3threads()
#> [1] 5
old <- setHW3threads(2)
getHW3threads()
#> [1] 2
setHW3threads(0)  # reset to default
#> [1] 5
```

</div>

</div>

</div>
