# Max Number of Concurrent Threads

Wrapper of
[std::thread::hardware_concurrency()](https://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency.html).
As stated by [cppreference](https://en.cppreference.com/w/), the
returned value should be considered only a hint.

## Usage

``` r
stdThreadMax()
```

## Value

An integer representing the number of concurrent threads supported by
the user implementation. If the value cannot be determined, `1L` is
returned.

## See also

[`detectCores`](https://rdrr.io/r/parallel/detectCores.html)

## Examples

``` r
stdThreadMax()
#> [1] 4
```
