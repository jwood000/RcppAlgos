# Sample the Cartesian Product

- Generate a specific (lexicographically) or random sample of the
  Cartesian product of the input vectors.

- Produce results in parallel using the `nThreads` arguments.

- GMP support allows for exploration where the number of results is
  large.

## Usage

``` r
expandGridSample(
    ..., n = NULL, sampleVec = NULL, seed = NULL,
    nThreads = NULL, namedSample = FALSE, return_df = FALSE
)
```

## Arguments

- ...:

  vectors, factors or a list containing these. (See
  [`?expand.grid`](https://rdrr.io/r/base/expand.grid.html)).

- n:

  Number of results to return. The default is `NULL`.

- sampleVec:

  A vector of numbers representing the lexicographical partition of
  groups to return. Accepts vectors of class `bigz` as well as vectors
  of characters

- seed:

  Random seed initialization. The default is `NULL`. N.B. If the gmp
  library is needed, this parameter must be set in order to have
  reproducible results (*E.g*
  [`set.seed()`](https://rdrr.io/r/base/Random.html) has no effect in
  these cases).

- nThreads:

  Specific number of threads to be used. The default is `NULL`.

- namedSample:

  Logical flag. If `TRUE`, `rownames` corresponding to the
  lexicographical result, will be added to the returned matrix. The
  default is `FALSE`.

- return_df:

  Logical flag to force the output to be a `data.frame`. The default is
  `FALSE`.

## Details

These algorithms rely on efficiently generating the \\n^{th}\\
lexicographical result.

## Value

When all of the input is of the same type, by default `expandGridSample`
produces a `matrix` (a `data.frame` otherwise). This can be ignored by
setting the argument `return_df = TRUE`.

## References

[Lexicographical
order](https://en.wikipedia.org/wiki/Lexicographical_order)

## Author

Joseph Wood

## Examples

``` r
## input vectors
lst = list(factor(state.abb), euro, islands)

## generate 10 random products
expandGridSample(lst, n = 10, seed = 123)
#>    Var1       Var2 Var3
#> 1    OK 1936.27000  306
#> 2    OK   40.33990  306
#> 3    WV    6.55957   13
#> 4    CO 1936.27000   73
#> 5    AR    6.55957   30
#> 6    WI 1936.27000   29
#> 7    CT    5.94573   25
#> 8    MN   13.76030   33
#> 9    GA   13.76030   84
#> 10   IL   40.33990  227

## using sampleVec to generate specific results
expandGridSample(lst, sampleVec = c(1, 100, 1e3))
#>   Var1     Var2  Var3
#> 1   AL 13.76030 11506
#> 2   AL  1.95583  2968
#> 3   AK  2.20371    16

all.equal(expandGridSample(lst, sampleVec = 1:expandGridCount(lst)),
         expandGrid(lst))
#> [1] TRUE

## Examples with enormous number of total results
big_lst = Map(function(x, y) x:y, 8:33, 15:40)
num = expandGridCount(big_lst)
gmp::log2.bigz(num)
#> [1] 78
## [1] 78

first = gmp::urand.bigz(n = 1, size = 78, seed = 123)
#> Seed initialisation
mySamp = do.call(c, lapply(0:10, function(x) gmp::add.bigz(first, x)))

class(mySamp)
#> [1] "bigz"
## [1] "bigz"

## using the sampling function
cartSamp = expandGridSample(big_lst, sampleVec = mySamp)

## using the standard function
cartGeneral = expandGrid(big_lst,
                         lower = first,
                         upper = gmp::add.bigz(first, 10))

identical(cartSamp, cartGeneral)
#> [1] TRUE
## [1] TRUE
```
