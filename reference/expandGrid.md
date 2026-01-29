# Cartesian Product

Generate the [Cartesian
Product](https://en.wikipedia.org/wiki/Cartesian_product) of the input
vectors. It is very similar to
[`expand.grid`](https://rdrr.io/r/base/expand.grid.html) however there
are some notable differences:

- Produces lexicographic ordered output consistent with other functions
  in `RcppAlgos`. Compared to `expand.grid` where the first column
  varies the fastest, `expandGrid` varies the first column the slowest.

- When all of the input is of the same type, by default `expandGrid`
  produce a `matrix` (a `data.frame` otherwise). This can be ignored by
  setting the argument `return_df = TRUE`.

- No attributes are added nor are strings converted to factors. In
  `expand.grid` we would achieve this by setting
  `KEEP.OUT.ATTRS = FALSE` and `stringsAsFactors = FALSE`.

- If it is possible to return a matrix, we can utilize the argument
  `nThreads` in order to produce results in parallel for maximal
  efficiency.

## Usage

``` r
expandGrid(..., lower = NULL, upper = NULL,
               nThreads = NULL, return_df = FALSE)
```

## Arguments

- ...:

  vectors, factors or a list containing these. (See
  [`?expand.grid`](https://rdrr.io/r/base/expand.grid.html)).

- lower:

  The lower bound. Cartesian products are generated lexicographically,
  thus utilizing this argument will determine which specific product to
  start generating from (*e.g.* `expandGrid(1:5, 3:11, lower = 6)` is
  equivalent to
  `expandGrid(1:5, 3:11)[6:expandGridCount(1:5, 3:11), ]`). This
  argument along with `upper` is very useful for generating products in
  chunks allowing for easy parallelization.

- upper:

  The upper bound. Similar to `lower`, however this parameter allows the
  user to *stop* generation at a specific product (*e.g.*
  `expandGrid(1:5, 3:11, upper = 5)` is equivalent to
  `expandGrid(1:5, 3:11)[1:5, ]`)

- nThreads:

  Specific number of threads to be used. The default is `NULL`.

- return_df:

  Logical flag to force the output to be a `data.frame`. The default is
  `FALSE`.

## Value

When all of the input is of the same type, by default `expandGrid`
produces a `matrix` (a `data.frame` otherwise). This can be ignored by
setting the argument `return_df = TRUE`.

## Author

Joseph Wood

## See also

[`comboGrid`](https://jwood000.github.io/RcppAlgos/reference/comboGrid.md)

## Examples

``` r
## description example
lst = list(1:2, 1:2)

## Using base R
t = expand.grid(lst)

## vs using expandGrid. N.B. Output is a matrix
expandGrid(lst)
#>      Var1 Var2
#> [1,]    1    1
#> [2,]    1    2
#> [3,]    2    1
#> [4,]    2    2

## Force a data.frame to be returned
expandGrid(lst, return_df = TRUE)
#>   Var1 Var2
#> 1    1    1
#> 2    1    2
#> 3    2    1
#> 4    2    2

lst = Map(function(x, y) x:y, 8:14, 15:21)

## Use multiple threads for greater efficiency
system.time(expandGrid(lst, nThreads = 2))
#>    user  system elapsed 
#>   0.016   0.006   0.011 
```
