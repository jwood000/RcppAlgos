# Unordered Cartesian Product

Efficient version of `expand.grid` where order does not matter. This is
a combinatorial variant where groups of elements are treated as
equivalent regardless of order. For example, given: `{1, 2}, {1, 2}`,
the unordered Cartesian product is `{1, 1}, {1, 2}, {2, 2}`. It is
loosely equivalent to the following:

- `t = expand.grid(lst)`

- `t = t[do.call(order, t), ]`

- `key = apply(t, 1, function(x) paste0(sort(x), collapse = ""))`

- `t[!duplicated(key), ]`

## Usage

``` r
comboGrid(..., repetition = TRUE, return_df = FALSE)
```

## Arguments

- ...:

  vectors, factors or a list containing these. (See
  [`?expand.grid`](https://rdrr.io/r/base/expand.grid.html)).

- repetition:

  Logical value indicating whether results should be with or without
  repetition. The default is `TRUE`.

- return_df:

  Logical flag to force the output to be a `data.frame`. The default is
  `FALSE`.

## Value

When all of the input is of the same type, by default `comboGrid`
produce a `matrix` (a `data.frame` otherwise). This can be ignored by
setting the argument `return_df = TRUE`.

## Author

Joseph Wood

## See also

[`expandGrid`](https://jwood000.github.io/RcppAlgos/reference/expandGrid.md)

## Examples

``` r
## description example
lst = list(1:2, 1:2)

t = expand.grid(lst)
t = t[do.call(order, t), ]
key = apply(t, 1, function(x) paste0(sort(x), collapse = ""))
t[!duplicated(key), ]
#>   Var1 Var2
#> 1    1    1
#> 3    1    2
#> 4    2    2

## vs using comboGrid. N.B. Output is a matrix
comboGrid(lst)
#>      Var1 Var2
#> [1,]    1    1
#> [2,]    1    2
#> [3,]    2    2

## Force a data.frame to be returned
comboGrid(lst, return_df = TRUE)
#>   Var1 Var2
#> 1    1    1
#> 2    1    2
#> 3    2    2

## Input vectors are of different type, so a data.frame is returned
expGridNoOrder = comboGrid(1:5, 3:9, letters[1:5], letters[c(1,4,5,8)])
head(expGridNoOrder)
#>   Var1 Var2 Var3 Var4
#> 1    1    3    a    a
#> 2    1    3    a    d
#> 3    1    3    a    e
#> 4    1    3    a    h
#> 5    1    3    b    a
#> 6    1    3    b    d
tail(expGridNoOrder)
#>     Var1 Var2 Var3 Var4
#> 539    5    9    c    h
#> 540    5    9    d    d
#> 541    5    9    d    e
#> 542    5    9    d    h
#> 543    5    9    e    e
#> 544    5    9    e    h

expGridNoOrderNoRep = comboGrid(1:5, 3:9, letters[1:5],
                                letters[c(1,4,5,8)], repetition = FALSE)

head(expGridNoOrderNoRep)
#>   Var1 Var2 Var3 Var4
#> 1    1    3    a    d
#> 2    1    3    a    e
#> 3    1    3    a    h
#> 4    1    3    b    a
#> 5    1    3    b    d
#> 6    1    3    b    e
tail(expGridNoOrderNoRep)
#>     Var1 Var2 Var3 Var4
#> 401    5    9    c    d
#> 402    5    9    c    e
#> 403    5    9    c    h
#> 404    5    9    d    e
#> 405    5    9    d    h
#> 406    5    9    e    h
```
