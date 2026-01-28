# Count of the Cartesian Product

Calculate the number of Cartesian products of the input vectors.

## Usage

``` r
expandGridCount(...)
```

## Arguments

- ...:

  vectors, factors or a list containing these. (See
  [`?expand.grid`](https://rdrr.io/r/base/expand.grid.html)).

## Value

When the number of results exceeds \\2^{53} - 1\\, a number of class
`bigz` is returned.

## Author

Joseph Wood

## See also

[`expandGrid`](https://jwood000.github.io/RcppAlgos/reference/expandGrid.md)

## Examples

``` r
## description example
lst = list(1:2, 1:2)

## Using base R
t = expand.grid(lst)
nrow(t)
#> [1] 4

## vs calling expandGridCount directly
expandGridCount(lst)
#> [1] 4

## Call it just like you would if you were generating the results
expandGridCount(1:5, 3:9, letters[1:5], letters[c(1,4,5,8)])
#> [1] 700

## Same as
nrow(expand.grid(1:5, 3:9, letters[1:5], letters[c(1,4,5,8)]))
#> [1] 700

lst = Map(function(x, y) x:y, 8:33, 15:40)

## Return object of class 'bigz'
expandGridCount(lst)
#> Big Integer ('bigz') :
#> [1] 302231454903657293676544
```
