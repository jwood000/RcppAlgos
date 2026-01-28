# expandGrid Iterator

- Returns an iterator for iterating over the Cartesian product of the
  input vectors.

- Supports random access via the `[[` method.

- GMP support allows for exploration of cases where the number of
  products is large.

- Use the [`next`](https://rdrr.io/r/base/Control.html) methods to
  obtain results in lexicographical order.

## Usage

``` r
expandGridIter(..., nThreads = NULL, return_df = FALSE)
```

## Arguments

- ...:

  vectors, factors or a list containing these. (See
  [`?expand.grid`](https://rdrr.io/r/base/expand.grid.html)).

- nThreads:

  Specific number of threads to be used. The default is `NULL`.

- return_df:

  Logical flag to force the output to be a `data.frame`. The default is
  `FALSE`.

## Value

- If `nextIter` is called, a named vector is returned if a `matrix` can
  be returned in the general case. Otherwise, a `data.frame` is
  returned.

- When `nextNIter` and `nextRemaining` are called, a named `matrix` is
  returned when all of the input is of the same type and
  `return_df = FALSE`. Otherwise, a `data.frame` is returned.

## Details

Once you initialize a new iterator, the following methods are available:

- `nextIter`:

  Retrieve the **next** lexicographical result

- `nextNIter`:

  Pass an integer *n* to retrieve the **next** *n* lexicographical
  results

- `nextRemaining`:

  Retrieve all remaining lexicographical results

- `currIter`:

  Returns the current iteration

- `startOver`:

  Resets the iterator

- `sourceVector`:

  View the source input

- `summary`:

  Returns a list of summary information about the iterator

- `front`:

  Retrieve the **first** lexicographical result

- `back`:

  Retrieve the **last** lexicographical result

- `[[`:

  Random access method. Pass a single value or a vector of valid
  indices. If a single value is passed, the internal index of the
  iterator will be updated, however if a vector is passed the internal
  state will not change. GMP support allows for flexible indexing.

## Note

- If `nThreads` is utilized, it will only take effect if the number of
  elements requested is greater than some threshold (determined
  internally). *E.g*:

      serial   <- expandGridIter(Map(\(x, y) x:y, 1:10, 11:20))
      multi    <- expandGridIter(Map(\(x, y) x:y, 1:10, 11:20), nThreads = 4)
      fetch1e6 <- multi@nextNIter(1e6)  ## much faster than serial@nextNIter(1e6)
      fetch1e3 <- multi@nextNIter(1e3)  ## only one thread used... same as serial@nextNIter(1e3)library(microbenchmark)
      microbenchmark(multi@nextNIter(1e6), serial@nextNIter(1e6), times = 20)
      microbenchmark(multi@nextNIter(1e3), serial@nextNIter(1e3), times = 20)

- The maximum number of expandGrid that can be generated at one time is
  \\2^{31} - 1\\.

## See also

[expandGrid](https://jwood000.github.io/RcppAlgos/reference/expandGrid.md)

## Author

Joseph Wood

## Examples

``` r
a = expandGridIter(factor(state.abb), euro, islands)
a@nextIter()
#>   Var1    Var2  Var3
#> 1   AL 13.7603 11506
a@nextNIter(3)
#>   Var1    Var2  Var3
#> 1   AL 13.7603  5500
#> 2   AL 13.7603 16988
#> 3   AL 13.7603  2968
a@front()
#>   Var1    Var2  Var3
#> 1   AL 13.7603 11506
all_remaining = a@nextRemaining()
dim(all_remaining)
#> [1] 26399     3
a@summary()
#> $description
#> [1] "Cartesian Product of the source (see the sourceVector method for more info)"
#> 
#> $currentIndex
#> [1] 26401
#> 
#> $totalResults
#> [1] 26400
#> 
#> $totalRemaining
#> [1] -1
#> 
a@back()
#>   Var1    Var2 Var3
#> 1   WY 200.482   82
a[[5]]
#>   Var1    Var2 Var3
#> 1   AL 13.7603   16
a@summary()
#> $description
#> [1] "Cartesian Product of the source (see the sourceVector method for more info)"
#> 
#> $currentIndex
#> [1] 5
#> 
#> $totalResults
#> [1] 26400
#> 
#> $totalRemaining
#> [1] 26395
#> 
a[[c(1, 17, 3)]]
#>   Var1    Var2  Var3
#> 1   AL 13.7603 11506
#> 2   AL 13.7603    13
#> 3   AL 13.7603 16988
a@summary()
#> $description
#> [1] "Cartesian Product of the source (see the sourceVector method for more info)"
#> 
#> $currentIndex
#> [1] 5
#> 
#> $totalResults
#> [1] 26400
#> 
#> $totalRemaining
#> [1] 26395
#> 
```
