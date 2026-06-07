# Partition/Composition Iterator

- Returns an iterator for iterating over partitions/compositions of a
  numbers.

- Supports random access via the `[[` method.

- GMP support allows for exploration of cases where the number of
  partitions/compositions is large.

- Use the [`next`](https://rdrr.io/r/base/Control.html) methods to
  obtain results in lexicographical order.

## Usage

``` r
partitionsIter(v, m = NULL, ...)
compositionsIter(v, m = NULL, ...)

# Default S3 method
partitionsIter(v, m = NULL, repetition = FALSE,
               freqs = NULL, target = NULL,
               nThreads = NULL, tolerance = NULL, ...)

# Default S3 method
compositionsIter(v, m = NULL, repetition = FALSE, freqs = NULL,
                 target = NULL, weak = FALSE, nThreads = NULL,
                 tolerance = NULL, ...)

# S3 method for class 'table'
partitionsIter(
    v, m = NULL, target = NULL, nThreads = NULL, tolerance = NULL, ...
)
# S3 method for class 'table'
compositionsIter(
    v, m = NULL, target = NULL, weak = FALSE, nThreads = NULL, tolerance = NULL, ...
)
```

## Arguments

- v:

  Source vector. If `v` is a positive integer, it will be converted to
  the sequence `1:v`. If `v` is a negative integer, it will be converted
  to the sequence `v:-1`. Only integer and numeric vectors are accepted.

- m:

  Width of the partition. If `m = NULL`, the length will be determined
  by the partitioning case (*e.g.* When we are generating distinct
  partitions of \\n\\, the width will be equal to the smallest \\m\\
  such that `sum(1:m) >= n`).

- ...:

  Further arguments passed to methods.

- repetition:

  Logical value indicating whether partitions/compositions should be
  with or without repetition. The default is `FALSE`.

- freqs:

  A vector of frequencies used for producing all partitions of a
  multiset of `v`. Each element of `freqs` represents how many times
  each element of the source vector, `v`, is repeated. It is analogous
  to the `times` argument in [`rep`](https://rdrr.io/r/base/rep.html).
  The default value is `NULL`.

- target:

  Number to be partitioned. If `NULL`, `max(v)` will be used.

- weak:

  (Compositions only) Logical flag indicating whether to allow terms of
  the sequence to be zero.

- nThreads:

  Specific number of threads to be used. The default is `NULL`.

- tolerance:

  A numeric value greater than or equal to zero. This parameter is
  utilized when a constraint is applied on a numeric vector. The default
  value is 0 when it can be determined that whole values are being
  utilized, otherwise it is `sqrt(.Machine$double.eps)` which is
  approximately \\1.5e-8\\. N.B. If the input vector is of type integer,
  this parameter will be ignored and strict equality will be enforced.

## Value

- If `nextIter` is called, a vector is returned

- Otherwise, a matrix with \\m\\ columns

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

  View the source vector

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

      serial   <- partitionsIter(1000, 10)
      multi    <- partitionsIter(1000, 10, nThreads = 4)
      fetch1e6 <- multi@nextNIter(1e6)  ## much faster than serial@nextNIter(1e6)
      fetch1e3 <- multi@nextNIter(1e3)  ## only one thread used... same as serial@nextNIter(1e3)library(microbenchmark)
      microbenchmark(multi@nextNIter(1e6), serial@nextNIter(1e6))
      microbenchmark(multi@nextNIter(1e3), serial@nextNIter(1e3))

- `nThreads` will be ignored in the following cases (i.e. Generating the
  \\n^{th}\\ partition in these cases are currently unavailable):

  - With standard multisets. If zero is the only element with a
    non-trivial multiplicity, multithreading is possible.

  - If the source vector is not isomorphic to `1:length(v)`

- The maximum number of partitions/compositions that can be generated at
  one time is \\2^{31} - 1\\.

## See also

[`partitionsGeneral`](https://jwood000.github.io/RcppAlgos/reference/partitionsGeneral.md),
[`compositionsGeneral`](https://jwood000.github.io/RcppAlgos/reference/partitionsGeneral.md)

## References

- [Lexicographical
  Order](https://en.wikipedia.org/wiki/Lexicographical_order)

- [Subset Sum Problem](https://en.wikipedia.org/wiki/Subset_sum_problem)

- [Partition (number
  theory)](https://en.wikipedia.org/wiki/Partition_(number_theory))

- [Composition
  (combinatorics))](https://en.wikipedia.org/wiki/Composition_(combinatorics))

## Author

Joseph Wood

## Examples

``` r
a = partitionsIter(0:10, repetition = TRUE)
a@nextIter()
#>  [1]  0  0  0  0  0  0  0  0  0 10
a@nextNIter(3)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    0    0    0    0    0    0    0    0    1     9
#> [2,]    0    0    0    0    0    0    0    0    2     8
#> [3,]    0    0    0    0    0    0    0    0    3     7
a@front()
#>  [1]  0  0  0  0  0  0  0  0  0 10
a@nextRemaining()
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#>  [1,]    0    0    0    0    0    0    0    0    1     9
#>  [2,]    0    0    0    0    0    0    0    0    2     8
#>  [3,]    0    0    0    0    0    0    0    0    3     7
#>  [4,]    0    0    0    0    0    0    0    0    4     6
#>  [5,]    0    0    0    0    0    0    0    0    5     5
#>  [6,]    0    0    0    0    0    0    0    1    1     8
#>  [7,]    0    0    0    0    0    0    0    1    2     7
#>  [8,]    0    0    0    0    0    0    0    1    3     6
#>  [9,]    0    0    0    0    0    0    0    1    4     5
#> [10,]    0    0    0    0    0    0    0    2    2     6
#> [11,]    0    0    0    0    0    0    0    2    3     5
#> [12,]    0    0    0    0    0    0    0    2    4     4
#> [13,]    0    0    0    0    0    0    0    3    3     4
#> [14,]    0    0    0    0    0    0    1    1    1     7
#> [15,]    0    0    0    0    0    0    1    1    2     6
#> [16,]    0    0    0    0    0    0    1    1    3     5
#> [17,]    0    0    0    0    0    0    1    1    4     4
#> [18,]    0    0    0    0    0    0    1    2    2     5
#> [19,]    0    0    0    0    0    0    1    2    3     4
#> [20,]    0    0    0    0    0    0    1    3    3     3
#> [21,]    0    0    0    0    0    0    2    2    2     4
#> [22,]    0    0    0    0    0    0    2    2    3     3
#> [23,]    0    0    0    0    0    1    1    1    1     6
#> [24,]    0    0    0    0    0    1    1    1    2     5
#> [25,]    0    0    0    0    0    1    1    1    3     4
#> [26,]    0    0    0    0    0    1    1    2    2     4
#> [27,]    0    0    0    0    0    1    1    2    3     3
#> [28,]    0    0    0    0    0    1    2    2    2     3
#> [29,]    0    0    0    0    0    2    2    2    2     2
#> [30,]    0    0    0    0    1    1    1    1    1     5
#> [31,]    0    0    0    0    1    1    1    1    2     4
#> [32,]    0    0    0    0    1    1    1    1    3     3
#> [33,]    0    0    0    0    1    1    1    2    2     3
#> [34,]    0    0    0    0    1    1    2    2    2     2
#> [35,]    0    0    0    1    1    1    1    1    1     4
#> [36,]    0    0    0    1    1    1    1    1    2     3
#> [37,]    0    0    0    1    1    1    1    2    2     2
#> [38,]    0    0    1    1    1    1    1    1    1     3
#> [39,]    0    0    1    1    1    1    1    1    2     2
#> [40,]    0    1    1    1    1    1    1    1    1     2
#> [41,]    1    1    1    1    1    1    1    1    1     1
a@summary()
#> $description
#> [1] "Partitions with repetition of 10 into 10 parts"
#> 
#> $currentIndex
#> [1] 43
#> 
#> $totalResults
#> [1] 42
#> 
#> $totalRemaining
#> [1] -1
#> 
a@back()
#>  [1] 1 1 1 1 1 1 1 1 1 1
a[[5]]
#>  [1] 0 0 0 0 0 0 0 0 4 6
a@summary()
#> $description
#> [1] "Partitions with repetition of 10 into 10 parts"
#> 
#> $currentIndex
#> [1] 5
#> 
#> $totalResults
#> [1] 42
#> 
#> $totalRemaining
#> [1] 37
#> 
a[[c(1, 17, 3)]]
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    0    0    0    0    0    0    0    0    0    10
#> [2,]    0    0    0    0    0    0    1    1    3     5
#> [3,]    0    0    0    0    0    0    0    0    2     8
a@summary()
#> $description
#> [1] "Partitions with repetition of 10 into 10 parts"
#> 
#> $currentIndex
#> [1] 5
#> 
#> $totalResults
#> [1] 42
#> 
#> $totalRemaining
#> [1] 37
#> 

## Multisets... no random access
b = partitionsIter(40, 5, freqs = rep(1:4, 10), target = 80)
b@nextIter()
#> [1]  1  2  2 35 40
b@nextNIter(10)
#>       [,1] [,2] [,3] [,4] [,5]
#>  [1,]    1    2    2   36   39
#>  [2,]    1    2    2   37   38
#>  [3,]    1    2    3   34   40
#>  [4,]    1    2    3   35   39
#>  [5,]    1    2    3   36   38
#>  [6,]    1    2    4   33   40
#>  [7,]    1    2    4   34   39
#>  [8,]    1    2    4   35   38
#>  [9,]    1    2    4   36   37
#> [10,]    1    2    5   32   40
b@summary()
#> $description
#> [1] "Partitions of a multiset of 80 into 5 parts"
#> 
#> $currentIndex
#> [1] 11
#> 
#> $totalResults
#> [1] 10144
#> 
#> $totalRemaining
#> [1] 10133
#> 
b@nextIter()
#> [1]  1  2  5 33 39
b@currIter()
#> [1]  1  2  5 33 39
```
