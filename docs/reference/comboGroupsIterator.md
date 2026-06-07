# comboGroups Iterator

- Returns an iterator for iterating over partitions of a vector into
  groups.

- Supports random access via the `[[` method.

- GMP support allows for exploration of cases where the number of
  comboGroups is large.

- Use the [`next`](https://rdrr.io/r/base/Control.html) methods to
  obtain results in lexicographical order.

## Usage

``` r
comboGroupsIter(v, numGroups = NULL, grpSizes = NULL,
                retType = "matrix", Parallel = FALSE,
                nThreads = NULL)
```

## Arguments

- v:

  Source vector. If `v` is a positive integer, it will be converted to
  the sequence `1:v`. If `v` is a negative integer, it will be converted
  to the sequence `v:-1`. All atomic types are supported (See
  [`is.atomic`](https://rdrr.io/r/base/is.recursive.html)).

- numGroups:

  An Integer. The number of groups that the vector will be partitioned
  into. The default is `NULL`. If provided and `grpSize` is `NULL`, it
  must divide the length of v (if v is a vector) or v (if v is a
  scalar).

- grpSizes:

  A vector of whole numbers representing the size of each group. The
  default is `NULL`. If provided, the sum of the elements must total the
  length of v (if v is a vector) or v (if v is a scalar).

- retType:

  A string, "3Darray" or "matrix", that determines the shape of the
  output. The default is "matrix". Note, "3Darray" can only be used when
  the size of each group is uniform. When the size of each group varies,
  the return output will always be a matrix.

- Parallel:

  Logical value indicating whether results should be generated in
  parallel using \\n - 1\\ threads, where \\n\\ is the maximum number of
  threads. The default is `FALSE`. If `nThreads` is not `NULL`, it will
  be given preference (*e.g.* if user has 8 threads with
  `Parallel = TRUE` and `nThreads = 4`, only 4 threads will be spawned).
  If your system is single-threaded, the arguments `Parallel` and
  `nThreads` are ignored.

- nThreads:

  Specific number of threads to be used. The default is `NULL`. See
  `Parallel`.

## Value

- If `nextIter` is called, a named vector is returned if
  `retType = "matrix"`. If `retType = "3Darray"`, a named matrix is
  returned.

- Otherwise a named matrix is returned when `retType = "matrix"` and a
  named 3D array is returned when `retType = "3Darray"`.

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

      serial   <- comboGroupsIter(50, 10)
      multi    <- comboGroupsIter(50, 10, nThreads = 4)
      fetch1e6 <- multi@nextNIter(1e6)  ## much faster than serial@nextNIter(1e6)
      fetch1e3 <- multi@nextNIter(1e3)  ## only one thread used... same as serial@nextNIter(1e3)library(microbenchmark)
      microbenchmark(multi@nextNIter(1e6), serial@nextNIter(1e6), times = 20)
      microbenchmark(multi@nextNIter(1e3), serial@nextNIter(1e3), times = 20)

- The maximum number of comboGroups that can be generated at one time is
  \\2^{31} - 1\\.

## See also

[comboGroups](https://jwood000.github.io/RcppAlgos/reference/comboGroups.md)

## Author

Joseph Wood

## Examples

``` r
a = comboGroupsIter(12, 3)
a@nextIter()
#> Grp1 Grp1 Grp1 Grp1 Grp2 Grp2 Grp2 Grp2 Grp3 Grp3 Grp3 Grp3 
#>    1    2    3    4    5    6    7    8    9   10   11   12 
a@nextNIter(3)
#>      Grp1 Grp1 Grp1 Grp1 Grp2 Grp2 Grp2 Grp2 Grp3 Grp3 Grp3 Grp3
#> [1,]    1    2    3    4    5    6    7    9    8   10   11   12
#> [2,]    1    2    3    4    5    6    7   10    8    9   11   12
#> [3,]    1    2    3    4    5    6    7   11    8    9   10   12
a@front()
#> Grp1 Grp1 Grp1 Grp1 Grp2 Grp2 Grp2 Grp2 Grp3 Grp3 Grp3 Grp3 
#>    1    2    3    4    5    6    7    8    9   10   11   12 
all_remaining = a@nextRemaining()
dim(all_remaining)
#> [1] 5774   12
a@summary()
#> $description
#> [1] "Partition of v of length 12 into 3 uniform groups"
#> 
#> $currentIndex
#> [1] 5776
#> 
#> $totalResults
#> [1] 5775
#> 
#> $totalRemaining
#> [1] -1
#> 
a@back()
#> Grp1 Grp1 Grp1 Grp1 Grp2 Grp2 Grp2 Grp2 Grp3 Grp3 Grp3 Grp3 
#>    1   10   11   12    2    7    8    9    3    4    5    6 
a[[5]]
#> Grp1 Grp1 Grp1 Grp1 Grp2 Grp2 Grp2 Grp2 Grp3 Grp3 Grp3 Grp3 
#>    1    2    3    4    5    6    7   12    8    9   10   11 
a@summary()
#> $description
#> [1] "Partition of v of length 12 into 3 uniform groups"
#> 
#> $currentIndex
#> [1] 5
#> 
#> $totalResults
#> [1] 5775
#> 
#> $totalRemaining
#> [1] 5770
#> 
a[[c(1, 17, 3)]]
#>      Grp1 Grp1 Grp1 Grp1 Grp2 Grp2 Grp2 Grp2 Grp3 Grp3 Grp3 Grp3
#> [1,]    1    2    3    4    5    6    7    8    9   10   11   12
#> [2,]    1    2    3    4    5    7    8   10    6    9   11   12
#> [3,]    1    2    3    4    5    6    7   10    8    9   11   12
a@summary()
#> $description
#> [1] "Partition of v of length 12 into 3 uniform groups"
#> 
#> $currentIndex
#> [1] 5
#> 
#> $totalResults
#> [1] 5775
#> 
#> $totalRemaining
#> [1] 5770
#> 
```
