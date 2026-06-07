# Partition a Vector into Groups

- Generate partitions of a vector into groups. See [Create Combinations
  in R by Groups](https://stackoverflow.com/q/57732672/4408538) on
  https://stackoverflow.com for a direct use case of when the groups
  sizes are equal.

- Produce results in parallel using the `Parallel` or `nThreads`
  arguments.

- GMP support allows for exploration where the number of results is
  large.

- The output is in lexicographical order by groups.

## Usage

``` r
comboGroups(v, numGroups = NULL, grpSizes = NULL,
            retType = "matrix", lower = NULL, upper = NULL,
            Parallel = FALSE, nThreads = NULL)
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

- lower:

  The lower bound. Partitions of groups are generated lexicographically,
  thus utilizing this argument will determine which specific result to
  start generating from (*e.g.* `comboGroups(8, 2, lower = 30)` is
  equivalent to `comboGroups(8, 2)[30:comboGroupsCount(8, 2), ]`). This
  argument along with `upper` is very useful for generating results in
  chunks allowing for easy parallelization.

- upper:

  The upper bound. Similar to `lower`, however this parameter allows the
  user to *stop* generation at a specific result (*e.g.*
  `comboGroups(8, 2, upper = 5)` is equivalent to
  `comboGroups(8, 2)[1:5, ]`)

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

## Details

Conceptually, this problem can be viewed as generating all permutations
of the vector `v` and removing the within group permutations. To
illustrate this, let us consider the case of generating partitions of
`1:8` into 2 groups each of size 4.

- To begin, generate the permutations of `1:8` and group the first/last
  four elements of each row.

  |     |      |     |      |      |      |      |     |      |
  |-----|------|-----|------|------|------|------|-----|------|
  |     |      |     | Grp1 |      |      | Grp2 |     |      |
  |     | C1   | C2  | C3   | C4   | C5   | C6   | C7  | C8   |
  | R1  | \| 1 | 2   | 3    | 4 \| | \| 5 | 6    | 7   | 8 \| |
  | R2  | \| 1 | 2   | 3    | 4 \| | \| 5 | 6    | 8   | 7 \| |
  | R3  | \| 1 | 2   | 3    | 4 \| | \| 5 | 7    | 6   | 8 \| |
  | R4  | \| 1 | 2   | 3    | 4 \| | \| 5 | 7    | 8   | 6 \| |
  | R5  | \| 1 | 2   | 3    | 4 \| | \| 5 | 8    | 6   | 7 \| |
  | R6  | \| 1 | 2   | 3    | 4 \| | \| 5 | 8    | 7   | 6 \| |

- Note that the permutations above are equivalent partitions of 2 groups
  of size 4 as only the last four elements are permuted. If we look at
  at the \\25^{th}\\ lexicographical permutation, we observe our second
  distinct partition.

  |         |          |       |       |          |          |       |       |          |
  |---------|----------|-------|-------|----------|----------|-------|-------|----------|
  |         |          |       | Grp1  |          |          | Grp2  |       |          |
  |         | C1       | C2    | C3    | C4       | C5       | C6    | C7    | C8       |
  | R24     | \| 1     | 2     | 3     | 4 \|     | \| 8     | 7     | 6     | 5 \|     |
  | **R25** | **\| 1** | **2** | **3** | **5 \|** | **\| 4** | **6** | **7** | **8 \|** |
  | R26     | \| 1     | 2     | 3     | 5 \|     | \| 4     | 6     | 8     | 7 \|     |
  | R27     | \| 1     | 2     | 3     | 5 \|     | \| 4     | 7     | 6     | 8 \|     |
  | R28     | \| 1     | 2     | 3     | 5 \|     | \| 4     | 7     | 8     | 6 \|     |

- Continuing on, we will reach the \\3,457^{th}\\ lexicographical
  permutation, which represents the last result:

  |           |          |       |       |          |          |       |       |          |
  |-----------|----------|-------|-------|----------|----------|-------|-------|----------|
  |           |          |       | Grp1  |          |          | Grp2  |       |          |
  |           | C1       | C2    | C3    | C4       | C5       | C6    | C7    | C8       |
  | R3454     | \| 1     | 6     | 7     | 5 \|     | \|8      | 3     | 4     | 2 \|     |
  | R3455     | \| 1     | 6     | 7     | 5 \|     | \|8      | 4     | 2     | 3 \|     |
  | R3456     | \| 1     | 6     | 7     | 5 \|     | \|8      | 4     | 3     | 2 \|     |
  | **R3457** | **\| 1** | **6** | **7** | **8 \|** | **\| 2** | **3** | **4** | **5 \|** |
  | R3458     | \| 1     | 6     | 7     | 8 \|     | \|2      | 3     | 5     | 4 \|     |

- For this small example, the method above will not be that
  computationally expensive. In fact, there are only 35 total partitions
  of `1:8` into 2 groups of size 4 out of a possible
  `factorial(8) = 40320` permutations. However, just doubling the size
  of the vector will make this approach infeasible as there are over 10
  trillion permutations of `1:16`.

- The algorithm in `comboGroups` avoids these duplicate partitions of
  groups by utilizing an efficient algorithm analogous to the
  [std::next_permutation](https://en.cppreference.com/w/cpp/algorithm/next_permutation.html)
  found in the standard algorithm library in C++.

## Value

By default, a matrix is returned with column names corresponding to the
associated group. If `retType = "3Darray"`, a named 3D array is
returned.

## Author

Joseph Wood

## Note

- The maximum number of partitions of groups that can be generated at
  one time is \\2^{31} - 1\\. Utilizing `lower` and `upper` makes it
  possible to generate additional combinations/permutations.

- The length of `grpSizes` must equal `numGroups` if both `grpSize` and
  `numGroups` are provided.

## Examples

``` r
## return a matrix
comboGroups(8, 2)
#>       Grp1 Grp1 Grp1 Grp1 Grp2 Grp2 Grp2 Grp2
#>  [1,]    1    2    3    4    5    6    7    8
#>  [2,]    1    2    3    5    4    6    7    8
#>  [3,]    1    2    3    6    4    5    7    8
#>  [4,]    1    2    3    7    4    5    6    8
#>  [5,]    1    2    3    8    4    5    6    7
#>  [6,]    1    2    4    5    3    6    7    8
#>  [7,]    1    2    4    6    3    5    7    8
#>  [8,]    1    2    4    7    3    5    6    8
#>  [9,]    1    2    4    8    3    5    6    7
#> [10,]    1    2    5    6    3    4    7    8
#> [11,]    1    2    5    7    3    4    6    8
#> [12,]    1    2    5    8    3    4    6    7
#> [13,]    1    2    6    7    3    4    5    8
#> [14,]    1    2    6    8    3    4    5    7
#> [15,]    1    2    7    8    3    4    5    6
#> [16,]    1    3    4    5    2    6    7    8
#> [17,]    1    3    4    6    2    5    7    8
#> [18,]    1    3    4    7    2    5    6    8
#> [19,]    1    3    4    8    2    5    6    7
#> [20,]    1    3    5    6    2    4    7    8
#> [21,]    1    3    5    7    2    4    6    8
#> [22,]    1    3    5    8    2    4    6    7
#> [23,]    1    3    6    7    2    4    5    8
#> [24,]    1    3    6    8    2    4    5    7
#> [25,]    1    3    7    8    2    4    5    6
#> [26,]    1    4    5    6    2    3    7    8
#> [27,]    1    4    5    7    2    3    6    8
#> [28,]    1    4    5    8    2    3    6    7
#> [29,]    1    4    6    7    2    3    5    8
#> [30,]    1    4    6    8    2    3    5    7
#> [31,]    1    4    7    8    2    3    5    6
#> [32,]    1    5    6    7    2    3    4    8
#> [33,]    1    5    6    8    2    3    4    7
#> [34,]    1    5    7    8    2    3    4    6
#> [35,]    1    6    7    8    2    3    4    5

## or a 3 dimensional array
temp = comboGroups(8, 2, retType = "3Darray")

## view the first partition
temp[1, , ]
#>      Grp1 Grp2
#> [1,]    1    5
#> [2,]    2    6
#> [3,]    3    7
#> [4,]    4    8

## Example with groups of varying size
comboGroups(8, grpSizes = c(3, 5))
#>       Grp1 Grp1 Grp1 Grp2 Grp2 Grp2 Grp2 Grp2
#>  [1,]    1    2    3    4    5    6    7    8
#>  [2,]    1    2    4    3    5    6    7    8
#>  [3,]    1    2    5    3    4    6    7    8
#>  [4,]    1    2    6    3    4    5    7    8
#>  [5,]    1    2    7    3    4    5    6    8
#>  [6,]    1    2    8    3    4    5    6    7
#>  [7,]    1    3    4    2    5    6    7    8
#>  [8,]    1    3    5    2    4    6    7    8
#>  [9,]    1    3    6    2    4    5    7    8
#> [10,]    1    3    7    2    4    5    6    8
#> [11,]    1    3    8    2    4    5    6    7
#> [12,]    1    4    5    2    3    6    7    8
#> [13,]    1    4    6    2    3    5    7    8
#> [14,]    1    4    7    2    3    5    6    8
#> [15,]    1    4    8    2    3    5    6    7
#> [16,]    1    5    6    2    3    4    7    8
#> [17,]    1    5    7    2    3    4    6    8
#> [18,]    1    5    8    2    3    4    6    7
#> [19,]    1    6    7    2    3    4    5    8
#> [20,]    1    6    8    2    3    4    5    7
#> [21,]    1    7    8    2    3    4    5    6
#> [22,]    2    3    4    1    5    6    7    8
#> [23,]    2    3    5    1    4    6    7    8
#> [24,]    2    3    6    1    4    5    7    8
#> [25,]    2    3    7    1    4    5    6    8
#> [26,]    2    3    8    1    4    5    6    7
#> [27,]    2    4    5    1    3    6    7    8
#> [28,]    2    4    6    1    3    5    7    8
#> [29,]    2    4    7    1    3    5    6    8
#> [30,]    2    4    8    1    3    5    6    7
#> [31,]    2    5    6    1    3    4    7    8
#> [32,]    2    5    7    1    3    4    6    8
#> [33,]    2    5    8    1    3    4    6    7
#> [34,]    2    6    7    1    3    4    5    8
#> [35,]    2    6    8    1    3    4    5    7
#> [36,]    2    7    8    1    3    4    5    6
#> [37,]    3    4    5    1    2    6    7    8
#> [38,]    3    4    6    1    2    5    7    8
#> [39,]    3    4    7    1    2    5    6    8
#> [40,]    3    4    8    1    2    5    6    7
#> [41,]    3    5    6    1    2    4    7    8
#> [42,]    3    5    7    1    2    4    6    8
#> [43,]    3    5    8    1    2    4    6    7
#> [44,]    3    6    7    1    2    4    5    8
#> [45,]    3    6    8    1    2    4    5    7
#> [46,]    3    7    8    1    2    4    5    6
#> [47,]    4    5    6    1    2    3    7    8
#> [48,]    4    5    7    1    2    3    6    8
#> [49,]    4    5    8    1    2    3    6    7
#> [50,]    4    6    7    1    2    3    5    8
#> [51,]    4    6    8    1    2    3    5    7
#> [52,]    4    7    8    1    2    3    5    6
#> [53,]    5    6    7    1    2    3    4    8
#> [54,]    5    6    8    1    2    3    4    7
#> [55,]    5    7    8    1    2    3    4    6
#> [56,]    6    7    8    1    2    3    4    5

total = comboGroupsCount(11, grpSizes = c(3, 3, 5))

## Start generating from particular index
comboGroups(11, grpSizes = c(3, 3, 5), lower = total - 20)
#>       Grp1 Grp1 Grp1 Grp2 Grp2 Grp2 Grp3 Grp3 Grp3 Grp3 Grp3
#>  [1,]    5    9   10    6    7   11    1    2    3    4    8
#>  [2,]    5    9   10    6    8   11    1    2    3    4    7
#>  [3,]    5    9   10    7    8   11    1    2    3    4    6
#>  [4,]    5    9   11    6    7    8    1    2    3    4   10
#>  [5,]    5    9   11    6    7   10    1    2    3    4    8
#>  [6,]    5    9   11    6    8   10    1    2    3    4    7
#>  [7,]    5    9   11    7    8   10    1    2    3    4    6
#>  [8,]    5   10   11    6    7    8    1    2    3    4    9
#>  [9,]    5   10   11    6    7    9    1    2    3    4    8
#> [10,]    5   10   11    6    8    9    1    2    3    4    7
#> [11,]    5   10   11    7    8    9    1    2    3    4    6
#> [12,]    6    7    8    9   10   11    1    2    3    4    5
#> [13,]    6    7    9    8   10   11    1    2    3    4    5
#> [14,]    6    7   10    8    9   11    1    2    3    4    5
#> [15,]    6    7   11    8    9   10    1    2    3    4    5
#> [16,]    6    8    9    7   10   11    1    2    3    4    5
#> [17,]    6    8   10    7    9   11    1    2    3    4    5
#> [18,]    6    8   11    7    9   10    1    2    3    4    5
#> [19,]    6    9   10    7    8   11    1    2    3    4    5
#> [20,]    6    9   11    7    8   10    1    2    3    4    5
#> [21,]    6   10   11    7    8    9    1    2    3    4    5
```
