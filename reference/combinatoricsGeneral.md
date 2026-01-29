# Generate Combinations and Permutations of a Vector with/without Constraints

- Generate combinations or permutations of a vector with or without
  constraints.

- Produce results in parallel using the `Parallel` or `nThreads`
  arguments. You can also apply each of the five compiled functions
  given by the argument `constraintFun` in parallel.

- The arguments `lower` and `upper` make it possible to generate
  combinations/permutations in chunks allowing for parallelization via
  the `parallel-package`. This is convenient when you want to apply a
  custom function to the output in parallel.

- Attack integer partition and general subset sum problems.

- GMP support allows for exploration of combinations/permutations of
  vectors with many elements.

- The output is in lexicographical order.

## Usage

``` r
comboGeneral(v, m = NULL, ...)
permuteGeneral(v, m = NULL, ...)

# S3 method for class 'numeric'
comboGeneral(v, m = NULL, repetition = FALSE, freqs = NULL,
             lower = NULL, upper = NULL, constraintFun = NULL,
             comparisonFun = NULL, limitConstraints = NULL,
             keepResults = NULL, FUN = NULL, Parallel = FALSE,
             nThreads = NULL, tolerance = NULL, FUN.VALUE = NULL, ...)

# S3 method for class 'numeric'
permuteGeneral(v, m = NULL, repetition = FALSE, freqs = NULL,
               lower = NULL, upper = NULL, constraintFun = NULL,
               comparisonFun = NULL, limitConstraints = NULL,
               keepResults = NULL, FUN = NULL, Parallel = FALSE,
               nThreads = NULL, tolerance = NULL, FUN.VALUE = NULL, ...)

# S3 method for class 'factor'
comboGeneral(
    v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL, upper = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, FUN.VALUE = NULL, ...
)
# S3 method for class 'factor'
permuteGeneral(
    v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL, upper = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, FUN.VALUE = NULL, ...
)

# Default S3 method
comboGeneral(v, m = NULL, repetition = FALSE,
             freqs = NULL, lower = NULL, upper = NULL,
             FUN = NULL, FUN.VALUE = NULL, ...)
# Default S3 method
permuteGeneral(v, m = NULL, repetition = FALSE,
               freqs = NULL, lower = NULL, upper = NULL,
               FUN = NULL, FUN.VALUE = NULL, ...)

# S3 method for class 'table'
comboGeneral(
    v, m = NULL, lower = NULL, upper = NULL, constraintFun = NULL,
    comparisonFun = NULL, limitConstraints = NULL, keepResults = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, tolerance = NULL,
    FUN.VALUE = NULL, ...
)
# S3 method for class 'table'
permuteGeneral(
    v, m = NULL, lower = NULL, upper = NULL, constraintFun = NULL,
    comparisonFun = NULL, limitConstraints = NULL, keepResults = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, tolerance = NULL,
    FUN.VALUE = NULL, ...
)

# S3 method for class 'list'
comboGeneral(v, m = NULL, repetition = FALSE,
             freqs = NULL, lower = NULL, upper = NULL, ...)
# S3 method for class 'list'
permuteGeneral(v, m = NULL, repetition = FALSE,
               freqs = NULL, lower = NULL, upper = NULL, ...)
```

## Arguments

- v:

  Source vector. If `v` is a positive integer, it will be converted to
  the sequence `1:v`. If `v` is a negative integer, it will be converted
  to the sequence `v:-1`. All atomic types are supported (See
  [`is.atomic`](https://rdrr.io/r/base/is.recursive.html)).

- m:

  Number of elements to choose. If `repetition = TRUE` or `freqs` is
  utilized, `m` can exceed the length of `v`. If `m = NULL`, the length
  will default to `length(v)` or `sum(freqs)`.

- ...:

  Further arguments passed to methods.

- repetition:

  Logical value indicating whether combinations/permutations should be
  with or without repetition. The default is `FALSE`.

- freqs:

  A vector of frequencies used for producing all
  combinations/permutations of a multiset of `v`. Each element of
  `freqs` represents how many times each element of the source vector,
  `v`, is repeated. It is analogous to the `times` argument in
  [`rep`](https://rdrr.io/r/base/rep.html). The default value is `NULL`.

- lower:

  The lower bound. Combinations/permutations are generated
  lexicographically, thus utilizing this argument will determine which
  specific combination/permutation to start generating from (*e.g.*
  `comboGeneral(5, 3, lower = 6)` is equivalent to
  `comboGeneral(5, 3)[6:choose(5, 3), ]`). This argument along with
  `upper` is very useful for generating combinations/permutations in
  chunks allowing for easy parallelization.

- upper:

  The upper bound. Similar to `lower`, however this parameter allows the
  user to *stop* generation at a specific combination/permutation
  (*e.g.* `comboGeneral(5, 3, upper = 5)` is equivalent to
  `comboGeneral(5, 3)[1:5, ]`)

  If the output is constrained and `lower` isn't supplied, `upper`
  serves as a cap for how many results will be returned that meet the
  criteria (*e.g.* setting `upper = 100` alone will return the first 100
  results that meet the criteria, while setting `lower = 1` and
  `upper = 100` will test the first 100 results against the criteria).

  In addition to the benefits listed for `lower`, this parameter is
  useful when the total number of combinations/permutations without
  constraint is large and you expect/need a small number of
  combinations/permutations that meet a certain criteria. Using `upper`
  can improve run time if used judiciously as we call the member
  function
  [reserve](https://en.cppreference.com/w/cpp/container/vector/reserve.html)
  of
  [std::vector](https://en.cppreference.com/w/cpp/container/vector.html).
  See examples below.

- constraintFun:

  Function to be applied to the elements of `v` that should be passed as
  a string (*e.g.* `constraintFun = "sum"`). The possible constraint
  functions are: `"sum"`, `"prod"`, `"mean"`, `"max"`, & `"min"`. The
  default is `NULL`, meaning no function is applied.

- comparisonFun:

  Comparison operator that will be used to compare `limitConstraints`
  with the result of `constraintFun` applied to `v`. It should be passed
  as a string or a vector of two strings (*e.g.* `comparisonFun = "<="`
  or `comparisonFun = c(">","<")`). The possible comparison operators
  are: `"<"`, `">"`, `"<="`, `">="`, `"=="`. The default is `NULL`.

  When `comparisonFun` is a vector of two comparison strings, *e.g*
  `comparisonFun = c(comp1, comp2)`, and `limitConstraints` is a vector
  of two numerical values, *e.g* `limitConstraints = c(x1, x2)`, the
  combinations/permutations will be filtered in one of the following two
  ways:

  1.  When `comp1` is one of the 'greater-than' operators (*i.e.* "\>="
      or "\>"), `comp2` is one of the 'less-than' operators (*i.e.*
      "\<=" or "\<"), and `x1 < x2`, the combinations/permutations that
      are returned will have a value (after `constraintFun` has been
      applied) between `x1` and `x2`.

  2.  When `comp1` and `comp2` are defined as in \#1 and `x1 > x2`, the
      combinations/permutations that are returned will have a value
      outside the range of `x1` and `x2`. See the examples below.

  In other words, the first comparison operator is applied to the first
  limit and the second operator is applied to the second limit.

- limitConstraints:

  This is the value(s) that will be used for comparison. Can be passed
  as a single value or a vector of two numerical values. The default is
  `NULL`. See the definition of `comparisonFun` as well as the examples
  below for more information.

- keepResults:

  A logical flag indicating if the result of `constraintFun` applied to
  `v` should be displayed; if `TRUE`, an additional column of results
  will be added to the resulting matrix. The default is `FALSE`. If user
  is only applying `constraintFun`, `keepResults` will default to
  `TRUE`.

  *E.g*. The following are equivalent and will produce a \\4^{th}\\
  column of row sums:

  - `comboGeneral(5, 3 constraintFun = "sum", keepResults = TRUE)`

  - `comboGeneral(5, 3 constraintFun = "sum")`

- FUN:

  Function to be applied to each combination/permutation. The default is
  `NULL`.

- Parallel:

  Logical value indicating whether combinations/permutations should be
  generated in parallel using \\n - 1\\ threads, where \\n\\ is the
  maximum number of threads. The default is `FALSE`. If `nThreads` is
  not `NULL`, it will be given preference (*e.g.* if user has 8 threads
  with `Parallel = TRUE` and `nThreads = 4`, only 4 threads will be
  spawned). If your system is single-threaded, the arguments `Parallel`
  and `nThreads` are ignored.

- nThreads:

  Specific number of threads to be used. The default is `NULL`. See
  `Parallel`.

- tolerance:

  A numeric value greater than or equal to zero. This parameter is
  utilized when a constraint is applied on a numeric vector. The default
  value is 0 when it can be determined that whole values are being
  utilized, otherwise it is `sqrt(.Machine$double.eps)` which is
  approximately \\1.5e-8\\. N.B. If the input vector is of type integer,
  this parameter will be ignored and strict equality will be enforced.

- FUN.VALUE:

  A template for the return value from `FUN`. See 'Details' of
  [`vapply`](https://rdrr.io/r/base/lapply.html) for more information.

## Value

- In general, a matrix with \\m\\ or \\m + 1\\ columns, depending on the
  value of `keepResults`

- If `FUN` is utilized and `FUN.VALUE = NULL`, a list is returned

- When both `FUN` and `FUN.VALUE` are not `NULL`, the return is modeled
  after the return of `vapply`. See the 'Value' section of
  [`vapply`](https://rdrr.io/r/base/lapply.html).

## Details

For the general case, finding all combinations/permutations with
constraints is optimized by organizing them in such a way that when
`constraintFun` is applied, a *partially* monotonic sequence is
produced. Combinations/permutations are added successively, until a
particular combination exceeds the given constraint value for a given
constraint/comparison function combo. After this point, we can safely
skip several combinations knowing that they will exceed the given
constraint value.

There are special cases where more efficient algorithms are dyncamically
deployed. These cases center around the subject of integer partitions.
See
[`partitionsGeneral`](https://jwood000.github.io/RcppAlgos/reference/partitionsGeneral.md)
for more information.

When there are any negative values in `v` and `constraintFun = "prod"`,
producing a monotonic set is non-trivial for the general case. As a
result, performance will suffer as all combinations/permutations must be
tested against the constraint criteria.

## Note

- `Parallel` and `nThreads` will be ignored in the following cases:

  - When the output is constrained (except for most partitions cases)

  - If the class of the vector passed is `character`, `raw`, and
    `complex` (N.B. `Rcpp::CharacterMatrix` is not thread safe).
    Alternatively, you can generate an indexing matrix in parallel.

  - If `FUN` is utilized.

- If either `constraintFun`, `comparisonFun` or `limitConstraints` is
  `NULL` –or– if the class of the vector passed is `logical`,
  `character`, `raw`, `factor`, or `complex`, the constraint check will
  not be carried out. This is equivalent to simply finding all
  combinations/permutations of \\v\\ choose \\m\\.

- The maximum number of combinations/permutations that can be generated
  at one time is \\2^{31} - 1\\. Utilizing `lower` and `upper` makes it
  possible to generate additional combinations/permutations.

- Factor vectors are accepted. Class and level attributes are preserved
  except when `FUN` is used.

- Lexicographical ordering isn't guaranteed for permutations if `lower`
  isn't supplied and the output is constrained.

- If `lower` is supplied and the output is constrained, the
  combinations/permutations that will be tested will be in the
  lexicographical range `lower` to `upper` or up to the total possible
  number of results if `upper` is not given. See the second paragraph
  for the definition of `upper`.

- `FUN` will be ignored if the constraint check is satisfied.

## Author

Joseph Wood

## References

- [Passing user-supplied C++
  functions](https://gallery.rcpp.org/articles/passing-cpp-function-pointers/)

- [Monotonic Sequence](https://en.wikipedia.org/wiki/Monotonic_function)

- [Multiset](https://en.wikipedia.org/wiki/Multiset)

- [Lexicographical
  Order](https://en.wikipedia.org/wiki/Lexicographical_order)

- [Subset Sum Problem](https://en.wikipedia.org/wiki/Subset_sum_problem)

- [Partition (Number
  Theory)](https://en.wikipedia.org/wiki/Partition_(number_theory))

## Examples

``` r
comboGeneral(4, 3)
#>      [,1] [,2] [,3]
#> [1,]    1    2    3
#> [2,]    1    2    4
#> [3,]    1    3    4
#> [4,]    2    3    4
permuteGeneral(3)
#>      [,1] [,2] [,3]
#> [1,]    1    2    3
#> [2,]    1    3    2
#> [3,]    2    1    3
#> [4,]    2    3    1
#> [5,]    3    1    2
#> [6,]    3    2    1

permuteGeneral(factor(letters[1:3]), repetition = TRUE)
#>       [,1] [,2] [,3]
#>  [1,] a    a    a   
#>  [2,] a    a    b   
#>  [3,] a    a    c   
#>  [4,] a    b    a   
#>  [5,] a    b    b   
#>  [6,] a    b    c   
#>  [7,] a    c    a   
#>  [8,] a    c    b   
#>  [9,] a    c    c   
#> [10,] b    a    a   
#> [11,] b    a    b   
#> [12,] b    a    c   
#> [13,] b    b    a   
#> [14,] b    b    b   
#> [15,] b    b    c   
#> [16,] b    c    a   
#> [17,] b    c    b   
#> [18,] b    c    c   
#> [19,] c    a    a   
#> [20,] c    a    b   
#> [21,] c    a    c   
#> [22,] c    b    a   
#> [23,] c    b    b   
#> [24,] c    b    c   
#> [25,] c    c    a   
#> [26,] c    c    b   
#> [27,] c    c    c   
#> Levels: a b c

## permutations of the multiset :
## c(1,1,1,2,2,3)
permuteGeneral(table(c(1,1,1,2,2,3)))
#>       [,1] [,2] [,3] [,4] [,5] [,6]
#>  [1,]    1    1    1    2    2    3
#>  [2,]    1    1    1    2    3    2
#>  [3,]    1    1    1    3    2    2
#>  [4,]    1    1    2    1    2    3
#>  [5,]    1    1    2    1    3    2
#>  [6,]    1    1    2    2    1    3
#>  [7,]    1    1    2    2    3    1
#>  [8,]    1    1    2    3    1    2
#>  [9,]    1    1    2    3    2    1
#> [10,]    1    1    3    1    2    2
#> [11,]    1    1    3    2    1    2
#> [12,]    1    1    3    2    2    1
#> [13,]    1    2    1    1    2    3
#> [14,]    1    2    1    1    3    2
#> [15,]    1    2    1    2    1    3
#> [16,]    1    2    1    2    3    1
#> [17,]    1    2    1    3    1    2
#> [18,]    1    2    1    3    2    1
#> [19,]    1    2    2    1    1    3
#> [20,]    1    2    2    1    3    1
#> [21,]    1    2    2    3    1    1
#> [22,]    1    2    3    1    1    2
#> [23,]    1    2    3    1    2    1
#> [24,]    1    2    3    2    1    1
#> [25,]    1    3    1    1    2    2
#> [26,]    1    3    1    2    1    2
#> [27,]    1    3    1    2    2    1
#> [28,]    1    3    2    1    1    2
#> [29,]    1    3    2    1    2    1
#> [30,]    1    3    2    2    1    1
#> [31,]    2    1    1    1    2    3
#> [32,]    2    1    1    1    3    2
#> [33,]    2    1    1    2    1    3
#> [34,]    2    1    1    2    3    1
#> [35,]    2    1    1    3    1    2
#> [36,]    2    1    1    3    2    1
#> [37,]    2    1    2    1    1    3
#> [38,]    2    1    2    1    3    1
#> [39,]    2    1    2    3    1    1
#> [40,]    2    1    3    1    1    2
#> [41,]    2    1    3    1    2    1
#> [42,]    2    1    3    2    1    1
#> [43,]    2    2    1    1    1    3
#> [44,]    2    2    1    1    3    1
#> [45,]    2    2    1    3    1    1
#> [46,]    2    2    3    1    1    1
#> [47,]    2    3    1    1    1    2
#> [48,]    2    3    1    1    2    1
#> [49,]    2    3    1    2    1    1
#> [50,]    2    3    2    1    1    1
#> [51,]    3    1    1    1    2    2
#> [52,]    3    1    1    2    1    2
#> [53,]    3    1    1    2    2    1
#> [54,]    3    1    2    1    1    2
#> [55,]    3    1    2    1    2    1
#> [56,]    3    1    2    2    1    1
#> [57,]    3    2    1    1    1    2
#> [58,]    3    2    1    1    2    1
#> [59,]    3    2    1    2    1    1
#> [60,]    3    2    2    1    1    1

## Example with list
comboGeneral(
    v = list(
        p1 = matrix(1:10, ncol = 2),
        p2 = data.frame(a = letters, b = 1:26),
        p3 = as.complex(1:10)
    ),
    m = 2
)
#> [[1]]
#> [[1]]$p1
#>      [,1] [,2]
#> [1,]    1    6
#> [2,]    2    7
#> [3,]    3    8
#> [4,]    4    9
#> [5,]    5   10
#> 
#> [[1]]$p2
#>    a  b
#> 1  a  1
#> 2  b  2
#> 3  c  3
#> 4  d  4
#> 5  e  5
#> 6  f  6
#> 7  g  7
#> 8  h  8
#> 9  i  9
#> 10 j 10
#> 11 k 11
#> 12 l 12
#> 13 m 13
#> 14 n 14
#> 15 o 15
#> 16 p 16
#> 17 q 17
#> 18 r 18
#> 19 s 19
#> 20 t 20
#> 21 u 21
#> 22 v 22
#> 23 w 23
#> 24 x 24
#> 25 y 25
#> 26 z 26
#> 
#> 
#> [[2]]
#> [[2]]$p1
#>      [,1] [,2]
#> [1,]    1    6
#> [2,]    2    7
#> [3,]    3    8
#> [4,]    4    9
#> [5,]    5   10
#> 
#> [[2]]$p3
#>  [1]  1+0i  2+0i  3+0i  4+0i  5+0i  6+0i  7+0i  8+0i  9+0i 10+0i
#> 
#> 
#> [[3]]
#> [[3]]$p2
#>    a  b
#> 1  a  1
#> 2  b  2
#> 3  c  3
#> 4  d  4
#> 5  e  5
#> 6  f  6
#> 7  g  7
#> 8  h  8
#> 9  i  9
#> 10 j 10
#> 11 k 11
#> 12 l 12
#> 13 m 13
#> 14 n 14
#> 15 o 15
#> 16 p 16
#> 17 q 17
#> 18 r 18
#> 19 s 19
#> 20 t 20
#> 21 u 21
#> 22 v 22
#> 23 w 23
#> 24 x 24
#> 25 y 25
#> 26 z 26
#> 
#> [[3]]$p3
#>  [1]  1+0i  2+0i  3+0i  4+0i  5+0i  6+0i  7+0i  8+0i  9+0i 10+0i
#> 
#> 

#### Examples using "upper" and "lower":
## See specific range of permutations
permuteGeneral(75, 10, freqs = rep(1:3, 25),
               lower = 1e12, upper = 1e12 + 10)
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#>  [1,]    1    2    2   10   22   18   69   57   47    25
#>  [2,]    1    2    2   10   22   18   69   57   47    26
#>  [3,]    1    2    2   10   22   18   69   57   47    27
#>  [4,]    1    2    2   10   22   18   69   57   47    28
#>  [5,]    1    2    2   10   22   18   69   57   47    29
#>  [6,]    1    2    2   10   22   18   69   57   47    30
#>  [7,]    1    2    2   10   22   18   69   57   47    31
#>  [8,]    1    2    2   10   22   18   69   57   47    32
#>  [9,]    1    2    2   10   22   18   69   57   47    33
#> [10,]    1    2    2   10   22   18   69   57   47    34
#> [11,]    1    2    2   10   22   18   69   57   47    35

## Researcher only needs 10 7-tuples of mySamp
## such that the sum is greater than 7200.
## Generate some random data
set.seed(1009)
mySamp = rnorm(75, 997, 23)
comboGeneral(mySamp, 7, constraintFun = "sum",
             comparisonFun = ">", limitConstraints = 7200, upper = 10)
#>           [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]
#>  [1,] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1025.575
#>  [2,] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1024.763
#>  [3,] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1021.563
#>  [4,] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1019.610
#>  [5,] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1019.538
#>  [6,] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1018.049
#>  [7,] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1017.300
#>  [8,] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1014.765
#>  [9,] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1013.674
#> [10,] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1012.299

## Similarly, you can use "lower" to obtain the last rows.
## Generate the last 10 rows
comboGeneral(mySamp, 7, lower = choose(75, 7) - 9)
#>            [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]
#>  [1,]  982.8377 1001.800 976.6742 981.6386 980.9146 963.3604 985.2435
#>  [2,]  982.8377  985.326 976.6742 981.6386 980.9146 963.3604 985.2435
#>  [3,]  992.2474 1001.800 985.3260 976.6742 981.6386 980.9146 963.3604
#>  [4,]  992.2474 1001.800 985.3260 976.6742 981.6386 980.9146 985.2435
#>  [5,]  992.2474 1001.800 985.3260 976.6742 981.6386 963.3604 985.2435
#>  [6,]  992.2474 1001.800 985.3260 976.6742 980.9146 963.3604 985.2435
#>  [7,]  992.2474 1001.800 985.3260 981.6386 980.9146 963.3604 985.2435
#>  [8,]  992.2474 1001.800 976.6742 981.6386 980.9146 963.3604 985.2435
#>  [9,]  992.2474  985.326 976.6742 981.6386 980.9146 963.3604 985.2435
#> [10,] 1001.7997  985.326 976.6742 981.6386 980.9146 963.3604 985.2435

## Or if you would like to generate a specific chunk,
## use both "lower" and "upper". E.g. Generate one
## million combinations starting with the 900,000,001
## lexicographic combination.
t1 = comboGeneral(mySamp, 7,
                  lower = 9*10^8 + 1,
                  upper = 9*10^8 + 10^6)

## class of the source vector is preserved
class(comboGeneral(5,3)[1,]) == class(1:5)
#> [1] TRUE
class(comboGeneral(c(1,2:5),3)[1,]) == class(c(1,2:5))
#> [1] TRUE
class(comboGeneral(factor(month.name),3)[1,]) == class(factor(month.name))
#> [1] TRUE

## Using keepResults will add a column of results
comboGeneral(-3, 6, TRUE,
             constraintFun = "sum",
             comparisonFun = "==",
             limitConstraints = -8,
             keepResults = TRUE)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> [1,]   -3   -1   -1   -1   -1   -1   -8
#> [2,]   -2   -2   -1   -1   -1   -1   -8

## Using multiple constraints:

## Get combinations such that the product
## is between 3000 and 4000 inclusive
comboGeneral(5, 7, TRUE, constraintFun = "prod",
             comparisonFun = c(">=","<="),
             limitConstraints = c(3000, 4000),
             keepResults = TRUE)
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#>  [1,]    1    1    5    5    5    5    5 3125
#>  [2,]    1    2    3    4    5    5    5 3000
#>  [3,]    1    2    3    5    5    5    5 3750
#>  [4,]    1    2    4    4    4    5    5 3200
#>  [5,]    1    2    4    4    5    5    5 4000
#>  [6,]    1    3    3    3    5    5    5 3375
#>  [7,]    1    3    3    4    4    5    5 3600
#>  [8,]    1    3    4    4    4    4    4 3072
#>  [9,]    1    3    4    4    4    4    5 3840
#> [10,]    2    2    2    3    5    5    5 3000
#> [11,]    2    2    2    4    4    5    5 3200
#> [12,]    2    2    2    4    5    5    5 4000
#> [13,]    2    2    3    3    4    5    5 3600
#> [14,]    2    2    3    4    4    4    4 3072
#> [15,]    2    2    3    4    4    4    5 3840
#> [16,]    2    3    3    3    3    4    5 3240
#> [17,]    2    3    3    3    4    4    4 3456
#> [18,]    3    3    3    3    3    3    5 3645
#> [19,]    3    3    3    3    3    4    4 3888

## Or, get the combinations such that the
## product is less than or equal to 10 or
## greater than or equal to 40000
comboGeneral(5, 7, TRUE, constraintFun = "prod",
             comparisonFun = c("<=",">="),
             limitConstraints = c(10, 40000),
             keepResults = TRUE)
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7]  [,8]
#>  [1,]    1    1    1    1    1    1    1     1
#>  [2,]    1    1    1    1    1    1    2     2
#>  [3,]    1    1    1    1    1    1    3     3
#>  [4,]    1    1    1    1    1    1    4     4
#>  [5,]    1    1    1    1    1    1    5     5
#>  [6,]    1    1    1    1    1    2    2     4
#>  [7,]    1    1    1    1    1    2    3     6
#>  [8,]    1    1    1    1    1    2    4     8
#>  [9,]    1    1    1    1    1    2    5    10
#> [10,]    1    1    1    1    1    3    3     9
#> [11,]    1    1    1    1    2    2    2     8
#> [12,]    5    5    5    5    5    5    5 78125
#> [13,]    5    5    5    5    5    5    4 62500
#> [14,]    5    5    5    5    5    5    3 46875
#> [15,]    5    5    5    5    5    4    4 50000
#> [16,]    5    5    5    5    4    4    4 40000

#### General subset sum problem
set.seed(516781810)
comboGeneral(runif(100, 0, 42), 5, constraintFun = "mean",
             comparisonFun = "==", limitConstraints = 30,
             tolerance = 0.0000002)
#>          [,1]     [,2]    [,3]     [,4]     [,5]
#> [1,] 7.554598 31.01344 34.7431 37.39019 39.29867


#### Integer Partitions
comboGeneral(0:5, 5, TRUE, constraintFun = "sum",
             comparisonFun = "==", limitConstraints = 5)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    0    0    0    5
#> [2,]    0    0    0    1    4
#> [3,]    0    0    0    2    3
#> [4,]    0    0    1    1    3
#> [5,]    0    0    1    2    2
#> [6,]    0    1    1    1    2
#> [7,]    1    1    1    1    1


## Using FUN
comboGeneral(10000, 5, lower = 20, upper = 22,
             FUN = function(x) {
                 which(cummax(x) %% 2 == 1)
             })
#> [[1]]
#> [1] 1 3
#> 
#> [[2]]
#> [1] 1 3 5
#> 
#> [[3]]
#> [1] 1 3
#> 

if (FALSE) { # \dontrun{
## Parallel example generating more than 2^31 - 1 combinations.
library(parallel)
numCores = detectCores() - 1

## 10086780 evenly divides choose(35, 15) and is "small enough" to
## generate quickly in chunks.
system.time(mclapply(seq(1, comboCount(35, 15), 10086780), function(x) {
    a = comboGeneral(35, 15, lower = x, upper = x + 10086779)
    ## do something
    x
}, mc.cores = numCores))


## Find 13-tuple combinations of 1:25 such
## that the mean is less than 10
system.time(myComb <- comboGeneral(25, 13, FALSE,
                                   constraintFun = "mean",
                                   comparisonFun = "<",
                                   limitConstraints = 10))

## Alternatively, you must generate all combinations and subsequently
## subset to obtain the combinations that meet the criteria
system.time(myComb2 <- combn(25, 13))
system.time(myCols <- which(colMeans(myComb2) < 10))
system.time(myComb2 <- myComb2[, myCols])

## Any variation is much slower
system.time(myComb2 <- combn(25, 13)[,combn(25, 13, mean) < 10])

## Test equality with myComb above
all.equal(myComb, t(myComb2))

## Fun example... see stackoverflow:
## https://stackoverflow.com/q/22218640/4408538
system.time(permuteGeneral(seq(0L,100L,10L), 8, TRUE,
                           constraintFun = "sum",
                           comparisonFun = "==",
                           limitConstraints = 100))

## These are called weak integer compositions. Below, we call
## compositionsGeneral which gives the same output except it
## in lexicographical order. See 'Note' above
system.time(compositionsGeneral(seq(0L,100L,10L), 8, TRUE, weak = TRUE))
} # }
```
