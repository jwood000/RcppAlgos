# Combination and Permutation Iterator

- Returns an iterator for iterating over combinations or permutations of
  a vector with or without constraints.

- Supports random access via the `[[` method.

- GMP support allows for exploration of combinations/permutations of
  vectors with many elements.

- The output is in lexicographical order for the
  [`next`](https://rdrr.io/r/base/Control.html) methods and reverse
  lexicographical order for the `prev` methods.

- Learn more in `vignette("iterators")`.

## Usage

``` r
comboIter(v, m = NULL, ...)
permuteIter(v, m = NULL, ...)

# S3 method for class 'numeric'
comboIter(v, m = NULL, repetition = FALSE, freqs = NULL,
          constraintFun = NULL, comparisonFun = NULL,
          limitConstraints = NULL, keepResults = NULL,
          FUN = NULL, Parallel = FALSE, nThreads = NULL,
          tolerance = NULL, FUN.VALUE = NULL, ...)

# S3 method for class 'numeric'
permuteIter(v, m = NULL, repetition = FALSE, freqs = NULL,
            constraintFun = NULL, comparisonFun = NULL,
            limitConstraints = NULL, keepResults = NULL,
            FUN = NULL, Parallel = FALSE, nThreads = NULL,
            tolerance = NULL, FUN.VALUE = NULL, ...)

# S3 method for class 'factor'
comboIter(
    v, m = NULL, repetition = FALSE, freqs = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, FUN.VALUE = NULL, ...
)
# S3 method for class 'factor'
permuteIter(
    v, m = NULL, repetition = FALSE, freqs = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, FUN.VALUE = NULL, ...
)

# Default S3 method
comboIter(
    v, m = NULL, repetition = FALSE, freqs = NULL,
    FUN = NULL, FUN.VALUE = NULL, ...
)
# Default S3 method
permuteIter(
    v, m = NULL, repetition = FALSE, freqs = NULL,
    FUN = NULL, FUN.VALUE = NULL, ...
)

# S3 method for class 'table'
comboIter(
    v, m = NULL, constraintFun = NULL, comparisonFun = NULL,
    limitConstraints = NULL, keepResults = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, tolerance = NULL, FUN.VALUE = NULL, ...
)
# S3 method for class 'table'
permuteIter(
    v, m = NULL, constraintFun = NULL, comparisonFun = NULL,
    limitConstraints = NULL, keepResults = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, tolerance = NULL, FUN.VALUE = NULL, ...
)

# S3 method for class 'list'
comboIter(v, m = NULL, repetition = FALSE, freqs = NULL, ...)
# S3 method for class 'list'
permuteIter(v, m = NULL, repetition = FALSE, freqs = NULL, ...)
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

- If `nextIter` or `prevIter` is called, a vector is returned

- Otherwise, a matrix with \\m\\ or \\m + 1\\ columns, depending on the
  value of `keepResults`

- If `FUN` is utilized, `FUN.VALUE = NULL`, and either `nextIter` or
  `prevIter` is called, the result will be determined by `FUN`,
  otherwise a list is returned.

- When both `FUN` and `FUN.VALUE` are not `NULL`, the return is modeled
  after the return of `vapply`. See the 'Value' section of
  [`vapply`](https://rdrr.io/r/base/lapply.html).

## Details

Once you initialize a new iterator, the following methods are available
via `@` (*e.g.* `a@nextIter()`) or `$` (*e.g.* `a$nextIter()`). The
preferred practice is to use `@` as it is much more efficient (See
examples below). Also note that not all of the methods below are
available in all cases. See
[`Combo-class`](https://jwood000.github.io/RcppAlgos/reference/Combo-Class.md),
[`Constraints-class`](https://jwood000.github.io/RcppAlgos/reference/Constraints-Class.md),
and
[`Partitions-class`](https://jwood000.github.io/RcppAlgos/reference/Partitions-Class.md):

- `nextIter`:

  Retrieve the **next** lexicographical result

- `nextNIter`:

  Pass an integer *n* to retrieve the **next** *n* lexicographical
  results

- `nextRemaining`:

  Retrieve all remaining lexicographical results

- `currIter`:

  Returns the current iteration

- `prevIter`:

  Retrieve the **previous** lexicographical result (the **next**
  *reverse* lexicographical result)

- `prevNIter`:

  Pass an integer *n* to retrieve the **previous** *n* lexicographical
  results (the **next** *n* *reverse* lexicographical results)

- `prevRemaining`:

  Retrieve all remaining *reverse* lexicographical results

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
  at one time is \\2^{31} - 1\\.

- Factor vectors are accepted. Class and level attributes are preserved
  except when `FUN` is used.

- Lexicographical ordering isn't guaranteed for permutations if the
  output is constrained.

- `FUN` will be ignored if the constraint check is satisfied.

## See also

[`comboGeneral`](https://jwood000.github.io/RcppAlgos/reference/combinatoricsGeneral.md),
[`permuteGeneral`](https://jwood000.github.io/RcppAlgos/reference/combinatoricsGeneral.md)

## References

- [Lexicographical
  Order](https://en.wikipedia.org/wiki/Lexicographical_order)

- [Reverse Lexicographical
  Order](https://oeis.org/wiki/Orderings#Reverse_lexicographic_order)

## Author

Joseph Wood

## Examples

``` r
## Typical usage
a = permuteIter(unique(state.region))
a@nextIter()
#> [1] South         West          Northeast     North Central
#> Levels: Northeast South North Central West
a@nextNIter(3)
#>      [,1]  [,2]      [,3]          [,4]         
#> [1,] South West      North Central Northeast    
#> [2,] South Northeast West          North Central
#> [3,] South Northeast North Central West         
#> Levels: Northeast South North Central West
a@front()
#> [1] South         West          Northeast     North Central
#> Levels: Northeast South North Central West
a@nextRemaining()
#>       [,1]          [,2]          [,3]          [,4]         
#>  [1,] South         West          North Central Northeast    
#>  [2,] South         Northeast     West          North Central
#>  [3,] South         Northeast     North Central West         
#>  [4,] South         North Central West          Northeast    
#>  [5,] South         North Central Northeast     West         
#>  [6,] West          South         Northeast     North Central
#>  [7,] West          South         North Central Northeast    
#>  [8,] West          Northeast     South         North Central
#>  [9,] West          Northeast     North Central South        
#> [10,] West          North Central South         Northeast    
#> [11,] West          North Central Northeast     South        
#> [12,] Northeast     South         West          North Central
#> [13,] Northeast     South         North Central West         
#> [14,] Northeast     West          South         North Central
#> [15,] Northeast     West          North Central South        
#> [16,] Northeast     North Central South         West         
#> [17,] Northeast     North Central West          South        
#> [18,] North Central South         West          Northeast    
#> [19,] North Central South         Northeast     West         
#> [20,] North Central West          South         Northeast    
#> [21,] North Central West          Northeast     South        
#> [22,] North Central Northeast     South         West         
#> [23,] North Central Northeast     West          South        
#> Levels: Northeast South North Central West
a@prevIter()
#> [1] North Central Northeast     West          South        
#> Levels: Northeast South North Central West
a@prevNIter(15)
#>       [,1]          [,2]          [,3]          [,4]         
#>  [1,] North Central Northeast     South         West         
#>  [2,] North Central West          Northeast     South        
#>  [3,] North Central West          South         Northeast    
#>  [4,] North Central South         Northeast     West         
#>  [5,] North Central South         West          Northeast    
#>  [6,] Northeast     North Central West          South        
#>  [7,] Northeast     North Central South         West         
#>  [8,] Northeast     West          North Central South        
#>  [9,] Northeast     West          South         North Central
#> [10,] Northeast     South         North Central West         
#> [11,] Northeast     South         West          North Central
#> [12,] West          North Central Northeast     South        
#> [13,] West          North Central South         Northeast    
#> [14,] West          Northeast     North Central South        
#> [15,] West          Northeast     South         North Central
#> Levels: Northeast South North Central West
a@summary()
#> $description
#> [1] "Permutations of 4 choose 4"
#> 
#> $currentIndex
#> [1] 9
#> 
#> $totalResults
#> [1] 24
#> 
#> $totalRemaining
#> [1] 15
#> 
a@back()
#> [1] North Central Northeast     West          South        
#> Levels: Northeast South North Central West
a@prevRemaining()
#>       [,1]          [,2]          [,3]          [,4]         
#>  [1,] North Central Northeast     South         West         
#>  [2,] North Central West          Northeast     South        
#>  [3,] North Central West          South         Northeast    
#>  [4,] North Central South         Northeast     West         
#>  [5,] North Central South         West          Northeast    
#>  [6,] Northeast     North Central West          South        
#>  [7,] Northeast     North Central South         West         
#>  [8,] Northeast     West          North Central South        
#>  [9,] Northeast     West          South         North Central
#> [10,] Northeast     South         North Central West         
#> [11,] Northeast     South         West          North Central
#> [12,] West          North Central Northeast     South        
#> [13,] West          North Central South         Northeast    
#> [14,] West          Northeast     North Central South        
#> [15,] West          Northeast     South         North Central
#> [16,] West          South         North Central Northeast    
#> [17,] West          South         Northeast     North Central
#> [18,] South         North Central Northeast     West         
#> [19,] South         North Central West          Northeast    
#> [20,] South         Northeast     North Central West         
#> [21,] South         Northeast     West          North Central
#> [22,] South         West          North Central Northeast    
#> [23,] South         West          Northeast     North Central
#> Levels: Northeast South North Central West
a[[5]]
#> [1] South         North Central West          Northeast    
#> Levels: Northeast South North Central West
a@summary()
#> $description
#> [1] "Permutations of 4 choose 4"
#> 
#> $currentIndex
#> [1] 5
#> 
#> $totalResults
#> [1] 24
#> 
#> $totalRemaining
#> [1] 19
#> 
a[[c(1, 17, 3)]]
#>      [,1]      [,2]          [,3]      [,4]         
#> [1,] South     West          Northeast North Central
#> [2,] Northeast North Central South     West         
#> [3,] South     Northeast     West      North Central
#> Levels: Northeast South North Central West
a@summary()
#> $description
#> [1] "Permutations of 4 choose 4"
#> 
#> $currentIndex
#> [1] 5
#> 
#> $totalResults
#> [1] 24
#> 
#> $totalRemaining
#> [1] 19
#> 

## See examples for comboGeneral where lower and upper are used
set.seed(1009)
mySamp = sort(rnorm(75, 997, 23))

b = comboIter(mySamp, 7,
              constraintFun = "sum",
              comparisonFun = ">",
              limitConstraints = 7200)
b@nextIter()
#> [1] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1025.575
b@nextNIter(3)
#>          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]
#> [1,] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1024.763
#> [2,] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1021.563
#> [3,] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1019.610
b@summary()
#> $description
#> [1] "Combinations of 75 choose 7 where the sum is > 7200"
#> 
#> $currentIndex
#> [1] 4
#> 
#> $totalResults
#> [1] NA
#> 
#> $totalRemaining
#> [1] NA
#> 
b@currIter()
#> [1] 1056.087 1038.314 1036.531 1035.189 1029.416 1026.804 1019.610

if (FALSE) { # \dontrun{
## We don't have random access or previous methods
b@back()
#> Error: no slot of name "back" for this object of class "Constraints"
b@prevIter()
#> Error: no slot of name "prevIter" for this object of class "Constraints"
} # }
```
