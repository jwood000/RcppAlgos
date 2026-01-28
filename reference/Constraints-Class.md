# S4-class for Exposing C++ Constraints Class

The `Constraints` class is an S4-class that exposes C++ classes that
provide access to iterators and other useful methods.

## Slots

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

## Author

Joseph Wood

## See also

[`Combo-class`](https://jwood000.github.io/RcppAlgos/reference/Combo-Class.md),
[`Partitions-class`](https://jwood000.github.io/RcppAlgos/reference/Partitions-Class.md)

## Examples

``` r
  showClass("Constraints")
#> Class "Constraints" [package "RcppAlgos"]
#> 
#> Slots:
#>                                                                             
#> Name:            ptr     startOver      nextIter     nextNIter nextRemaining
#> Class:   externalptr      function      function      function      function
#>                                                 
#> Name:       currIter  sourceVector       summary
#> Class:      function      function      function
```
