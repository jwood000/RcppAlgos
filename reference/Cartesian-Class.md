# S4-class for Exposing C++ Cartesian Class

The `Cartesian` class is an S4-class that exposes C++ classes that
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

- `front`:

  Retrieve the **first** lexicographical result

- `back`:

  Retrieve the **last** lexicographical result

- `randomAccess`:

  Random access method. Pass a single value or a vector of valid
  indices. If a single value is passed, the internal index of the
  iterator will be updated, however if a vector is passed the internal
  state will not change. GMP support allows for flexible indexing.

## Author

Joseph Wood

## See also

[`ComboGroups-class`](https://jwood000.github.io/RcppAlgos/reference/ComboGroups-Class.md),
[`Combo-class`](https://jwood000.github.io/RcppAlgos/reference/Combo-Class.md),
[`Partitions-class`](https://jwood000.github.io/RcppAlgos/reference/Partitions-Class.md)

## Examples

``` r
  showClass("Cartesian")
#> Class "Cartesian" [package "RcppAlgos"]
#> 
#> Slots:
#>                                                                             
#> Name:            ptr     startOver      nextIter     nextNIter nextRemaining
#> Class:   externalptr      function      function      function      function
#>                                                                             
#> Name:       currIter  randomAccess  sourceVector         front          back
#> Class:      function      function      function      function      function
#>                     
#> Name:        summary
#> Class:      function
```
