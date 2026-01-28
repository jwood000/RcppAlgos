# S4-classes for Exposing C++ Combinatorial Classes

The `Combo` class family are S4-classes that expose C++ classes that
provide access to iterators and other useful methods.

## Slots

of `"Combo"` and all classes inheriting from it:

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

- `randomAccess`:

  Random access method. Pass a single value or a vector of valid
  indices. If a single value is passed, the internal index of the
  iterator will be updated, however if a vector is passed the internal
  state will not change. GMP support allows for flexible indexing.

## Author

Joseph Wood

## See also

[`Partitions-class`](https://jwood000.github.io/RcppAlgos/reference/Partitions-Class.md),
[`Constraints-class`](https://jwood000.github.io/RcppAlgos/reference/Constraints-Class.md)

## Examples

``` r
  showClass("Combo")
#> Class "Combo" [package "RcppAlgos"]
#> 
#> Slots:
#>                                                                             
#> Name:            ptr     startOver      nextIter     nextNIter nextRemaining
#> Class:   externalptr      function      function      function      function
#>                                                                             
#> Name:       prevIter     prevNIter prevRemaining      currIter  randomAccess
#> Class:      function      function      function      function      function
#>                                                               
#> Name:   sourceVector         front          back       summary
#> Class:      function      function      function      function
#> 
#> Known Subclasses: "ComboApply", "ComboRes"
```
