reprex::reprex({
    #'
    #' This document introduces integer partitions in `RcppAlgos`. An integer partition is a way of writing a number as a sum of positive integers where order does not matter. Although partitions can be expressed as constrained combinations, their additional structure allows us to employ dedicated algorithms that are far more efficient than the more general algorithms discussed in [Constraints in RcppAlgos: Constraint-Driven Combinatorial Enumeration](<https://jwood000.github.io/RcppAlgos/articles/CombPermConstraints.html>)
    #'
    #' The focus here is on partition-specific functionality, including classic integer partitions with repetition, distinct partitions, restricted lengths, specific multiplicities, and parallel performance. Notably, `RcppAlgos` includes a specialized algorithm for generating partitions of multisets, a case that appears to be largely absent from existing software and published literature.
    #'
    #' ***
    #'
    #' # Integer Partitions
    #'
    #' Specialized algorithms are employed when it can be determined that we are looking for [integer partitions](<https://en.wikipedia.org/wiki/Partition_(number_theory)>).
    #'
    #' As of version `2.5.0`, we now have added `partitionsGeneral` which is similar to `comboGeneral` with `constraintFun = "sum"` and `comparisonFun = "=="`. Instead of using the very general `limitConstraints` parameter, we use `target` with a default of `max(v)` as it seems more fitting for partitions.
    #'
    #' ## Standard Partitions
    #'
    #' When most folks are introduced to the subject of integer partitions, they are typically presented with the idea of breaking _N_ down into parts in an unrestricted fashion. That is, parts are allowed to repeat, and the specific length is not of concern. After a few examples, partitions of a specific length are studied, followed by partitions with distinct parts, odd parts, and so on. Indeed, many textbooks present the subject in this way.
    #'
    #' In this section, we will illustrate how to attack partitions where the parts are allowed to repeat.
    #'
    #' ### Case 1: All Integer Partitions of _N_
    #'
    #' We need `v = 0:N`, `repetition = TRUE`. When we leave `m = NULL`, `m` is internally set to the length of the longest non-zero combination (this is true for all cases below).
    #'

    library(RcppAlgos)
    options(width = 90)

    packageVersion("RcppAlgos")

    cat(paste(capture.output(sessionInfo())[1:3], collapse = "\n"))

    partitionsGeneral(0:5, repetition = TRUE)

    ## Note that we could also use comboGeneral:
    ## comboGeneral(0:5, repetition = TRUE,
    ##              constraintFun = "sum",
    ##              comparisonFun = "==", limitConstraints = 5)
    ##
    ## The same goes for any of the examples below

    #'
    #' ### Case 2: Integer Partitions of _N_ of Length _m_
    #'
    #' Everything is the same as above except for explicitly setting the desired length and deciding whether to include zero or not.
    #'

    ## Including zero
    partitionsGeneral(0:5, 3, repetition = TRUE)

    ## Zero not included
    partitionsGeneral(5, 3, repetition = TRUE)

    #'
    #' ## Distinct Partitions
    #'
    #' Partitions where the parts are distinct have many interesting properties. Of note, Leonhard Euler proved that the number of partitions with distinct parts is equivalent to the number of partitions with only odd parts (See [Odd parts and distinct parts](<https://en.wikipedia.org/wiki/Integer_partition#Odd_parts_and_distinct_parts>)). We will see this in action in the section discussing [`freqs`](#using-freqs-to-refine-length).
    #'
    #' In this section, we will look at partitions where each non-zero part cannot occur more than one time.
    #'
    #' ### Case 3: Integer Partitions of _N_ into Distinct Parts
    #'
    #' Same as `Case 1 & 2` except now we have `repetition = FALSE`.
    #'

    partitionsGeneral(0:10)

    ## Zero not included and restrict the length
    partitionsGeneral(10, 3)

    ## Include zero and restrict the length
    partitionsGeneral(0:10, 3)

    ## partitions of 10 into distinct parts of every length
    lapply(1:4, function(x) {
        partitionsGeneral(10, x)
    })

    #'
    #' #### Using `freqs` to Refine Length
    #'
    #' We can utilize the `freqs` argument to obtain more distinct partitions by allowing for repeated zeros. This super optimized algorithm will only be carried out if zero is included and the number of repetitions for every number except zero is one.
    #'
    #' For example, given `v = 0:N` and `J >= 1`, if `freqs = c(J, rep(1, N))`, then the super optimized algorithm will be used, however if `freqs = c(J, 2, rep(1, N - 1))`, the general algorithm will be used. It should be noted that the general algorithms are still highly optimized so one should not fear using them.
    #'
    #' A pattern that is guaranteed to retrieve all distinct partitions of _N_ is to set `v = 0:N` and `freqs = c(N, rep(1, N))` (any unused zeros beyond the maximum partition width will be omitted). For example, when `N = 10`, the first partition is returned as `0 0 0 10`, not `0 0 0 0 0 0 0 0 0 0 10`, because the maximum width is 4 in this case.
    #'

    ## Obtain all distinct partitions of 10
    partitionsGeneral(0:10, freqs = c(10, rep(1, 10)))    ## Same as c(3, rep(1, 10))

    #'
    #' #### Euler's Theorem in Action (Odd Parts and Distinct Parts)
    #'
    #' Consider $N = 100$. The number of partitions with odd parts only is given by:
    #'

    sum(
        sapply(1:100, function(width) {
            partitionsCount(
                seq.int(1L, 99L, by = 2L), width,
                repetition = TRUE, target = 100
            )
        })
    )

    #'
    #' And now, the number of distinct parts of 100:
    #'

    partitionsCount(0:100, freqs = c(100, rep(1, 100)))

    #'
    #' *Et Voilà!* Note that, when generating the odd parts, we passed the vector `seq.int(1L, 99L, by = 2L)`. One might believe that, since the source vector is not of the standard form (i.e. `1:n` or `0:n`), standard partition algorithms cannot be applied. In `RcppAlgos`, we sanitize and transform the data into an isomorphic standard partition case if possible. Indeed, the non-exported function `partitionsDesign()` can be used to inspect the isomorphic standard representation used internally.
    #'

    ## Inspecting the partitions of 100 into odd parts of length 10
    RcppAlgos:::partitionsDesign(
        seq.int(1L, 99L, by = 2L), 10, TRUE, target = 100, showDesign = TRUE
    )

    #'
    #' #### Caveats Using `freqs`
    #'
    #' As noted in `Case 1`, if `m = NULL`, the length of the output will be determined by the longest non-zero combination that sums to _N_. On the other hand, if the user provides `m`, that value determines the output length.
    #'

    ## m is NOT NULL and output has at most 2 zeros
    partitionsGeneral(0:10, 3, freqs = c(2, rep(1, 10)))

    ## m is NULL and output has at most 2 zeros
    partitionsGeneral(0:10, freqs = c(2, rep(1, 10)))

    ## m is NOT NULL, m is > maximum width, and
    ## more zeros than maximum width are provided
    partitionsGeneral(0:10, 7, freqs = c(6, rep(1, 10)))

    ## Similar to above, however we are missing the partition
    ## with part 10 as we don't have enough zeros
    partitionsGeneral(0:10, 8, freqs = c(6, rep(1, 10)))

    ## m is simply too big... we don't have enough zeros, so no results
    partitionsGeneral(0:10, 7, freqs = c(2, rep(1, 10)))

    #'
    #' ## Partitions of Multisets
    #'
    #' In `RcppAlgos`, we've covered combinations and permutations of multisets thoroughly. This is not surprising, as these structures are well known in both academia and industry alike. The same can be said for standard and distinct partitions, which we covered above. To our knowledge, algorithms for generating partitions where each part, $p_i$, may repeat up to $r_i$ times are not readily available in combinatorics software/literature.
    #'
    #' Put simply, finding all partitions of _N_ under part-specific multiplicity constraints has proven elusive. In `RcppAlgos`, however, this is no problem.
    #'
    #' ### Case 4: Integer Partitions of _N_ into Parts of Varying Multiplicity
    #'

    ## partitions of 12 into 4 parts where each part can
    ## be used a specific number of times (e.g. 2 or 3)
    partitionsGeneral(12, 4, freqs = rep(2:3, 6))


    #'
    #' ## Using the `table` S3 Method
    #'
    #' For both the multiset case and the distinct case with repeating zeros, we can take advantage of the `table()` S3 method. This is a natural fit for multisets, since problems are often framed in terms of multiplicities:
    #'
    #' > _Given the multiset $[p_0^{r_0}, p_1^{r_1}, \ldots, p_k^{r_k}]$, where each $p_i$ repeats $r_i$ times, find ..._
    #'
    #' Without a dedicated interface, setting up `v` and `freqs` can be cumbersome, often requiring some combination of `table()`, `unique()`, `sort()`, and/or `rle()`.
    #'
    #' As of version `2.8.3`, `RcppAlgos` provides an S3 method for `table()` objects, so a multiset can be passed directly.
    #'

    ms = c(1,1,1,2,2,2,2,3,3,3,4,4,4,4,5,5,5,6)
    tab = table(ms)
    tab

    ## Find all partitions of 30 of length 10 under the multiplicities in `ms`
    head(partitionsGeneral(tab, 10, target = 30))
    tail(partitionsGeneral(tab, 10, target = 30))

    ## Find all distinct partitions of 10 (zeros supply padding capacity)
    partitionsGeneral(table(c(0L, 0L, 0L, 1:10)))

    #'
    #' ## The Role of `target`
    #'
    #' The `target` argument specifies the integer being partitioned. That is, we are finding all ways to write `target` as a sum of elements from `v`, subject to the specified multiplicity rules.
    #'
    #' When `v` is of the standard form `0:N` or `1:N`, `target` defaults to `max(v)`, which corresponds to the classical problem of partitioning _N_. However, `target` may be set independently of `v`, allowing for more general capped or restricted partition problems.
    #'
    #' Observe:
    #'

    ## Here we partition 30 using only distinct parts up to 10
    partitionsGeneral(10, 5, target = 30)

    ## Here we partition 15 using only parts up to 4 with zero included
    partitionsGeneral(0:4, 5, repetition = TRUE, target = 15)

    ## Here we partition 22 using only parts up to 8 with zero(s) included
    partitionsGeneral(0:8, 6, freqs = c(4, rep(1, 8)), target = 22)

    ## Same as above, just making use of the table method
    partitionsGeneral(table(c(rep(0L, 4), 1:8)), 6, target = 22)

    #'
    #' ## Efficiency Generating Partitions
    #'
    #' Note, as of version `2.5.0`, one can generate partitions in parallel using the `nThreads` argument.
    #'

    ## partitions of 60
    partitionsCount(0:60, repetition = TRUE)

    ## Single threaded
    system.time(partitionsGeneral(0:60, repetition = TRUE))

    ## Using nThreads
    system.time(partitionsGeneral(0:60, repetition = TRUE, nThreads = 4))

    ## partitions of 120 into distinct parts
    partitionsCount(0:120, freqs = c(120, rep(1, 120)))

    system.time(partitionsGeneral(0:120, freqs = c(120, rep(1, 120))))

    system.time(partitionsGeneral(0:120, freqs = c(120, rep(1, 120)), nThreads = 4))

    ## partitions of 100 into parts of 15 with specific multiplicity
    partitionsCount(100, 15, freqs = rep(4:8, 20))

    ## Over 6 million almost instantly!
    system.time(partitionsGeneral(100, 15, freqs = rep(4:8, 20)))

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
