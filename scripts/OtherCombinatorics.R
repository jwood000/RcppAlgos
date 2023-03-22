reprex::reprex({
    #'
    #' In combinatorics, sometimes it can be difficult to figure out exactly which tool/method we need to attack our problem. Do we need combinations? permutations? Cartesian product? partitions? compositions? What about repetition or multiplicity? The list goes on. Oftentimes the solution ends up being overly complicated and prone to error, or more commonly, a simple brute force solution is employed. The latter is okay in some situations, but in many real world problems, this approach becomes untenable very quickly. The two types of problems addressed below fall into this category.
    #'
    #' ## Cartesian Product where Order does not Matter
    #'
    #' Given a list of vectors, _v<sub>1</sub>_, _v<sub>2</sub>_, ... , _v<sub>n</sub>_, where the intersection of two or more vectors in non-empty, find all unique combinations (order does not matter) of elements of the Cartesian product of all of the vectors.
    #'
    #' For example, lets say we have: `v1 = 1:4` and `v2 = 2:5`. The Cartesian product is given by `expand.grid(v1, v2)` (We continue to use the `ht` function defined in the [Combination and Permutation Basics](<https://jwood000.github.io/RcppAlgos/articles/GeneralCombinatorics.html>) vignette):
    #'

    ht <- function(d, m = 5, n = m) {
      ## print the head and tail together
      cat("head -->\n")
      print(head(d, m))
      cat("--------\n")
      cat("tail -->\n")
      print(tail(d, n))
    }

    expand.grid(1:4, 2:5)

    #'
    #' If we don't care about order, the following row pairs are considered equal and can therefore be pruned to obtain our desired results:
    #'
    #' * (r3, r6)
    #' * (r4, r10)
    #' * (r8, r11)
    #'
    #' With `comboGrid` no duplicates are generated:
    #'

    library(RcppAlgos)
    comboGrid(1:4, 2:5)

    #'
    #' Note that the order of `expand.grid` and `comboGrid` differ. The order of `comboGrid` is lexicographical meaning that the last column will vary the fastest whereas with `expand.grid`, the first column will vary the fastest.
    #'
    #' You will also note that the output of `expand.grid` is a `data.frame` whereas with `comboGrid`, we get a `matrix`. With `comboGrid`, we only get a `data.frame` when the classes of each vector are different as generally speaking, working with matrices is preferable.
    #'
    #' With the small example above, we only had to filter out 3 out of 16 total results (less than 20%). That isn't that bad. If this was the general case, we might as well just stick with `expand.grid` as it is very efficient. Unfortunately, this is not the general case and as the number of vectors with overlap increases, filtering will become impractical.
    #'
    #' Consider the following example:
    #'

    pools = list(c(1, 10, 14, 6),
                 c(7, 2, 4, 8, 3, 11, 12),
                 c(11, 3, 13, 4, 15, 8, 6, 5),
                 c(10, 1, 3, 2, 9, 5,  7),
                 c(1, 5, 10, 3, 8, 14),
                 c(15, 3, 7, 10, 4, 5, 8, 6),
                 c(14, 9, 11, 15),
                 c(7, 6, 13, 14, 10, 11, 9, 4),
                 c(6,  3,  2, 14,  7, 12,  9),
                 c(6, 11,  2,  5, 15,  7))

    ## If we used expand.grid, we would have to filter
    ## more than 100 million results
    prod(lengths(pools))

    ## With comboGrid, this is no problem
    system.time(myCombs <- comboGrid(pools))

    print(object.size(myCombs), unit = "Mb")

    ht(myCombs)

    ## This is just the time to create the cartesian product
    ## Generating keys, then filtering will take much more time
    system.time(cartProd <- expand.grid(pools))

    ## Creates huge object
    print(object.size(cartProd), unit = "Mb")

    ## What if we want results with unique values...
    ## Simply set repetition = FALSE
    system.time(myCombsNoRep <- comboGrid(pools, repetition = FALSE))

    ht(myCombsNoRep)

    #'
    #' The function `comboGrid` was highly inspired by the following question on stackoverflow:
    #'
    #' * [Picking unordered combinations from pools with overlap](<https://stackoverflow.com/q/51834467/4408538>)
    #'
    #' Currenlty, the underlying algorithm is not the gold standard. By that, we mean that results are not generated one by one. Efforts are underway to achieve this, but up until this point it has proven quite difficult (See the comprehensive answer by Tim Peters (yes, that [Tim Peters](<https://stackoverflow.com/users/2705542/tim-peters>))).
    #'
    #' The algorithm in `comboGrid` leverages [The Fundamental Theorem of Arithmetic](<https://en.wikipedia.org/wiki/Fundamental_theorem_of_arithmetic>) to efficiently generate keys that will be used in a hash function to determine if a particular combination of elements have been encountered. For greater efficiency, we make use of deduplication as [user2357112](<https://stackoverflow.com/a/51886857/4408538>) suggests.
    #'
    #' ### In the Wild
    #'
    #' * [R - Expand Grid Without Duplicates](<https://stackoverflow.com/q/68047141/4408538>)
    #' * [Non-redundant version of expand.grid](<https://stackoverflow.com/a/68050873/4408538>)
    #'
    #' ## Partitions of Groups of Varying Sizes with `comboGroups`
    #'
    #' Given a vector of length _n_ and _k_ groups, where _k_ divides _n_, each group is comprised of a combination of the vector chosen _g = n / k_ at a time. As is stated in the documentation (see `?comboGroups`), these can be constructed by first generating all permutations of the vector and subsequently removing entries with permuted groups. Let us consider the following example. Given `v = 1:12`, generate all partitions `v` into 3 groups each of size 4.
    #'

    funBruteGrp <- function(myLow = 1, myUp) {
        mat <- do.call(rbind, permuteGeneral(12, lower = myLow, upper = myUp,
            FUN = function(x) {
            sapply(seq(0, 8, 4), function(y) {
                 paste0(c("(", x[(y + 1):(y + 4)], ")"), collapse = " ")
            })
        }))
        colnames(mat) <- paste0("Grp", 1:3)
        rownames(mat) <- myLow:myUp
        mat
    }

    ## All of these are the same as only the 3rd group is being permuted
    funBruteGrp(myUp = 6)

    ## We found our second distinct partition
    funBruteGrp(myLow = 23, myUp = 26)

    funBruteGrp(myLow = 48, myUp = 50)

    #'
    #' We are starting to see a pattern. Each new partition is exactly 24 spots away. This makes sense as there are `factorial(4) = 24` permutations of size 4. Now, this is an oversimplification as if we simply generate every _24<sup>th</sup>_ permutation, we will still get duplication as they start to carry over to the other groups. Observe:
    #'

    do.call(rbind, lapply(seq(1, 169, 24), function(x) {
        funBruteGrp(myLow = x, myUp = x)
    }))

    #'
    #' This only gets more muddled as the number of groups increases. It is also very inefficient, however this exercise hopefully serves to better illustrate these structures.
    #'
    #' The algorithm in `comboGroups` avoids all of this duplication by implementing a novel algorithm akin to [std::next_permutation](https://en.cppreference.com/w/cpp/algorithm/next_permutation) from the algorithm library in `C++`.
    #'

    system.time(comboGroups(12, numGroups = 3))

    ht(comboGroups(12, numGroups = 3))

    #'
    #' Just as in `{combo|permute}General`, we can utilize the arguments `lower`, `upper`, `Parallel`, and `nThreads`.
    #'

    comboGroupsCount(30, 6)

    system.time(a1 <- comboGroups(30, numGroups = 6,
                                  lower = "123378675000000000",
                                  upper = "123378675005000000"))

    ## Use specific number of threads
    system.time(a2 <- comboGroups(30, numGroups = 6,
                                  lower = "123378675000000000",
                                  upper = "123378675005000000", nThreads = 4))

    ## Use n - 1 number of threads (in this case, there are 7)
    system.time(a3 <- comboGroups(30, numGroups = 6,
                                  lower = "123378675000000000",
                                  upper = "123378675005000000", Parallel = TRUE))

    identical(a1, a2)

    identical(a1, a3)

    #'
    #' As of `2.8.+` we can generate partitions of groups of varying sizes. For example, say we want to generate all partitions of the vector `v = 1:14` into 2 groups of 3 and 2 groups of 4:
    #'

    system.time(a4 <- comboGroups(14, grpSizes = c(3, 3, 4, 4)))

    ht(a4)

    #'
    #' There is one additional argument (i.e. `retType`) not present in the other two general functions that allows the user to specify the type of object returned. The user can select between `"matrix"` (the default) and `"3Darray"`. This structure has a natural connection to 3D space when the size of each group is uniform. We have a particular result (_1<sup>st</sup>_ dimension) broken down into groups (_2<sup>nd</sup>_ dimension) of a certain size (_3<sup>rd</sup>_ dimension).
    #'

    my3D <- comboGroups(factor(month.abb), 4, retType = "3Darray")
    my3D[1, , ]

    comboGroupsCount(12, 4)

    my3D[15400, , ]

    #' ### Relevant Posts on Stackoverflow as well as OEIS.
    #'
    #' * [Iterating through combinations of groups of 4 within a group of 16](<https://stackoverflow.com/a/51754958/4408538>)
    #' * [Create Combinations in R by Groups](<https://stackoverflow.com/q/57732672/4408538>)
    #' * [Algorithm that can create all combinations and all groups of those combinations](<https://stackoverflow.com/q/39126712/4408538>)
    #' * [R expand.grid for repeated combinations of a vector in groups](https://stackoverflow.com/q/74160916/4408538)
    #' * https://oeis.org/A025035 (See also sequences A025036-A025042)

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
