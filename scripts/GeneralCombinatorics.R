reprex::reprex({
    #'
    #' This document serves as an overview for attacking _common_ combinatorial problems in `R`. One of the goals of `RcppAlgos` is to provide a comprehensive and accessible suite of functionality so that users can easily get to the heart of their problem. As a bonus, the functions in `RcppAlgos` are extremely efficient and are constantly being improved with every release.
    #'
    #' It should be noted that this document only covers common problems. For more information on other combinatorial problems addressed by `RcppAlgos`, see the following vignettes:
    #'
    #' * [Combinatorial Sampling](<https://jwood000.github.io/RcppAlgos/articles/CombinatorialSampling.html>)
    #' * [Constraints, Partitions, and Compositions](<https://jwood000.github.io/RcppAlgos/articles/CombPermConstraints.html>)
    #' * [Attacking Problems Related to the Subset Sum Problem](<https://jwood000.github.io/RcppAlgos/articles/SubsetSum.html>)
    #' * [Combinatorial Iterators in RcppAlgos](<https://jwood000.github.io/RcppAlgos/articles/CombinatoricsIterators.html>)
    #'
    #' For much of the output below, we will be using the following function obtained here [combining head and tail methods in R
    #' ](https://stackoverflow.com/a/11601162/4408538) (credit to user @flodel)
    #'

    ht <- function(d, m = 5, n = m) {
      ## print the head and tail together
      cat("head -->\n")
      print(head(d, m))
      cat("--------\n")
      cat("tail -->\n")
      print(tail(d, n))
    }

    #'
    #' ***
    #'
    #' ## Introducing `comboGeneral` and `permuteGeneral`
    #'
    #' Easily executed with a very simple interface. The output is in [lexicographical order](<https://en.wikipedia.org/wiki/Lexicographical_order>).
    #'
    #' We first look at getting results without repetition. You can pass an integer *n* and it will be converted to the sequence `1:n`, or you can pass any vector with an atomic type (i.e. `logical`, `integer`, `numeric`, `complex`, `character`, and `raw`).
    #'

    library(RcppAlgos)
    options(width = 90)

    ## combn output for reference
    combn(4, 3)

    ## This is the same as combn expect the output is transposed
    comboGeneral(4, 3)

    ## Find all 3-tuple permutations without
    ## repetition of the numbers c(1, 2, 3, 4).
    head(permuteGeneral(4, 3))

    ## If you don't specify m, the length of v (if v is a vector) or v (if v is a
    ## scalar (see the examples above)) will be used
    v <- c(2, 3, 5, 7, 11, 13)
    comboGeneral(v)

    head(permuteGeneral(v))

    ## They are very efficient...
    system.time(comboGeneral(25, 12))

    comboCount(25, 12)

    ht(comboGeneral(25, 12))

    ## And for permutations... over 8 million instantly
    system.time(permuteGeneral(13, 7))

    permuteCount(13, 7)

    ht(permuteGeneral(13, 7))

    ## Factors are preserved
    permuteGeneral(factor(c("low", "med", "high"),
                   levels = c("low", "med", "high"),
                   ordered = TRUE))

    #'
    #' ## Combinations/Permutations with Repetition
    #'
    #' There are many problems in combinatorics which require finding combinations/permutations with repetition. This is easily achieved by setting `repetition` to `TRUE`.
    #'

    fourDays <- weekdays(as.Date("2019-10-09") + 0:3, TRUE)
    ht(comboGeneral(fourDays, repetition = TRUE))

    ## When repetition = TRUE, m can exceed length(v)
    ht(comboGeneral(fourDays, 8, TRUE))

    fibonacci <- c(1L, 2L, 3L, 5L, 8L, 13L, 21L, 34L)
    permsFib <- permuteGeneral(fibonacci, 5, TRUE)

    ht(permsFib)

    ## N.B. class is preserved
    class(fibonacci)

    class(permsFib[1, ])

    ## Binary representation of all numbers from 0 to 1023
    ht(permuteGeneral(0:1, 10, T))

    #'
    #' ## Working with Multisets
    #'
    #' Sometimes, the standard combination/permutation functions don't quite get us to our desired goal. For
    #' example, one may need all permutations of a vector with some of the elements repeated a specific
    #' number of times (i.e. a multiset). Consider the following vector `a <- c(1,1,1,1,2,2,2,7,7,7,7,7)` and one
    #' would like to find permutations of `a` of length 6. Using traditional methods, we would need to generate all
    #' permutations, then eliminate duplicate values. Even considering that `permuteGeneral` is very efficient,
    #' this approach is clunky and not as fast as it could be. Observe:
    #'

    getPermsWithSpecificRepetition <- function(z, n) {
        b <- permuteGeneral(z, n)
        myDupes <- duplicated(b)
        b[!myDupes, ]
    }

    a <- as.integer(c(1, 1, 1, 1, 2, 2, 2, 7, 7, 7, 7, 7))

    system.time(test <- getPermsWithSpecificRepetition(a, 6))

    #'
    #' ### Enter `freqs`
    #'
    #' Situations like this call for the use of the `freqs` argument. Simply, enter the number
    #' of times each unique element is repeated and Voila!
    #'

    ## Using the S3 method for class 'table'
    system.time(test2 <- permuteGeneral(table(a), 6))

    identical(test, test2)

    #'
    #' Here are some more general examples with multisets:
    #'

    ## Generate all permutations of a vector with specific
    ## length of repetition for each element (i.e. multiset)
    ht(permuteGeneral(3, freqs = c(1,2,2)))

    ## or combinations of a certain length
    comboGeneral(3, 2, freqs = c(1,2,2))

    #'
    #' ## Parallel Computing
    #'
    #' Using the parameter `Parallel` or `nThreads`, we can generate combinations/permutations with greater efficiency.
    #'

    library(microbenchmark)

    ## RcppAlgos uses the "number of threads available minus one" when Parallel is TRUE
    RcppAlgos::stdThreadMax()

    comboCount(26, 13)

    ## Compared to combn using 4 threads
    microbenchmark(combn = combn(26, 13),
                   serAlgos = comboGeneral(26, 13),
                   parAlgos = comboGeneral(26, 13, nThreads = 4),
                   times = 10,
                   unit = "relative")

    ## Using 7 cores w/ Parallel = TRUE
    microbenchmark(
        serial = comboGeneral(20, 10, freqs = rep(1:4, 5)),
        parallel = comboGeneral(20, 10, freqs = rep(1:4, 5), Parallel = TRUE),
        unit = "relative"
    )

    #'
    #' ### Using arguments `lower` and `upper`
    #'
    #' There are arguments `lower` and `upper` that can be utilized to generate chunks of combinations/permutations without having to generate all of them followed by subsetting.  As the output is in lexicographical order, these arguments specify where to start and stop generating. For example, `comboGeneral(5, 3)` outputs 10 combinations of the vector `1:5` chosen 3 at a time. We can set `lower` to 5 in order to start generation from the _5<sup>th</sup>_ lexicographical combination. Similarly, we can set `upper` to 4 in order to only generate the first 4 combinations. We can also use them together to produce only a certain chunk of combinations. For example, setting `lower` to 4 and `upper` to 6 only produces the _4<sup>th</sup>_, _5<sup>th</sup>_, and _6<sup>th</sup>_ lexicographical combinations. Observe:
    #'

    comboGeneral(5, 3, lower = 4, upper = 6)

    ## is equivalent to the following:
    comboGeneral(5, 3)[4:6, ]

    #'
    #' ### Generating Results Beyond `.Machine$integer.max`
    #'
    #' In addition to being useful by avoiding the unnecessary overhead of generating all combination/permutations followed by subsetting just to see a few specific results, lower and upper can be utilized to generate large number of combinations/permutations in parallel (see this [stackoverflow post](<https://stackoverflow.com/a/51595866/4408538>) for a real use case). Observe:

    ## Over 3 billion results
    comboCount(35, 15)

    ## 10086780 evenly divides 3247943160, otherwise you need to ensure that
    ## upper does not exceed the total number of results (E.g. see below, we
    ## would have "if ((x + foo) > 3247943160) {myUpper = 3247943160}" where
    ## foo is the size of the increment you choose to use in seq()).

    system.time(lapply(seq(1, 3247943160, 10086780), function(x) {
         temp <- comboGeneral(35, 15, lower = x, upper = x + 10086779)
         ## do something
         x
    }))

    ## Enter parallel
    library(parallel)
    system.time(mclapply(seq(1, 3247943160, 10086780), function(x) {
         temp <- comboGeneral(35, 15, lower = x, upper = x + 10086779)
         ## do something
         x
    }, mc.cores = 6))

    #'
    #' ## GMP Support
    #'
    #' The arguments `lower` and `upper` are also useful when one needs to explore combinations/permutations where the number of results is large:
    #'

    set.seed(222)
    myVec <- rnorm(1000)

    ## HUGE number of combinations
    comboCount(myVec, 50, repetition = TRUE)

    ## Let's look at one hundred thousand combinations in the range (1e15 + 1, 1e15 + 1e5)
    system.time(b <- comboGeneral(myVec, 50, TRUE,
                                  lower = 1e15 + 1,
                                  upper = 1e15 + 1e5))

    b[1:5, 45:50]

    #'
    #' ## User Defined Functions
    #'
    #' You can also pass user defined functions by utilizing the argument `FUN`. This feature's main purpose is for convenience, however it is somewhat more efficient than generating all combinations/permutations and then using a function from the `apply` family (N.B. the argument `Parallel` has no effect when `FUN` is employed).
    #'

    funCustomComb = function(n, r) {
        combs = comboGeneral(n, r)
        lapply(1:nrow(combs), function(x) cumprod(combs[x,]))
    }

    identical(funCustomComb(15, 8), comboGeneral(15, 8, FUN = cumprod))

    microbenchmark(f1 = funCustomComb(15, 8),
                   f2 = comboGeneral(15, 8, FUN = cumprod), unit = "relative")

    comboGeneral(15, 8, FUN = cumprod, upper = 3)

    ## An example involving the powerset... Note, we could
    ## have used the FUN.VALUE parameter here instead of
    ## calling unlist. See the next section.
    unlist(comboGeneral(c("", letters[1:3]), 3,
                        freqs = c(2, rep(1, 3)),
                        FUN = function(x) paste(x, collapse = "")))

    #'
    #' ### Using `FUN.VALUE`
    #'
    #' As of version `2.5.0`, we can make use of `FUN.VALUE` which serves as a template for the return value from `FUN`. The behavior is nearly identical to `vapply`:
    #'

    ## Example from earlier involving the power set
    comboGeneral(c("", letters[1:3]), 3, freqs = c(2, rep(1, 3)),
                 FUN = function(x) paste(x, collapse = ""), FUN.VALUE = "a")

    comboGeneral(15, 8, FUN = cumprod, upper = 3, FUN.VALUE = as.numeric(1:8))

    ## Fun example with binary representations... consider the following:
    permuteGeneral(0:1, 3, TRUE)

    permuteGeneral(c(FALSE, TRUE), 3, TRUE, FUN.VALUE = 1,
                   FUN = function(x) sum(2^(which(rev(x)) - 1)))

    #'
    #' ### Passing additional arguments with `...`
    #'
    #' As of version `2.8.3`, we have added the ability to pass further arguments to `FUN` via `...`.
    #'

    ## Again, same example with the power set only this time we
    ## conveniently pass the additional arguments to paste via '...'
    comboGeneral(c("", letters[1:3]), 3, freqs = c(2, rep(1, 3)),
                 FUN = paste, collapse = "", FUN.VALUE = "a")

    #'
    #' This concludes our discussion around user defined functions. There are several nice features that allow the user to more easily get the desired output with fewer function calls as well as fewer keystrokes. This was most clearly seen in our example above with the power set.
    #'
    #' We started with wrapping our call to `comboGeneral` with `unlist`, which was alleviated by the parameter `FUN.VALUE`. We then further simplified our usage of `FUN` by allowing additional arguments to be passed via `...`.
    #'
    #' ## S3 methods
    #'
    #' As of version `2.8.3`, we have added several S3 methods for convenience.
    #'
    #' Take our earlier example where we were talking about multisets.
    #'

    a <- as.integer(c(1, 1, 1, 1, 2, 2, 2, 7, 7, 7, 7, 7))

    ## Explicitly utilizing the freqs argument and determining the unique
    ## values for v... Still works, but clunky
    t1 <- permuteGeneral(rle(a)$values, 6, freqs = rle(a)$lengths)

    ## Now using the table method... much cleaner
    t2 <- permuteGeneral(table(a), 6)

    identical(t1, t2)

    #'
    #' There are other S3 methods defined that simplify the interface. Take for example the case when we want to pass a character vector. We know underneath the hood, character vectors are not thread safe so the `Parallel` and `nThreads` argument are ignored. We also know that the constraints parameters are only applicable to numeric vectors. For these reason, our default method's interface is greatly simplified:
    #'
    #' <p align="center"> <img src='default_method.png' width="400px" /> </p>
    #'
    #' We see only the necessary options. With numeric types, the options are more numerous:
    #'
    #' <p align="center"> <img src='numeric_method.png' width="400px" /> </p>
    #'
    #' There is also a `list` method that allows one to find combinations or permutations of lists:
    #'

    comboGeneral(
        list(
            numbers   = rnorm(4),
            states    = state.abb[1:5],
            some_data = data.frame(a = c('a', 'b'), b = c(10, 100))
        ),
        m = 2
    )

    #'
    #' This feature was inspired by [ggrothendieck](<https://github.com/ggrothendieck>) here: [Issue 20](<https://github.com/jwood000/RcppAlgos/issues/20>).
    #'

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")

