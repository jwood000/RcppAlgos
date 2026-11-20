test_that("check parallel iterators", {

    check_iterators = function(a, b, parallel_limit = 20000L) {
        a@startOver()
        nr = as.integer(nrow(b))
        res = vector(mode = "logical", length = 3L)
        gap = as.integer(round(nr / 3))
        # Keep every forward batch above the production worker cutoff, even
        # though these fixtures enumerate much smaller spaces than before.
        expect_gte(gap, parallel_limit)
        expect_gte(nr - 2L * gap, parallel_limit)
        a1 = a@nextNIter(n = gap)
        res[1L] = identical(a1, b[1:gap, ])
        a1 = a@nextNIter(n = gap)
        res[2L] = identical(a1, b[(gap + 1):(2 * gap), ])
        a1 = a@nextRemaining()
        res[3L] = identical(a1, b[(2 * gap + 1):nr, ])

        ## Test traversal state after exhaustion via nextRemaining
        if ("prevIter" %in% slotNames(a)) {
            res = c(res, identical(b[nr, ], a@prevIter()))
            res = c(res, identical(b[nr - 1L, ], a@prevIter()))

            res = c(
                res, identical(
                    b[(nr - 2L):(nr - gap - 1L), ],
                    a@prevNIter(n = gap)
                )
            )

            res = c(
                res, identical(
                    b[(nr - gap):(nr - 1L), ],
                    a@nextNIter(n = gap)
                )
            )
        }

        a[[2 * gap]]
        a1  = a@nextNIter(n = 2 * gap)
        res = c(res, identical(a1, b[(2 * gap + 1):nr, ]))
        msg <- capture.output(noMore <- a@currIter())
        res = c(res, is.null(noMore))
        res = c(res, grepl("No more results", msg[1]))

        ## Test traversal state after exhaustion via an oversized
        ## nextNIter request
        if ("prevIter" %in% slotNames(a)) {
            res = c(res, identical(b[nr, ], a@prevIter()))
            res = c(res, identical(b[nr - 1L, ], a@prevIter()))

            res = c(
                res, identical(
                    b[(nr - 2L):(nr - gap - 1L), ],
                    a@prevNIter(n = gap)
                )
            )

            res = c(
                res, identical(
                    b[(nr - gap):(nr - 1L), ],
                    a@nextNIter(n = gap)
                )
            )
        }

        return(all(res))
    }

    ###### ************************ Combinations ************************ ######

    ## Combinations Distinct
    a = comboIter(19, 9, nThreads = 2)
    b = comboGeneral(19, 9)
    expect_true(check_iterators(a, b))

    ### Serial
    a = comboIter(19, 9)
    expect_true(check_iterators(a, b))

    rm(a, b)

    ## Combinations w/ Reps
    a = comboIter(12, 8, TRUE, nThreads = 2)
    b = comboGeneral(12, 8, TRUE)
    expect_true(check_iterators(a, b))

    ### Serial
    a = comboIter(12, 8, TRUE)
    expect_true(check_iterators(a, b))

    rm(a, b)

    ## Combinations of Multisets
    a = comboIter(14, 8, freqs = rep(1:4, length.out = 14), nThreads = 2)
    b = comboGeneral(14, 8, freqs = rep(1:4, length.out = 14))
    expect_true(check_iterators(a, b))

    ### Serial
    a = comboIter(14, 8, freqs = rep(1:4, length.out = 14))
    expect_true(check_iterators(a, b))

    rm(a, b)

    ###### ************************ Permutations ************************ ######

    ## Permutations Distinct
    a = permuteIter(12, 5, nThreads = 2)
    b = permuteGeneral(12, 5)
    expect_true(check_iterators(a, b))

    ### Serial
    a = permuteIter(12, 5)
    expect_true(check_iterators(a, b))

    rm(a, b)

    ## Permutations w/ Rep
    a = permuteIter(10, 5, TRUE, nThreads = 2)
    b = permuteGeneral(10, 5, TRUE)
    expect_true(check_iterators(a, b))

    ### Serial
    a = permuteIter(10, 5, TRUE)
    expect_true(check_iterators(a, b))

    rm(a, b)

    ## Permutations of Multisets
    a = permuteIter(10, 5, freqs = c(1:5, 1:5), nThreads = 2)
    b = permuteGeneral(10, 5, freqs = c(1:5, 1:5))
    expect_true(check_iterators(a, b))

    ### Serial
    a = permuteIter(10, 5, freqs = c(1:5, 1:5))
    expect_true(check_iterators(a, b))

    rm(a, b)

    ###### ************************ ComboGroups ************************* ######

    ## comboGroups Same Grp Size
    a = comboGroupsIter(15, 3, nThreads = 2)
    b = comboGroups(15, 3)
    expect_true(check_iterators(a, b))

    ### Serial
    a = comboGroupsIter(15, 3)
    expect_true(check_iterators(a, b))

    rm(a, b)

    ## comboGroups Distinct Grp Size
    a = comboGroupsIter(13, grpSizes = c(3, 4, 6), nThreads = 2)
    b = comboGroups(13, grpSizes = c(3, 4, 6))
    expect_true(check_iterators(a, b))

    ### Serial
    a = comboGroupsIter(13, grpSizes = c(3, 4, 6))
    expect_true(check_iterators(a, b))

    rm(a, b)

    ## comboGroups General Grp Sizes
    a = comboGroupsIter(15, grpSizes = c(2, 2, 2, 9), nThreads = 2)
    b = comboGroups(15, grpSizes = c(2, 2, 2, 9))
    expect_true(check_iterators(a, b))

    ### Serial
    a = comboGroupsIter(15, grpSizes = c(2, 2, 2, 9))
    expect_true(check_iterators(a, b))

    rm(a, b)

    ###### ************************* Cartesian ************************** ######

    ## Cartesian Product
    set.seed(12345678)
    v = rep(list(sample(1e6, 16)), 4)
    a = expandGridIter(v, nThreads = 2)
    b = expandGrid(v)
    expect_true(check_iterators(a, b))

    ### Serial
    a = expandGridIter(v)
    expect_true(check_iterators(a, b))

    rm(a, b)

    ###### ********************* Integer Partitions ********************* ######

    ## Partitions w/ Rep
    a = partitionsIter(v = 70, m = 10, repetition = TRUE, nThreads = 2)
    b = partitionsGeneral(v = 70, m = 10, repetition = TRUE)
    expect_true(check_iterators(a, b, parallel_limit = 40000L))

    ### Serial
    a = partitionsIter(v = 70, m = 10, repetition = TRUE)
    expect_true(check_iterators(a, b, parallel_limit = 40000L))

    rm(a, b)

    ## Partitions Distinct
    a = partitionsIter(v = 102, m = 7, nThreads = 2)
    b = partitionsGeneral(v = 102, m = 7)
    expect_true(check_iterators(a, b, parallel_limit = 40000L))

    ### Serial
    a = partitionsIter(v = 102, m = 7)
    expect_true(check_iterators(a, b, parallel_limit = 40000L))

    rm(a, b)

    ###### ******************** Integer Compositions ******************** ######

    ## Compositions w/ Rep
    a = compositionsIter(v = 31, m = 6, repetition = TRUE, nThreads = 2)
    b = compositionsGeneral(v = 31, m = 6, repetition = TRUE)
    expect_true(check_iterators(a, b, parallel_limit = 40000L))

    ### Serial
    a = compositionsIter(v = 31, m = 6, repetition = TRUE)
    expect_true(check_iterators(a, b, parallel_limit = 40000L))

    rm(a, b)

    ## Compositions Distinct
    a = compositionsIter(v = 49, m = 5, nThreads = 2)
    b = compositionsGeneral(v = 49, m = 5)
    expect_true(check_iterators(a, b, parallel_limit = 40000L))

    ### Serial
    a = compositionsIter(v = 49, m = 5)
    expect_true(check_iterators(a, b, parallel_limit = 40000L))
})
