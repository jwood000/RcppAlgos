test_that("check parallel iterators", {

    check_iterators = function(a, b) {
        a@startOver()
        nr = as.integer(nrow(b))
        res = vector(mode = "logical", length = 3L)
        gap = as.integer(round(nr / 3))
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
    a = comboIter(22, 10, nThreads = 2)
    b = comboGeneral(22, 10)
    expect_true(check_iterators(a, b))

    ### Serial
    a = comboIter(22, 10)
    expect_true(check_iterators(a, b))

    rm(a, b)
    gc()

    ## Combinations w/ Reps
    a = comboIter(17, 8, TRUE, nThreads = 2)
    b = comboGeneral(17, 8, TRUE)
    expect_true(check_iterators(a, b))

    ### Serial
    a = comboIter(17, 8, TRUE)
    expect_true(check_iterators(a, b))

    rm(a, b)
    gc()

    ## Combinations of Multisets
    a = comboIter(18, 8, freqs = rep(1:6, 3), nThreads = 2)
    b = comboGeneral(18, 8, freqs = rep(1:6, 3))
    expect_true(check_iterators(a, b))

    ### Serial
    a = comboIter(18, 8, freqs = rep(1:6, 3))
    expect_true(check_iterators(a, b))

    rm(a, b)
    gc()

    ###### ************************ Permutations ************************ ######

    ## Permutations Distinct
    a = permuteIter(12, 6, nThreads = 2)
    b = permuteGeneral(12, 6)
    expect_true(check_iterators(a, b))

    ### Serial
    a = permuteIter(12, 6)
    expect_true(check_iterators(a, b))

    rm(a, b)
    gc()

    ## Permutations w/ Rep
    a = permuteIter(10, 6, TRUE, nThreads = 2)
    b = permuteGeneral(10, 6, TRUE)
    expect_true(check_iterators(a, b))

    ### Serial
    a = permuteIter(10, 6, TRUE)
    expect_true(check_iterators(a, b))

    rm(a, b)
    gc()

    ## Permutations of Multisets
    a = permuteIter(10, 6, freqs = c(1:5, 1:5), nThreads = 2)
    b = permuteGeneral(10, 6, freqs = c(1:5, 1:5))
    expect_true(check_iterators(a, b))

    ### Serial
    a = permuteIter(10, 6, freqs = c(1:5, 1:5))
    expect_true(check_iterators(a, b))

    rm(a, b)
    gc()

    ###### ************************ ComboGroups ************************* ######

    ## comboGroups Same Grp Size
    a = comboGroupsIter(15, 5, nThreads = 2)
    b = comboGroups(15, 5)
    expect_true(check_iterators(a, b))

    ### Serial
    a = comboGroupsIter(15, 5)
    expect_true(check_iterators(a, b))

    rm(a, b)
    gc()

    ## comboGroups Distinct Grp Size
    a = comboGroupsIter(14, grpSizes = c(2, 3, 4, 5), nThreads = 2)
    b = comboGroups(14, grpSizes = c(2, 3, 4, 5))
    expect_true(check_iterators(a, b))

    ### Serial
    a = comboGroupsIter(14, grpSizes = c(2, 3, 4, 5))
    expect_true(check_iterators(a, b))

    rm(a, b)
    gc()

    ## comboGroups General Grp Sizes
    a = comboGroupsIter(14, grpSizes = c(3, 3, 3, 5), nThreads = 2)
    b = comboGroups(14, grpSizes = c(3, 3, 3, 5))
    expect_true(check_iterators(a, b))

    ### Serial
    a = comboGroupsIter(14, grpSizes = c(3, 3, 3, 5))
    expect_true(check_iterators(a, b))

    rm(a, b)
    gc()

    ###### ************************* Cartesian ************************** ######

    ## Cartesian Product
    set.seed(12345678)
    v = rep(list(sample(1e6, 45)), 4)
    a = expandGridIter(v, nThreads = 2)
    b = expandGrid(v)
    expect_true(check_iterators(a, b))

    ### Serial
    a = expandGridIter(v)
    expect_true(check_iterators(a, b))

    rm(a, b)
    gc()

    ###### ********************* Integer Partitions ********************* ######

    ## Partitions w/ Rep
    a = partitionsIter(v = 100, m = 10, repetition = TRUE, nThreads = 2)
    b = partitionsGeneral(v = 100, m = 10, repetition = TRUE)
    expect_true(check_iterators(a, b))

    ### Serial
    a = partitionsIter(v = 100, m = 10, repetition = TRUE)
    expect_true(check_iterators(a, b))

    rm(a, b)
    gc()

    ## Partitions Distinct
    a = partitionsIter(v = 130, m = 7, nThreads = 2)
    b = partitionsGeneral(v = 130, m = 7)
    expect_true(check_iterators(a, b))

    ### Serial
    a = partitionsIter(v = 130, m = 7)
    expect_true(check_iterators(a, b))

    rm(a, b)
    gc()

    ###### ******************** Integer Compositions ******************** ######

    ## Compositions w/ Rep
    a = compositionsIter(v = 40, m = 7, repetition = TRUE, nThreads = 2)
    b = compositionsGeneral(v = 40, m = 7, repetition = TRUE)
    expect_true(check_iterators(a, b))

    ### Serial
    a = compositionsIter(v = 40, m = 7, repetition = TRUE)
    expect_true(check_iterators(a, b))

    rm(a, b)
    gc()

    ## Compositions Distinct
    a = compositionsIter(v = 60, m = 5, nThreads = 2)
    b = compositionsGeneral(v = 60, m = 5)
    expect_true(check_iterators(a, b))

    ### Serial
    a = compositionsIter(v = 60, m = 5)
    expect_true(check_iterators(a, b))
})
