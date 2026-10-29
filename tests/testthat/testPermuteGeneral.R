context("testing permuteGeneral")

test_that("permuteGeneral produces correct results with no constraints and no repetition", {
    expect_equal(nrow(permuteGeneral(3, 3)), 6)
    expect_equal(as.vector(permuteGeneral(1,1)), 1)

    expect_equal(permuteGeneral(10, 5)[500:600, ], permuteGeneral(10, 5,
                                                              lower = 500,
                                                              upper = 600))

    expect_equal(permuteGeneral(month.abb[1:7], 7)[4000:5000, ],
                 permuteGeneral(month.abb[1:7], 7, lower = 4000, upper = 5000))

    expect_equal(as.vector(permuteGeneral(100,1)), 1:100)

    ## Constraint should not be carried out if no comparisonFun is given
    expect_equal(permuteGeneral(3, 3, constraintFun = "sum",
                                limitConstraints = 100), permuteGeneral(1:3, 3))

    expect_equal(nrow(permuteGeneral(8, 8)), factorial(8))

    set.seed(11)
    myNums <- rnorm(5)
    expect_equal(permuteGeneral(myNums, 3),
                 permuteGeneral(myNums, 3, freqs = rep(1, 5)))

    expect_equal(ncol(permuteGeneral(5, 3)), 3)
    expect_equal(ncol(permuteGeneral(5, 3, FALSE, constraintFun = "prod",
                                     keepResults = TRUE)), 4)
    expect_equal(nrow(permuteGeneral(5, 3, upper = 20)), 20)
    expect_equal(nrow(permuteGeneral(5, 3, FALSE, constraintFun = "prod",
                                     keepResults = TRUE, upper = 10L)), 10)
})

test_that("permuteGeneral produces correct results with no constraints and has repetition", {
    expect_equal(as.vector(permuteGeneral(1,1,TRUE)), 1)
    expect_equal(as.vector(permuteGeneral(1,5,TRUE)), rep(1, 5))
    expect_equal(permuteGeneral(letters[1:9], 5, TRUE)[59000:permuteCount(9,5,T), ],
                 permuteGeneral(letters[1:9], 5, TRUE, lower = 59000L))
    expect_equal(ncol(permuteGeneral(5, 3, TRUE)), 3)
    expect_equal(nrow(permuteGeneral(2, 2, TRUE)), 4)
    expect_equal(nrow(permuteGeneral(5, 3, TRUE, constraintFun = "prod", upper = 10)), 10)

    set.seed(111)
    myNums <- rnorm(5)
    expect_equal(permuteGeneral(myNums, 3, TRUE),
                 permuteGeneral(myNums, 3, freqs = rep(3, 5)))

    expect_true(all(permuteGeneral(3, 3, TRUE) ==
                        as.matrix(expand.grid(1:3, 1:3, 1:3))[,3:1]))

    expect_equal(nrow(permuteGeneral(5, 3, TRUE, upper = 10)), 10)
    expect_equal(ncol(permuteGeneral(5, 3, TRUE, constraintFun = "prod", keepResults = TRUE)), 4)

    ## In older versions the test below would fail b/c it would produce NaNs during prep
    expect_equal(nrow(permuteGeneral(2, 180, freqs = c(180, 2))),
                 permuteCount(2, 180, freqs = c(180, 2)))
})

test_that("permuteGeneral produces correct results with no constraints for multisets", {
    expect_equal(nrow(permuteGeneral(5, 5, freqs = 1:5, upper = 10)), 10)
    expect_equal(ncol(permuteGeneral(5, 3, FALSE,
                                     constraintFun = "prod", freqs = c(1,2,1,2,4),
                                     keepResults = TRUE)), 4)
    expect_equal(as.vector(permuteGeneral(1, 2, freqs = 2)), c(1, 1))

    expect_equal(permuteGeneral(LETTERS[1:5], 3),
                 permuteGeneral(LETTERS[1:5], 3, freqs = rep(1, 5)))

    expect_equal(permuteGeneral(month.name[1:5], 3, TRUE),
                 permuteGeneral(month.name[1:5], 3, freqs = rep(3, 5)))

    all.equal(t(as.matrix(partitions::multiset(rep(1:4, times = c(1:3, 6))))),
              permuteGeneral(4, freqs = c(1:3, 6), nThreads = 2))

    expect_equal(permuteGeneral(letters[1:3], freqs = 1:3),
                 matrix(letters[1:3][t(partitions::multiset(rep(1:3, times = 1:3)))], ncol = 6))

    expect_equal(permuteGeneral(3, lower = 3), permuteGeneral(3)[3:6, ])

    myNums2 <- 1:10 / 3
    expect_equal(permuteGeneral(myNums2, 5, freqs = rep(2, 10))[80000:90000, ],
                 permuteGeneral(myNums2, 5, freqs = rep(2, 10),
                                lower = 80000, upper = 90000))

    expect_equal(sum(permuteGeneral(3, 3, freqs = c(1, 1, 1))),
                 sum(permuteGeneral(3, 3)))

    expect_equal(permuteGeneral(factor(1:5, ordered = TRUE), 5, freqs = rep(3, 5)),
                 permuteSample(factor(1:5, ordered = TRUE), 5, freqs = rep(3, 5),
                               sampleVec = 1:permuteCount(5, 5, freqs = rep(3, 5))))

    expect_equal(permuteGeneral(5, 5),
                 permuteGeneral(5, 5, freqs = rep(1, 5)))

    expect_equal(permuteCount(30, freqs = rep(1:2, 15)),
                 permuteCount(30, 45, freqs = rep(1:2, 15)))

    expect_equal(do.call(rbind, lapply(seq(1L, 1680, 168), function(x) {
        permuteGeneral(4, freqs = c(2,1,2,3), lower = x, upper = x+167)
    })), permuteGeneral(4, 20, freqs = c(2,1,2,3)))
})

test_that("permuteGeneral produces correct results with constraints", {

    perms <- permuteGeneral(10, 5, TRUE, constraintFun = "sum",
                            comparisonFun = "==", limitConstraints = 10)
    comps <- compositionsGeneral(10, 5, TRUE)
    expect_equal(comps, perms[do.call(order, as.data.frame(perms)), ])
    ## weak does nothing here
    expect_equal(compositionsGeneral(10, 5, TRUE, weak = TRUE), comps)

    perms <- permuteGeneral(0:10, 5, TRUE, constraintFun = "sum",
                            comparisonFun = "==", limitConstraints = 10)
    comps <- compositionsGeneral(0:10, 5, TRUE, weak = TRUE)
    expect_equal(comps, perms[do.call(order, as.data.frame(perms)), ])

    expect_equal(
        nrow(
            permuteGeneral(
                15, 7, comparisonFun = "==", constraintFun = "sum",
                limitConstraints = 80, upper = 100
            )
        ),
        100
    )

    expect_equal(nrow(permuteGeneral(15, 7, TRUE,
                                     comparisonFun = "==", constraintFun = "sum",
                                     limitConstraints = 80, upper = 200)), 200)

    expect_equal(nrow(permuteGeneral(15, 7, freqs = rep(1:5, 3),
                                     comparisonFun = "==", constraintFun = "sum",
                                     limitConstraints = 80, upper = 220)), 220)

    expect_equal(nrow(permuteGeneral(3, 3, FALSE, constraintFun = "sum",
                                     comparisonFun = "==", limitConstraints = 6)), 6)

    ## NA should be removed when constraint check is carried out
    expect_equal(unique(permuteGeneral(c(NA,1:5), 5, TRUE,
                                       constraintFun = "sum",
                                       comparisonFun = "==", limitConstraints = 9,
                                       keepResults = TRUE)[,6]), 9)

    expect_true(all(permuteGeneral(5, 5L, TRUE,
                                   constraintFun = "min",
                                   comparisonFun = "<", limitConstraints = 3,
                                   keepResults = TRUE)[,6] < 3))

    expect_true(all(permuteGeneral(5, 5, TRUE,
                                   constraintFun = "prod",
                                   comparisonFun = ">", limitConstraints = 100,
                                   keepResults = TRUE)[,6] > 100))

    expect_true(all(permuteGeneral(5, 3, FALSE,
                                   constraintFun = "max",
                                   comparisonFun = "=<", limitConstraints = 4,
                                   keepResults = TRUE)[,4] <= 4))

    expect_true(all(permuteGeneral(3, 5, TRUE,
                                   constraintFun = "mean",
                                   comparisonFun = ">=", limitConstraints = 2,
                                   keepResults = TRUE)[,6] >= 2))

    expect_true(all(permuteGeneral(5, 5, FALSE, constraintFun = "sum",
                                   comparisonFun = ">", limitConstraints = 18,
                                   freqs = c(1,2,1,2,4),
                                   keepResults = TRUE)[,6] > 18))

    expect_equal(sum(permuteGeneral(4, 6, freqs = c(1,3,2,2),
                                   constraintFun = "sum",
                                   keepResults = TRUE)[,7] < 16),
                sum(apply(permuteGeneral(c(3,1,4,2), 6, freqs = c(2,1,2,3),
                                         constraintFun = "sum",
                                         comparisonFun = "<",
                                         limitConstraints = 16), 1, sum) < 16))

    expect_true(min(permuteGeneral(15, 5, constraintFun = "prod",
                                   comparisonFun = ">",
                                   limitConstraints = 3*10^5,
                                   keepResults = TRUE)[,6]) > 3*10^5)

    a <- permuteGeneral(8, 5, freqs = rep(3, 8))
    b <- apply(a, 1, min)
    expect_equal(permuteGeneral(8, 5, freqs = rep(3, 8),
                              constraintFun = "min",
                              comparisonFun = "==",
                              limitConstraints = 3,
                              lower = 17900, upper = 18500), a[(17900:18500)[b[17900:18500] == 3], ])

    a <- permuteGeneral(c(-1L,1:5), 7, T)
    b <- apply(a, 1, prod)
    expect_equal(nrow(permuteGeneral(c(-1L,1:5), 7, TRUE, constraintFun = "prod",
                              comparisonFun = c(">=","<="),
                              limitConstraints = c(2000, 5000))),
                 nrow(a[which(b >= 2000 & b <= 5000), ]))
})

test_that("permuteGeneral produces correct results with exotic constraints", {

    comp1 <- c("<", "<=")
    comp2 <- c(">", ">=")
    allPerms <- permuteGeneral(c(-6:(-1),1:2), 5, freqs = c(rep(1:3, 2), 2:3),
                              constraintFun = "prod", keepResults = TRUE)

    theSum <- allPerms[, 6]
    allPerms <- allPerms[, 1:5]
    q <- quantile(theSum)

    for (i in 1:2) {

        if (i == 1) {
            a <- comp1
            b <- comp2
        } else {
            a <- comp2
            b <- comp1
        }

        for (j in a) {
            for (k in b) {
                myComp <- c(j, k)
                myTest <- permuteGeneral(c(-6:(-1),1:2), 5,
                                         freqs = c(rep(1:3, 2), 2:3),
                                         constraintFun = "prod",
                                         comparisonFun = myComp,
                                         limitConstraints = c(q[2], q[4]))
                fun1 <- match.fun(j)
                fun2 <- match.fun(k)

                if (i == 1) {
                    temp <- allPerms[fun1(theSum, q[2]) | fun2(theSum, q[4]),]
                } else {
                    temp <- allPerms[fun1(theSum, q[2]) & fun2(theSum, q[4]),]
                }

                expect_equal(temp, myTest, info = c(myComp, q[2], q[4]))
            }
        }
    }

    allPerms <- permuteGeneral(c(-6:(-1),1:2), 5, TRUE,
                              constraintFun = "prod", keepResults = TRUE)

    theSum <- allPerms[, 6]
    allPerms <- allPerms[, 1:5]
    q <- quantile(theSum)

    for (i in 1:2) {

        if (i == 1) {
            a <- comp1
            b <- comp2
        } else {
            a <- comp2
            b <- comp1
        }

        for (j in a) {
            for (k in b) {
                myComp <- c(j, k)
                myTest <- permuteGeneral(c(-6:(-1),1:2), 5, TRUE,
                                       constraintFun = "prod", comparisonFun = myComp,
                                       limitConstraints = c(q[2], q[4]))
                fun1 <- match.fun(j)
                fun2 <- match.fun(k)

                if (i == 1) {
                    temp <- allPerms[fun1(theSum, q[2]) | fun2(theSum, q[4]),]
                } else {
                    temp <- allPerms[fun1(theSum, q[2]) & fun2(theSum, q[4]),]
                }

                expect_equal(temp, myTest)
            }
        }
    }

    allPerms <- permuteGeneral(c(-6:(-1),1:4), 5,
                              constraintFun = "prod", keepResults = TRUE)

    theSum <- allPerms[, 6]
    allPerms <- allPerms[, 1:5]
    q <- quantile(theSum)

    for (i in 1:2) {

        if (i == 1) {
            a <- comp1
            b <- comp2
        } else {
            a <- comp2
            b <- comp1
        }

        for (j in a) {
            for (k in b) {
                myComp <- c(j, k)
                myTest <- permuteGeneral(c(-6:(-1),1:4), 5,
                                       constraintFun = "prod", comparisonFun = myComp,
                                       limitConstraints = c(q[2], q[4]))
                fun1 <- match.fun(j)
                fun2 <- match.fun(k)

                if (i == 1) {
                    temp <- allPerms[fun1(theSum, q[2]) | fun2(theSum, q[4]),]
                } else {
                    temp <- allPerms[fun1(theSum, q[2]) & fun2(theSum, q[4]),]
                }
                expect_equal(temp, myTest)
            }
        }
    }
})

test_that("permuteGeneral produces correct results with use of FUN", {
    test <- permuteGeneral(6, 6, constraintFun = "mean")[, 7]
    expect_equal(as.vector(test), unlist(permuteGeneral(6, 6, FUN = mean)))

    expect_equal(sum(unlist(permuteGeneral(as.complex(c(1, -1, -1i, 1i)), 3,
                                           FUN = function(x) sum(Re(x))))), 0)

    expect_equal(class(permuteGeneral(as.complex(1:5 + 1i),
                                      3, TRUE, FUN = prod, FUN.VALUE = 1i)),
                 "complex")
    expect_equal(permuteGeneral(c("A", "B", "C"), 2,
                                FUN = function(x) paste(x, collapse = ""),
                                FUN.VALUE = "A"),
                 c("AB", "AC", "BA", "BC", "CA", "CB"))
    expect_equal(apply(expand.grid(rep(list(as.complex(1:3 + 1i)), 5)), 1,
                         function(x) sum(rev(x) / (as.complex(1:5 + 1i)))),
                 permuteGeneral(as.complex(1:3 + 1i), 5, TRUE,
                                FUN = function(x) sum(x / (as.complex(1:5 + 1i))),
                                FUN.VALUE = as.complex(1)))
    expect_equal(permuteGeneral(as.complex(1:3), FUN = function(x) as.numeric(mean(x)),
                                FUN.VALUE = 2), rep(2, permuteCount(3)))
    expect_equal(class(permuteGeneral(c(TRUE, FALSE), 5, TRUE,
                                      FUN = any, FUN.VALUE = 1L)), "integer")

    test <- permuteGeneral(6, 6, lower = 100, constraintFun = "prod")[, 7]
    expect_equal(as.vector(test), unlist(permuteGeneral(6, 6, lower = 100, FUN = prod)))

    test <- permuteGeneral(10, 5, constraintFun = "sum", keepResults = TRUE)
    expect_equal(as.vector(test[,6]), unlist(permuteGeneral(10, 5, FUN = sum)))

    test <- permuteGeneral(10, 4, TRUE)
    testFun <- apply(test, 1, function(x) mean(x) * 2)
    expect_equal(testFun, unlist(permuteGeneral(10, 4, T,
                                                FUN = function(x) {mean(x) * 2})))

    test <- permuteGeneral(8, 4, freqs = rep(1:4, 2))
    testFun <- lapply(1:nrow(test), function(x) cumsum(test[x, ]))
    expect_equal(testFun, permuteGeneral(8, 4, freqs = rep(1:4, 2), FUN = cumsum))

    test <- apply(permuteGeneral(4, 8, freqs = c(1,3,1,3)), 1, {
        function(x) paste0(cumprod(x), collapse = "")
    })

    testFun <- unlist(permuteGeneral(4, 8, freqs = c(1,3,1,3), FUN = function(x) {
        paste0(cumprod(x), collapse = "")
    }))

    expect_equal(test, testFun)

    expect_equal(permuteGeneral(4, 8, freqs = c(1,3,1,3),
                                constraintFun = "sum", nThreads = 2)[, 9],
                 rowSums(permuteGeneral(4, 8, freqs = c(1,3,1,3), nThreads = 2)))
})

test_that("permuteGeneral produces correct results with very large results", {
    ##******** BIG TESTS *********##

    ## NO REPETITION
    numR <- permuteCount(1000, 10)
    n1 <- gmp::sub.bigz(numR, 99)

    ## accepts raw values
    expect_equal(nrow(permuteGeneral(1000, 10, lower = n1)), 100)
    ## accepts characters
    expect_equal(nrow(permuteGeneral(1000, 10, lower = as.character(n1))), 100)
    expect_equal(as.vector(permuteGeneral(1000, 10, lower = numR)),
                 1000:991)

    ## WITH REPETITION
    numR <- permuteCount(1000, 10, TRUE)
    n1 <- gmp::sub.bigz(numR, 99)

    expect_equal(nrow(permuteGeneral(1000, 10, TRUE, lower = n1)), 100)
    expect_equal(nrow(permuteGeneral(1000, 10, TRUE, lower = as.character(n1))), 100)
    expect_equal(as.vector(permuteGeneral(1000, 10, TRUE, lower = numR)),
                 rep(1000, 10))

    ## MULTISETS
    numR <- permuteCount(100, 10, freqs = rep(1:4, 25))
    n1 <- gmp::sub.bigz(numR, 99)

    expect_equal(nrow(permuteGeneral(100, 10, freqs = rep(1:4, 25), lower = n1)), 100)
    expect_equal(nrow(permuteGeneral(100, 10, freqs = rep(1:4, 25), lower = as.character(n1))), 100)
    expect_equal(as.vector(permuteGeneral(100, 10, freqs = rep(1:4, 25), lower = numR)),
                 rep(100:97, times = 4:1))
})

test_that("permuteCount produces correct results under partition constraints", {
    expect_equal(
        permuteCount(
            10, 5, TRUE, constraintFun = "sum",
            comparisonFun = "==", limitConstraints = 10
        ),
        compositionsCount(10, 5, TRUE)
    )

    ## See issue #72
    expect_true(
        class(permuteCount(v = 2, m = 2000, freqs = rep(1000, 2))) == "bigz"
    )

    expect_identical(
        permuteCount(v = 2, m = 2000, freqs = rep(1000, 2)),
        gmp::chooseZ(2000, 1000)
    )

    expect_equal(compositionsCount(0, 7, freqs = 8), 1)
    expect_equal(compositionsCount(0, 7, freqs = 8, weak = TRUE), 1)

    expect_equal(
        permuteCount(
            100, 15, TRUE, constraintFun = "sum",
            comparisonFun = "==", limitConstraints = 100
        ),
        compositionsCount(100, 15, TRUE)
    )

    expect_equal(
        permuteCount(
            0:100, 10, constraintFun = "sum",
            comparisonFun = "==", limitConstraints = 100
        ),
        compositionsCount(0:100, 10, weak = TRUE)
    )

    expect_equal(
        permuteCount(
            0:100, 10, freqs = c(4, rep(1, 100)),
            constraintFun = "sum", comparisonFun = "==",
            limitConstraints = 100
        ),
        compositionsCount(0:100, 10, freqs = c(4, rep(1, 100)), weak = TRUE)
    )

    expect_equal(
        permuteCount(
            100, 10, constraintFun = "sum",
            comparisonFun = "==", limitConstraints = 100
        ),
        compositionsCount(100, 10)
    )

    expect_equal(
        permuteCount(
            table(c(rep(0L, 5), 1:100)), constraintFun = "sum",
            comparisonFun = "==", limitConstraints = 100
        ),
        compositionsCount(table(c(rep(0L, 5), 1:100)), weak = TRUE)
    )

    expect_equal(
        permuteCount(
            1:15, 8, freqs = rep(5:1, each = 3), constraintFun = "sum",
            comparisonFun = "==", limitConstraints = 40
        ),
        compositionsCount(1:15, 8, freqs = rep(5:1, each = 3), target = 40)
    )
})


test_that("Multiset permutation boundaries have independent counts", {
    freqs <- c(3L, 1L, 1L, 1L)
    for (width in c(4L, 6L)) {
        expected <- if (width == 4L) 72L else 120L
        rows <- permuteGeneral(0:3, width, freqs = freqs)
        expect_identical(dim(rows), c(expected, width))
        expect_identical(nrow(unique(rows)), expected)
    }

    rows <- permuteGeneral(
        0:3, 4L, freqs = freqs, lower = 1L,
        constraintFun = "sum", comparisonFun = "==", limitConstraints = 5L
    )
    # Exactly 4! / 2! arrangements of two zeros, a two, and a three.
    expect_identical(dim(rows), c(12L, 4L))
    expect_identical(nrow(unique(rows)), 12L)
    expect_true(all(apply(rows, 1L, function(row) {
        identical(sort(row), c(0L, 0L, 2L, 3L))
    })))
})

test_that("Permutation sum bounds filter ordinary candidate ranks", {
    expect_identical(
        permutePartsDesign(0:3, 4, freqs = c(3, 1, 1, 1),
                           target = 5)$partition_type,
        "PrmDstPrtCapMZ"
    )
    cases <- list(
        list(args = list(v = 1:8, m = 3, repetition = TRUE), target = 8),
        list(args = list(v = 0:8, m = 3, repetition = TRUE), target = 8),
        list(args = list(v = 1:5, m = 3, repetition = TRUE), target = 9),
        list(args = list(v = 1:8, m = 3), target = 8),
        list(args = list(v = 0:8, m = 3), target = 8),
        list(args = list(v = 0:8, m = 3, freqs = c(2, rep(1, 8))),
             target = 8),
        list(args = list(v = 1:5, m = 3), target = 9),
        list(args = list(v = 0:5, m = 3, freqs = c(2, rep(1, 5))),
             target = 9),
        list(args = list(v = 0:5, m = 3, freqs = rep(1:3, 2)), target = 8),
        # Four distinct values, but six occurrences: width four is partial.
        list(args = list(v = 0:3, m = 4, freqs = c(3, 1, 1, 1)), target = 5),
        list(args = list(v = table(c(0L, 0L, 0L, 1:3)), m = 4), target = 5),
        list(args = list(v = 0:3, m = 6, freqs = c(3, 1, 1, 1)), target = 6),
        list(args = list(v = (1:5) * 2e9 + 7, m = 3,
                         freqs = c(1, 2, 3, 2, 1)), target = 8 * 2e9 + 21),
        list(args = list(v = c(0L, 1L, 3L, 7L), m = 3,
                         freqs = c(2, 1, 2, 1)), target = 8),
        list(args = list(v = 1:3, freqs = c(2, 1, 1)), target = 7),
        list(args = list(v = 1:5, m = 1), target = 3),
        list(args = list(v = 1:5, m = 3), target = 100)
    )
    for (case in cases) {
        candidates <- do.call(permuteGeneral, case$args)
        sums <- rowSums(candidates)
        total <- nrow(candidates)
        args <- c(case$args, list(constraintFun = "sum", comparisonFun = "==",
                                 limitConstraints = case$target))
        for (bounds in list(c(1L, min(25L, total)),
                            c(2L, total), c(total, total))) {
            rows <- seq.int(bounds[1L], bounds[2L])
            expected <- candidates[rows[sums[rows] == case$target], ,
                                   drop = FALSE]
            for (keep in c(FALSE, TRUE)) {
                actual <- do.call(permuteGeneral, c(args, list(
                    lower = bounds[1L], upper = bounds[2L], keepResults = keep
                )))
                bench <- if (keep) {
                    unname(cbind(expected, rep(case$target, nrow(expected))))
                } else expected
                expect_equal(actual, bench)
            }
        }
        expect_equal(
            do.call(permuteGeneral, c(args, list(lower = 2))),
            candidates[which(sums == case$target & seq_len(total) >= 2L), ,
                       drop = FALSE]
        )
        # Omitting lower preserves the existing partition-based output order.
        if (any(sums == case$target)) {
            all_matches <- do.call(permuteGeneral, args)
            size <- min(3L, nrow(all_matches))
            expect_equal(
                do.call(permuteGeneral, c(args, list(upper = size))),
                all_matches[seq_len(size), , drop = FALSE]
            )
        }
    }
})

test_that("Permutation candidate bounds preserve exact GMP ranks", {
    lower <- gmp::as.bigz(2)^53 + 1L
    upper <- lower + 40L
    for (mode in list(list(), list(repetition = TRUE),
                      list(freqs = rep(2L, 30)))) {
        args <- c(list(v = 1:30, m = 12), mode)
        candidates <- do.call(permuteSample, c(args, list(
            sampleVec = lower + 0:40
        )))
        sums <- rowSums(candidates)
        target <- sums[21L]
        expected <- candidates[sums == target, , drop = FALSE]
        expect_gt(nrow(expected), 0L)
        expect_lt(nrow(expected), nrow(candidates))
        for (bound in list(lower, as.character(lower))) {
            expect_equal(
                do.call(permuteGeneral, c(args, list(
                    constraintFun = "sum", comparisonFun = "==",
                    limitConstraints = target, lower = bound, upper = upper
                ))),
                expected
            )
        }
    }
})

test_that("Permutation candidate ranges agree across one and two threads", {
    # The special-constraint worker threshold is 20,000 candidates.
    args <- list(v = 1:10, m = 5, freqs = rep(1:5, 2))
    lower <- 17L
    upper <- lower + 20000L
    candidates <- do.call(permuteGeneral, c(args, list(
        lower = lower, upper = upper
    )))
    sums <- rowSums(candidates)
    target <- sums[10001L]
    rows <- sums == target
    expected <- unname(cbind(candidates[rows, , drop = FALSE], sums[rows]))
    for (threads in 1:2) {
        expect_equal(
            do.call(permuteGeneral, c(args, list(
                constraintFun = "sum", comparisonFun = "==",
                limitConstraints = target, lower = lower, upper = upper,
                keepResults = TRUE, nThreads = threads
            ))),
            expected
        )
    }
})


test_that("Distinct permutation caps stop within the final partition block", {
    for (width in 3:4) {
        target <- if (width == 3L) 8L else 14L
        for (first in 0:1) {
            v <- first:target
            # Independently form each partition's lexicographic permutations.
            parts <- comboGeneral(v, width)
            parts <- parts[rowSums(parts) == target, , drop = FALSE]
            blocks <- lapply(seq_len(nrow(parts)), function(i) {
                permuteGeneral(parts[i, ], width)
            })
            expected <- do.call(rbind, blocks)
            block_size <- factorial(width)
            total <- nrow(expected)
            expect_gt(total, block_size)
            expect_identical(
                permuteGeneral(v, width, constraintFun = "sum",
                               comparisonFun = "==", limitConstraints = target),
                expected
            )
            for (upper in unique(c(1L, 3L, block_size - 1L, block_size,
                                   block_size + 1L, total - 1L, total))) {
                prefix <- expected[seq_len(upper), , drop = FALSE]
                for (keep in c(FALSE, TRUE)) {
                    actual <- permuteGeneral(
                        v, width, constraintFun = "sum", comparisonFun = "==",
                        limitConstraints = target, upper = upper,
                        keepResults = keep
                    )
                    bench <- if (keep) {
                        unname(cbind(prefix, rep(target, upper)))
                    } else prefix
                    expect_identical(actual, bench)
                }
            }
        }
    }
})
