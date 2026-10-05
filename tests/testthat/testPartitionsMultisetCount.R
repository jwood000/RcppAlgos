context("testing partitions multiset count")

test_that("Double DP correctness against brute-force enumeration", {

    check_count <- function(x) {
        expect_identical(x$brute, x$algo)
        expect_true(x$check)
    }

    ## Multiset partitions
    check_count(
        partitionsMultisetCount(
            20, 5, freqs = rep(3:1, length.out = 20)
        )
    )

    check_count(
        partitionsMultisetCount(
            30, 7, freqs = rep(c(1, 3, 2, 4), length.out = 30)
        )
    )

    ## Ordinary multiset compositions (no zero)
    check_count(
        compositionsMultisetCount(
            20, 5, freqs = rep(3:1, length.out = 20)
        )
    )

    check_count(
        compositionsMultisetCount(
            40, 7, freqs = rep(5:1, 8)
        )
    )

    ## Non-weak multiset compositions with zero padding
    check_count(
        compositionsMultisetCount(
            0:12, 5,
            freqs = c(3, rep(c(2, 1), 6))
        )
    )

    check_count(
        compositionsMultisetCount(
            0:15, 6,
            freqs = c(4, rep(c(1, 3, 2), 5))
        )
    )

    ## Weak multiset compositions
    check_count(
        compositionsMultisetCount(
            0:10, 5,
            freqs = c(3, rep(c(2, 1), 5)),
            weak = TRUE
        )
    )

    check_count(
        compositionsMultisetCount(
            0:12, 6,
            freqs = c(4, rep(c(1, 2, 3), 4)),
            weak = TRUE
        )
    )

    check_count(
        partitionsMultisetCount(
            1:12, 5,
            freqs = c(3, 2, 1, 4, 2, 1, 3, 1, 2, 1, 2, 1),
            target = 24
        )
    )

    check_count(
        partitionsMultisetCount(
            2:14, 6,
            freqs = c(2, 1, 3, 2, 1, 2, 4, 1, 1, 2, 1, 3, 1),
            target = 38
        )
    )

    check_count(
        compositionsMultisetCount(
            1:10, 5,
            freqs = c(3, 2, 1, 4, 2, 1, 3, 1, 2, 1),
            target = 23
        )
    )

    check_count(
        compositionsMultisetCount(
            1:12, 6,
            freqs = c(2, 3, 1, 2, 4, 1, 2, 1, 3, 1, 2, 1),
            target = 31
        )
    )

    check_count(
        compositionsMultisetCount(
            0:10, 5,
            freqs = c(4, 2, 1, 3, 2, 1, 2, 1, 3, 1, 2),
            target = 18
        )
    )

    check_count(
        compositionsMultisetCount(
            0:12, 6,
            freqs = c(5, 1, 3, 2, 1, 2, 4, 1, 2, 1, 3, 1, 2),
            target = 24
        )
    )

    check_count(
        compositionsMultisetCount(
            0:9, 5,
            freqs = c(3, 2, 1, 3, 2, 1, 2, 1, 3, 1),
            target = 16,
            weak = TRUE
        )
    )

    check_count(
        compositionsMultisetCount(
            0:11, 6,
            freqs = c(4, 1, 3, 2, 1, 2, 3, 1, 2, 1, 2, 1),
            target = 21,
            weak = TRUE
        )
    )

    expect_error(
        compositionsMultisetCount(
            1:8, 4,
            freqs = c(2, 1, 2, 1, 2, 1, 2, 1),
            target = 100
        ),
        "Unexpected PartitionType in PartitionsMultisetCount: NoSolution"
    )

    expect_error(
        compositionsMultisetCount(
            1:8, 4,
            freqs = c(2, 1, 2, 1, 2, 1, 2, 1),
            target = 100,
            checkGmp = TRUE
        ),
        "GMP checker does not support this partition type."
    )
})

test_that("GMP correctness below the double exact-integer boundary", {

    check_gmp <- function(x) {
        expect_lt(x$dbl_algo, 2^53)
        expect_identical(x$gmp_algo, x$dbl_algo)
        expect_true(x$check)
    }

    ## Multiset partitions
    check_gmp(
        partitionsMultisetCount(
            25, 6,
            freqs = rep(c(3, 1, 2), length.out = 25),
            checkGmp = TRUE
        )
    )

    check_gmp(
        partitionsMultisetCount(
            1:15, 6,
            freqs = c(4, 1, 3, 2, 1, 2, 3, 1, 2, 1, 3, 1, 2, 1, 2),
            target = 32,
            checkGmp = TRUE
        )
    )

    ## Ordinary multiset compositions
    check_gmp(
        compositionsMultisetCount(
            30, 6,
            freqs = rep(c(4, 2, 1, 3), length.out = 30),
            checkGmp = TRUE
        )
    )

    check_gmp(
        compositionsMultisetCount(
            1:14, 6,
            freqs = c(3, 2, 1, 4, 2, 1, 3, 1, 2, 1, 3, 2, 1, 2),
            target = 30,
            checkGmp = TRUE
        )
    )

    ## Non-weak zero-padded multiset compositions
    check_gmp(
        compositionsMultisetCount(
            0:12, 6,
            freqs = c(5, 2, 1, 3, 2, 1, 4, 1, 2, 1, 3, 1, 2),
            target = 24,
            checkGmp = TRUE
        )
    )

    ## Weak multiset compositions
    check_gmp(
        compositionsMultisetCount(
            0:11, 6,
            freqs = c(4, 1, 3, 2, 1, 2, 3, 1, 2, 1, 2, 1),
            target = 21,
            weak = TRUE,
            checkGmp = TRUE
        )
    )

    check_gmp(
        partitionsMultisetCount(
            30, 8,
            freqs = c(8, 7, 6, 5, 4, 3, 2, 1, rep(1, 22)),
            checkGmp = TRUE
        )
    )

    check_gmp(
        compositionsMultisetCount(
            35, 8,
            freqs = rep(c(5, 4, 3, 2, 1), 7),
            checkGmp = TRUE
        )
    )

    ##### ***** Large exact counts near the double integer boundary ***** #####
    check_gmp(
        partitionsMultisetCount(
            392, 25,
            freqs = c(20, 7, rep(10:1, 39)),
            checkGmp = TRUE
        )
    )

    check_gmp(
        partitionsMultisetCount(
            300, 20,
            freqs = rep(10:1, 30),
            target = 439,
            checkGmp = TRUE
        )
    )

    check_gmp(
        partitionsMultisetCount(
            0:300, 20,
            freqs = c(5, rep(10:1, 30)),
            target = 418,
            checkGmp = TRUE
        )
    )

    set.seed(123456789)
    fr <- sample(9, 300, TRUE)

    check_gmp(
        partitionsMultisetCount(
            300, 15,
            freqs = fr,
            target = 550,
            checkGmp = TRUE
        )
    )

    check_gmp(
        compositionsMultisetCount(
            150, 11, freqs = fr[1:150],
            checkGmp = TRUE
        )
    )

    check_gmp(
        compositionsMultisetCount(
            0:150, 11, freqs = fr[1:151],
            checkGmp = TRUE
        )
    )

    check_gmp(
        compositionsMultisetCount(
            0:150, 11, freqs = fr[1:151],
            weak = TRUE,
            checkGmp = TRUE
        )
    )
})

test_that("GMP-only multiset counts satisfy exact counting invariants", {

    ## Multiplicity bounds cannot bind when every frequency >= m.
    ## Therefore bounded multiset counts must equal repetition counts.

    expect_identical(
        partitionsCount(
            1000, 50,
            freqs = rep(50, 1000)
        ),
        partitionsCount(
            1000, 50,
            repetition = TRUE
        )
    )

    expect_identical(
        compositionsCount(
            500, 20,
            freqs = rep(20, 500)
        ),
        compositionsCount(
            500, 20,
            repetition = TRUE
        )
    )

    expect_identical(
        compositionsCount(
            0:500, 14,
            freqs = rep(20, 501)
        ),
        compositionsCount(
            0:500, 14,
            repetition = TRUE
        )
    )

    expect_identical(
        compositionsCount(
            0:500, 14,
            freqs = rep(20, 501),
            weak = TRUE
        ),
        compositionsCount(
            0:500, 14,
            repetition = TRUE,
            weak = TRUE
        )
    )

    expect_identical(
        partitionsCount(
            1:500, 30,
            freqs = rep(30, 500),
            target = 900
        ),
        partitionsCount(
            1:500, 30,
            repetition = TRUE,
            target = 900
        )
    )

    expect_identical(
        partitionsCount(
            0:500, 30,
            freqs = rep(30, 501),
            target = 900
        ),
        partitionsCount(
            0:500, 30,
            repetition = TRUE,
            target = 900
        )
    )

    expect_identical(
        compositionsCount(
            1:500, 30,
            freqs = rep(30, 500),
            target = 900
        ),
        compositionsCount(
            1:500, 30,
            repetition = TRUE,
            target = 900
        )
    )

    expect_identical(
        compositionsCount(
            0:500, 30,
            freqs = rep(30, 501),
            target = 900
        ),
        compositionsCount(
            0:500, 30,
            repetition = TRUE,
            target = 900
        )
    )

    expect_identical(
        compositionsCount(
            0:500, 30,
            freqs = rep(30, 501),
            target = 900,
            weak = TRUE
        ),
        compositionsCount(
            0:500, 30,
            repetition = TRUE,
            target = 900,
            weak = TRUE
        )
    )

    ## monotonicity in the multiplicities
    x <- partitionsCount(
        1000, 50,
        freqs = rep(10, 1000)
    )

    y <- partitionsCount(
        1000, 50,
        freqs = rep(20, 1000)
    )

    z <- partitionsCount(
        1000, 50,
        repetition = TRUE
    )

    expect_s3_class(x, "bigz")
    expect_lt(x, y)
    expect_lt(y, z)

    x <- compositionsCount(
        500, 20,
        freqs = rep(5, 500)
    )

    y <- compositionsCount(
        500, 20,
        freqs = rep(10, 500)
    )

    z <- compositionsCount(
        500, 20,
        repetition = TRUE
    )

    expect_s3_class(x, "bigz")
    expect_lt(x, y)
    expect_lt(y, z)
})

test_that("Zero-semantics tests", {

    check_count <- function(x) {
        expect_identical(x$brute, x$algo)
        expect_true(x$check)
    }

    ## No zero in the input: all multiplicities belong to positive values.
    check_count(
        compositionsMultisetCount(
            1:10, 5,
            freqs = c(3, 2, 1, 4, 2, 1, 3, 1, 2, 1)
        )
    )

    ## Zero present, non-weak: zero is padding and must not participate
    ## in the positive-part DP count.
    check_count(
        compositionsMultisetCount(
            0:10, 5,
            freqs = c(4, 2, 1, 3, 2, 1, 2, 1, 3, 1, 2)
        )
    )

    check_count(
        compositionsMultisetCount(
            0:12, 6,
            freqs = c(5, 1, 3, 2, 1, 2, 4, 1, 2, 1, 3, 1, 2),
            target = 24
        )
    )

    ## Zero present, weak: zero is an actual part and its multiplicity
    ## must remain part of the counting problem.
    check_count(
        compositionsMultisetCount(
            0:10, 5,
            freqs = c(4, 2, 1, 3, 2, 1, 2, 1, 3, 1, 2),
            weak = TRUE
        )
    )

    check_count(
        compositionsMultisetCount(
            0:12, 6,
            freqs = c(5, 1, 3, 2, 1, 2, 4, 1, 2, 1, 3, 1, 2),
            target = 24,
            weak = TRUE
        )
    )

    nonweak <- compositionsCount(
        0:8, 4,
        freqs = c(3, 2, 2, 1, 2, 1, 2, 1, 2),
        target = 10
    )

    weak <- compositionsCount(
        0:8, 4,
        freqs = c(3, 2, 2, 1, 2, 1, 2, 1, 2),
        target = 10,
        weak = TRUE
    )

    expect_true(nonweak != weak)

    expect_identical(
        compositionsDesign(
            0:10, 5,
            freqs = rep(1, 11),
            weak = TRUE
        )$partition_type,
        "CompDistinctWeak"
    )

    expect_identical(
        compositionsDesign(
            0:10, 5,
            freqs = rep(1, 11)
        )$partition_type,
        "CompDistinctZero"
    )

    x <- compositionsCount(
        0:10, 6,
        freqs = c(6, rep(2, 10)),
        target = 18
    )

    y <- compositionsCount(
        0:10, 6,
        freqs = c(20, rep(2, 10)),
        target = 18
    )

    expect_identical(x, y)

    x <- compositionsCount(
        0:10, 6,
        freqs = c(1, rep(2, 10)),
        target = 18,
        weak = TRUE
    )

    y <- compositionsCount(
        0:10, 6,
        freqs = c(6, rep(2, 10)),
        target = 18,
        weak = TRUE
    )

    expect_lt(x, y)

    x <- compositionsCount(
        0:10, 6,
        freqs = c(4, rep(2, 10)),
        target = 18
    )

    y <- compositionsCount(
        0:10, 6,
        freqs = c(20, rep(2, 10)),
        target = 18
    )

    expect_identical(x, y)

    z <- compositionsCount(
        0:10, 6,
        freqs = c(3, rep(2, 10)),
        target = 18
    )

    expect_lt(z, x)
})

test_that("Length-summing behavior for non-weak compositions", {

    ## Non-weak zero-padded compositions may contain fewer than m
    ## positive parts. The count should therefore equal the sum of
    ## the exact positive-length counts over the admissible range.

    zero_pad_check <- function(v, m, fr, tar) {

        stopifnot(v[1L] == 0)

        x <- compositionsCount(
            v, m, freqs = fr, target = tar
        )

        ## This will be the same result as compositions
        first <- partitionsGeneral(
            v, m, freqs = fr, target = tar, upper = 1
        )[1, ]

        strtLen <- m - which(first > 0)[1L] + 1L

        y <- sum(
            do.call(
                c,
                lapply(seq.int(strtLen, m), \(k) {
                    res <- compositionsCount(
                        v[-1L], k,
                        freqs = fr[-1L],
                        target = tar
                    )

                    if (inherits(x, "bigz")) {
                        gmp::as.bigz(res)
                    } else {
                        as.numeric(res)
                    }
                })
            )
        )

        expect_equal(x, y)
    }

    zero_pad_check(0:10, 6, c(6, rep(3, 10)), 18)

    zero_pad_check(
        0:100, 10, c(10, rep(c(5, rep(2, 4), 1, rep(3, 4)), 10)), 220
    )

    zero_pad_check(
        0:100, 10, c(10, rep(c(5, rep(2, 4), 1, rep(3, 4)), 10)), 420
    )

    zero_pad_check(0:18, 7, c(3, rep(3:1, 6)), 35)

    ## strtLen = 1
    zero_pad_check(
        0:20, 8,
        c(8, rep(8, 20)),
        20
    )

    ## strtLen = 2
    zero_pad_check(
        0:20, 8,
        c(8, rep(8, 20)),
        39
    )

    ## strtLen = 4
    zero_pad_check(
        0:20, 8,
        c(8, rep(8, 20)),
        75
    )

    ## strtLen = m -- only one DP row contributes
    zero_pad_check(
        0:20, 8,
        c(8, rep(8, 20)),
        155
    )

    zero_pad_check(
        0:10, 6,
        c(
            6,
            rep(2, 7),  # values 1:7
            1, 1, 1     # values 8:10
        ),
        28
    )

    zero_pad_check(
        0:20, 8,
        c(2, rep(5, 20)),
        40
    )

    zero_pad_check(
        0:20, 8,
        c(8, rep(5, 20)),
        40
    )

    zero_pad_check(
        0:300, 20,
        c(20, rep(10:1, 30)),
        418
    )
})

test_that("Mapped-domain / allowed behavior", {

    fr <- c(5, 1, 3, 2, 6, 1, 4, 2, 1, 5)

    ## Positive partition mapping:
    ## shifting every value by c shifts the target by c * m.
    expect_identical(
        partitionsCount(
            1:10, 5,
            freqs = fr,
            target = 28
        ),
        partitionsCount(
            11:20, 5,
            freqs = fr,
            target = 28 + 10 * 5
        )
    )

    ## Same invariant for compositions.
    expect_identical(
        compositionsCount(
            1:10, 5,
            freqs = fr,
            target = 28
        ),
        compositionsCount(
            11:20, 5,
            freqs = fr,
            target = 28 + 10 * 5
        )
    )

    check_count <- function(x) {
        expect_identical(x$brute, x$algo)
        expect_true(x$check)
    }

    check_count(
        partitionsMultisetCount(
            10:24, 6,
            freqs = c(7, 1, 4, 2, 1, 5, 2, 6, 1, 3, 2, 1, 4, 2, 3),
            target = 93
        )
    )

    check_count(
        compositionsMultisetCount(
            8:20, 6,
            freqs = c(6, 1, 3, 2, 5, 1, 4, 2, 1, 6, 2, 3, 1),
            target = 79
        )
    )

    fr1 <- c(6, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    fr2 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 6)

    x <- compositionsCount(
        1:10, 5,
        freqs = fr1,
        target = 15
    )

    y <- compositionsCount(
        1:10, 5,
        freqs = fr2,
        target = 15
    )

    expect_false(identical(x, y))

    check_count(
        compositionsMultisetCount(
            0:12, 6,
            freqs = c(
                5,                 # zero
                7, 1, 4, 2, 1, 5, 2, 6, 1, 3, 2, 1
            ),
            target = 31
        )
    )

    fr <- rep(c(7, 2, 5, 1, 4, 3), length.out = 100)

    x <- compositionsCount(
        1:100, 15,
        freqs = fr,
        target = 600
    )

    y <- compositionsCount(
        101:200, 15,
        freqs = fr,
        target = 600 + 100 * 15
    )

    expect_identical(x, y)
})

test_that("Boundary and degenerate counts", {

    ## Impossible target
    expect_identical(
        partitionsCount(
            1:10, 5,
            freqs = rep(3, 10),
            target = 1000
        ),
        0L
    )

    expect_identical(
        compositionsCount(
            1:10, 5,
            freqs = rep(3, 10),
            target = 1000
        ),
        0L
    )

    ## m = 1: exactly one usable value hits the target
    expect_identical(
        partitionsCount(
            1:10, 1,
            freqs = rep(3, 10),
            target = 7
        ),
        1L
    )

    expect_identical(
        compositionsCount(
            1:10, 1,
            freqs = rep(3, 10),
            target = 7
        ),
        1L
    )

    ## m = 1: target not present
    expect_identical(
        partitionsCount(
            c(2, 4, 6, 8), 1,
            freqs = c(2, 3, 1, 4),
            target = 7
        ),
        0L
    )

    ## Only one partition is possible
    expect_identical(
        partitionsCount(
            1:10, 5,
            freqs = rep(5, 10),
            target = 5
        ),
        1L
    )

    ## Same values for compositions: only (1,1,1,1,1)
    expect_identical(
        compositionsCount(
            1:10, 5,
            freqs = rep(5, 10),
            target = 5
        ),
        1L
    )

    ## Insufficient multiplicity kills an otherwise obvious solution
    expect_identical(
        partitionsCount(
            1:10, 5,
            freqs = c(4, rep(5, 9)),
            target = 5
        ),
        0L
    )

    expect_identical(
        compositionsCount(
            1:10, 5,
            freqs = c(4, rep(5, 9)),
            target = 5
        ),
        0L
    )

    ## Only one way to hit the maximum possible sum
    expect_identical(
        partitionsCount(
            1:10, 5,
            freqs = rep(5, 10),
            target = 50
        ),
        1L
    )

    expect_identical(
        compositionsCount(
            1:10, 5,
            freqs = rep(5, 10),
            target = 50
        ),
        1L
    )

    ## Not enough total elements to form width m
    expect_error(
        partitionsCount(
            1:4, 6,
            freqs = c(1, 1, 1, 1),
            target = 10
        ),
        "m must be less than or equal to the length of v"
    )

    expect_identical(
        partitionsCount(
            1:4, 7,
            freqs = c(2, 1, 2, 1),
            target = 10
        ),
        0L
    )

    expect_identical(
        compositionsCount(
            0:10, 6,
            freqs = c(6, rep(2, 10)),
            target = 0
        ),
        1L
    )

    expect_identical(
        compositionsCount(
            0:10, 6,
            freqs = c(5, rep(2, 10)),
            target = 0
        ),
        0L
    )
})

test_that("Design routing / PartitionType selection", {

    ## Ordinary multiset partition
    expect_identical(
        partitionsDesign(
            1:10, 5,
            freqs = c(3, 2, 1, 4, 2, 1, 3, 1, 2, 1)
        )$partition_type,
        "Multiset"
    )

    ## Multiset composition, no zero
    expect_identical(
        compositionsDesign(
            1:10, 5,
            freqs = c(3, 2, 1, 4, 2, 1, 3, 1, 2, 1)
        )$partition_type,
        "CompMultiset"
    )

    ## Zero present, non-weak:
    ## zero is internal padding rather than an actual weak-composition part.
    expect_identical(
        compositionsDesign(
            0:10, 5,
            freqs = c(4, 2, 1, 3, 2, 1, 2, 1, 3, 1, 2)
        )$partition_type,
        "CompMultisetZero"
    )

    expect_identical(
        compositionsDesign(
            0:10, 5,
            freqs = rep(1, 11),
            target = 10
        )$partition_type,
        "CompDistinctZero"
    )

    expect_identical(
        compositionsDesign(
            0:10, 5,
            freqs = rep(1, 11),
            target = 10,
            weak = TRUE
        )$partition_type,
        "CompDistinctWeak"
    )

    ## Zero present, weak:
    ## zero is an actual returned part.
    expect_identical(
        compositionsDesign(
            0:10, 5,
            freqs = c(4, 2, 1, 3, 2, 1, 2, 1, 3, 1, 2),
            weak = TRUE
        )$partition_type,
        "CompMultisetWeak"
    )

    expect_identical(
        partitionsDesign(
            0:10, 5,
            freqs = rep(1, 11),
            target = 10
        )$partition_type,
        "DistinctOneZero"
    )

    ## All multiplicities one route to distinct algorithms,
    ## not the multiset counter.

    expect_identical(
        partitionsDesign(
            1:10, 5,
            freqs = rep(1, 10)
        )$partition_type,
        "NoSolution"
    )

    expect_identical(
        partitionsDesign(
            1:10, 5,
            freqs = rep(1, 10),
            target = 15
        )$partition_type,
        "DistinctNoZero"
    )

    expect_identical(
        compositionsDesign(
            1:10, 5,
            freqs = rep(1, 10)
        )$partition_type,
        "NoSolution"
    )

    expect_identical(
        compositionsDesign(
            1:10, 5,
            freqs = rep(1, 10),
            target = 15
        )$partition_type,
        "CompDistinctNoZero"
    )

    expect_identical(
        compositionsDesign(
            0:10, 5,
            freqs = rep(1, 11)
        )$partition_type,
        "CompDistinctZero"
    )

    expect_identical(
        compositionsDesign(
            0:10, 5,
            freqs = rep(1, 11),
            weak = TRUE
        )$partition_type,
        "CompDistinctWeak"
    )

    nonweak <- compositionsDesign(
        0:10, 5,
        freqs = c(3, 2, 1, 3, 2, 1, 2, 1, 3, 1, 2)
    )

    weak <- compositionsDesign(
        0:10, 5,
        freqs = c(3, 2, 1, 3, 2, 1, 2, 1, 3, 1, 2),
        weak = TRUE
    )

    expect_identical(nonweak$partition_type, "CompMultisetZero")
    expect_identical(weak$partition_type, "CompMultisetWeak")

    expect_identical(
        partitionsDesign(
            1:100, 10,
            freqs = rep(10, 100)
        )$partition_type,
        "Multiset"
    )

    expect_identical(
        partitionsDesign(
            1:100, 10,
            repetition = TRUE
        )$partition_type,
        "RepNoZero"
    )

    expect_identical(
        compositionsDesign(
            1:100, 10,
            freqs = rep(10, 100)
        )$partition_type,
        "CompMultiset"
    )

    expect_identical(
        permutePartsDesign(
            1:10, 5,
            freqs = c(3, 2, 1, 4, 2, 1, 3, 1, 2, 1)
        )$partition_type,
        "PrmMultiset"
    )

    expect_identical(
        partitionsDesign(
            1:8, 3,
            freqs = rep(1, 8)
        )$partition_type,
        "DistinctNoZero"
    )

    x <- partitionsDesign(
        1:8, 4,
        freqs = c(2, rep(1, 7))
    )

    y <- partitionsDesign(
        0:4, 4,
        freqs = c(2, rep(1, 4))
    )

    expect_identical(x$partition_type, y$partition_type)
    expect_identical(x$mapped_target, y$mapped_target)
    expect_identical(x$first_index_vector, y$first_index_vector)
})

test_that("Generation count consistency", {

    check_partition_generation <- function(v, m, fr, tar) {

        cnt <- partitionsCount(
            v, m,
            freqs = fr,
            target = tar
        )

        gen <- partitionsGeneral(
            v, m,
            freqs = fr,
            target = tar
        )

        expect_identical(nrow(gen), as.integer(cnt))
    }

    ## Ordinary genuine multiset
    check_partition_generation(
        1:12, 5,
        c(3, 2, 1, 4, 2, 1, 3, 1, 2, 1, 2, 1),
        24
    )

    ## Uneven multiplicities / larger target
    check_partition_generation(
        1:20, 7,
        rep(c(4, 1, 3, 2), length.out = 20),
        52
    )

    ## Mapping-heavy case
    check_partition_generation(
        10:25, 6,
        c(5, 1, 3, 2, 4, 1, 2, 5, 1, 3, 2, 1, 4, 2, 3, 1),
        97
    )

    ## Interesting routing case: maps to DistinctMZ
    check_partition_generation(
        1:8, 4,
        c(2, rep(1, 7)),
        8
    )

    ## Tight multiplicity constraints
    check_partition_generation(
        1:10, 6,
        c(1, 2, 1, 3, 1, 2, 1, 2, 1, 3),
        27
    )

    check_composition_generation <- function(v, m, fr, tar, weak = FALSE) {

        cnt <- compositionsCount(
            v, m,
            freqs = fr,
            target = tar,
            weak = weak
        )

        perms_parts <- if (weak || !(0 %in% v)) {
            ## Generate valid multiset partitions first, then permute
            ## each one using the existing permutation machinery.
            nrow(
                permuteGeneral(
                    v, m,
                    freqs = fr,
                    constraintFun = "sum",
                    comparisonFun = "==",
                    limitConstraints = tar
                )
            )
        } else {
            first <- partitionsGeneral(
                v, m,
                freqs = fr,
                target = tar,
                upper = 1
            )[1, ]

            strtLen <- m - which(first > 0)[1L] + 1L

            sum(
                sapply(
                    seq.int(strtLen, m), \(k) {
                        nrow(
                            permuteGeneral(
                                v[-1L], k,
                                freqs = fr[-1L],
                                constraintFun = "sum",
                                comparisonFun = "==",
                                limitConstraints = tar
                            )
                        )
                    }
                )
            )
        }

        expect_equal(perms_parts, as.integer(cnt))
    }

    ## Ordinary multiset compositions
    check_composition_generation(
        1:10, 5,
        c(3, 2, 1, 4, 2, 1, 3, 1, 2, 1),
        23
    )

    ## More uneven multiplicities
    check_composition_generation(
        1:12, 6,
        c(4, 1, 3, 2, 1, 5, 2, 1, 3, 2, 1, 4),
        31
    )

    ## Zero present, non-weak
    check_composition_generation(
        0:10, 6,
        c(4, 2, 1, 3, 2, 1, 2, 1, 3, 1, 2),
        18
    )

    ## Same underlying problem, but weak
    check_composition_generation(
        0:10, 6,
        c(4, 2, 1, 3, 2, 1, 2, 1, 3, 1, 2),
        18,
        weak = TRUE
    )

    ## Larger target with several admissible positive lengths
    check_composition_generation(
        0:15, 7,
        c(5, rep(c(3, 1, 2), 5)),
        38
    )

    ## Positive-only case with many repeated-part possibilities
    check_composition_generation(
        1:15, 7,
        rep(c(4, 2, 3, 1, 2), 3),
        42
    )

    ## easy to verify by hand
    check_partition_generation(
        1:5, 3,
        c(2, 1, 2, 1, 2),
        7
    )

    check_composition_generation(
        1:5, 3,
        c(2, 1, 2, 1, 2),
        7
    )
})

test_that("Limited-output generation", {

    check_limited <- function(v, m, fr, tar, upper) {

        full <- partitionsGeneral(
            v, m,
            freqs = fr,
            target = tar
        )

        limited <- partitionsGeneral(
            v, m,
            freqs = fr,
            target = tar,
            upper = upper
        )

        n_expected <- min(upper, nrow(full))

        expect_equal(nrow(limited), n_expected)

        if (n_expected > 0L) {
            expect_equal(
                limited,
                full[seq_len(n_expected), , drop = FALSE]
            )
        }
    }

    ## upper = 1
    check_limited(
        1:12, 5,
        c(3, 2, 1, 4, 2, 1, 3, 1, 2, 1, 2, 1),
        24,
        1
    )

    ## Small prefix from a reasonably populated result set
    check_limited(
        1:20, 7,
        rep(c(4, 1, 3, 2), length.out = 20),
        52,
        3
    )

    check_limited(
        1:20, 7,
        rep(c(4, 1, 3, 2), length.out = 20),
        52,
        7
    )

    ## Mapping-heavy input
    check_limited(
        10:25, 6,
        c(5, 1, 3, 2, 4, 1, 2, 5, 1, 3, 2, 1, 4, 2, 3, 1),
        97,
        4
    )

    ## DistinctMZ mapped path
    check_limited(
        1:8, 4,
        c(2, rep(1, 7)),
        14,
        1
    )

    ## Tight multiplicities
    check_limited(
        1:10, 6,
        c(1, 2, 1, 3, 1, 2, 1, 2, 1, 3),
        27,
        2
    )

    v <- 1:12
    m <- 5
    fr <- c(3, 2, 1, 4, 2, 1, 3, 1, 2, 1, 2, 1)
    tar <- 24

    cnt <- partitionsCount(
        v, m,
        freqs = fr,
        target = tar
    )

    full <- partitionsGeneral(
        v, m,
        freqs = fr,
        target = tar
    )

    limited <- partitionsGeneral(
        v, m,
        freqs = fr,
        target = tar,
        upper = cnt
    )

    expect_identical(limited, full)

    expect_error(
        partitionsGeneral(
            v, m,
            freqs = fr,
            target = tar,
            upper = cnt + 10L
        ),
        "bounds cannot exceed the maximum number of possible results"
    )
})
test_that("Permutation/multiset generation", {

    check_prm_multiset <- function(v, m, fr, tar) {

        cnt <- compositionsCount(
            v, m,
            freqs = fr,
            target = tar,
            weak = TRUE
        )

        gen <- permuteGeneral(
            v, m,
            freqs = fr,
            constraintFun = "sum",
            comparisonFun = "==",
            limitConstraints = tar
        )

        expect_identical(nrow(gen), as.integer(cnt))
        expect_identical(nrow(unique(gen)), nrow(gen))
        expect_true(all(rowSums(gen) == tar))
    }

    ## Tiny, human-auditable case
    check_prm_multiset(
        1:5, 3,
        c(2, 1, 2, 1, 2),
        7
    )

    ## Ordinary multiset composition
    check_prm_multiset(
        1:10, 5,
        c(3, 2, 1, 4, 2, 1, 3, 1, 2, 1),
        23
    )

    ## Uneven multiplicities
    check_prm_multiset(
        1:12, 6,
        c(5, 1, 3, 2, 1, 4, 2, 1, 3, 2, 1, 5),
        31
    )

    ## More repeated-part opportunities
    check_prm_multiset(
        1:15, 7,
        rep(c(4, 2, 3, 1, 2), 3),
        42
    )

    ## Weak. Non-weak checked in 'Generation count consistency' section
    check_prm_multiset(
        0:10, 6,
        c(4, 2, 1, 3, 2, 1, 2, 1, 3, 1, 2),
        18
    )
})
test_that("Unsupported ranking/sampling behavior", {

    expect_multiset_unsupported <- function(expr, ranking = FALSE) {
        myCase <- if (ranking) "Ranking" else "Sampling"
        expect_error(
            expr,
            paste(myCase, "not available for this case.")
        )
    }

    fr <- c(3, 2, 1, 4, 2, 1, 3, 1, 2, 1)

    ## Genuine multiset partition
    expect_identical(
        partitionsDesign(
            1:10, 5,
            freqs = fr,
            target = 20
        )$partition_type,
        "Multiset"
    )

    expect_multiset_unsupported(
        partitionsSample(
            1:10, 5,
            freqs = fr,
            target = 20,
            n = 5
        )
    )

    x <- partitionsGeneral(
        1:10, 5,
        freqs = fr,
        target = 20,
        upper = 1
    )

    expect_multiset_unsupported(
        partitionsRank(
            x,
            v = 1:10,
            freqs = fr,
            target = 20,
            n = 5
        ),
        ranking = TRUE
    )

    ## Positive multiset composition
    expect_identical(
        compositionsDesign(
            1:10, 5,
            freqs = fr,
            target = 20
        )$partition_type,
        "CompMultiset"
    )

    expect_multiset_unsupported(
        compositionsSample(
            1:10, 5,
            freqs = fr,
            target = 20,
            n = 5
        )
    )

    expect_multiset_unsupported(
        compositionsRank(
            x,
            v = 1:10,
            freqs = fr,
            target = 20,
            n = 5
        ),
        ranking = TRUE
    )

    ## Zero-padded non-weak multiset composition
    fr_zero <- c(4, 2, 1, 3, 2, 1, 2, 1, 3, 1, 2)

    expect_identical(
        compositionsDesign(
            0:10, 6,
            freqs = fr_zero,
            target = 18
        )$partition_type,
        "CompMultisetZero"
    )

    x <- partitionsGeneral(
        0:10, 6,
        freqs = fr_zero,
        target = 18,
        upper = 1
    )

    expect_multiset_unsupported(
        compositionsSample(
            0:10, 6,
            freqs = fr_zero,
            target = 18,
            n = 5
        )
    )

    expect_multiset_unsupported(
        compositionsRank(
            x,
            v = 0:10,
            freqs = fr_zero,
            target = 18,
            n = 5
        ),
        ranking = TRUE
    )

    ## Weak multiset composition
    expect_identical(
        compositionsDesign(
            0:10, 6,
            freqs = fr_zero,
            target = 18,
            weak = TRUE
        )$partition_type,
        "CompMultisetWeak"
    )

    expect_multiset_unsupported(
        compositionsSample(
            0:10, 6,
            freqs = fr_zero,
            target = 18,
            weak = TRUE,
            n = 5
        )
    )

    expect_multiset_unsupported(
        compositionsRank(
            x,
            v = 0:10,
            freqs = fr_zero,
            target = 18,
            weak = TRUE,
            n = 5
        ),
        ranking = TRUE
    )
})
