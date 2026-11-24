test_that("GMP rank vectors reject zeros in any position", {
    for (indices in list(c(0, 1), c(1, 0), c(1, 0, 0))) {
        expect_error(
            comboSample(100, 20, sampleVec = gmp::as.bigz(indices)),
            "Each element in sampleVec must be a positive number"
        )
    }
})

test_that("GMP scalar bounds and seeds reject zero", {
    zero <- gmp::as.bigz(0)
    expect_error(comboGeneral(100, 20, lower = zero, upper = 1),
                 "lower must be a positive number")
    expect_error(comboGeneral(100, 20, upper = zero),
                 "upper must be a positive number")
    # The public wrapper also passes seed through base::set.seed().
    expect_error(comboSample(100, 20, n = 1, seed = 0),
                 "seed must be a positive number")
})

test_that("Positive GMP ranks, bounds and seeds retain exact conversions", {
    indices <- c("1", "9007199254740993", "9007199254740994")
    expect_identical(
        comboSample(100, 20, sampleVec = gmp::as.bigz(indices)),
        comboSample(100, 20, sampleVec = indices)
    )
    expect_identical(
        comboGeneral(100, 20, lower = gmp::as.bigz(indices[2]),
                     upper = gmp::as.bigz(indices[3])),
        comboGeneral(100, 20, lower = indices[2], upper = indices[3])
    )
    expect_identical(
        comboSample(100, 20, n = 3, seed = 42),
        comboSample(100, 20, n = 3, seed = "42")
    )
})

test_that("GMP decoding advances past leading zero magnitude words", {
    # Encode ranks 3 and 4, with an extra leading zero word in rank 3.
    # Its encoded size exceeds the size of the normalized GMP value.
    padded <- writeBin(c(2L, 2L, 1L, 0L, 3L, 1L, 1L, 4L), raw(), size = 4L)
    scalar <- writeBin(c(1L, 2L, 1L, 0L, 3L), raw(), size = 4L)

    for (case in list(c(5L, 3L), c(100L, 20L))) {
        expected <- comboGeneral(case[1], case[2], lower = 3L, upper = 4L)
        expect_identical(comboSample(case[1], case[2], sampleVec = padded),
                         expected)
        iterator <- comboIter(case[1], case[2])
        expect_identical(iterator[[padded]], expected)
        expect_identical(comboGeneral(case[1], case[2], lower = scalar,
                                      upper = 4L), expected)
    }
})
