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
