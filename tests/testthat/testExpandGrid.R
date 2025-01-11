context("testing expandGrid")

test_that("expandGrid generates correct output", {
    baseRCheck <- function(myList) {
        testOut <- expandGrid(myList)

        ## Match the original order of the passed vectors
        get_order <- function(df, orig_args) {
            temp <- as.data.frame(
                mapply(function(v, z) match(v, z), df, orig_args)
            )

            do.call(order, temp)
        }

        t <- expand.grid(myList, stringsAsFactors = FALSE)
        t <- t[get_order(t, myList), ]

        if (length(unique(sapply(t, class))) == 1) {
            t <- as.matrix(t)
            dimnames(t) <- dimnames(testOut)
        } else {
            attributes(t) <- attributes(testOut)
        }

        all.equal(t, testOut)
    }

    set.seed(123456789)
    expect_equal(expandGrid(), expand.grid())

    myList <- list(1:5, 2:6, 3:7)
    expect_true(baseRCheck(myList))

    ## shuffle the order. We will follow this pattern below
    expect_true(baseRCheck(lapply(myList, sample)))

    myList <- list(1:5, 2:6, letters[1:4], letters[2:4])
    expect_true(baseRCheck(myList))
    expect_true(baseRCheck(lapply(myList, sample)))

    myList <- list("name1" = 1:5, 2:6)
    expect_true(baseRCheck(myList))
    expect_true(baseRCheck(lapply(myList, sample)))

    myList <- list(1:5, factor(2:6))
    expect_true(baseRCheck(myList))
    expect_true(baseRCheck(lapply(myList, sample)))

    myList <- list(1:5 + 0.1, 2:6 + 0.1)
    expect_true(baseRCheck(myList))
    expect_true(baseRCheck(lapply(myList, sample)))

    myList <- list(1:5 + 0.1, factor(2:6 + 0.1))
    expect_true(baseRCheck(myList))
    expect_true(baseRCheck(lapply(myList, sample)))

    myList <- list(letters[1:5], letters[2:6], letters[3:7])
    expect_true(baseRCheck(myList))
    expect_true(baseRCheck(lapply(myList, sample)))

    myList <- rep(list(c(TRUE, FALSE)), 10)
    expect_true(baseRCheck(myList))
    expect_true(baseRCheck(lapply(myList, sample)))
})

test_that("expandGrid and expandGridSample generates
          correct output with indexing", {

    reset_row_number <- function(df) {
        row.names(df) <- 1:nrow(df)
        df
    }

    set.seed(987654321)

    ## INTEGERS
    myList <- list(1:5, 2:6, 3:7)
    total  <- expandGridCount(myList)
    expect_equal(total, nrow(expand.grid(myList)))

    idx <- sample(total, 5)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        expandGrid(myList)[idx, ]
    )

    ## shuffle the order. We will follow this pattern below
    myList <- lapply(myList, sample)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        expandGrid(myList)[idx, ]
    )

    ## testing upper and lower
    expect_equal(
        expandGrid(myList, lower = min(idx), upper = max(idx)),
        expandGrid(myList)[min(idx):max(idx), ]
    )

    ## RAW
    myList <- list(as.raw(1:5), as.raw(2:6), as.raw(3:7))
    total  <- expandGridCount(myList)
    expect_equal(total, nrow(expand.grid(myList)))

    idx <- sample(total, 5)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        expandGrid(myList)[idx, ]
    )

    expect_equal(
        expandGrid(myList, lower = min(idx), upper = max(idx)),
        expandGrid(myList)[min(idx):max(idx), ]
    )

    myList <- lapply(myList, sample)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        expandGrid(myList)[idx, ]
    )

    ## COMPLEX
    myList <- list(as.complex(1:5 + 1i), as.complex(2:6 + 1i),
                   as.complex(3:7 + 1i))
    total  <- expandGridCount(myList)
    expect_equal(total, nrow(expand.grid(myList)))

    idx <- sample(total, 5)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        expandGrid(myList)[idx, ]
    )

    expect_equal(
        expandGrid(myList, lower = min(idx), upper = max(idx)),
        expandGrid(myList)[min(idx):max(idx), ]
    )

    myList <- lapply(myList, sample)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        expandGrid(myList)[idx, ]
    )

    ## DATA.FRAME
    myList <- list(1:5, as.complex(2:6 + 1i), letters[1:4], letters[2:4])
    total  <- expandGridCount(myList)
    expect_equal(total, nrow(expand.grid(myList)))

    idx <- sample(total, 5)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        reset_row_number(expandGrid(myList)[idx, ])
    )

    expect_equal(
        expandGrid(myList, lower = min(idx), upper = max(idx)),
        reset_row_number(expandGrid(myList)[min(idx):max(idx), ])
    )

    myList <- lapply(myList, sample)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        reset_row_number(expandGrid(myList)[idx, ])
    )

    ## INTEGERS w/ custom names
    myList <- list("name1" = 1:5, 2:6)
    total  <- expandGridCount(myList)
    expect_equal(total, nrow(expand.grid(myList)))

    idx <- sample(total, 5)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        expandGrid(myList)[idx, ]
    )

    expect_equal(
        expandGrid(myList, lower = min(idx), upper = max(idx)),
        expandGrid(myList)[min(idx):max(idx), ]
    )

    myList <- lapply(myList, sample)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        expandGrid(myList)[idx, ]
    )

    ## FACTORS with different levels. Should return a data.frame.
    myList <- list("v1" = factor(1:5), "v2" = factor(2:6))
    total  <- expandGridCount(myList)
    expect_equal(total, nrow(expand.grid(myList)))
    expect_true(is.data.frame(expandGrid(myList)))

    idx  <- sample(total, 3)
    samp <- expandGridSample(myList, sampleVec = idx)
    expect_true(is.data.frame(samp))
    expect_equal(samp, reset_row_number(expandGrid(myList)[idx, ]))

    ## Levels should be the same as the original
    expect_equal(lapply(samp, levels), lapply(myList, levels))

    expect_equal(
        expandGrid(myList, lower = min(idx), upper = max(idx)),
        reset_row_number(expandGrid(myList)[min(idx):max(idx), ])
    )

    myList <- lapply(myList, sample)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        reset_row_number(expandGrid(myList)[idx, ])
    )

    ## FACTORS with same levels. Should return a matrix.
    myList <- list("v1" = factor(1:5, levels = 1:10),
                   "v2" = factor(2:6, levels = 1:10))
    expect_true(is.matrix(expandGrid(myList)))
    expect_true(is.data.frame(expandGrid(myList, return_df = TRUE)))

    expect_equal(levels(expandGrid(myList)), levels(myList$v1))
    expect_true(
        Reduce(
            identical,
            lapply(expandGrid(myList, return_df = TRUE), levels)
        )
    )

    expect_equal(
        Reduce(
            union,
            lapply(expandGrid(myList, return_df = TRUE), levels)
        ),
        levels(myList$v1)
    )

    ## DATA.FRAME w/ integers and factors
    myList <- list(1:5, factor(2:6), as.raw(sample(10, 4)), c(TRUE, FALSE))
    total  <- expandGridCount(myList)
    expect_equal(total, nrow(expand.grid(myList)))

    idx <- sample(total, 5)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        reset_row_number(expandGrid(myList)[idx, ])
    )

    expect_equal(
        expandGrid(myList, lower = min(idx), upper = max(idx)),
        reset_row_number(expandGrid(myList)[min(idx):max(idx), ])
    )

    myList <- lapply(myList, sample)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        reset_row_number(expandGrid(myList)[idx, ])
    )

    ## DOUBLE
    myList <- list(1:5 + 0.1, 2:6 + 0.1)
    total  <- expandGridCount(myList)
    expect_equal(total, nrow(expand.grid(myList)))

    idx <- sample(total, 5)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        expandGrid(myList)[idx, ]
    )

    expect_equal(
        expandGrid(myList, lower = min(idx), upper = max(idx)),
        expandGrid(myList)[min(idx):max(idx), ]
    )

    myList <- lapply(myList, sample)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        expandGrid(myList)[idx, ]
    )

    ## DATA.FRAME w/ doubles and factors
    myList <- list(1:5 + 0.1, factor(2:6 + 0.1))
    total  <- expandGridCount(myList)
    expect_equal(total, nrow(expand.grid(myList)))

    idx <- sample(total, 5)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        reset_row_number(expandGrid(myList)[idx, ])
    )

    expect_equal(
        expandGrid(myList, lower = min(idx), upper = max(idx)),
        reset_row_number(expandGrid(myList)[min(idx):max(idx), ])
    )

    myList <- lapply(myList, sample)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        reset_row_number(expandGrid(myList)[idx, ])
    )

    ## CHARACTER
    myList <- list(letters[1:5], letters[2:6], letters[3:7])
    total  <- expandGridCount(myList)
    expect_equal(total, nrow(expand.grid(myList)))

    idx <- sample(total, 5)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        expandGrid(myList)[idx, ]
    )

    expect_equal(
        expandGrid(myList, lower = min(idx), upper = max(idx)),
        expandGrid(myList)[min(idx):max(idx), ]
    )

    myList <- lapply(myList, sample)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        expandGrid(myList)[idx, ]
    )

    ## BOOLEAN
    myList <- rep(list(c(TRUE, FALSE)), 10)
    total  <- expandGridCount(myList)
    expect_equal(total, nrow(expand.grid(myList)))

    idx <- sample(total, 5)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        expandGrid(myList)[idx, ]
    )

    expect_equal(
        expandGrid(myList, lower = min(idx), upper = max(idx)),
        expandGrid(myList)[min(idx):max(idx), ]
    )

    myList <- lapply(myList, sample)
    expect_equal(
        expandGridSample(myList, sampleVec = idx),
        expandGrid(myList)[idx, ]
    )
})

test_that("expandGrid and expandGridSample generates correct
          output with multiple threads and IsGmp", {

    set.seed(591827364)

    ## INTEGERS
    myList  <- Map(\(x, y) x:y, 10:30, 30:50)
    total   <- expandGridCount(myList)
    last100 <- total - 99

    expect_equal(
        expandGrid(myList, lower = last100),
        expandGridSample(myList, sampleVec = do.call(
            c, lapply(0:99, \(x) last100 + x)
        ))
    )

    capture.output(
        rand100 <- gmp::urand.bigz(
            100, size = gmp::log2.bigz(total), seed = 11111
        )
    )

    expect_equal(
        expandGridSample(myList, sampleVec = rand100),
        expandGridSample(myList, sampleVec = rand100, nThreads = 2)
    )

    ## DOUBLES
    myList  <- Map(\(x, y) x:y + rnorm(1), 10:30, 30:50)
    total   <- expandGridCount(myList)
    last100 <- total - 99

    expect_equal(
        expandGrid(myList, lower = last100),
        expandGridSample(myList, sampleVec = do.call(
            c, lapply(0:99, \(x) last100 + x)
        ))
    )

    capture.output(
        rand100 <- gmp::urand.bigz(
            100, size = gmp::log2.bigz(total), seed = 22222
        )
    )

    expect_equal(
        expandGridSample(myList, sampleVec = rand100),
        expandGridSample(myList, sampleVec = rand100, nThreads = 2)
    )

    ## Non-GMP
    myList  <- Map(\(x, y) x:y + rnorm(1), 10:15, 16:21)
    total   <- expandGridCount(myList)
    last100 <- total - 99

    expect_equal(
        expandGrid(myList, lower = last100),
        expandGridSample(myList, sampleVec = last100 + 0:99)
    )

    ## Test nThreads
    expect_equal(
        expandGrid(myList),
        expandGrid(myList, nThreads = 2)
    )

    rand100 <- sample(total, 100)

    expect_equal(
        expandGridSample(myList, sampleVec = rand100),
        expandGridSample(myList, sampleVec = rand100, nThreads = 2)
    )

    ## BOOLEANS
    myList  <- rep(list(c(TRUE, FALSE)), 60)
    total   <- expandGridCount(myList)
    last100 <- total - 99

    expect_equal(
        expandGrid(myList, lower = last100),
        expandGridSample(myList, sampleVec = do.call(
            c, lapply(0:99, \(x) last100 + x)
        ))
    )

    expect_equal(
        expandGridSample(myList, sampleVec = do.call(
            c, lapply(0:99, \(x) last100 + x)
        )),
        expandGridSample(myList, sampleVec = do.call(
            c, lapply(0:99, \(x) last100 + x)
        ), nThreads = 2)
    )

    ## CHARACTERS
    myList  <- rep(list(letters, LETTERS, state.abb), 4)
    total   <- expandGridCount(myList)
    last100 <- total - 99

    expect_equal(
        expandGrid(myList, lower = last100),
        expandGridSample(myList, sampleVec = do.call(
            c, lapply(0:99, \(x) last100 + x)
        ))
    )

    ## DATA.FRAME
    myList <- rep(
        list(LETTERS, unique(state.region), c(TRUE, FALSE),
        as.raw(1:10), as.complex(10:1 + 1i), 1:100), 4
    )
    total   <- expandGridCount(myList)
    last100 <- total - 99

    expect_equal(
        expandGrid(myList, lower = last100),
        expandGridSample(myList, sampleVec = do.call(
            c, lapply(0:99, \(x) last100 + x)
        ))
    )
})
