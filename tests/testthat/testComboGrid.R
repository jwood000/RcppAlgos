context("testing comboGrid")

test_that("comboGrid generates correct output", {
    bruteCheck <- function(myList, rep = TRUE) {
        testOut <- comboGrid(myList, repetition = rep)

        t  <- expand.grid(myList, stringsAsFactors = FALSE)
        t  <- t[do.call(order, t), ]
        t1 <- apply(t, 1, function(x) paste0(sort(x), collapse = ""))
        t2 <- t[!duplicated(t1), ]

        if (!rep) {
            t2 <- t2[!apply(t2, 1, function(x) any(duplicated(x))), ]
        }

        if (length(unique(sapply(t2, class))) == 1) {
            t2 <- as.matrix(t2)
            dimnames(t2) <- dimnames(testOut)
        } else {
            attributes(t2) <- attributes(testOut)
        }

        all.equal(t2, testOut)
    }

    expect_equal(comboGrid(), expand.grid())

    myList <- list(1:5, 2:6, 3:7)
    expect_true(bruteCheck(myList, rep = TRUE))
    expect_true(bruteCheck(myList, rep = FALSE))

    myList <- list(1:5, 2:6, letters[1:4], letters[2:4])
    expect_true(bruteCheck(myList, rep = TRUE))
    expect_true(bruteCheck(myList, rep = FALSE))

    myList <- list("name1" = 1:5, 2:6)
    expect_true(bruteCheck(myList, rep = TRUE))

    myList <- list(1:5, factor(2:6))
    expect_true(bruteCheck(myList, rep = TRUE))

    myList <- list(1:5 + 0.1, 2:6 + 0.1)
    expect_true(bruteCheck(myList, rep = TRUE))

    myList <- list(1:5 + 0.1, factor(2:6 + 0.1))
    expect_true(bruteCheck(myList, rep = TRUE))

    myList <- list(letters[1:5], letters[2:6], letters[3:7])
    expect_true(bruteCheck(myList, rep = TRUE))
    expect_true(bruteCheck(myList, rep = FALSE))

    myList <- rep(list(c(T, F)), 10)
    expect_true(bruteCheck(myList, rep = TRUE))

    ## Huge test... This will trigger mpz_t keys in ComboCartesian.cpp
    pools <- list(c(1, 10, 14, 6),
                  c(7, 2, 4, 8, 3, 11, 12),
                  c(11, 3, 13, 4, 15, 8, 6, 5),
                  c(10, 1, 3, 2, 9, 5,  7),
                  c(1, 5, 10, 3, 8, 14),
                  c(15, 3, 7, 10, 4, 5, 8, 6),
                  c(14, 9, 11, 15),
                  c(7, 6, 13, 14, 10, 11, 9, 4),
                  c(6,  3,  2, 14,  7, 12,  9),
                  c(6, 11,  2,  5, 15,  7),
                  16:19, 20:23, 24:28,
                  c(16, 20, 28, 1))

    hugeTest = comboGrid(pools, repetition = FALSE)
    expect_equal(ncol(hugeTest), length(pools))

    ## If NULL is an input, we return an empty named data.frame
    expect_equal(dim(comboGrid(NA, NA, 1:10, NA, 1:5, NULL)),
                 c(0, 5))

    ## If NA is an input, we tack on a column of NA's
    df <- comboGrid(NA, 1:10)
    expect_equal(df$Var1, rep(NA, 10))

    ## Test factors where there are more levels than elements
    df <- comboGrid(factor(letters[1:5], levels = letters),
                    factor(letters[c(16, 22)]),
                    letters[4:9],
                    c(7, 17))

    expect_equal(levels(df$Var1), letters)
    expect_equal(levels(df$Var2), letters[c(16, 22)])
    expect_false(is.factor(df$Var3))
    expect_false(is.factor(df$Var4))
})
