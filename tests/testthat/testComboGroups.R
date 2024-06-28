context("testing comboGroup")

test_that("comboGroups produces correct results when groups are uniform", {

    ## See the excellent answer by user @Cole here:
    ##       https://stackoverflow.com/a/57847300/4408538
    group <- function(input, step) {
        len <- length(input)
        combination[1, step] <<- input[1]

        for (i1 in 2:(len - 1)) {
            combination[2, step] <<- input[i1]

            for (i2 in (i1 + 1):(len - 0)) {
                combination[3, step] <<- input[i2]

                if (step == m) {
                    print(z);
                    result[z, ,] <<- combination
                    z <<- z + 1
                } else {
                    rest <- setdiff(input, input[c(i1, i2, 1)])

                    #recursive if there are still additional possibilities
                    group(rest, step + 1)
                }
            }
        }
    }

    group_N <- function(input, k = 2) {
        N = length(input)
        m = N/k
        combos <- factorial(N) / (factorial(k)^m * factorial(m))

        result <- array(NA_integer_, dim = c(combos, m, k))
        combination = matrix(NA_integer_, nrow = k, ncol = m)

        z = 1
        group_f_start = 'group <- function(input, step){\n len <- length(input) \n combination[1,  step] <<- input[1] \n '
        i_s <- paste0('i', seq_len(k-1))

        group_f_fors = paste0(
            'for (', i_s, ' in ',
            c('2', if (length(i_s) != 1) {paste0('(', i_s[-length(i_s)], '+1)')}),
            ':(len-', rev(seq_len(k)[-k]) - 1, ')) { \n combination[', seq_len(k)[-1],
            ', step] <<- input[', i_s, '] \n',
            collapse = '\n '
        )

        group_f_inner = paste0(
            'if (step == m) { \n result[z, ,] <<- combination \n z <<- z+1 \n } else { \n rest <- setdiff(input, input[c(',
            paste0(i_s, collapse = ','),
            ', 1)]) \n group(rest, step +1) \n }'
        )

        eval(parse(text = paste0(
            group_f_start, group_f_fors, group_f_inner,
            paste0(rep('}', times = k), collapse = ' \n ')
        )))

        group(input, 1)
        return(result)
    }

    base_R <- function(v, numGroups) {
        v <- RcppAlgos:::GetV(v)
        n <- length(v)
        size_grp <- as.integer(round(n / numGroups))
        temp <- group_N(v, size_grp)
        mat <- t(apply(temp, 1, as.vector))
        colnames(mat) <- rep(paste0("Grp", 1:numGroups), each = size_grp)
        mat
    }

    set.seed(42)
    cmplx_v <- complex(real = runif(8), imaginary = runif(8))

    for (i in c(1, 2, 4)) {
        expect_equal(base_R(8, i), comboGroups(8, i))
        expect_equal(base_R(letters[1:8], i), comboGroups(letters[1:8], i))
        expect_equal(base_R(cmplx_v, i), comboGroups(cmplx_v, i))
    }

    for (i in c(1, 2, 5)) {
        expect_equal(base_R(10, i), comboGroups(10, i))
        expect_equal(base_R(letters[1:10], i), comboGroups(letters[1:10], i))
    }

    expect_equal(base_R(12, 4), comboGroups(12, 4))
    expect_equal(base_R(12, 3), comboGroups(12, 3))
    expect_equal(base_R(12, 4), comboGroups(12, 4))


    expect_equal(nrow(comboGroups(letters[1:12], 4, retType = "3Darray")),
                 comboGroupsCount(12, 4))
    expect_equal(nrow(comboGroups(factor(letters[1:12]), 3,
                                  retType = "3Darray")),
                 comboGroupsCount(12, 3))
    expect_equal(nrow(comboGroups(12, 1, retType = "3Darray")),
                 comboGroupsCount(12, 1))
    expect_equal(nrow(comboGroups(12, 12, retType = "3Darray")),
                 comboGroupsCount(12, 12))

    # comboGroupsCount(50, 5)
    # Big Integer ('bigz') :
    # [1] 402789797982510165934296910320
    expect_equal(dim(comboGroups(50, 5, retType = "3Darray",
                                 lower = "402789797982510165934296910301")),
                 c(20, 10, 5))
    expect_equal(comboGroups(10, 5, lower = 201, upper = 300),
                 comboGroups(10, 5)[201:300, ])

    expect_equal(comboGroups(12, 4),
                 comboGroupsSample(12, 4, sampleVec = 1:comboGroupsCount(12, 4)))

    expect_equal(comboGroups(12, 4),
                 comboGroups(12, 4, nThreads = 2))
    expect_equal(
        comboGroups(12, 4),
        comboGroupsSample(12, 4,
                          sampleVec = 1:comboGroupsCount(12, 4), nThreads = 2)
    )

    expect_equal(dim(comboGroups(200, 5, retType = "3Darray", upper = 100)),
                 c(100, 40, 5))

    expect_equal(comboGroups(30, 15, lower = 5e15, upper = 5e15 + 1e5),
                 comboGroupsSample(30, 15, nThreads = 2,
                                   sampleVec = (5e15):(5e15 + 1e5)))

    expect_equal(comboGroups(14, 7, nThreads = 2), comboGroups(14, 7))
    expect_equal(sum(comboGroups(c(T, F), 2, retType = "3Darray")), 1)

    expect_equal(rownames(comboGroupsSample(9, 3, retType = "3Darray",
                                            sampleVec = c(67, 15, 248),
                                            namedSample = TRUE)),
                 as.character(c(67, 15, 248)))

    expect_equal(dim(comboGroups(as.raw(1:4), 2, retType = "3Darray")),
                 dim(comboGroups(as.complex(c(1, -1, 1i, -1i)), 2,
                                 retType = "3Darray")))

    ## test class preservations
    expect_equal(class(comboGroups(as.raw(1:4), 2)[1,]), "raw")
    expect_equal(class(comboGroups(as.character(1:4), 2)[1,]), "character")
    expect_equal(class(comboGroups(1:4, 2)[1,]), "integer")
    expect_equal(class(comboGroups(1:4 + 0.1, 2)[1,]), "numeric")
    expect_equal(class(comboGroups(c(T, F), 2)[1, ]), "logical")
    expect_equal(class(comboGroups(as.complex(c(1, -1, 1i, -1i)), 2)[1,]), "complex")

    ## comboGroupsCount(9, 3)
    ## [1] 280
    expect_equal(comboGroups(as.raw(1:9), 3)[c(1, 100, 280), ],
                 comboGroupsSample(as.raw(1:9), 3, sampleVec = c(1, 100, 280)))
    expect_equal(comboGroups(LETTERS[1:9], 3)[c(1, 100, 280), ],
                 comboGroupsSample(LETTERS[1:9], 3, sampleVec = c(1, 100, 280)))
    cmp_v = c(1, -1, 1i, -1i, 2, -2, 2i, -2i, -3i)
    expect_equal(comboGroups(cmp_v, 3)[c(1, 100, 280), ],
                 comboGroupsSample(cmp_v, 3, sampleVec = c(1, 100, 280)))

    ## The logical case is bit strange... with only two values, we can only have
    ## result with numGroups = 1 or 2
    expect_equal(comboGroups(c(T, F), 1),
                 comboGroupsSample(c(T, F), 1, n = 1))
    expect_equal(comboGroups(c(T, F), 2),
                 comboGroupsSample(c(T, F), 2, n = 1))

    expect_equal(rownames(comboGroupsSample(30, 5, n = 2,
                                            seed = 1, namedSample = TRUE)),
                 c("7162662695786451", "3525427663529072"))

    expect_equal(
        comboGroups(1:4 + 0.1,       2, retType = "3Darray"),
        comboGroupsSample(1:4 + 0.1, 2, retType = "3Darray", sampleVec = 1:3)
    )
})

test_that("comboGroups produces correct results when groups vary", {

    retIdxMat <- function(x, g) {
        while (is.list(x[[1]])) {
            x <- unlist(x, recursive = FALSE)
        }

        res <- do.call(rbind, x)
        s <- cumsum(c(0, g[-length(g)])) + 1
        r <- rle(g)

        g_idx <- 1L

        for (k in r$lengths) {
            if (k == 1) {
                g_idx <- g_idx + 1L
                next
            }

            for (j in 1:(k - 1L)) {
                for (q in j:(k - 1L)) {
                    test <- res
                    c1 <- s[g_idx + (j - 1L)]
                    c2 <- s[g_idx + q]

                    idx    <- which(res[, c1] > res[, c2])
                    g_size <- g[g_idx]

                    rng1 <- c1:(c1 + g_size - 1L)
                    rng2 <- c2:(c2 + g_size - 1L)

                    for (i in idx) {
                        temp <- test[i, rng1]
                        test[i, rng1] <- test[i, rng2]
                        test[i, rng2] <- temp
                    }

                    b   <- apply(test, 1, paste, collapse = "_")
                    res <- res[!duplicated(b), ]
                }
            }

            g_idx <- g_idx + k
        }

        return(res)
    }

    create_txt_fun <- function(g) {
        first <- paste("combn(x,", g[1], ", function(a) {")

        lst <- lapply(seq_len(length(g) - 1), function(i) {
            v <- if (i == 1) letters[1] else paste(letters[1:i],
                                                   collapse = ", ")
            paste("combn(setdiff(x, c(", v, ")),", g[i + 1],
                  ", function(", letters[i + 1] , ") {")
        })

        eqn <- do.call(paste, append(lst, first, 0))
        eqn <- paste(eqn, "c(", paste(letters[seq_along(g)], collapse = ", "),
                     ")}, simplify = FALSE)")
        return(paste(eqn, paste(rep("}, simplify = FALSE)",
                                length(g) - 1), collapse = "")))
    }

    gen_base_R <- function(v, g) {
        isFactor <- is.factor(v)
        eqn <- create_txt_fun(g)
        v   <- RcppAlgos:::GetV(v)
        x   <- seq_along(v)
        lst <- eval(parse(text = eqn))
        res <- v[retIdxMat(lst, g)]
        dim(res) <- c(length(res) / length(x), length(x))
        colnames(res) <- paste0("Grp", rep(seq_along(g), times = g))
        return(res)
    }

    expect_equal(gen_base_R(10, 1:4), comboGroups(10, grpSizes = 1:4))
    expect_equal(gen_base_R(letters[1:10], c(2, 2, 3, 3)),
                 comboGroups(letters[1:10], grpSizes = c(2, 2, 3, 3)))
    expect_equal(gen_base_R(factor(letters[1:10]), c(1, 3, 3, 3)),
                 comboGroups(factor(letters[1:10]), grpSizes = c(1, 3, 3, 3)))

    set.seed(42)
    cmplx_v <- complex(real = runif(10), imaginary = runif(10))
    expect_equal(gen_base_R(cmplx_v, c(2, 2, 2, 4)),
                 comboGroups(cmplx_v, grpSizes = c(2, 2, 2, 4)))

    # comboGroupsCount(30, grpSizes = rep(1:4, 3))
    # Big Integer ('bigz') :
    # [1] 8567963547433290000000
    expect_equal(dim(comboGroups(30, grpSizes = rep(1:4, 3),
                                 lower = "8567963547433289999981")), c(20, 30))
    expect_equal(comboGroups(10, grpSizes = 1:4, lower = 201, upper = 300),
                 comboGroups(10, grpSizes = 1:4)[201:300, ])

    # comboGroupsCount(v = sum(2:8), grpSizes = 2:8)
    # Big Integer ('bigz') :
    #     [1] 2043503369875369321920000
    v <- seq_len(sum(2:8))
    r <- c()

    for (i in 8:2) {
        r <- c(v[1:i], r)
        v <- setdiff(v, v[1:i])
    }

    t <- comboGroups(v = sum(2:8), grpSizes = 2:8,
                     lower = "2043503369875369321920000")[1, ]

    expect_equal(unname(t), r)

    ## See issue 53
    ##  ComboGroups for various group size seems missing some combination
    sol <- matrix(c("a", "a", "a", "b", "b", "c",
                    "b", "c", "d", "c", "d", "d",
                    "c", "b", "b", "a", "a", "a",
                    "d", "d", "c", "d", "c", "b"), ncol = 4)
    expect_equal(unname(comboGroups(letters[1:4], grpSizes = c(1, 1, 2))),
                 sol)

    expect_equal(comboGroups(3, grpSizes = 1:2),
                 comboGroupsSample(
                     3, grpSizes = 1:2,
                     sampleVec = 1:comboGroupsCount(3, grpSizes = 1:2)
                 ))

    expect_equal(comboGroups(5, grpSizes = c(1, 2, 2)),
                 comboGroupsSample(
                     5, grpSizes = c(1, 2, 2),
                     sampleVec = 1:comboGroupsCount(5, grpSizes = c(1, 2, 2))
                 ))

    expect_equal(comboGroups(4, grpSizes = c(1, 1, 2)),
                 comboGroupsSample(
                     4, grpSizes = c(1, 1, 2),
                     sampleVec = 1:comboGroupsCount(4, grpSizes = c(1, 1, 2))
                 ))

    expect_equal(comboGroups(5, grpSizes = c(1, 2, 2)),
                 comboGroupsSample(
                     5, grpSizes = c(1, 2, 2),
                     sampleVec = 1:comboGroupsCount(5, grpSizes = c(1, 2, 2))
                 ))

    expect_equal(comboGroups(5, grpSizes = c(1, 1, 1, 2)),
                 comboGroupsSample(
                     5, grpSizes = c(1, 1, 1, 2),
                     sampleVec = 1:comboGroupsCount(5, grpSizes = c(1, 1, 1, 2))
                 ))

    expect_equal(comboGroups(7, grpSizes = c(1, 2, 2, 2)),
                 comboGroupsSample(
                     7, grpSizes = c(1, 2, 2, 2),
                     sampleVec = 1:comboGroupsCount(7, grpSizes = c(1, 2, 2, 2))
                 ))

    bigCount <- comboGroupsCount(12, grpSizes = c(1, 1, 1, 1, 2, 2, 2, 2))
    expect_equal(comboGroups(12, grpSizes = c(1, 1, 1, 1, 2, 2, 2, 2)),
                 comboGroupsSample(
                     12, grpSizes = c(1, 1, 1, 1, 2, 2, 2, 2),
                     sampleVec = 1:bigCount
                 ))

    expect_equal(
        comboGroups(12, grpSizes = c(1, 1, 1, 1, 2, 2, 2, 2)),
        comboGroupsSample(12, grpSizes = c(1, 1, 1, 1, 2, 2, 2, 2),
                          sampleVec = 1:bigCount, nThreads = 2)
    )

    expect_equal(dim(comboGroups(200, grpSizes = rep(5:8, times = 6:9),
                                 upper = 100)),
                 c(100, 200))

    expect_equal(comboGroups(30, grpSizes = 2 * (4:6),
                             lower = 5e11, upper = 5e11 + 1e5),
                 comboGroupsSample(30, grpSizes = 2 * (4:6), nThreads = 2,
                                   sampleVec = (5e11):(5e11 + 1e5)))

    expect_equal(comboGroups(30, grpSizes = rep(4:6, each = 2),
                             lower = 5e15, upper = 5e15 + 1e5),
                 comboGroupsSample(30, grpSizes = rep(4:6, each = 2),
                                   nThreads = 2,
                                   sampleVec = (5e15):(5e15 + 1e5)))

    expect_equal(comboGroups(14, grpSizes = c(4, 5, 5), nThreads = 2),
                 comboGroups(14, grpSizes = c(4, 5, 5)))

    expect_equal(comboGroups(14, grpSizes = c(3, 5, 6), nThreads = 2),
                 comboGroups(14, grpSizes = c(3, 5, 6)))

    expect_equal(rownames(comboGroupsSample(30, grpSizes = c(5, 10, 15),
                                            sampleVec = c(337806368805,
                                                          448334922606,
                                                          433232041762),
                                            namedSample = TRUE)),
                 as.character(c(337806368805,
                                448334922606,
                                433232041762)))

    expect_equal(rownames(comboGroupsSample(30, grpSizes = c(5, 10, 15), n = 2,
                                            seed = 1, namedSample = TRUE)),
                 c("866248189", "372205330702"))

    expect_equal(
        rownames(comboGroupsSample(
            80, grpSizes = c(5, 5, 5, 10, 10, 15, 15, 15), n = 2,
            seed = 1, namedSample = TRUE
        )),
        c("1424186085660521051631145763878676451483834974622463043287356",
          "17562722460072173386395565529258043289650490250407423764080023")
    )
})

