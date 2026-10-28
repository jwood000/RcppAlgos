test_that("Multiset ranking handles reservations and exhausted buckets", {
    expect_identical(
        partitionsRank(rbind(c(1L, 3L, 4L), c(2L, 3L, 3L)),
                       v = 1:4, freqs = c(1, 1, 2, 1), target = 8),
        1:2
    )
    expect_identical(
        partitionsRank(c(1, 2, 3, 3), v = 1:4,
                       freqs = c(1, 1, 2, 1), target = 9),
        1L
    )
    expect_identical(
        partitionsRank(rbind(c(1L, 3L), c(2L, 2L)), v = 1:3,
                       freqs = c(1, 2, 1), target = 4),
        1:2
    )
    expect_identical(
        partitionsRank(1:3, v = 5, freqs = 1:5, target = 6),
        1L
    )
})

test_that("Multiset ranking preserves mapped values and zero-valued parts", {
    freqs <- c(3, rep(2:1, 5))
    rows <- partitionsGeneral(0:10, 5, freqs = freqs, target = 20)
    expect_identical(nrow(rows), 65L)

    for (mapping in list(c(1, 0), c(1, 1), c(2, -7), c(2e9, 0))) {
        scale <- mapping[1L]
        shift <- mapping[2L]
        v <- (0:10) * scale + shift
        mapped <- rows * scale + shift
        target <- 20 * scale + 5 * shift
        for (threads in 1:2) {
            expect_identical(
                partitionsRank(mapped, v = v, freqs = freqs,
                               target = target, nThreads = threads),
                seq_len(nrow(rows))
            )
        }
    }

    # Row order, duplicates, and order within a partition must not affect rank.
    selected <- c(65L, 1L, 32L, 32L)
    expect_identical(
        partitionsRank(rows[selected, 5:1], v = 0:10, freqs = freqs,
                       target = 20, nThreads = 2),
        selected
    )
})

test_that("Multiset double ranks above the integer limit round-trip", {
    v <- 1:50
    freqs <- rep(2:4, length.out = length(v))
    total <- partitionsCount(v, 12, freqs = freqs, target = 240)
    expect_type(total, "double")
    expect_gt(total, .Machine$integer.max)
    expect_lt(total, 2^53)
    indices <- c(1, 2, 2^31 - 1, 2^31, 2^31 + 1,
                 floor(total / 2), total - 1, total)
    rows <- partitionsSample(v, 12, freqs = freqs, target = 240,
                             sampleVec = indices, nThreads = 2)
    for (threads in 1:2) {
        expect_identical(
            partitionsRank(rows, v = v, freqs = freqs,
                           target = 240, nThreads = threads),
            indices
        )
    }
})

test_that("Multiset ranking still validates multiplicities", {
    expect_error(
        partitionsRank(c(2, 2, 4), v = 1:4,
                       freqs = c(1, 1, 2, 1), target = 8),
        "Input frequencies do not match supplied freqs"
    )
})

test_that("GMP multiset ranks are exact at large ranks and block boundaries", {
    v <- 1:100
    freqs <- rep(2:4, length.out = length(v))
    total <- partitionsCount(v, 20, freqs = freqs, target = 500)
    expect_s3_class(total, "bigz")
    limit <- gmp::as.bigz("9007199254740992")

    # Count the entire first block independently: fix the first part at 1.
    remaining <- freqs
    remaining[1L] <- remaining[1L] - 1L
    boundary <- partitionsCount(v, 19, freqs = remaining, target = 499)
    expect_true(boundary > limit)
    expect_true(boundary < total)
    indices <- c(gmp::as.bigz(1:2), limit - 1, limit, limit + 1,
                 boundary - 1, boundary, boundary + 1, total - 1, total)
    rows <- partitionsSample(v, 20, freqs = freqs, target = 500,
                             sampleVec = indices, nThreads = 2)
    expect_identical(rows[7L, 1L], 1L)
    expect_gt(rows[8L, 1L], 1L)
    for (threads in 1:2) {
        ranks <- partitionsRank(rows, v = v, freqs = freqs,
                                target = 500, nThreads = threads)
        expect_s3_class(ranks, "bigz")
        expect_identical(as.character(ranks), as.character(indices))
    }

    expect_identical(
        as.character(partitionsRank(rows - 1L, v = v - 1L, freqs = freqs,
                                    target = 480, nThreads = 2)),
        as.character(indices)
    )
    expect_identical(
        as.character(partitionsRank(rows * 2e9, v = v * 2e9, freqs = freqs,
                                    target = 1e12, nThreads = 2)),
        as.character(indices)
    )

    # Sequential enumeration independently checks the first ranks.
    first <- partitionsGeneral(v, 20, freqs = freqs, target = 500, upper = 12)
    expect_identical(
        as.character(partitionsRank(first, v = v, freqs = freqs,
                                    target = 500, nThreads = 2)),
        as.character(1:12)
    )
})

test_that("GMP multiset ranking agrees with unrestricted repetition", {
    v <- 1:100
    width <- 20L
    freqs <- rep(width, length(v))
    total <- partitionsCount(v, width, repetition = TRUE, target = 500)
    expect_s3_class(total, "bigz")
    indices <- c(gmp::as.bigz(1), gmp::as.bigz("9007199254740993"),
                 total %/% 2, total - 1, total)
    rows <- partitionsSample(v, width, repetition = TRUE, target = 500,
                             sampleVec = indices, nThreads = 2)
    expect_identical(
        as.character(partitionsRank(rows, v = v, freqs = freqs,
                                    target = 500, nThreads = 2)),
        as.character(partitionsRank(rows, v = v, repetition = TRUE,
                                    target = 500, nThreads = 2))
    )
})
