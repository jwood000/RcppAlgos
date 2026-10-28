# A small independent oracle: choose values in nondecreasing order and
# consume their multiplicities. No package counting/unranking code is used.
enumerate_multiset_parts <- function(freqs, width) {
    rows <- list()
    visit <- function(row, first, remaining) {
        if (length(row) == width) {
            rows[[length(rows) + 1L]] <<- row
            return(invisible(NULL))
        }
        for (value in seq.int(first, length(freqs))) {
            if (remaining[value] == 0L) next
            next_reps <- remaining
            next_reps[value] <- next_reps[value] - 1L
            visit(c(row, value), value, next_reps)
        }
    }
    visit(integer(), 1L, freqs)
    if (!length(rows)) return(matrix(integer(), 0L, width))
    do.call(rbind, rows)
}

test_that("multiset public counts and samples match exhaustive enumeration", {
    for (cap in 1:5) {
        frequencies <- expand.grid(rep(list(1:3), cap))
        for (case in seq_len(nrow(frequencies))) {
            freqs <- as.integer(frequencies[case, ])
            # The public API caps m at the total available multiplicity.
            for (width in seq_len(min(6L, sum(freqs)))) {
                rows <- enumerate_multiset_parts(freqs, width)
                sums <- rowSums(rows)
                targets <- seq.int(width, cap * width)
                expected_counts <- vapply(targets, function(target) {
                    sum(sums == target)
                }, integer(1))
                actual_counts <- vapply(targets, function(target) {
                    as.numeric(partitionsCount(
                        seq_len(cap), width, freqs = freqs, target = target
                    ))
                }, numeric(1))
                inputs <- paste0("freqs=", paste(freqs, collapse = ","),
                                 "; width=", width)
                expect_equal(actual_counts, expected_counts, info = inputs)

                for (target in targets[expected_counts > 0L]) {
                    expected <- rows[sums == target, , drop = FALSE]
                    actual <- partitionsSample(
                        seq_len(cap), width, freqs = freqs, target = target,
                        sampleVec = seq_len(nrow(expected)), nThreads = 1L
                    )
                    expect_identical(
                        actual, expected,
                        info = paste0(inputs, "; target=", target)
                    )
                }
            }
        }
    }
})
