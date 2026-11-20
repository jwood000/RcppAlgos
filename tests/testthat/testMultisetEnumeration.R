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

test_that("multiset counts, samples and ranks match exhaustive enumeration", {
    for (cap in 1:5) {
        frequencies <- expand.grid(rep(list(1:3), cap))
        for (case in seq_len(nrow(frequencies))) {
            freqs <- as.integer(frequencies[case, ])
            # The public API caps m at the total available multiplicity.
            widths <- seq_len(min(6L, sum(freqs)))
            expected_counts <- actual_counts <-
                expected_samples <- actual_samples <-
                expected_ranks <- actual_ranks <-
                setNames(vector("list", length(widths)), paste0("width=", widths))

            for (width in widths) {
                rows <- enumerate_multiset_parts(freqs, width)
                sums <- rowSums(rows)
                targets <- seq.int(width, cap * width)
                counts <- vapply(targets, function(target) {
                    sum(sums == target)
                }, integer(1))
                expected_counts[[width]] <- counts
                actual_counts[[width]] <- vapply(targets, function(target) {
                    as.numeric(partitionsCount(
                        seq_len(cap), width, freqs = freqs, target = target
                    ))
                }, numeric(1))

                nonempty <- targets[counts > 0L]
                names(nonempty) <- paste0("target=", nonempty)
                expected_samples[[width]] <- lapply(nonempty, function(target) {
                    rows[sums == target, , drop = FALSE]
                })
                actual_samples[[width]] <- lapply(nonempty, function(target) {
                    partitionsSample(
                        seq_len(cap), width, freqs = freqs, target = target,
                        sampleVec = seq_len(counts[target - width + 1L]),
                        nThreads = 1L
                    )
                })
                expected_ranks[[width]] <- lapply(expected_samples[[width]], function(rows) {
                    seq_len(nrow(rows))
                })
                actual_ranks[[width]] <- lapply(nonempty, function(target) {
                    partitionsRank(
                        rows[sums == target, , drop = FALSE], v = seq_len(cap),
                        freqs = freqs, target = target, nThreads = 1L
                    )
                })
            }

            # Check every original width, target and rank, while avoiding tens
            # of thousands of individual testthat expectations. Named list
            # entries retain the width/target context if a comparison fails.
            inputs <- paste0("freqs=", paste(freqs, collapse = ","))
            expect_equal(actual_counts, expected_counts, info = inputs)
            expect_identical(actual_samples, expected_samples, info = inputs)
            expect_identical(actual_ranks, expected_ranks, info = inputs)
        }
    }
})
