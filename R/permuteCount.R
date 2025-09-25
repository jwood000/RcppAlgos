permuteCount <- function(v, m = NULL, ...) {
    UseMethod("permuteCount")
}

permuteCount.default <- function(v, m = NULL, repetition = FALSE,
                                 freqs = NULL, ...) {

    lst <- PermuteArgs(...)

    if (lst$is_part) {
        return(.Call(`_RcppAlgos_PartitionsCount`, lst$target,
                     v, m, repetition, freqs, FALSE, "==", NULL,
                     NULL, FALSE, FALSE, FALSE, FALSE))
    } else {
        ComboPermuteCount(v, m, repetition, freqs, FALSE)
    }
}

permuteCount.table <- function(v, m = NULL, ...) {

    clean <- ResolveVFreqs(v)
    lst <- PermuteArgs(...)

    if (lst$is_part) {
        return(.Call(`_RcppAlgos_PartitionsCount`, lst$target,
                     clean$v, m, FALSE, clean$freqs, FALSE,
                     "==", NULL, NULL, FALSE, FALSE, FALSE, FALSE))
    } else {
        ComboPermuteCount(clean$v, m, FALSE, clean$freqs, FALSE)
    }
}

permuteCount.list <- function(v, m = NULL, repetition = FALSE,
                              freqs = NULL, ...) {
    ComboPermuteCount(length(v), m, repetition, freqs, FALSE)
}
