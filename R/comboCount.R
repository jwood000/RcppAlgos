comboCount <- function(v, m = NULL, ...) {
    UseMethod("comboCount")
}

comboCount.default <- function(v, m = NULL, repetition = FALSE,
                               freqs = NULL, ...) {
    ComboPermuteCount(v, m, repetition, freqs, TRUE)
}

comboCount.table <- function(v, m = NULL, ...) {
    clean <- ResolveVFreqs(v)
    ComboPermuteCount(clean$v, m, FALSE, clean$freqs, TRUE)
}

comboCount.list <- function(v, m = NULL, repetition = FALSE,
                            freqs = NULL, ...) {
    ComboPermuteCount(length(v), m, repetition, freqs, TRUE)
}
