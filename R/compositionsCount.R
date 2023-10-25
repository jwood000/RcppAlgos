compositionsCount <- function(v, m = NULL, repetition = FALSE,
                              freqs = NULL, target = NULL, weak = FALSE) {
    stopifnot(is.numeric(v))
    UseMethod("compositionsCount")
}

compositionsCount.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL, weak = FALSE
) {
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target),
                 v, m, repetition, freqs, "==", NULL, NULL, FALSE,
                 FALSE, TRUE, weak))
}

compositionsCount.table <- function(
        v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL, weak = FALSE
) {
    clean <- ResolveVFreqs(v, freqs)
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(clean$v, target),
                 clean$v, m, repetition, clean$freqs, "==", NULL, NULL,
                 FALSE, FALSE, TRUE, weak))
}
