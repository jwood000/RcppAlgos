partitionsCount <- function(v, m = NULL, repetition = FALSE,
                            freqs = NULL, target = NULL) {

    stopifnot(is.numeric(v))
    UseMethod("partitionsCount")
}

partitionsCount.default <- function(v, m = NULL, repetition = FALSE,
                                    freqs = NULL, target = NULL) {
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target),
                 v, m, repetition, freqs, "==", NULL, NULL, FALSE,
                 FALSE, FALSE, FALSE))
}

partitionsCount.table <- function(v, m = NULL, repetition = FALSE,
                                  freqs = NULL, target = NULL) {
    clean <- ResolveVFreqs(v, freqs)
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(clean$v, target),
                 v, m, repetition, clean$freqs, "==", NULL, NULL, FALSE,
                 FALSE, FALSE, FALSE))
}
