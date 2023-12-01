partitionsCount <- function(v, m = NULL, ...) {
    stopifnot(is.numeric(v))
    UseMethod("partitionsCount")
}

partitionsCount.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL, ...
) {
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target),
                 v, m, repetition, freqs, "==", NULL, NULL, FALSE,
                 FALSE, FALSE, FALSE))
}

partitionsCount.table <- function(v, m = NULL, target = NULL, ...) {
    clean <- ResolveVFreqs(v)
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(clean$v, target),
                 clean$v, m, FALSE, clean$freqs, "==", NULL, NULL, FALSE,
                 FALSE, FALSE, FALSE))
}
