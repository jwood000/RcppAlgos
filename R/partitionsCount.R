partitionsCount <- function(v, m = NULL, ...) {
    stopifnot(is.numeric(v))
    UseMethod("partitionsCount")
}

partitionsCount.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL, ...
) {
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target),
                 v, m, repetition, freqs, TRUE, "==", NULL, NULL,
                 FALSE, FALSE, FALSE, FALSE))
}

partitionsCount.table <- function(v, m = NULL, target = NULL, ...) {
    clean <- ResolveVFreqs(v)
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(clean$v, target),
                 clean$v, m, FALSE, clean$freqs, TRUE, "==", NULL, NULL,
                 FALSE, FALSE, FALSE, FALSE))
}

partitionsMultisetCount <- function(v, m = NULL, ...) {
    stopifnot(is.numeric(v))
    UseMethod("partitionsMultisetCount")
}

partitionsMultisetCount.default <- function(
    v, m = NULL, freqs = NULL, target = NULL, checkGmp = FALSE, ...
) {
    return(.Call(
        `_RcppAlgos_PartitionsMultisetCount`, GetTarget(v, target),
        v, m, freqs, TRUE, "==", NULL, NULL, FALSE, FALSE, checkGmp
    ))
}

partitionsMultisetCount.table <- function(
    v, m = NULL, target = NULL, checkGmp = FALSE, ...
) {
    clean <- ResolveVFreqs(v)
    return(.Call(
        `_RcppAlgos_PartitionsMultisetCount`, GetTarget(clean$v, target),
        clean$v, m, clean$freqs, TRUE, "==", NULL, NULL, FALSE, FALSE, checkGmp
    ))
}
