compositionsCount <- function(v, m = NULL, ...) {
    stopifnot(is.numeric(v))
    UseMethod("compositionsCount")
}

compositionsCount.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL,
    target = NULL, weak = FALSE, ...
) {
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target),
                 v, m, repetition, freqs, FALSE, "==", NULL, NULL,
                 FALSE, FALSE, TRUE, weak))
}

compositionsCount.table <- function(v, m = NULL, target = NULL,
                                    weak = FALSE, ...) {
    clean <- ResolveVFreqs(v)
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(clean$v, target),
                 clean$v, m, FALSE, clean$freqs, FALSE, "==", NULL, NULL,
                 FALSE, FALSE, TRUE, weak))
}

compositionsMultisetCount <- function(v, m = NULL, ...) {
    stopifnot(is.numeric(v))
    UseMethod("compositionsMultisetCount")
}

compositionsMultisetCount.default <- function(
    v, m = NULL, freqs = NULL, target = NULL,
    weak = FALSE, checkGmp = FALSE, ...
) {
    return(.Call(
        `_RcppAlgos_PartitionsMultisetCount`, GetTarget(v, target),
        v, m, freqs, FALSE, "==", NULL, NULL, TRUE, weak, checkGmp
    ))
}

compositionsMultisetCount.table <- function(
    v, m = NULL, target = NULL, weak = FALSE, checkGmp = FALSE, ...
) {
    clean <- ResolveVFreqs(v)
    return(.Call(
        `_RcppAlgos_PartitionsMultisetCount`, GetTarget(clean$v, target),
        clean$v, m, clean$freqs, FALSE, "==", NULL, NULL, TRUE, weak, checkGmp
    ))
}
