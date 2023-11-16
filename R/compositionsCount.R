compositionsCount <- function(v, m = NULL, ...) {
    stopifnot(is.numeric(v))
    UseMethod("compositionsCount")
}

compositionsCount.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL,
    target = NULL, weak = FALSE, ...
) {
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target),
                 v, m, repetition, freqs, "==", NULL, NULL, FALSE,
                 FALSE, TRUE, weak))
}

compositionsCount.table <- function(v, m = NULL, target = NULL,
                                    weak = FALSE, ...) {
    clean <- ResolveVFreqs(v)
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(clean$v, target),
                 clean$v, m, FALSE, clean$freqs, "==", NULL, NULL, FALSE,
                 FALSE, TRUE, weak))
}
