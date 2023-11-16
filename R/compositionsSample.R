compositionsSample <- function(v, m = NULL, ...) {
    stopifnot(is.numeric(v))
    UseMethod("compositionsSample")
}

compositionsSample.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL,
    weak = FALSE, n = NULL, sampleVec = NULL, seed = NULL,
    nThreads = NULL, namedSample = FALSE, ...
) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    return(.Call(
        `_RcppAlgos_SamplePartitions`, v, m, repetition, freqs, sampleVec,
        seed, n, sample, nThreads, pkgEnv$nThreads, namedSample, "==",
        GetTarget(v, target), new.env(), TRUE, weak
    ))
}

compositionsSample.table <- function(
    v, m = NULL, target = NULL, weak = FALSE, n = NULL, sampleVec = NULL,
    seed = NULL, nThreads = NULL, namedSample = FALSE, ...
) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    clean <- ResolveVFreqs(v)
    return(.Call(
        `_RcppAlgos_SamplePartitions`, clean$v, m, FALSE, clean$freqs,
        sampleVec, seed, n, sample, nThreads, pkgEnv$nThreads, namedSample,
        "==", GetTarget(clean$v, target), new.env(), TRUE, weak
    ))
}
