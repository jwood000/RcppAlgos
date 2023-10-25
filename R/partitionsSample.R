partitionsSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                             target = NULL, n = NULL, sampleVec = NULL,
                             seed = NULL, nThreads = NULL,
                             namedSample = FALSE) {

    stopifnot(is.numeric(v))

    if (!is.null(seed)) {
        set.seed(seed)
    }

    UseMethod("partitionsSample")
}

partitionsSample.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL, n = NULL,
    sampleVec = NULL, seed = NULL, nThreads = NULL, namedSample = FALSE
) {
    return(.Call(`_RcppAlgos_SamplePartitions`, v, m, repetition, freqs,
                 sampleVec, seed, n, sample, FALSE, nThreads,
                 pkgEnv$nThreads, namedSample, "==",
                 GetTarget(v, target), NULL, new.env(), FALSE, FALSE))
}

partitionsSample.table <- function(
        v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL, n = NULL,
        sampleVec = NULL, seed = NULL, nThreads = NULL, namedSample = FALSE
) {
    clean <- ResolveVFreqs(v, freqs)
    return(.Call(`_RcppAlgos_SamplePartitions`, clean$v, m, repetition,
                 clean$freqs, sampleVec, seed, n, sample, FALSE, nThreads,
                 pkgEnv$nThreads, namedSample, "==",
                 GetTarget(clean$v, target), NULL, new.env(), FALSE, FALSE))
}
