comboSample <- function(v, m = NULL, ...) {
    UseMethod("comboSample")
}

comboSample.integer <-
comboSample.numeric <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
    sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
    nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL, ...
) {
    ComboPermuteSample(
        v, m, repetition, freqs, n, sampleVec, seed, FUN,
        Parallel, nThreads, namedSample, FUN.VALUE, TRUE
    )
}

comboSample.factor <-
comboSample.logical <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
    sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
    nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL, ...
) {
    ComboPermuteSample(
        v, m, repetition, freqs, n, sampleVec, seed, FUN,
        Parallel, nThreads, namedSample, FUN.VALUE, TRUE
    )
}

comboSample.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL, sampleVec = NULL,
    seed = NULL, FUN = NULL, namedSample = FALSE, FUN.VALUE = NULL, ...
) {
    ComboPermuteSample(
        v, m, repetition, freqs, n, sampleVec, seed, FUN,
        FALSE, NULL, namedSample, FUN.VALUE, TRUE
    )
}

comboSample.table <- function(
    v, m = NULL, n = NULL, sampleVec = NULL, seed = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL, ...
) {
    clean <- ResolveVFreqs(v)
    ComboPermuteSample(
        clean$v, m, FALSE, clean$freqs, n, sampleVec, seed, FUN,
        Parallel, nThreads, namedSample, FUN.VALUE, TRUE
    )
}

comboSample.list <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
    sampleVec = NULL, seed = NULL, namedSample = FALSE, ...
) {
    ComboPermuteSample(
        seq_along(v), m, repetition, freqs, n, sampleVec, seed,
        FUN = function(x) v[x], FALSE, NULL, namedSample, NULL, TRUE
    )
}
