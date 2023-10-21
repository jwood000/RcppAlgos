comboSample <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
    sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
    nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL
) {
    UseMethod("comboSample")
}

comboSample.default <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
    sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
    nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL
) {
    ComboPermSample(
        v, m, repetition, freqs, n, sampleVec, seed, FUN,
        Parallel, nThreads, namedSample, FUN.VALUE, TRUE
    )
}

comboSample.table <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
    sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
    nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL
) {
    clean <- ResolveVFreqs(v, freqs)
    ComboPermSample(
        clean$v, m, repetition, clean$freqs, n, sampleVec, seed, FUN,
        Parallel, nThreads, namedSample, FUN.VALUE, TRUE
    )
}

comboSample.list <- function(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
    sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
    nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL
) {
    ComboPermSample(
        seq_along(v), m, repetition, freqs, n, sampleVec, seed,
        FUN = function(x) v[x], Parallel, nThreads, namedSample,
        FUN.VALUE = NULL, TRUE
    )
}
