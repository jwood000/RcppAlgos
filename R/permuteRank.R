permuteRank <- function(..., v, repetition = FALSE,
                        freqs = NULL, nThreads = NULL) {
    GetRank(..., v = v, repetition = repetition,
            freqs = freqs, IsComb = FALSE, nThreads = nThreads)
}
