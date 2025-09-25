compositionsRank <- function(..., v, repetition = FALSE, freqs = NULL,
                             target = NULL, weak = FALSE, nThreads = NULL) {
    GetRankPart(..., v = v, repetition = repetition, freqs = freqs,
                target = target, IsComposition = TRUE, weak = weak,
                nThreads = nThreads)
}
