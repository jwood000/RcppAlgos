compositionsRank <- function(..., v, repetition = FALSE,
                             freqs = NULL, target = NULL, weak = FALSE) {
    GetRankPart(..., v = v, repetition = repetition, freqs = freqs,
                target = target, IsComposition = TRUE, weak = weak)
}
