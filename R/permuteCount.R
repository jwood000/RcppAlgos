permuteCount <-  function(v, m = NULL, ...) {
    UseMethod("permuteCount")
}

permuteCount.default <-  function(v, m = NULL, repetition = FALSE,
                                  freqs = NULL, ...) {
    ComboPermuteCount(v, m, repetition, freqs, FALSE)
}

permuteCount.table <-  function(v, m = NULL, ...) {
    clean <- ResolveVFreqs(v)
    ComboPermuteCount(clean$v, m, FALSE, clean$freqs, FALSE)
}

permuteCount.list <-  function(v, m = NULL, repetition = FALSE,
                               freqs = NULL, ...) {
    ComboPermuteCount(length(v), m, repetition, freqs, FALSE)
}
