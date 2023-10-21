permuteCount <-  function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    UseMethod("permuteCount")
}

permuteCount.default <-  function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    ComboPermCount(v, m, repetition, freqs, FALSE)
}

permuteCount.table <-  function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    clean <- ResolveVFreqs(v, freqs)
    ComboPermCount(clean$v, m, repetition, clean$freqs, FALSE)
}

permuteCount.list <-  function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    ComboPermCount(length(v), m, repetition, freqs, FALSE)
}
