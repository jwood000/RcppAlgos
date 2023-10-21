comboCount <-  function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    UseMethod("comboCount")
}

comboCount.default <-  function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    ComboPermCount(v, m, repetition, freqs, TRUE)
}

comboCount.table <-  function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    clean <- ResolveVFreqs(v, freqs)
    ComboPermCount(clean$v, m, repetition, clean$freqs, TRUE)
}

comboCount.list <-  function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    ComboPermCount(length(v), m, repetition, freqs, TRUE)
}
