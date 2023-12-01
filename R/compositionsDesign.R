compositionsDesign <- function(v, m = NULL, repetition = FALSE,
                               freqs = NULL, target = NULL,
                               weak = FALSE, showDesign = FALSE) {

    stopifnot(is.numeric(v))
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target),
                 v, m, repetition, freqs, "==", NULL, NULL, TRUE,
                 showDesign, TRUE, weak))
}
