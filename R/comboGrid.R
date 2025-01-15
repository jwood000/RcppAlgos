comboGrid <- function(..., repetition = TRUE, return_df = FALSE) {

    lst <- GridInputs(...)

    if (lst$early_return) {
        return(lst$res)
    }

    pools <- lst$p
    nmc   <- lst$nmc
    iArgs <- lst$iArgs

    if (length(lst$i)) {
        if (length(lst$i) == lst$n) {
            return(as.data.frame(pools))
        }

        pools <- pools[-lst$i]
    }

    pools <- lapply(pools, function(x) {
        sort(unique(x), na.last = FALSE)
    })

    res <- .Call(`_RcppAlgos_ComboGridCpp`, pools, repetition, return_df)

    if (length(lst$i)) {
        res <- as.data.frame(res)
        names(res) <- nmc[setdiff(iArgs, lst$i)]

        for (idx in lst$i) {
            res[nmc[idx]] <- NA
        }

        res <- res[, nmc]
    } else if (is.matrix(res)) {
        colnames(res) <- nmc
    }

    return(res)
}
