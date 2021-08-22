comboGrid <- function(..., repetition = TRUE) {
    ## This is from expand.grid
    nargs <- length(args <- list(...))
    if (!nargs)
        return(as.data.frame(list()))
    if (nargs == 1L && is.list(a1 <- args[[1L]]))
        nargs <- length(args <- a1)
    if (nargs == 0L)
        return(as.data.frame(list()))

    iArgs <- seq_len(nargs)
    nmc <- paste0("Var", iArgs)
    nm <- names(args)
    if (is.null(nm))
        nm <- nmc
    else if (any(ng0 <- nzchar(nm)))
        nmc[ng0] <- nm[ng0]

    pools <- args
    names(pools) <- nmc

    numChars <- sum(sapply(pools, class) == "character")
    convertCharToFac <- numChars < length(pools) && length(pools) > 1

    pools <- lapply(pools, function(x) {
        t <- sort(unique(x))

        if (convertCharToFac && class(t) == "character") {
            return(factor(t, levels = t))
        } else {
            return(t)
        }
    })

    .Call(ComboGridCpp, pools, repetition, PACKAGE = "RcppAlgos")
}