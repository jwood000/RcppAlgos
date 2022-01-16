# GetIndices <- function(obj) {
#
#     if (class(obj))
#
# }
#

comboRank <- function(..., v, repetition = FALSE, freqs = NULL) {
    ## This varies somewhat from expand.grid
    nargs <- length(args <- list(...))

    if (!nargs)
        return(integer(0))

    if (nargs == 1L) {
        a1 <- args[[1L]]

        if (is.list(a1)) {
            mat <- sapply(a1, identity)

            if (is.list(mat)) {
                stop("malformed arguments... all vectors should have the same length")
            }
        }

        if (is.matrix(a1)) {
            mat <- t(a1)
        } else if (is.vector(a1)) {
            mat <- matrix(a1, ncol = 1)
        } else {
            stop("input not supported")
        }
    }

    if (nargs == 0L)
        return(integer(0))

    iArgs <- seq_len(nargs)
    nmc <- paste0("Var", iArgs)
    nm <- names(args)

    if (is.null(nm))  {
        nm <- nmc
    } else if (any(ng0 <- nzchar(nm))) {
        nmc[ng0] <- nm[ng0]
    }


}

cleanVec <- function(v) {
    v <- sort(v)
    myRle <- rle(v)
    list(x = myRle$values, freqs = myRle$lengths)
}