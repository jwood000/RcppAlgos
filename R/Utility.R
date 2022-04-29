stdThreadMax <- function() {
    nThreads <- .Call(`_RcppAlgos_cpp11GetNumThreads`)
    if (nThreads < 1L || is.na(nThreads) || is.null(nThreads)) nThreads = 1L
    return(nThreads)
}

GetV <- function(v) {
    if (!is.atomic(v) || (is.raw(v) && !is.null(attributes(v)))) {
        stop("Only atomic types are supported for v")
    }

    if (is.numeric(v) && length(v) == 1L) {
        return(if (v < 0) v:-1L else if (v == 0) 0 else 1:v)
    }

    return(v)
}

GetTarget <- function(v, target) {
    if (is.null(target)) {target = max(v, na.rm = TRUE)}
    return(target)
}

GetRank <- function(..., v, repetition = FALSE,
                    freqs = NULL, IsComb = TRUE) {

    n_args <- length(arg_s <- list(...))

    if (!n_args) {
        return(integer(0))
    } else if (n_args > 1L) {
        arg_s <- list(arg_s)
    }

    input <- arg_s[[1L]]
    msg   <- "Inputs must be a subset of v"
    v     <- GetV(v)

    if (is.list(input)) {
        return(
            Map(function(obj) {
                if (!is.atomic(obj)) stop("Inputs must be atomic")
                idx <- match(if (is.matrix(obj)) t(input) else input, v)
                if (any(is.na(idx))) stop(msg)
                .Call(`_RcppAlgos_RankCombPerm`, idx, v, repetition, freqs,
                      if (is.matrix(obj)) ncol(obj) else length(obj), IsComb)
            }, input)
        )
    } else if (is.matrix(input)) {
        idx <- match(t(input), v)
        if (any(is.na(idx))) stop(msg)
        return(.Call(`_RcppAlgos_RankCombPerm`, idx, v,
                     repetition, freqs, ncol(input), IsComb));
    } else if (is.atomic(input)) {
        idx <- match(input, v)
        if (any(is.na(idx))) stop(msg)
        return(.Call(`_RcppAlgos_RankCombPerm`, idx, v,
                     repetition, freqs, length(input), IsComb));
    } else {
        stop("Input not supported")
    }
}
