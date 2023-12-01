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

## Currently only used when an object of class table is passed
ResolveVFreqs <- function(v) {

    nms <- names(v)

    suppressWarnings({
        ## round-trip conversion methodology.
        ##
        ## N.B. We are unable to get back factors.
        ##
        ## N.B. We don't test for raws as table doesn't support raws
        ## Error in order(y) : unimplemented type 'raw' in 'orderVector1'
        conv_v <- if (identical(nms, as.character(as.integer(nms)))) {
            as.integer(nms)
        } else if (identical(nms, as.character(as.logical(nms)))) {
            as.logical(nms)
        } else if (isTRUE(all.equal(nms, as.character(as.numeric(nms))))) {
            as.numeric(nms)
        } else if (isTRUE(all.equal(nms, as.character(as.complex(nms))))) {
            as.complex(nms)
        } else {
            nms
        }
    })

    ## Get the tabulated frequencies
    attributes(v) <- NULL
    return(list(v = conv_v, freqs = v))
}

GetTarget <- function(v, target) {
    if (is.null(target)) {target <- max(v, na.rm = TRUE)}
    return(target)
}

GetRank <- function(..., v, repetition = FALSE,
                    freqs = NULL, IsComb = TRUE) {

    n_args <- length(arg_s <- list(...))

    if (!n_args) {
        return(integer(0))
    } else if (n_args > 1L) {
        arg_s <- list(arg_s)
    } else if (!is.atomic(arg_s[[1]])) {
        stop("Input must be atomic")
    }

    input <- arg_s[[1L]]
    msg   <- "Inputs must be a subset of v"
    v     <- GetV(v)

    if (is.list(input)) {
        return(
            Map(function(obj) {
                if (!is.atomic(obj)) stop("Inputs must be atomic")
                idx <- match(if (is.matrix(obj)) t(obj) else obj, v)
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
    } else {
        idx <- match(input, v)
        if (any(is.na(idx))) stop(msg)
        return(.Call(`_RcppAlgos_RankCombPerm`, idx, v,
                     repetition, freqs, length(input), IsComb));
    }
}

GetRankPart <- function(..., v, repetition = FALSE, freqs = NULL,
                        target = NULL, IsComposition = FALSE, weak = FALSE) {

    n_args   <- length(arg_s <- list(...))
    target   <- GetTarget(v, target)
    msg_sub  <- "Inputs must be a subset of v"
    msg_part <- paste("Inputs must be a partition of", target)
    msg_cls  <- "Inputs must be of class numeric or integer"

    if (!n_args) {
        return(integer(0))
    } else if (n_args > 1L) {
        arg_s <- list(arg_s)
    } else if (!is.numeric(arg_s[[1]])) {
        stop(msg_cls)
    }

    input <- arg_s[[1L]]
    v     <- GetV(v)

    if (!is.numeric(v)) {
        stop("v must be of class numeric or integer")
    }

    if (is.list(input)) {
        return(
            Map(function(obj) {
                if (!is.numeric(obj)) stop(msg_cls)
                if ((is.matrix(obj) && any(rowSums(obj) != target)) ||
                    sum(obj) != target) stop(msg_part)
                idx <- match(if (is.matrix(obj)) t(obj) else obj, v)
                if (any(is.na(idx))) stop(msg_sub)
                .Call(`_RcppAlgos_RankPartitionMain`, idx, v, repetition, freqs,
                      if (is.matrix(obj)) ncol(obj) else length(obj),
                      "==", target, NULL, IsComposition, weak)
            }, input)
        )
    } else if (!is.numeric(input)) {
        stop(msg_cls)
    } else if (is.matrix(input)) {
        if (any(rowSums(input) != target)) stop(msg_part)
        idx <- match(t(input), v)
        if (any(is.na(idx))) stop(msg_sub)
        return(.Call(`_RcppAlgos_RankPartitionMain`, idx, v,
                     repetition, freqs, ncol(input), "==", target,
                     NULL, IsComposition, weak));
    } else {
        if (sum(input) != target) stop(msg_part)
        idx <- match(input, v)
        if (any(is.na(idx))) stop(msg_sub)
        return(.Call(`_RcppAlgos_RankPartitionMain`, idx, v,
                     repetition, freqs, length(input), "==", target,
                     NULL, IsComposition, weak));
    }
}
