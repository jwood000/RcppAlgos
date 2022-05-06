partitionsGeneral <- function(v, m = NULL, repetition = FALSE,
                              freqs = NULL, target = NULL, lower = NULL,
                              upper = NULL, nThreads = NULL,
                              tolerance = NULL) {

    stopifnot(is.numeric(v))
    return(.Call(`_RcppAlgos_CombinatoricsCnstrt`, v, m, repetition, freqs,
                 lower, upper, "sum", "==", GetTarget(v, target), TRUE,
                 FALSE, FALSE, nThreads, pkgEnv$nThreads, tolerance))
}

partitionsCount <- function(v, m = NULL, repetition = FALSE,
                            freqs = NULL, target = NULL) {

    stopifnot(is.numeric(v))
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target), v, m,
                 repetition, freqs, "==", NULL, NULL, FALSE, FALSE))
}

partitionsDesign <- function(v, m = NULL, repetition = FALSE,
                             freqs = NULL, target = NULL,
                             showDesign = FALSE) {

    stopifnot(is.numeric(v))
    return(.Call(`_RcppAlgos_PartitionsCount`, GetTarget(v, target), v, m,
                 repetition, freqs, "==", NULL, NULL, TRUE, showDesign))
}

partitionsRank <- function(..., v, repetition = FALSE,
                           freqs = NULL, target = NULL) {

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
                      "==", target, NULL)
            }, input)
        )
    } else if (!is.numeric(input)) {
        stop(msg_cls)
    } else if (is.matrix(input)) {
        if (any(rowSums(input) != target)) stop(msg_part)
        idx <- match(t(input), v)
        if (any(is.na(idx))) stop(msg_sub)
        return(.Call(`_RcppAlgos_RankPartitionMain`, idx, v,
                     repetition, freqs, ncol(input), "==", target, NULL));
    } else {
        if (sum(input) != target) stop(msg_part)
        idx <- match(input, v)
        if (any(is.na(idx))) stop(msg_sub)
        return(.Call(`_RcppAlgos_RankPartitionMain`, idx, v,
                     repetition, freqs, length(input), "==", target, NULL));
    }
}

partitionsSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL,
                             target = NULL, n = NULL, sampleVec = NULL,
                             seed = NULL, nThreads = NULL,
                             namedSample = FALSE) {

    stopifnot(is.numeric(v))

    if (!is.null(seed)) {
        set.seed(seed)
    }

    return(.Call(`_RcppAlgos_SamplePartitions`, v, m, repetition, freqs,
                 sampleVec, seed, n, sample, FALSE, nThreads,
                 pkgEnv$nThreads, namedSample, "==",
                 GetTarget(v, target), NULL, new.env()))
}

partitionsIter <- function(v, m = NULL, repetition = FALSE,
                           freqs = NULL, target = NULL,
                           nThreads = NULL, tolerance = NULL) {

    stopifnot(is.numeric(v))
    InitVals <- .Call(`_RcppAlgos_GetClassVals`, v, m, repetition, freqs,
                      TRUE, NULL, nThreads, pkgEnv$nThreads, TRUE)

    new("Partitions", InitVals, FALSE, "sum", "==",
        GetTarget(v, target), FALSE, tolerance, is.null(m))
}
