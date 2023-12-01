comboGrid <- function(..., repetition = TRUE) {
    ## This is from expand.grid
    n_args <- length(arg_s <- list(...))

    if (any(sapply(arg_s, is.null))) {
        return(expand.grid(arg_s))
    }

    if (!n_args) {
        return(as.data.frame(list()))
    }

    if (n_args == 1L && is.list(a1 <- arg_s[[1L]])) {
        n_args <- length(arg_s <- a1)
    }

    if (n_args == 0L) {
        return(as.data.frame(list()))
    }

    iArgs <- seq_len(n_args)
    nmc   <- paste0("Var", iArgs)
    nm    <- names(arg_s)

    if (is.null(nm)) {
        nm <- nmc
    } else if (any(ng0 <- nzchar(nm))) {
        nmc[ng0] <- nm[ng0]
    }

    idx_nas <- which(sapply(arg_s, function(x) all(is.na(x))))

    pools <- arg_s
    names(pools) <- nmc

    if (length(idx_nas)) {
        if (length(idx_nas) == n_args) {
            return(as.data.frame(pools))
        }

        pools <- pools[-idx_nas]
    }

    numChars <- sum(sapply(pools, is.character))
    convertCharToFac <- numChars < length(pools) && length(pools) > 1

    pools <- lapply(pools, function(x) {
        t <- sort(unique(x), na.last = FALSE)

        if (convertCharToFac && is.character(t)) {
            return(factor(t, levels = t))
        } else {
            return(t)
        }
    })

    res <- .Call(`_RcppAlgos_ComboGridCpp`, pools, repetition)

    if (length(idx_nas)) {
        res <- as.data.frame(res)
        names(res) <- nmc[setdiff(iArgs, idx_nas)]

        for (idx in idx_nas) {
            res[nmc[idx]] <- NA
        }

        res <- res[, nmc]
    } else if (is.matrix(res)) {
        colnames(res) <- nmc
    }

    return(res)
}
