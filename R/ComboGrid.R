comboGrid <- function(..., repetition = TRUE) {
    ## This is from expand.grid
    nargs <- length(args <- list(...))

    if (any(sapply(args, is.null))) {
        return(expand.grid(args))
    }

    if (!nargs) {
        return(as.data.frame(list()))
    }

    if (nargs == 1L && is.list(a1 <- args[[1L]])) {
        nargs <- length(args <- a1)
    }

    if (nargs == 0L) {
        return(as.data.frame(list()))
    }

    iArgs <- seq_len(nargs)
    nmc   <- paste0("Var", iArgs)
    nm    <- names(args)

    if (is.null(nm)) {
        nm <- nmc
    } else if (any(ng0 <- nzchar(nm))) {
        nmc[ng0] <- nm[ng0]
    }

    idx_nas <- which(sapply(args, function(x) all(is.na(x))))

    pools <- args
    names(pools) <- nmc

    if (length(idx_nas)) {
        if (length(idx_nas) == nargs) {
            return(as.data.frame(pools))
        }

        pools <- pools[-idx_nas]
    }

    numChars <- sum(sapply(pools, class) == "character")
    convertCharToFac <- numChars < length(pools) && length(pools) > 1

    pools <- lapply(pools, function(x) {
        t <- sort(unique(x), na.last = FALSE)

        if (convertCharToFac && class(t) == "character") {
            return(factor(t, levels = t))
        } else {
            return(t)
        }
    })

    res <- .Call(Algos_ComboGridCpp, pools, repetition)

    if (length(idx_nas)) {
        res <- as.data.frame(res)

        for (idx in idx_nas) {
            res[, nmc[idx]] <- NA
        }

        res <- res[, nmc]
    }

    return(res)
}