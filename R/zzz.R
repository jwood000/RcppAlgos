# nocov start

pkgEnv <- new.env(parent = emptyenv())
pkgEnv$nCores <- NULL
pkgEnv$maxThreads <- NULL

physicalCoreCount <- function() {
    ## Estimate physical cores using R's default detection path. This value is
    ## used for cache-size heuristics, not as the authoritative thread limit.
    ## If physical core detection is unavailable, fall back to logical cores
    ## and finally to 1.
    n <- parallel::detectCores(logical = FALSE)

    if (is.na(n) || n < 1L) {
        n <- parallel::detectCores(logical = TRUE)
    }

    if (is.na(n) || n < 1L) {
        1L
    } else {
        as.integer(n)
    }
}

## This will set the maximum number of cores
## and number of threads on a given machine
## when the package is loaded
.onLoad <- function(libname, pkgname) {
    CheckLinkedVersion(pkgname)
    pkgEnv$nCores <- physicalCoreCount()
    tempThreads <- stdThreadMax()

    pkgEnv$maxThreads <- if (is.na(tempThreads) || tempThreads < 1L) {
        1L
    } else {
        as.integer(tempThreads)
    }

    invisible()
}

# nocov end
