# nocov start

pkgEnv <- new.env(parent = emptyenv())
pkgEnv$nCores <- NULL
pkgEnv$maxThreads <- NULL

physicalCoreCount <- function() {
    ## Use R's default core-detection path only. This is a conservative
    ## package fallback, not an attempt to determine exact CPU topology.
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
