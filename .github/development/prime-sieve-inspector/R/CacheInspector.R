# Internal developer helper; see .github/development/cache-inspector.md.
primeSieveInspect <- function(bound1, bound2 = NULL, nThreads = NULL,
                              segmentBytes = NULL) {
    result <- .Call(`_RcppAlgos_PrimeSieveInspectCpp`, bound1, bound2,
                    nThreads, pkgEnv$nCores, pkgEnv$maxThreads, segmentBytes)
    result$host <- list(
        system = Sys.info()[c("sysname", "release", "machine")],
        R = R.version.string,
        physical_cores = pkgEnv$nCores,
        max_threads = pkgEnv$maxThreads
    )
    result
}
