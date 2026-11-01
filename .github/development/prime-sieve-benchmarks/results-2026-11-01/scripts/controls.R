# Compare ordinary and inspected calls across all kernels, using normal policy.
# Run separately from the main suite, with no concurrent benchmark processes.
# Rscript controls.R LIBRARY OUTPUT
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
library(RcppAlgos, lib.loc = normalizePath(args[1]))
out <- args[2]
inspect <- RcppAlgos:::primeSieveInspect
set.seed(20261102)
cases <- list(small_dense = c(1, 1e8),
              small_interior = c(1e8 + 17, 1e8 + 1e6 + 17),
              medium_short = c(1e12 + 17, 1e12 + 1e6 + 17),
              medium_wide = c(1e12 + 17, 1e12 + 1e8 + 17),
              big_short = c(4e15 + 17, 4e15 + 1e6 + 17),
              big_wide = c(4e15 + 17, 4e15 + 1e8 + 17))
rows <- list(); phases <- list(); index <- 0L
for (name in names(cases)) for (threads in c(1L, 4L)) {
    bounds <- cases[[name]]
    reference <- primeSieve(bounds[1], bounds[2], nThreads = 1L)
    sample_result <- inspect(bounds[1], bounds[2], nThreads = threads)
    stopifnot(identical(sample_result$primes, reference))
    phases[[length(phases) + 1L]] <- data.frame(
        case = name, requested = threads,
        actual = sample_result$summary$effective_threads,
        kernel = paste(unique(sample_result$kernels$kernel), collapse = "+"),
        sample_result$phases)
    rm(sample_result); gc()
    for (round in 1:15) for (mode in sample(c("ordinary", "inspected"))) {
        # Sum separately timed calls; comparisons and disposal are outside.
        # Batching reduces the relative effect of the ~1 ms R timer resolution.
        elapsed <- 0
        for (iteration in 1:10) {
            start <- proc.time()[["elapsed"]]
            result <- if (mode == "ordinary") {
                primeSieve(bounds[1], bounds[2], nThreads = threads)
            } else {
                inspect(bounds[1], bounds[2], nThreads = threads)
            }
            elapsed <- elapsed + proc.time()[["elapsed"]] - start
            actual <- if (mode == "ordinary") result else result$primes
            stopifnot(identical(actual, reference))
            rm(result, actual)
        }
        index <- index + 1L
        rows[[index]] <- data.frame(case = name, lower = bounds[1], upper = bounds[2],
            requested = threads, mode, round, calls = 10L, per_call_ms = elapsed * 100)
    }
    cat(name, threads, "complete\n"); flush.console()
    write.csv(do.call(rbind, rows), file.path(out, "controls.csv"), row.names = FALSE)
}
write.csv(do.call(rbind, phases), file.path(out, "control-phases.csv"), row.names = FALSE)
writeLines(c(capture.output(sessionInfo()),
    paste("Policy cores:", get("pkgEnv", asNamespace("RcppAlgos"))$nCores),
    "15 paired randomized rounds; 10 calls per mode per round; seed 20261102.",
    "Timing includes R allocation, garbage collection when triggered, and inspector report assembly."),
    file.path(out, "controls-metadata.txt"))
