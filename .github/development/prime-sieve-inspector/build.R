# Rscript build.R OUTPUT_DIRECTORY [GIT_REVISION]
# Omitting the revision uses current working-tree sources, including local edits.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || length(args) > 2L) {
    stop("Usage: Rscript build.R OUTPUT_DIRECTORY [GIT_REVISION]")
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
toolkit <- dirname(normalizePath(script, mustWork = TRUE))
root <- normalizePath(file.path(toolkit, "../../.."), mustWork = TRUE)
output <- file.path(normalizePath(dirname(args[1L]), mustWork = TRUE), basename(args[1L]))
if (file.exists(output)) stop("Choose a new output directory; existing builds are never overwritten")
if (startsWith(paste0(output, "/"), paste0(root, "/"))) {
    stop("Place the build outside the repository")
}

run <- function(command, arguments) {
    value <- system2(command, shQuote(arguments), stdout = TRUE, stderr = TRUE)
    status <- attr(value, "status")
    if (!is.null(status) && status != 0L) stop(paste(value, collapse = "\n"))
    value
}
git <- function(arguments) run("git", c("-C", root, arguments))
revision <- git(c("rev-parse", "--verify", if (length(args) == 2L) {
    paste0(args[2L], "^{commit}")
} else "HEAD"))
study <- file.path(toolkit, "studies", "2026-11-01")
pinned <- length(args) == 2L &&
    revision == readLines(file.path(study, "base-revision.txt"), warn = FALSE)
assets <- if (pinned) study else toolkit
dir.create(output)
package <- file.path(output, "package")
library <- file.path(output, "library")
dir.create(package); dir.create(library)
metadata <- list(source_revision = revision,
    mode = if (length(args) == 2L) "Git revision" else "working tree",
    toolkit = toolkit, package = package, library = library, status = "preparing",
    instrumentation_bundle = if (pinned) "studies/2026-11-01" else "current",
    patch_md5 = unname(tools::md5sum(file.path(assets, "instrumentation.patch"))))
saveRDS(metadata, file.path(output, "build.rds"))

if (length(args) == 2L) {
    archive <- file.path(output, "base.tar")
    invisible(git(c("archive", "--format=tar", paste0("--output=", archive), revision)))
    input <- file.path(output, "base")
    dir.create(input)
    utils::untar(archive, exdir = input)
    files <- list.files(input, recursive = TRUE, all.files = TRUE, no.. = TRUE)
} else {
    input <- root
    tracked <- git(c("-c", "core.quotePath=false", "ls-files", "--cached"))
    untracked <- git(c("-c", "core.quotePath=false", "ls-files", "--others", "--exclude-standard"))
    # Include new package sources, but not unrelated untracked root artifacts.
    untracked <- untracked[grepl("^(R|src|inst|tests|man|data|vignettes|exec|tools)/", untracked)]
    files <- unique(c(tracked, untracked))
}
rules <- readLines(file.path(input, ".Rbuildignore"), warn = FALSE)
rules <- c(rules[nzchar(rules)], "(^|/)\\.git(/|$)", "\\.(o|so|dll|dylib)$")
for (rule in rules) files <- files[!grepl(rule, files, perl = TRUE, ignore.case = TRUE)]
files <- files[file.exists(file.path(input, files))]
for (name in files) {
    target <- file.path(package, name)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(input, name), target)) stop("Could not copy ", name)
}
write.csv(data.frame(path = files, md5 = unname(tools::md5sum(file.path(package, files)))),
          file.path(output, "base-files.csv"), row.names = FALSE)

patch <- file.path(assets, "instrumentation.patch")
# Always check first; no fuzzy/partial application, no changes to the real checkout.
invisible(run("git", c("-C", package, "apply", "--check", patch)))
invisible(run("git", c("-C", package, "apply", patch)))
helpers <- c("src/CacheInspector.cpp", "inst/include/NumbersUtils/CacheInspector.h",
             "R/CacheInspector.R", "tests/testthat/testCacheInspector.R")
for (name in helpers) {
    target <- file.path(package, name)
    if (file.exists(target)) stop("Inspector helper already exists in source: ", name)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(assets, name), target)) stop("Could not copy ", name)
}

instrumented <- c("inst/include/NumbersUtils/Eratosthenes.h", "src/Eratosthenes.cpp",
                  "src/PrimeSieve.cpp", "src/cpp11.cpp", helpers)
actual <- data.frame(path = instrumented,
    md5 = unname(tools::md5sum(file.path(package, instrumented))))
write.csv(actual, file.path(output, "instrumented-files.csv"), row.names = FALSE)
if (pinned) {
    expected <- read.csv(file.path(assets, "instrumented-files.csv"), stringsAsFactors = FALSE)
    if (!identical(actual, expected)) stop("Pinned reconstruction did not match the preserved source hashes")
}

metadata$status <- "installing"
saveRDS(metadata, file.path(output, "build.rds"))
log <- file.path(output, "install.log")
cat("Installing instrumented copy; build log:", log, "\n")
status <- system2(file.path(R.home("bin"), "R"), shQuote(c("CMD", "INSTALL", "--preclean",
    "--no-multiarch", paste0("--library=", library), package)), stdout = log, stderr = log)
if (status != 0L) stop("Installation failed; see ", log)

metadata$status <- "testing"
saveRDS(metadata, file.path(output, "build.rds"))
validation <- file.path(output, "validate.R")
writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "library(RcppAlgos, lib.loc = args[1L])",
    "stopifnot(is.function(RcppAlgos:::primeSieveInspect))",
    "for (name in c('testCacheInspector.R', 'testPrimeSieve.R')) {",
    "    testthat::test_file(file.path(args[2L], 'tests/testthat', name),",
    "        reporter = 'summary', stop_on_failure = TRUE)",
    "}"), validation)
test_log <- file.path(output, "tests.log")
cat("Running focused tests; test log:", test_log, "\n")
status <- system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(validation, library, package)), stdout = test_log, stderr = test_log)
if (status != 0L) stop("Validation failed; see ", test_log)
metadata$status <- "validated"
metadata$completed <- format(Sys.time(), usetz = TRUE)
saveRDS(metadata, file.path(output, "build.rds"))
cat("Validated developer build.\nSource:", package, "\nLibrary:", library, "\n")
