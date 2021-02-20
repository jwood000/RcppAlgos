# comboGroups <- function(v, numGroups, retType = "matrix", lower = NULL,
#                         upper = NULL, Parallel = FALSE, nThreads = NULL) {
#     
#     ComboGroupsRcpp(v, numGroups, retType, lower, upper, is.factor(v),
#                     Parallel, nThreads, pkgEnv$nThreads, FALSE, NULL, 
#                     NULL, NULL, sample, FALSE)
# }
# 
# comboGroupsSample <- function(v, numGroups, retType = "matrix", n = NULL,
#                               sampleVec = NULL, seed = NULL, Parallel = FALSE, 
#                               nThreads = NULL, namedSample = FALSE) {
#     
#     if (!is.null(seed)) {set.seed(seed)}
#     ComboGroupsRcpp(v, numGroups, retType, NULL, NULL, is.factor(v),
#                     Parallel, nThreads, pkgEnv$nThreads, TRUE, sampleVec,
#                     seed, n, sample, namedSample)
# }
# 
# comboGroupsCount <- function(v, numGroups) {
#     ComboGroupsCountCpp(v, numGroups)
# }