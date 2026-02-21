test_that("linked DLL version matches package version", {
    ns <- asNamespace("RcppAlgos")
    linked_ver_ptr <- ns[[paste0("_RcppAlgos_linked_version")]]
    linked_ver <- do.call(".Call", list(linked_ver_ptr))
    expect_identical(linked_ver, as.character(packageVersion("RcppAlgos")))
})
