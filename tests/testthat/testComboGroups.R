context("testing comboGroup")

test_that("comboGroups produces correct results", {
    expect_equal(nrow(comboGroups(12, 4, "matrix")), comboGroupsCount(12, 4))
    expect_equal(nrow(comboGroups(12, 6, "matrix")), comboGroupsCount(12, 6))
    expect_equal(nrow(comboGroups(12, 2, "matrix")), comboGroupsCount(12, 2))
    expect_equal(nrow(comboGroups(12, 3, "matrix")), comboGroupsCount(12, 3))
    expect_equal(nrow(comboGroups(12, 1, "matrix")), comboGroupsCount(12, 1))
    expect_equal(nrow(comboGroups(12, 12, "matrix")), comboGroupsCount(12, 12))
})

