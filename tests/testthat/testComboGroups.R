context("testing comboGroup")

test_that("comboGroups produces correct results", {
    expect_equal(nrow(comboGroups(letters[1:12], 4, "matrix")), comboGroupsCount(12, 4))
    expect_equal(nrow(comboGroups(factor(letters[1:12]), 3, "matrix")),
                 comboGroupsCount(12, 3))
    expect_equal(nrow(comboGroups(12, 1, "matrix")), comboGroupsCount(12, 1))
    expect_equal(nrow(comboGroups(12, 12, "matrix")), comboGroupsCount(12, 12))
    
    # comboGroupsCount(50, 5)
    # Big Integer ('bigz') :
    # [1] 402789797982510165934296910320
    expect_equal(dim(comboGroups(50, 5, lower = "402789797982510165934296910301")),
                 c(20, 10, 5))
    expect_equal(comboGroups(10, 5, lower = 201, upper = 300),
                 comboGroups(10, 5)[201:300, ,])
    
    expect_equal(comboGroups(12, 4), 
                 comboGroupsSample(12, 4, sampleVec = 1:comboGroupsCount(12, 4)))
    
    expect_equal(comboGroups(12, 4), 
                 comboGroups(12, 4, nThreads = 2))
    expect_equal(comboGroups(12, 4), 
                 comboGroupsSample(12, 4, sampleVec = 1:comboGroupsCount(12, 4), nThreads = 2))
})

