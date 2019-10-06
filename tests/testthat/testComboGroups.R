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
    
    expect_equal(dim(comboGroups(200, 5, upper = 100)),
                 c(100, 40, 5))
    
    expect_equal(comboGroups(14, 7, nThreads = 2), comboGroups(14, 7))
    expect_equal(sum(comboGroups(c(T, F), 2, "matrix")), 1)
    
    expect_equal(rownames(comboGroupsSample(9, 3, "matrix",
                                            sampleVec = c(67, 15, 248),
                                            namedSample = TRUE)),
                 as.character(c(67, 15, 248)))
    
    expect_equal(dim(comboGroups(as.raw(1:4), 2, "matrix")),
                 dim(comboGroups(as.complex(c(1, -1, 1i, -1i)), 2, "matrix")))
    
    expect_equal(rownames(comboGroupsSample(30, 5, n = 2,
                                            seed = 1, namedSample = TRUE)),
                 c("7162662695786451", "3525427663529072"))
    
    expect_equal(comboGroups(1:4 + 0.1, 2, "matrix"), 
                 comboGroupsSample(1:4 + 0.1, 2, "matrix", sampleVec = 1:3))
})

