context("testing comboGrid")

test_that("divisorsSieve generates correct numbers", {
    bruteCheck <- function(myList, rep = TRUE) {
        testOut <- comboGrid(myList, repetition = rep)
        
        t <- expand.grid(myList)
        t <- t[do.call(order, t), ]
        t1 <- apply(t, 1, function(x) paste0(sort(x), collapse = ""))
        t2 <- t[!duplicated(t1), ]
        
        if (!rep) {
            t2 <- t2[!apply(t2, 1, function(x) any(duplicated(x))), ]
        }
        
        if (length(unique(sapply(t2, class))) == 1) {
            t2 <- as.matrix(t2)
            dimnames(t2) <- dimnames(testOut)
        } else {
            attributes(t2) <- attributes(testOut)
        }
        
        all.equal(t2, testOut)
    }
    
    expect_equal(comboGrid(), expand.grid())
    
    myList <- list(1:5, 2:6, 3:7)
    expect_true(bruteCheck(myList, rep = TRUE))
    expect_true(bruteCheck(myList, rep = FALSE))
    
    myList <- list(1:5, 2:6, letters[1:4], letters[2:4])
    expect_true(bruteCheck(myList, rep = TRUE))
    expect_true(bruteCheck(myList, rep = FALSE))
    
    myList <- list("name1" = 1:5, 2:6)
    expect_true(bruteCheck(myList, rep = TRUE))
    
    myList <- list(1:5, factor(2:6))
    expect_true(bruteCheck(myList, rep = TRUE))
    
    myList <- list(1:5 + 0.1, 2:6 + 0.1)
    expect_true(bruteCheck(myList, rep = TRUE))
    
    myList <- list(1:5 + 0.1, factor(2:6 + 0.1))
    expect_true(bruteCheck(myList, rep = TRUE))
    
    myList <- list(letters[1:5], letters[2:6], letters[3:7])
    expect_true(bruteCheck(myList, rep = TRUE))
    expect_true(bruteCheck(myList, rep = FALSE))
    
    myList <- rep(list(c(T, F)), 10)
    expect_true(bruteCheck(myList, rep = TRUE))
})
