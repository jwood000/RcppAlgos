context("testing comboGrid")

test_that("divisorsSieve generates correct numbers", {
    brute <- function(myList, rep = TRUE, testOut) {
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
        
        return(t2)
    }
    
    myL1 <- list(1:5, 2:6, 3:7)
    myL2 <- list(1:5, 2:6, letters[1:4], letters[2:4])
    
    myOut1 <- comboGrid(myL1)
    myOut2 <- comboGrid(myL2)
    
    expect_equal(brute(myL1, rep = TRUE, myOut1), myOut1)
    expect_equal(brute(myL2, rep = TRUE, myOut2), myOut2)
    
    myOut3 <- comboGrid(myL1, repetition = FALSE)
    myOut4 <- comboGrid(myL2, repetition = FALSE)
    
    expect_equal(brute(myL1, rep = FALSE, myOut3), myOut3)
    expect_equal(brute(myL2, rep = FALSE, myOut4), myOut4)
})
