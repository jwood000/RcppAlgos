context("testing errors in all functions")

## All error testing has been consolidated to one file for ease
## of testing r-devel-ubsan-clang on docker. Once we are ready
## to test, we simply remove this one file instead of tediously
## commenting out every error check in every test file.

test_that("comboGeneral produces appropriate error messages", {
    expect_error(comboGeneral(9,4,TRUE,NULL,NULL,NULL,"summ","<",10), "prod, sum, mean, max, or min")
    expect_error(comboGeneral(9,4,TRUE,NULL,NULL,NULL,"sum","=<>",10), "'>', '>=', '<', '<=', or '=='")
    expect_error(comboGeneral(9,4,TRUE,NULL,NULL,NULL,"sum",60,10), "must be passed as a character")
    expect_error(comboGeneral(9,4,FALSE,NULL,NULL,NULL,sum,"<",10), "must be passed as a character")
    expect_error(comboGeneral(9,4,TRUE,NULL,NULL,-1,"sum","<",10), "upper must be a positive whole number")
    expect_error(comboGeneral(170,7,FALSE,NULL,NULL,10^10,"sum","<",100), "number of rows cannot exceed")
    
    expect_error(comboGeneral(50, 5, lower = -100), "lower must be a positive whole number")
    expect_error(comboGeneral(50, 5, upper = -100), "upper must be a positive whole number")
    
    expect_error(comboGeneral(5, 50), "m must be less than or equal to the length of v")
    
    expect_error(comboGeneral(5, 3, upper = 100), "bounds cannot exceed the maximum number of possible results")
    expect_error(comboGeneral(5, 3, lower = 100), "bounds cannot exceed the maximum number of possible results")
    
    expect_error(comboGeneral(50, 3, lower = 100, upper = 10), "The number of rows must be positive")
    
    expect_error(comboGeneral(10,7,FALSE,NULL,NULL,NULL,"sum","<",c(20,30,40)), 
                 "there cannot be more than 2 limitConstraints")
    expect_error(comboGeneral(10,7,FALSE,NULL,NULL,NULL,"sum","<",c(20,NA)), 
                 "limitConstraints cannot be NA")
    
    expect_error(comboGeneral(10,7,FALSE,NULL,NULL,NULL,"sum",c("<",">","=="),c(20,30)), 
                 "there cannot be more than 2 comparison operator")
    
    expect_error(comboGeneral(10,7,FALSE,NULL,NULL,NULL,"sum",c("<","=="),c(20,30)), 
                 "If comparing against two limitConstraints, the")
    
    expect_error(comboGeneral(10,7,FALSE,NULL,NULL,NULL,"sum",c("<","=<"),c(20,30)), 
                 "Cannot have two 'less than' comparisonFuns or two")
    
    expect_error(comboGeneral(10,7,FALSE,NULL,NULL,NULL,"sum",c("<=","<"),c(20,20)), 
                 "The limitConstraints must be different")
    
    expect_error(comboGeneral(-100:100, 25, constraintFun = "prod",
                              comparisonFun = "<=", limitConstraints = 10, upper = 100), 
                 "The number of rows cannot exceed")
    
    expect_error(comboGeneral(5, 3, TRUE, constraintFun = "product"), 
                 "contraintFun must be one of the following:")
    expect_error(comboGeneral((1:5)+.01, 3, TRUE, constraintFun = "product", keepResults = TRUE), 
                 "contraintFun must be one of the following:")
    
    expect_error(comboGeneral(5,3,freqs = c(1,2,3,-2,1)), "in freqs must be a positive")
    expect_error(comboGeneral(5,1000,freqs = rep(5000, 5)), "number of rows cannot exceed")
    expect_error(comboGeneral(5,freqs = rep(1,6)), "m must be less than or equal to the length of v")
    
    numR = comboCount(1000, 10, TRUE)
    nextNum = gmp::add.bigz(numR, 1)
    expect_error(comboGeneral(1000, 10, TRUE, lower = nextNum),
                 "bounds cannot exceed the maximum number of possible results")
    expect_error(comboGeneral(1000, 10, TRUE, lower = numR, upper = nextNum),
                 "bounds cannot exceed the maximum number of possible results")
    expect_error(comboGeneral(1000, 10, TRUE, lower = -100),
                 "lower must be a positive number")
    expect_error(comboGeneral(1000, 10, TRUE, upper = -100),
                 "upper must be a positive number")
    expect_error(comboGeneral(1000, 10, TRUE, lower = 10, upper = 9),
                 "The number of rows must be positive")
    expect_error(comboGeneral(1000, 10, freqs = rep(1:4, 250)),
                 "number of rows cannot exceed")
    
    expect_error(comboGeneral(5, 3, FUN = 2), "FUN must be a function")
    
    expect_error(comboGeneral(5, 3.3), "must be a whole number")
    expect_error(comboGeneral(gmp::as.bigz(1:5), 3), "Only atomic types are supported for v")
})

test_that("comboGroups related functions produces appropriate error messages", {
    expect_error(comboGroups(10, 4), 
                 "The length of v \\(if v is a vector\\) or v \\(if v is a scalar\\) must be divisible by numGroups")
    expect_error(comboGroups(100, 4), "The number of rows cannot exceed 2\\^31 - 1.")
    expect_error(comboGroups(10, 2, retType = "crazyObj"), "retType must be '3Darray' or 'matrix'")
    expect_error(comboGroupsSample(10, 5, n = 3, namedSample = "TRUE"), 
                 "Only logical values are supported for namedSample")
    expect_error(comboGroupsSample(10, 5), 
                 "n and sampleVec cannot both be NULL")
})

test_that("divisorsRcpp produces appropriate error messages", {
    expect_error(divisorsRcpp(2^53), "The abs value of each element in v must be less than")
    expect_error(divisorsRcpp(-2^53), "The abs value of each element in v must be less than")
    expect_error(divisorsRcpp("10"), "must be of type numeric or integer")
    expect_error(divisorsRcpp(c(-2^53, 1:100)), "The abs value of each element in v must be less than")
    expect_error(divisorsRcpp(100, namedList = "TRUE"), 
                 "Only logical values are supported for namedList")
    expect_error(divisorsRcpp(100:200, nThreads = "9"), "must be of type numeric or integer")
})

test_that("divisorsSieve produces appropriate error messages", {
    expect_error(divisorsSieve(-1), "bound1 must be a positive whole number")
    expect_error(divisorsSieve(0), "bound1 must be a positive whole number")
    expect_error(divisorsSieve(2^53), "bound1 must be less than")
    expect_error(divisorsSieve(2^53, 1), "bound1 must be less than")
    expect_error(divisorsSieve(1, 2^53), "must be less than")
    expect_error(divisorsSieve("10"), "must be of type numeric or integer")
    expect_error(divisorsSieve(2, "10"), "must be of type numeric or integer")
    expect_error(divisorsSieve(100, namedList = "TRUE"), 
                 "Only logical values are supported for namedList")
    expect_error(divisorsSieve(100000, nThreads = "9"), "must be of type numeric or integer")
})

test_that("eulerPhiSieve produces appropriate error messages", {
    expect_error(eulerPhiSieve(-1), "bound1 must be a positive whole number")
    expect_error(eulerPhiSieve(0), "bound1 must be a positive whole number")
    expect_error(eulerPhiSieve(2^53), "bound1 must be less than")
    expect_error(eulerPhiSieve(2^53, 1), "must be less than")
    expect_error(eulerPhiSieve(1, 2^53), "must be less than")
    expect_error(eulerPhiSieve("10"), "must be of type numeric or integer")
    expect_error(eulerPhiSieve(2, "10"), "must be of type numeric or integer")
    expect_error(eulerPhiSieve(100, namedVector = "TRUE"), 
                 "Only logical values are supported for namedVector")
})

test_that("isPrimeRcpp produces appropriate error messages", {
    expect_error(isPrimeRcpp(0), "Each element in v must be a positive number")
    expect_error(isPrimeRcpp(-1), "Each element in v must be a positive number")
    expect_error(isPrimeRcpp(2^53), "The abs value of each element in v must be less than")
    expect_error(isPrimeRcpp("100000"), "must be of type numeric or integer")
})

test_that("numDivisorSieve produces appropriate error messages", {
    expect_error(numDivisorSieve(-1), "bound1 must be a positive whole number")
    expect_error(numDivisorSieve(0), "bound1 must be a positive whole number")
    expect_error(numDivisorSieve(2^53), "bound1 must be less than")
    expect_error(numDivisorSieve(2^53, 1), "must be less than")
    expect_error(numDivisorSieve(1, 2^53), "must be less than")
    expect_error(numDivisorSieve("10"), "must be of type numeric or integer")
    expect_error(numDivisorSieve(2, "10"), "must be of type numeric or integer")
    expect_error(numDivisorSieve(100, namedVector = "TRUE"), "Only logical values are supported for namedVector")
})

test_that("combo/permuteGeneral produces correct error messages with Parallel", {
    expect_error(permuteGeneral(10, 5, Parallel = "TRUE"), 
                 "Only logical values are supported for Parallel")
})

test_that("permuteGeneral produces appropriate error messages", {
    expect_error(permuteGeneral(100, 15, upper = 2^32), 
                 "The number of rows cannot exceed")
    expect_error(permuteGeneral(100, 15, lower = 100, upper = 2^32), 
                 "The number of rows cannot exceed")
    expect_error(permuteGeneral(100, 15, lower = 2^32), 
                 "The number of rows cannot exceed")
    expect_error(permuteGeneral(9,4,TRUE,constraintFun = "summ",
                                comparisonFun = "<",limitConstraints = 10),
                 "prod, sum, mean, max, or min")
    expect_error(permuteGeneral(9,4,TRUE,constraintFun = "sum",
                                comparisonFun = "=<>",limitConstraints = 10),
                 "'>', '>=', '<', '<=', or '=='")
    expect_error(permuteGeneral(9,4,TRUE,constraintFun = "sum",
                                comparisonFun = 60,limitConstraints = 10),
                 "must be passed as a character")
    expect_error(permuteGeneral(9,4,FALSE,constraintFun = sum,
                                comparisonFun = "<",limitConstraints = 10),
                 "must be passed as a character")
    expect_error(permuteGeneral(9,4,TRUE,constraintFun = "sum",
                                comparisonFun = "<",limitConstraints = 10,upper = -1),
                 "upper must be a positive whole number")
    expect_error(permuteGeneral(170,7,FALSE,constraintFun = "sum",
                                comparisonFun = "<",
                                limitConstraints = 100,
                                upper = 10^10), "number of rows cannot exceed")
    expect_error(permuteGeneral(5, 1:5), "length of m must be 1")
    expect_error(permuteGeneral(5, -5), "m must be a positive whole number")
    expect_error(permuteCount(5, 5, "TRUE"), "Only logical values are supported for repetition")
    
    expect_error(permuteGeneral(5,3,freqs = c(1,2,3,-2,1)), "in freqs must be a positive")
    expect_error(permuteGeneral(5,15,freqs = c(5,5,5,5,5)), "number of rows cannot exceed")
    expect_error(permuteGeneral(5,freqs = c(5,5,5,5,5)), "number of rows cannot exceed")
    expect_error(permuteGeneral(5,freqs = rep(1,6)), "m must be less than or equal to the length of v")
    
    numR <- permuteCount(1000, 10, TRUE)
    nextNum <- gmp::add.bigz(numR, 1)
    expect_error(permuteGeneral(1000, 10, TRUE, lower = nextNum),
                 "bounds cannot exceed the maximum number of possible results")
    expect_error(permuteGeneral(1000, 10, TRUE, lower = numR, upper = nextNum),
                 "bounds cannot exceed the maximum number of possible results")
    expect_error(permuteGeneral(1000, 10, TRUE, lower = -100),
                 "lower must be a positive number")
    expect_error(permuteGeneral(1000, 10, TRUE, upper = -100),
                 "upper must be a positive number")
    expect_error(permuteGeneral(1000, 10, TRUE, lower = 10, upper = 9),
                 "The number of rows must be positive")
    expect_error(permuteGeneral(1000, 10, freqs = rep(1:4, 250)),
                 "number of rows cannot exceed")
    
    expect_error(permuteGeneral(5, 3, FUN = 2), "FUN must be a function")
})

test_that("primeCount produces appropriate error messages", {
    expect_error(primeCount(0), "must be a positive")
    expect_error(primeCount(-1), "must be a positive")
    expect_error(primeCount(2^53), "n must be less than")
    expect_error(primeCount("100000"), "must be of type numeric or integer")
})

test_that("primeFactorize produces appropriate error messages", {
    expect_error(primeFactorize(2^53), "The abs value of each element in v must be less than")
    expect_error(primeFactorize(-2^53), "The abs value of each element in v must be less than")
    expect_error(primeFactorize(c(-2^53, 1:100)), "The abs value of each element in v must be less than")
    expect_error(primeFactorize("10"), "must be of type numeric or integer")
    expect_error(primeFactorize(100, namedList = "TRUE"), "Only logical values are supported for namedList")
})

test_that("primeFactorizeSieve produces appropriate error messages", {
    expect_error(primeFactorizeSieve(-1), "must be a positive whole number")
    expect_error(primeFactorizeSieve(2^53), "must be less than")
    expect_error(primeFactorizeSieve(2^53, 1), "must be less than")
    expect_error(primeFactorizeSieve(1, 2^53), "must be less than")
    expect_error(primeFactorizeSieve("10"), "must be of type numeric or integer")
    expect_error(primeFactorizeSieve(2, "10"), "must be of type numeric or integer")
    expect_error(primeFactorizeSieve(100, namedList = "TRUE"), "Only logical values are supported for namedList")
})

test_that("primeSieve produces appropriate error messages", {
    expect_error(primeSieve(-1), "must be a positive")
    expect_error(primeSieve(1,-1), "must be a positive whole number")
    expect_error(primeSieve(1,2^53), "must be less than")
    expect_error(primeSieve(2^53), "must be less than")
    expect_error(primeSieve(2^53, 1), "must be less than")
    expect_error(primeSieve(2^4, "1"), "must be of type numeric or integer")
    expect_error(primeSieve("500"), "must be of type numeric or integer")
})

test_that("comboSample produces appropriate error messages", {
    expect_error(comboSample(5, 3), "n and sampleVec cannot both be NULL")
    expect_error(comboSample(5,3,freqs = c(1,2,3,-2,1)), "in freqs must be a positive")
    expect_error(comboSample(5,3, n = 100), "n exceeds the maximum number of possible results")
    expect_error(comboSample(5,3, comboSample(5,3, sampleVec = 1:200)), "exceeds the maximum number of possible results")
    expect_error(comboSample(5,freqs = rep(1,6)), "m must be less than or equal to the length of v")
    expect_error(comboSample(5,3, n = 5, FUN = "sum"), "FUN must be a function!")
    expect_error(comboSample(5,3, n = "5"), "n must be of type numeric or integer")
    expect_error(comboSample(5,3, n = 1:5), "length of n must be 1")
    expect_error(comboSample(100000, 10, sampleVec = -1L),
                 "sampleVec must be a positive number")
    expect_error(comboSample(100000, 10, sampleVec = NaN),
                 "sampleVec cannot be NA or NaN")
    expect_error(comboSample(100000, 10, sampleVec = 2^1000),
                 "Number is too large for double precision. Consider using gmp::as.bigz or as.character for sampleVec")
    expect_error(comboSample(100000, 10, sampleVec = as.complex(1)),
                 "This type is not supported! No conversion possible for sampleVec")
    
    expect_error(comboSample(1000, 100, sampleVec = gmp::as.bigz(c(1, -1))),
                 "Each element in sampleVec must be a positive number")
    expect_error(comboSample(1000, 100, sampleVec = NA_character_), "sampleVec cannot be NA or NaN")
    expect_error(comboSample(1000, 100, sampleVec = 2^53), 
                 "Number is too large for double precision. Consider using gmp::as.bigz or as.character for sampleVec")
    expect_error(comboSample(1000, 15, sampleVec = gmp::as.bigz(c(NA, 1))),
                 "Each element in sampleVec cannot be NA or NaN")
})

test_that("permuteSample produces appropriate error messages", {
    expect_error(permuteSample(5, 6, n = 3), "m must be less than or equal to the length of v")
    expect_error(permuteSample(5, 1:5), "length of m must be 1")
    expect_error(permuteSample(5, -4), "m must be a positive whole number")
    expect_error(permuteSample(5, 3), "n and sampleVec cannot both be NULL")
    expect_error(permuteSample(5,3,freqs = c(1,2,3,-2,1)), "in freqs must be a positive")
    expect_error(permuteSample(5,3, n = 100), 
                 "n exceeds the maximum number of possible results")
    expect_error(permuteSample(5,3, n = "5"), "n must be of type numeric or integer")
    expect_error(permuteSample(5,3, n = 1:5), "length of n must be 1")
    expect_error(permuteSample(5,3, permuteSample(5,3, sampleVec = 1:200)), 
                 "exceeds the maximum number of possible results")
    expect_error(permuteSample(5,freqs = rep(1,6)), "m must be less than or equal to the length of v")
    expect_error(permuteSample(5, 4, sampleVec = "adfs"), 
                 "sampleVec must be a positive whole number")
    expect_error(permuteSample(5, 4, sampleVec = 1.1), 
                 "Each element in sampleVec must be a whole number")
    expect_error(permuteSample(5000, 40, sampleVec = 1.1),
                 "sampleVec must be a whole number")
    expect_error(permuteSample(100, 10, n = 2^50), "n must be less than or equal to 2147483647")
    expect_error(permuteSample(100, 10, n = -200), "n must be a positive whole number")
    expect_error(permuteSample(100, 10, sampleVec = c("62815650955529472001")),
                 "One or more of the requested values in sampleVec")
    expect_error(permuteSample(3.3, 2, n = 1), 
                 "If v is not a character and of length 1, it must be a whole number")
    expect_error(permuteSample(3, 2.2, n = 1), "m must be a whole number")
    expect_error(permuteSample(3, 2, n = 1.1), "n must be a whole number")
    expect_error(permuteSample(3, 2, sampleVec = as.complex(1)), 
                 "This type is not supported! No conversion possible for sampleVec")
    expect_error(permuteSample(3, 2, Parallel = "TRUE"), 
                 "Only logical values are supported for Parallel")
    expect_error(permuteSample(5, 3, n = 5, nThreads = "5"), 
                 "nThreads must be of type numeric or integer")
    expect_error(permuteSample(5, 3, n = 5, nThreads = 3.2), "nThreads must be a whole number")
    expect_error(permuteSample(5000, 10, n = 5.5), "n must be a whole number")
    expect_error(permuteSample(5000, 10, n = NA_integer_), "n cannot be NA or NaN")
    expect_error(permuteSample(5000, 10, n = NA), 
                 "This type is not supported! No conversion possible for n")
    expect_error(permuteSample(5000, 10, sampleVec = -1), 
                 "sampleVec must be a positive number")
    expect_error(permuteSample(5000, 10, sampleVec = c(1e9, -1)), 
                 "Each element in sampleVec must be a positive number")
    expect_error(permuteSample(NA_integer_, 3, n = 5), 
                 "If v is not a character and of length 1, it cannot be NA or NaN")
    expect_error(permuteSample(1000, 20, sampleVec = c(NA, "1234567890")),
                 "Each element in sampleVec cannot be NA or NaN")
})

test_that("combo/permuteGeneral produces appropriate error messages for subset sum", {
    expect_error(comboGeneral(0:130, 130, TRUE,
                              constraintFun = "sum",
                              comparisonFun = "==",
                              limitConstraints = 130), "The number of rows cannot exceed")
    
    expect_error(permuteGeneral(35, 18, TRUE,
                                constraintFun = "sum",
                                comparisonFun = "==",
                                limitConstraints = 35), "The number of rows cannot exceed")
})
