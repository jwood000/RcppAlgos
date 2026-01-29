\name{comboGeneral}
\alias{comboGeneral}
\alias{comboGeneral.integer}
\alias{comboGeneral.numeric}
\alias{comboGeneral.factor}
\alias{comboGeneral.logical}
\alias{comboGeneral.default}
\alias{comboGeneral.table}
\alias{comboGeneral.list}
\alias{permuteGeneral}
\alias{permuteGeneral.integer}
\alias{permuteGeneral.numeric}
\alias{permuteGeneral.factor}
\alias{permuteGeneral.logical}
\alias{permuteGeneral.default}
\alias{permuteGeneral.table}
\alias{permuteGeneral.list}
\title{
Generate Combinations and Permutations of a Vector with/without Constraints
}
\description{
    \itemize{
        \item Generate combinations or permutations of a vector with or without constraints.
        \item Produce results in parallel using the \code{Parallel} or \code{nThreads} arguments. You can also apply each of the five compiled functions given by the argument \code{constraintFun} in parallel.
        \item The arguments \code{lower} and \code{upper} make it possible to generate combinations/permutations in chunks allowing for parallelization via the \code{\link{parallel-package}}. This is convenient when you want to apply a custom function to the output in parallel.
        \item Attack integer partition and general subset sum problems.
        \item GMP support allows for exploration of combinations/permutations of vectors with many elements.
        \item The output is in lexicographical order.
    }
}
\usage{
comboGeneral(v, m = NULL, ...)
permuteGeneral(v, m = NULL, ...)

\method{comboGeneral}{numeric}(v, m = NULL, repetition = FALSE, freqs = NULL,
             lower = NULL, upper = NULL, constraintFun = NULL,
             comparisonFun = NULL, limitConstraints = NULL,
             keepResults = NULL, FUN = NULL, Parallel = FALSE,
             nThreads = NULL, tolerance = NULL, FUN.VALUE = NULL, ...)

\method{permuteGeneral}{numeric}(v, m = NULL, repetition = FALSE, freqs = NULL,
               lower = NULL, upper = NULL, constraintFun = NULL,
               comparisonFun = NULL, limitConstraints = NULL,
               keepResults = NULL, FUN = NULL, Parallel = FALSE,
               nThreads = NULL, tolerance = NULL, FUN.VALUE = NULL, ...)

\method{comboGeneral}{factor}(
    v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL, upper = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, FUN.VALUE = NULL, ...
)
\method{permuteGeneral}{factor}(
    v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL, upper = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, FUN.VALUE = NULL, ...
)

\method{comboGeneral}{default}(v, m = NULL, repetition = FALSE,
             freqs = NULL, lower = NULL, upper = NULL,
             FUN = NULL, FUN.VALUE = NULL, ...)
\method{permuteGeneral}{default}(v, m = NULL, repetition = FALSE,
               freqs = NULL, lower = NULL, upper = NULL,
               FUN = NULL, FUN.VALUE = NULL, ...)

\method{comboGeneral}{table}(
    v, m = NULL, lower = NULL, upper = NULL, constraintFun = NULL,
    comparisonFun = NULL, limitConstraints = NULL, keepResults = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, tolerance = NULL,
    FUN.VALUE = NULL, ...
)
\method{permuteGeneral}{table}(
    v, m = NULL, lower = NULL, upper = NULL, constraintFun = NULL,
    comparisonFun = NULL, limitConstraints = NULL, keepResults = NULL,
    FUN = NULL, Parallel = FALSE, nThreads = NULL, tolerance = NULL,
    FUN.VALUE = NULL, ...
)

\method{comboGeneral}{list}(v, m = NULL, repetition = FALSE,
             freqs = NULL, lower = NULL, upper = NULL, ...)
\method{permuteGeneral}{list}(v, m = NULL, repetition = FALSE,
               freqs = NULL, lower = NULL, upper = NULL, ...)
}
\arguments{
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. All atomic types are supported (See \code{\link{is.atomic}}).}
  \item{m}{Number of elements to choose.  If \code{repetition = TRUE} or \code{freqs} is utilized, \code{m} can exceed the length of \code{v}. If \code{m = NULL}, the length will default to \code{length(v)} or \code{sum(freqs)}.}
  \item{\dots}{Further arguments passed to methods.}
  \item{repetition}{Logical value indicating whether combinations/permutations should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all combinations/permutations of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
  \item{lower}{The lower bound. Combinations/permutations are generated lexicographically, thus utilizing this argument will determine which specific combination/permutation to start generating from (\emph{e.g.} \code{comboGeneral(5, 3, lower = 6)} is equivalent to \code{comboGeneral(5, 3)[6:choose(5, 3), ]}). This argument along with \code{upper} is very useful for generating combinations/permutations in chunks allowing for easy parallelization.}
  \item{upper}{The upper bound. Similar to \code{lower}, however this parameter allows the user to \emph{stop} generation at a specific combination/permutation (\emph{e.g.} \code{comboGeneral(5, 3, upper = 5)} is equivalent to \code{comboGeneral(5, 3)[1:5, ]})

  If the output is constrained and \code{lower} isn't supplied, \code{upper} serves as a cap for how many results will be returned that meet the criteria (\emph{e.g.} setting \code{upper = 100} alone will return the first 100 results that meet the criteria, while setting \code{lower = 1} and \code{upper = 100} will test the first 100 results against the criteria).

  In addition to the benefits listed for \code{lower}, this parameter is useful when the total number of combinations/permutations without constraint is large and you expect/need a small number of combinations/permutations that meet a certain criteria. Using \code{upper} can improve run time if used judiciously as we call the member function \href{https://en.cppreference.com/w/cpp/container/vector/reserve.html}{reserve} of \href{https://en.cppreference.com/w/cpp/container/vector.html}{std::vector}. See examples below.}
  \item{constraintFun}{Function to be applied to the elements of \code{v} that should be passed as a string (\emph{e.g.} \code{constraintFun = "sum"}). The possible constraint functions are: \code{"sum"}, \code{"prod"}, \code{"mean"}, \code{"max"}, & \code{"min"}. The default is \code{NULL}, meaning no function is applied.}
  \item{comparisonFun}{Comparison operator that will be used to compare \code{limitConstraints} with the result of \code{constraintFun} applied to \code{v}. It should be passed as a string or a vector of two strings (\emph{e.g.} \code{comparisonFun = "<="} or  \code{comparisonFun = c(">","<")}). The possible comparison operators are: \code{"<"}, \code{">"}, \code{"<="}, \code{">="}, \code{"=="}. The default is \code{NULL}.

  When \code{comparisonFun} is a vector of two comparison strings, \emph{e.g} \code{comparisonFun = c(comp1, comp2)}, and \code{limitConstraints} is a vector of two numerical values, \emph{e.g} \code{limitConstraints = c(x1, x2)}, the combinations/permutations will be filtered in one of the following two ways:
  \enumerate{
      \item{When \code{comp1} is one of the 'greater-than' operators (\emph{i.e.} ">=" or ">"), \code{comp2} is one of the 'less-than' operators (\emph{i.e.} "<=" or "<"), and \code{x1 < x2}, the combinations/permutations that are returned will have a value (after \code{constraintFun} has been applied) between \code{x1} and \code{x2}.}
      \item{When \code{comp1} and \code{comp2} are defined as in #1 and \code{x1 > x2}, the combinations/permutations that are returned will have a value outside the range of \code{x1} and \code{x2}. See the examples below.}
        }
  In other words, the first comparison operator is applied to the first limit and the second operator is applied to the second limit.}
  \item{limitConstraints}{This is the value(s) that will be used for comparison. Can be passed as a single value or a vector of two numerical values. The default is \code{NULL}. See the definition of \code{comparisonFun} as well as the examples below for more information.}
  \item{keepResults}{A logical flag indicating if the result of \code{constraintFun} applied to \code{v} should be displayed; if \code{TRUE}, an additional column of results will be added to the resulting matrix. The default is \code{FALSE}. If user is only applying \code{constraintFun}, \code{keepResults} will default to \code{TRUE}.

      \emph{E.g}. The following are equivalent and will produce a \eqn{4^{th}}{4th} column of row sums:
      \itemize{
            \item \code{comboGeneral(5, 3 constraintFun = "sum", keepResults = TRUE)}
            \item \code{comboGeneral(5, 3 constraintFun = "sum")}
        }
    }
  \item{FUN}{Function to be applied to each combination/permutation. The default is \code{NULL}.}
  \item{Parallel}{Logical value indicating whether combinations/permutations should be generated in parallel using \eqn{n - 1} threads, where \eqn{n} is the maximum number of threads. The default is \code{FALSE}. If \code{nThreads} is not \code{NULL}, it will be given preference (\emph{e.g.} if user has 8 threads with \code{Parallel = TRUE} and \code{nThreads = 4}, only 4 threads will be spawned). If your system is single-threaded, the arguments \code{Parallel} and \code{nThreads} are ignored.}
  \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}. See \code{Parallel}.}
  \item{tolerance}{A numeric value greater than or equal to zero. This parameter is utilized when a constraint is applied on a numeric vector. The default value is 0 when it can be determined that whole values are being utilized, otherwise it is \code{sqrt(.Machine$double.eps)} which is approximately \eqn{1.5e-8}. N.B. If the input vector is of type integer, this parameter will be ignored and strict equality will be enforced.}
  \item{FUN.VALUE}{A template for the return value from \code{FUN}. See 'Details' of \code{\link[base]{vapply}} for more information.}
}
\value{
    \itemize{
        \item In general, a matrix with \eqn{m} or \eqn{m + 1} columns, depending on the value of \code{keepResults}
        \item If \code{FUN} is utilized and \code{FUN.VALUE = NULL}, a list is returned
        \item When both \code{FUN} and \code{FUN.VALUE} are not \code{NULL}, the return is modeled after the return of \code{vapply}. See the 'Value' section of \code{\link[base]{vapply}}.
    }
}
\details{
For the general case, finding all combinations/permutations with constraints is optimized by organizing them in such a way that when \code{constraintFun} is applied, a \emph{partially} monotonic sequence is produced. Combinations/permutations are added successively, until a particular combination exceeds the given constraint value for a given constraint/comparison function combo. After this point, we can safely skip several combinations knowing that they will exceed the given constraint value.

There are special cases where more efficient algorithms are dyncamically deployed. These cases center around the subject of integer partitions. See \code{\link{partitionsGeneral}} for more information.

When there are any negative values in \code{v} and \code{constraintFun = "prod"}, producing a monotonic set is non-trivial for the general case. As a result, performance will suffer as all combinations/permutations must be tested against the constraint criteria.
}
\note{
    \itemize{
        \item {\code{Parallel} and \code{nThreads} will be ignored in the following cases:
            \itemize{
                \item When the output is constrained (except for most partitions cases)
                \item If the class of the vector passed is \code{character}, \code{raw}, and \code{complex} (N.B. \code{Rcpp::CharacterMatrix} is not thread safe). Alternatively, you can generate an indexing matrix in parallel.
                \item If \code{FUN} is utilized.
            }
        }
        \item If either \code{constraintFun}, \code{comparisonFun} or \code{limitConstraints} is \code{NULL} --or-- if the class of the vector passed is \code{logical}, \code{character}, \code{raw}, \code{factor}, or \code{complex}, the constraint check will not be carried out. This is equivalent to simply finding all combinations/permutations of \eqn{v} choose \eqn{m}.
        \item The maximum number of combinations/permutations that can be generated at one time is \eqn{2^{31} - 1}{2^31 - 1}. Utilizing \code{lower} and \code{upper} makes it possible to generate additional combinations/permutations.
        \item Factor vectors are accepted. Class and level attributes are preserved except when \code{FUN} is used.
        \item Lexicographical ordering isn't guaranteed for permutations if \code{lower} isn't supplied and the output is constrained.
        \item If \code{lower} is supplied and the output is constrained, the combinations/permutations that will be tested will be in the lexicographical range \code{lower} to \code{upper} or up to the total possible number of results if \code{upper} is not given. See the second paragraph for the definition of \code{upper}.
        \item \code{FUN} will be ignored if the constraint check is satisfied.
    }
}
\author{
Joseph Wood
}
\references{
    \itemize{
        \item{\href{https://gallery.rcpp.org/articles/passing-cpp-function-pointers/}{Passing user-supplied C++ functions}}
        \item{\href{https://en.wikipedia.org/wiki/Monotonic_function}{Monotonic Sequence}}
        \item{\href{https://en.wikipedia.org/wiki/Multiset}{Multiset}}
        \item{\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical Order}}
        \item{\href{https://en.wikipedia.org/wiki/Subset_sum_problem}{Subset Sum Problem}}
        \item{\href{https://en.wikipedia.org/wiki/Partition_(number_theory)}{Partition (Number Theory)}}
    }
}
\examples{
comboGeneral(4, 3)
permuteGeneral(3)

permuteGeneral(factor(letters[1:3]), repetition = TRUE)

## permutations of the multiset :
## c(1,1,1,2,2,3)
permuteGeneral(table(c(1,1,1,2,2,3)))

## Example with list
comboGeneral(
    v = list(
        p1 = matrix(1:10, ncol = 2),
        p2 = data.frame(a = letters, b = 1:26),
        p3 = as.complex(1:10)
    ),
    m = 2
)

#### Examples using "upper" and "lower":
## See specific range of permutations
permuteGeneral(75, 10, freqs = rep(1:3, 25),
               lower = 1e12, upper = 1e12 + 10)

## Researcher only needs 10 7-tuples of mySamp
## such that the sum is greater than 7200.
## Generate some random data
set.seed(1009)
mySamp = rnorm(75, 997, 23)
comboGeneral(mySamp, 7, constraintFun = "sum",
             comparisonFun = ">", limitConstraints = 7200, upper = 10)

## Similarly, you can use "lower" to obtain the last rows.
## Generate the last 10 rows
comboGeneral(mySamp, 7, lower = choose(75, 7) - 9)

## Or if you would like to generate a specific chunk,
## use both "lower" and "upper". E.g. Generate one
## million combinations starting with the 900,000,001
## lexicographic combination.
t1 = comboGeneral(mySamp, 7,
                  lower = 9*10^8 + 1,
                  upper = 9*10^8 + 10^6)

## class of the source vector is preserved
class(comboGeneral(5,3)[1,]) == class(1:5)
class(comboGeneral(c(1,2:5),3)[1,]) == class(c(1,2:5))
class(comboGeneral(factor(month.name),3)[1,]) == class(factor(month.name))

## Using keepResults will add a column of results
comboGeneral(-3, 6, TRUE,
             constraintFun = "sum",
             comparisonFun = "==",
             limitConstraints = -8,
             keepResults = TRUE)

## Using multiple constraints:

## Get combinations such that the product
## is between 3000 and 4000 inclusive
comboGeneral(5, 7, TRUE, constraintFun = "prod",
             comparisonFun = c(">=","<="),
             limitConstraints = c(3000, 4000),
             keepResults = TRUE)

## Or, get the combinations such that the
## product is less than or equal to 10 or
## greater than or equal to 40000
comboGeneral(5, 7, TRUE, constraintFun = "prod",
             comparisonFun = c("<=",">="),
             limitConstraints = c(10, 40000),
             keepResults = TRUE)

#### General subset sum problem
set.seed(516781810)
comboGeneral(runif(100, 0, 42), 5, constraintFun = "mean",
             comparisonFun = "==", limitConstraints = 30,
             tolerance = 0.0000002)


#### Integer Partitions
comboGeneral(0:5, 5, TRUE, constraintFun = "sum",
             comparisonFun = "==", limitConstraints = 5)


## Using FUN
comboGeneral(10000, 5, lower = 20, upper = 22,
             FUN = function(x) {
                 which(cummax(x) \%\% 2 == 1)
             })

\dontrun{
## Parallel example generating more than 2^31 - 1 combinations.
library(parallel)
numCores = detectCores() - 1

## 10086780 evenly divides choose(35, 15) and is "small enough" to
## generate quickly in chunks.
system.time(mclapply(seq(1, comboCount(35, 15), 10086780), function(x) {
    a = comboGeneral(35, 15, lower = x, upper = x + 10086779)
    ## do something
    x
}, mc.cores = numCores))


## Find 13-tuple combinations of 1:25 such
## that the mean is less than 10
system.time(myComb <- comboGeneral(25, 13, FALSE,
                                   constraintFun = "mean",
                                   comparisonFun = "<",
                                   limitConstraints = 10))

## Alternatively, you must generate all combinations and subsequently
## subset to obtain the combinations that meet the criteria
system.time(myComb2 <- combn(25, 13))
system.time(myCols <- which(colMeans(myComb2) < 10))
system.time(myComb2 <- myComb2[, myCols])

## Any variation is much slower
system.time(myComb2 <- combn(25, 13)[,combn(25, 13, mean) < 10])

## Test equality with myComb above
all.equal(myComb, t(myComb2))

## Fun example... see stackoverflow:
## https://stackoverflow.com/q/22218640/4408538
system.time(permuteGeneral(seq(0L,100L,10L), 8, TRUE,
                           constraintFun = "sum",
                           comparisonFun = "==",
                           limitConstraints = 100))

## These are called weak integer compositions. Below, we call
## compositionsGeneral which gives the same output except it
## in lexicographical order. See 'Note' above
system.time(compositionsGeneral(seq(0L,100L,10L), 8, TRUE, weak = TRUE))
}
}
\keyword{combinations}
\keyword{permutations}
\keyword{combinatorics}
