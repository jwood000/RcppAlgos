\name{comboGeneral}
\alias{comboGeneral}
\alias{permuteGeneral}
\title{
Generate all Combinations/Permutations of a Vector with/without Constraints
}
\description{
Quickly generate all combinations or permutations of a vector, chosen \eqn{m} at a time, with or without constraints using Rcpp.
}
\usage{
comboGeneral(v, m = NULL, repetition = FALSE, constraintFun = NULL,
                comparisonFun = NULL, limitConstraints = NULL,
                rowCap = NULL, keepResults = FALSE, freqs = NULL)

permuteGeneral(v, m = NULL, repetition = FALSE, constraintFun = NULL,
                  comparisonFun = NULL, limitConstraints = NULL,
                  rowCap = NULL, keepResults = FALSE, freqs = NULL)
}
\arguments{
  \item{v}{Source vector. If \code{v} is an integer, it will be converted to the sequence \code{1:v}.}
  \item{m}{Number of elements to choose. If \code{repetition = TRUE}, \code{m} can exceed the length of \code{v}. The default is \code{NULL}. This is to allow for combinations/permutations of multisets utilizing the \code{freqs} argument (See below).}
  \item{repetition}{Logical value indicating whether combinations should be with or without repetition. The default is \code{FALSE}.}
  \item{constraintFun}{Function to be applied to the elements of \code{v} that should be passed as a string (\emph{E.g.} \code{constraintFun = "sum"}). The possible contraint functions are: \code{"sum"}, \code{"prod"}, \code{"mean"}, \code{"max"}, & \code{"min"}. The default is \code{NULL}, meaning no function is applied.}
  \item{comparisonFun}{Comparison operator that will be used to compare \code{limitConstraints} with the result of \code{contraintFun} applied to \code{v}. It should be passed as a string (\emph{E.g.} \code{comparisonFun = "<="}). The possible comparison operators are: \code{"<"}, \code{">"}, \code{"<="}, \code{">="}, \code{"=="}. The default is \code{NULL}.}
  \item{limitConstraints}{This is the value that will be used for comparison. The default is \code{NULL}.}
  \item{rowCap}{The maximal number of expected results when a contraint is applied. Can also be used if you only need a specific number of combinations/permutations. This is useful when the total number of combinations/permutations without contraint is large and you expect/need a small number of combinations/permutations. Using \code{rowCap} can drastically improve run time and avoid unnecessary crashes due to lack of memory. See examples below.}
  \item{keepResults}{A logical flag indicating if the result of \code{constraintFun} applied to \code{v} should be displayed; if \code{TRUE}, an additional column of results will be added to the resulting matrix. The default is \code{FALSE}.}
\item{freqs}{A vector of frequencies used for producing all permutations/combinations of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.
}
}
\value{
In general, a matrix is returned with each row containing a vector of length \eqn{m} or \eqn{m + 1} depending on the value of \code{keepResults}. If \code{m} isn't supplied and \code{freqs} is given, a matrix is returned with each row containing a vector of length \code{sum(freqs)}.
}
\details{
Finding all combinations/permutations with constraints is optimized by organizing them in such a way that when \code{constraintFun} is applied, a monotonic sequence is produced. Combinations/permutations are added successively, until a particular combination exceeds the given constraint value for a given constraint/comparison function combo. After this point, we can safely skip several combinations knowing that they will exceed the given constraint value.

When there are any negative values in \code{v} and \code{constraintFun = "prod"}, producing a monotonic set is non-trivial for the general case. As a result, performance will suffer as all combinations/permutations must be tested against the constraint criteria. Additionally, \code{rowCap} will not have its normal effectiveness (\emph{i.e.} it will only limit the number of rows after producing all combinations/permutations).
}
\note{
\itemize{
\item If either \code{constraintFun}, \code{comparisonFun} or \code{limitConstraints} is \code{NULL}, the contraint check will not be carried out. This is equivalent to simply finding all combinations/permutations of \eqn{v} choose \eqn{m}.
\item The maximum number of rows is \eqn{2^{31} - 1}{2^31 - 1}.
\item Factor vectors are accepted. Class and level attributes are preserved.
}
}
\author{
Joseph Wood
}
\references{
\itemize{
\item{\href{http://gallery.rcpp.org/articles/passing-cpp-function-pointers/}{Passing user-supplied C++ functions}}
\item{\href{https://en.wikipedia.org/wiki/Heap\%27s_algorithm}{Heap's Algorithm}}
\item{\href{https://en.wikipedia.org/wiki/Monotonic_function}{Monotonic Sequence}}
\item{\href{https://en.wikipedia.org/wiki/Multiset}{Multiset}}
}
}
\examples{
system.time(comboGeneral(17,8))
system.time(permuteGeneral(13,6))

system.time(comboGeneral(13,10,repetition = TRUE))
system.time(permuteGeneral(factor(letters[1:9]),6,TRUE))

## permutations of the multiset (with or w/o setting m) :
## c(1,1,1,1,2,2,2,3,3,4,4,4,4,4)
system.time(permuteGeneral(4, freqs = c(4,3,2,5)))

permuteGeneral(4, m = 2, freqs = c(4,3,2,5))

## or combinations
comboGeneral(4, m = 2, freqs = c(4,3,2,5))

## Generate some random data
set.seed(1009)
mySamp <- rnorm(75, 997, 23)

## How to use rowCap example:
## Researcher only needs 1000 7-tuples of mySamp
## such that the sum is greater than 7200.
system.time(comboGeneral(mySamp,7,FALSE,"sum",">",7200,1000))

## If you leave rowCap NULL in examples like the above,
## it can take much longer when the total number of
## possible combinations is large (still fast enough
## most of the time). In some cases,it can crash your
## computer as the underlying code allocates enough
## space to account for every combination (e.g. In our
## example, there are choose(75, 7) = 1984829850 rows,
## 7 columns, with each cell occupying 8 bytes. This
## gives a total over 100 GB (i.e. choose(75,7)*7*8/(2^30))

## You can use rowCap with any combinatorial function
## regardless of the existence of constraints. E.g.
comboGeneral(1000, 20, rowCap = 10)

## class of the source vector is preserved
class(comboGeneral(5,3)[1,]) == class(1:5)
class(comboGeneral(c(1,2:5),3)[1,]) == class(c(1,2:5))

\dontrun{
## Using keepResults will add a column of results
permuteGeneral(-3,6,TRUE,"prod",keepResults = TRUE)
comboGeneral(-3,6,TRUE,"sum","==",-8,keepResults = TRUE)

## Find 13-tuple combinations of 1:25 such
## that the mean is less than 10
system.time(myComb <- comboGeneral(25,13,FALSE,"mean","<",10))

## Alternatively, you must generate all combinations and subsequently
## subset to obtain the combinations that meet the criteria
system.time(myComb2 <- combn(25,13))
system.time(myCols <- which(colMeans(myComb2) < 10))
system.time(myComb2 <- myComb2[,myCols])

## Any variation is much slower
system.time(myComb2 <- combn(25,13)[,combn(25,13,mean) < 10])

## Test equality with myComb above
all.equal(myComb,t(myComb2))

## Fun example... see stackoverflow:
## https://stackoverflow.com/q/22218640/4408538
system.time(permuteGeneral(seq(0L,100L,10L),8,TRUE,"sum","==",100,10^5))
}
}