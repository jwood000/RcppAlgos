\name{comboGeneral}
\alias{comboGeneral}
\title{
Generate all Combinations of a Vector with/without Constraints
}
\description{
Quickly generate all combinations of a vector, chosen \eqn{m} at a time, with or without constraints using Rcpp.
}
\usage{
comboGeneral(v, m, repetition = FALSE, constraintFun = NULL,
                comparisonFun = NULL, limitConstraints = NULL, rowCap = NULL)
}
\arguments{
  \item{v}{Source vector. If \code{v} is an integer, it will be converted to the sequence \code{1:v}.}
  \item{m}{Number of elements to choose.}
  \item{repetition}{Logical value indicating whether combinations should be with or without repetition. The default is \code{FALSE}.}
  \item{constraintFun}{Function to be applied to the elements of \code{v} that should be passed as a string (\emph{E.g.} \code{constraintFun = "sum"}). The possible contraint functions are: \code{"sum"}, \code{"prod"}, \code{"mean"}, \code{"max"}, & \code{"min"}. The default is \code{NULL}, meaning no function is applied.}
  \item{comparisonFun}{Comparison operator that will be used to compare \code{limitConstraints} with the result of \code{contraintFun} applied to \code{v}. Again, it should be passed as a string (\emph{E.g.} \code{comparisonFun = "<="}). The possible comparison operators are: \code{"<"}, \code{">"}, \code{"<="}, \code{">="}, \code{"=="}. The default is \code{NULL}.}
  \item{limitConstraints}{This is the value that will be used for comparison. The default is \code{NULL}.}
  \item{rowCap}{The maximal number of expected results when a contraint is applied. Can also be used if you only need a specific number of combinations. This is useful when the total number of combinations without contraint is large and you expect/need a small number of combinations that meet the criteria. Using rowCap can drastically improve run time and avoid unnecessary crashes due to lack of memory. See examples below.}
}
\references{
\url{http://gallery.rcpp.org/articles/passing-cpp-function-pointers/}
}
\value{
Returns a matrix where each row contains a vector of length \code{m}.
}
\author{
Joseph Wood
}
\note{
If either \code{constraintFun}, \code{comparisonFun} or \code{limitConstraints} is \code{NULL}, the contraint check will not be carried out. This is equivalent to simply finding all combinations of \eqn{v} choose \eqn{m}.
}
\examples{
## Quickly generate all 13-tuples of 1:25
## w/o repetition. Over 5 million combinations!
system.time(comboGeneral(25,13))

## Same as above, but this time we only want
## combinations such that the mean is less than 10
system.time(myComb <- comboGeneral(25,13,FALSE,"mean","<",10))

## Quickly generate all 9-tuples of 1:18
## w/ repetition. Over 13 million combinations!
system.time(comboGeneral(18,9,repetition = TRUE))

## Generate some random data
set.seed(1009)
mySamp <- rnorm(100, 997, 23)

## Using rowCap example:
## Researcher only needs 1000 7-tuples of mySamp
## such that the sum is greater than 7200.
system.time(comboGeneral(mySamp,7,FALSE,"sum",">",7200,1000))

## If you leave rowCap as NULL, it can take much longer
## (still fast enough in most cases) and in some case crach
## your computer as the underlying code allocates
## enough space to account for every combination
## (e.g. In our example, there are choose(100, 7)
## = 16007560800 rows, 7 columns, with each cell occupying
## 8 bytes. This gives a total over 830 GB 
## (i.e. choose(100, 7)*7*8/(2^30)).

\dontrun{
## You must generate all combinations and subsequently
## subset to obtain the combinations that meet the criteria
system.time(myComb2 <- combn(25,13))
ystem.time(myCols <- which(apply(myComb2, 2, mean) < 10))
system.time(myComb2 <- myComb2[,myCols])

## Test equality with generated object above
all.equal(myComb,t(myComb2))

## Any variation is much slower
system.time(myComb2 <- combn(25,13)[,combn(25,13,mean) < 10])
}
}