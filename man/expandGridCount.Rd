\name{expandGridCount}
\alias{expandGridCount}
\title{
    Count of the Cartesian Product
}
\description{
    Calculate the number of Cartesian products of the input vectors.
}
\usage{
    expandGridCount(...)
}
\arguments{
  \item{...}{vectors, factors or a list containing these. (See \code{?expand.grid}).}
}
\value{
When the number of results exceeds \eqn{2^{53} - 1}{2^53 - 1}, a number of class \code{bigz} is returned.
}
\author{
Joseph Wood
}
\seealso{
\code{\link{expandGrid}}
}
\examples{
## description example
lst = list(1:2, 1:2)

## Using base R
t = expand.grid(lst)
nrow(t)

## vs calling expandGridCount directly
expandGridCount(lst)

## Call it just like you would if you were generating the results
expandGridCount(1:5, 3:9, letters[1:5], letters[c(1,4,5,8)])

## Same as
nrow(expand.grid(1:5, 3:9, letters[1:5], letters[c(1,4,5,8)]))

lst = Map(function(x, y) x:y, 8:33, 15:40)

## Return object of class 'bigz'
expandGridCount(lst)
}
\keyword{expand}
\keyword{product}
\keyword{cartesian}
\keyword{expandGrid}
\keyword{expand.grid}
\keyword{combinatorics}
