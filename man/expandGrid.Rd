\name{expandGrid}
\alias{expandGrid}
\title{
    Cartesian Product
}
\description{
    Generate the \href{https://en.wikipedia.org/wiki/Cartesian_product}{Cartesian Product} of the input vectors. It is very similar to \code{\link{expand.grid}} however there are some notable differences:
    \itemize{
        \item{Produces lexicographic ordered output consistent with other functions in \code{RcppAlgos}. Compared to \code{expand.grid} where the first column varies the fastest, \code{expandGrid} varies the first column the slowest.}
        \item{When all of the input is of the same type, by default \code{expandGrid} produce a \code{matrix} (a \code{data.frame} otherwise). This can be ignored by setting the argument \code{return_df = TRUE}.}
        \item{No attributes are added nor are strings converted to factors. In \code{expand.grid} we would achieve this by setting \code{KEEP.OUT.ATTRS = FALSE} and \code{stringsAsFactors = FALSE}.}
        \item{If it is possible to return a matrix, we can utilize the argument \code{nThreads} in order to produce results in parallel for maximal efficiency.}
      }
}
\usage{
    expandGrid(..., lower = NULL, upper = NULL,
               nThreads = NULL, return_df = FALSE)
}
\arguments{
  \item{...}{vectors, factors or a list containing these. (See \code{?expand.grid}).}
  \item{lower}{The lower bound. Cartesian products are generated lexicographically, thus utilizing this argument will determine which specific product to start generating from (\emph{e.g.} \code{expandGrid(1:5, 3:11, lower = 6)} is equivalent to \code{expandGrid(1:5, 3:11)[6:expandGridCount(1:5, 3:11), ]}). This argument along with \code{upper} is very useful for generating products in chunks allowing for easy parallelization.}
  \item{upper}{The upper bound. Similar to \code{lower}, however this parameter allows the user to \emph{stop} generation at a specific product (\emph{e.g.} \code{expandGrid(1:5, 3:11, upper = 5)} is equivalent to \code{expandGrid(1:5, 3:11)[1:5, ]})}
  \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
  \item{return_df}{Logical flag to force the output to be a \code{data.frame}. The default is \code{FALSE}.}
}
\value{
When all of the input is of the same type, by default \code{expandGrid} produces a \code{matrix} (a \code{data.frame} otherwise). This can be ignored by setting the argument \code{return_df = TRUE}.
}
\author{
Joseph Wood
}
\seealso{
\code{\link{comboGrid}}
}
\examples{
## description example
lst = list(1:2, 1:2)

## Using base R
t = expand.grid(lst)

## vs using expandGrid. N.B. Output is a matrix
expandGrid(lst)

## Force a data.frame to be returned
expandGrid(lst, return_df = TRUE)

lst = Map(function(x, y) x:y, 8:14, 15:21)

## Use multiple threads for greater efficiency
system.time(expandGrid(lst, nThreads = 2))
}
\keyword{expand}
\keyword{product}
\keyword{cartesian}
\keyword{expandGrid}
\keyword{expand.grid}
\keyword{combinatorics}
