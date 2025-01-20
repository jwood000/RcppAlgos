\name{comboGrid}
\alias{comboGrid}
\title{
    Unordered Cartesian Product
}
\description{
    Efficient version of \code{expand.grid} where order does not matter. This is a combinatorial variant where groups of elements are treated as equivalent regardless of order. For example, given: \code{{1, 2}, {1, 2}}, the unordered Cartesian product is \code{{1, 1}, {1, 2}, {2, 2}}. It is loosely equivalent to the following:
    \itemize{
        \item{\code{t = expand.grid(lst)}}
        \item{\code{t = t[do.call(order, t), ]}}
        \item{\code{key = apply(t, 1, function(x) paste0(sort(x), collapse = ""))}}
        \item{\code{t[!duplicated(key), ]}}
      }
}
\usage{
    comboGrid(..., repetition = TRUE, return_df = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{vectors, factors or a list containing these. (See \code{?expand.grid}).}
  \item{repetition}{Logical value indicating whether results should be with or without repetition. The default is \code{TRUE}.}
  \item{return_df}{Logical flag to force the output to be a \code{data.frame}. The default is \code{FALSE}.}
}
\value{
When all of the input is of the same type, by default \code{comboGrid} produce a \code{matrix} (a \code{data.frame} otherwise). This can be ignored by setting the argument \code{return_df = TRUE}.
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

t = expand.grid(lst)
t = t[do.call(order, t), ]
key = apply(t, 1, function(x) paste0(sort(x), collapse = ""))
t[!duplicated(key), ]

## vs using comboGrid. N.B. Output is a matrix
comboGrid(lst)

## Force a data.frame to be returned
comboGrid(lst, return_df = TRUE)

## Input vectors are of different type, so a data.frame is returned
expGridNoOrder = comboGrid(1:5, 3:9, letters[1:5], letters[c(1,4,5,8)])
head(expGridNoOrder)
tail(expGridNoOrder)

expGridNoOrderNoRep = comboGrid(1:5, 3:9, letters[1:5],
                                letters[c(1,4,5,8)], repetition = FALSE)

head(expGridNoOrderNoRep)
tail(expGridNoOrderNoRep)
}
\keyword{expand}
\keyword{comboGrid}
\keyword{expandGrid}
\keyword{cartesian}
\keyword{expand.grid}
\keyword{combinatorics}
