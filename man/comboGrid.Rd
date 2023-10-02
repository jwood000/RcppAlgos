\name{comboGrid}
\alias{comboGrid}
\title{
Efficient Version of \code{expand.grid} Where order Does Not Matter
}
\description{
    This function efficiently generates Cartesian-product-like output where order does not matter. It is loosely equivalent to the following:
    \itemize{
        \item{\code{t = expand.grid(list)}}
        \item{\code{t = t[do.call(order, t), ]}}
        \item{\code{key = apply(t, 1, function(x) paste0(sort(x), collapse = ""))}}
        \item{\code{t[!duplicated(key), ]}}
      }
}
\usage{
comboGrid(..., repetition = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{vectors, factors or a list containing these. (See \code{?expand.grid}).}
  \item{repetition}{Logical value indicating whether results should be with or without repetition. The default is \code{TRUE}.}
}
\value{
If items with different classes are passed, a data frame will be returned, otherwise a matrix will be returned.
}
\author{
Joseph Wood
}
\examples{
## return a matrix
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
\keyword{cartesian}
\keyword{expand.grid}
\keyword{combinatorics}
