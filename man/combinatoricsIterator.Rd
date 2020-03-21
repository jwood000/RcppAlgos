\name{comboIter}
\alias{comboIter}
\alias{permuteIter}
\title{
Combination and Permutation Iterator
}
\description{
    \itemize{
        \item Returns an iterator for iterating over combinations or permutations of a vector.
        \item Supports random access via the \code{[[} method.
        \item GMP support allows for exploration of combinations/permutations of vectors with many elements.
        \item The output is in lexicographical order for the \code{next} methods and reverse lexicographical order for the \code{prev} methods.
    }
}
\usage{
comboIter(v, m = NULL, repetition = FALSE, freqs = NULL,
          constraintFun = NULL, comparisonFun = NULL,
          limitConstraints = NULL, FUN = NULL, tolerance = NULL)

permuteIter(v, m = NULL, repetition = FALSE, freqs = NULL,
            constraintFun = NULL, comparisonFun = NULL,
            limitConstraints = NULL, FUN = NULL, tolerance = NULL)
}
\arguments{
  \item{v}{Source vector. If \code{v} is an integer (including nonpositive integers), it will be converted to the sequence \code{1:v}. All atomic types are supported (See \code{?is.atomic}).}
  \item{m}{Number of elements to choose. If \code{repetition = TRUE} or \code{freqs} is utilized, \code{m} can exceed the length of \code{v}. If \code{m = NULL}, the length will default to \code{length(v)} or \code{sum(freqs)}.}
  \item{repetition}{Logical value indicating whether combinations/permutations should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all combinations/permutations of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
  \item{constraintFun}{This feature will be available in future releases. The default is \code{NULL}.}
  \item{comparisonFun}{This feature will be available in future releases. The default is \code{NULL}.}
  \item{limitConstraints}{This feature will be available in future releases. The default is \code{NULL}.}
  \item{FUN}{Function to be applied to each combination/permutation. The default is \code{NULL}.}
  \item{tolerance}{This feature will be available in future releases. The default is \code{NULL}.}
}

\value{
A vector of length \emph{m} is returned if \code{nextIter} or \code{prevIter} is called, otherwise a matrix is returned with each row containing a vector of length \eqn{m}. If \code{FUN} is utilized and either \code{nextIter} or \code{prevIter} is called, the result will be determined by \code{FUN}, otherwise a list is returned.
}
\details{
Once you initialize a new iterator, the following methods are available:
  \describe{
    \item{\code{nextIter}}{Retrieve the \bold{next} lexicographical result}
    \item{\code{nextNIter}}{Pass an integer \emph{n} to retrieve the \bold{next} \emph{n} lexicographical results}
    \item{\code{nextRemaining}}{Retrieve all remaining lexicographical results}
    \item{\code{currIter}}{Returns the current iteration}
    \item{\code{prevIter}}{Retrieve the \bold{previous} lexicographical result (the \bold{next} \emph{reverse} lexicographical result)}
    \item{\code{prevNIter}}{Pass an integer \emph{n} to retrieve the \bold{previous} \emph{n} lexicographical results (the \bold{next} \emph{n} \emph{reverse} lexicographical results)}
    \item{\code{prevRemaining}}{Retrieve all remaining \emph{reverse} lexicographical results}
    \item{\code{startOver}}{Resets the iterator}
    \item{\code{sourceVector}}{View the source vector}
    \item{\code{summary}}{Returns a list of summary information about the iterator}
    \item{\code{front}}{Retrieve the \bold{first} lexicographical result}
    \item{\code{back}}{Retrieve the \bold{last} lexicographical result}
    \item{\code{[[}}{Random access method. Pass a single value or a vector of valid indices. If a single value is passed, the internal index of the iterator will be updated, however if a vector is passed the internal state will not change. GMP support allows for flexible indexing.}
  }
}
\note{
    \itemize{
        \item The maximum number of combinations/permutations that can be generated at one time is \eqn{2^{31} - 1}{2^31 - 1}.
        \item Factor vectors are accepted. Class and level attributes are preserved.
        \item Currently, constrained iterators are not supported (Will be available in future releases).
    }
}
\seealso{
\code{\link{comboGeneral}}, \code{\link{permuteGeneral}}
}
\references{
    \itemize{
        \item{\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical Order}}
        \item{\href{https://oeis.org/wiki/Orderings#Reverse_lexicographic_order}{Reverse Lexicographical Order}}
    }
}
\author{
Joseph Wood
}
\examples{
a = permuteIter(unique(state.region))
a$nextIter()
a$nextNIter(3)
a$front()
a$nextRemaining()
a$prevIter()
a$prevNIter(15)
a$summary()
a$back()
a$prevRemaining()
a[[5]]
a$summary()
a[[c(1, 17, 3)]]
a$summary()
}