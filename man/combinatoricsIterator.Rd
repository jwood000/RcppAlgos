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
          limitConstraints = NULL, keepResults = NULL,
          FUN = NULL, Parallel = FALSE, nThreads = NULL,
          tolerance = NULL, FUN.VALUE = NULL)

permuteIter(v, m = NULL, repetition = FALSE, freqs = NULL,
            constraintFun = NULL, comparisonFun = NULL,
            limitConstraints = NULL, keepResults = NULL,
            FUN = NULL, Parallel = FALSE, nThreads = NULL,
            tolerance = NULL, FUN.VALUE = NULL)
}
\arguments{
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. All atomic types are supported (See \code{?is.atomic}).}
  \item{m}{Number of elements to choose.  If \code{repetition = TRUE} or \code{freqs} is utilized, \code{m} can exceed the length of \code{v}. If \code{m = NULL}, the length will default to \code{length(v)} or \code{sum(freqs)}.}
  \item{repetition}{Logical value indicating whether combinations/permutations should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all combinations/permutations of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
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
a@nextIter()
a@nextNIter(3)
a@front()
a@nextRemaining()
a@prevIter()
a@prevNIter(15)
a@summary()
a@back()
a@prevRemaining()
a[[5]]
a@summary()
a[[c(1, 17, 3)]]
a@summary()
}