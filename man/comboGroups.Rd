\name{comboGroups}
\alias{comboGroups}
\title{
Partition a Vector into Groups of Equal Size
}
\description{
    \itemize{
        \item Generate partitions of a vector into groups of equal size. See \href{https://stackoverflow.com/q/57732672/4408538}{Create Combinations in R by Groups} on http://stackoverflow.com for a direct use case.
        \item Produce results in parallel using the \code{Parallel} or \code{nThreads} arguments.
        \item GMP support allows for exploration where the number of results is large.
        \item The output is in lexicographical order by groups.
    }
}
\usage{
comboGroups(v, numGroups, retType = "matrix", 
            lower = NULL, upper = NULL, Parallel = FALSE,
            nThreads = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{Source vector. If \code{v} is an integer (including nonpositive integers), it will be converted to the sequence \code{1:v}. All atomic types are supported (See \code{?is.atomic}).}
  \item{numGroups}{An Integer. The number of groups that the vector will be partitioned into. Must divide the length of v (if v is a vector) or v (if v is a scalar).
}
  \item{retType}{A string, "3Darray" or "matrix", that determines the shape of the output. The default is "matrix".
}
  \item{lower}{The lower bound. Partitions of groups are generated lexicographically, thus utilizing this argument will determine which specific result to start generating from (\emph{e.g.} \code{comboGroups(8, 2, lower = 30)} is equivalent to \code{comboGroups(8, 2)[30:comboGroupsCount(8, 2), ]}). This argument along with \code{upper} is very useful for generating results in chunks allowing for easy parallelization.}
  \item{upper}{The upper bound. Similar to \code{lower}, however this parameter allows the user to \emph{stop} generation at a specific result (\emph{e.g.} \code{comboGroups(8, 2, upper = 5)} is equivalent to \code{comboGroups(8, 2)[1:5, ]})}
  \item{Parallel}{Logical value indicating whether results should be generated in parallel using \eqn{n - 1} threads, where \eqn{n} is the maximum number of threads. The default is \code{FALSE}. If \code{nThreads} is not \code{NULL}, it will be given preference (\emph{e.g.} if user has 8 threads with \code{Parallel = TRUE} and \code{nThreads = 4}, only 4 threads will be spawned). If your system is single-threaded, the arguments \code{Parallel} and \code{nThreads} are ignored.}
  \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}. See \code{Parallel}.}
}
\details{
Conceptually, this problem can be viewed as generating all permutations of the vector \code{v} and removing the within group permutations. To illustrate this, let us consider the case of generating partitions of \code{1:8} into 2 groups of size 4.
  \itemize{
    \item{To begin, generate the permutations of \code{1:8} and group the first/last four elements of each row.
      \tabular{ccccccccc}{
      \tab \tab \tab Grp1 \tab \tab \tab Grp2 \tab \tab \cr
      \tab C1\tab C2\tab C3\tab C4\tab C5\tab C6\tab C7\tab C8\cr
      R1 \tab  |1 \tab 2 \tab 3 \tab 4| \tab |5 \tab 6 \tab 7 \tab 8|\cr
      R2 \tab  |1 \tab 2 \tab 3 \tab 4| \tab |5 \tab 6 \tab 8 \tab 7|\cr
      R3 \tab  |1 \tab 2 \tab 3 \tab 4| \tab |5 \tab 7 \tab 6 \tab 8|\cr
      R4 \tab  |1 \tab 2 \tab 3 \tab 4| \tab |5 \tab 7 \tab 8 \tab 6|\cr
      R5 \tab  |1 \tab 2 \tab 3 \tab 4| \tab |5 \tab 8 \tab 6 \tab 7|\cr
      R6 \tab  |1 \tab 2 \tab 3 \tab 4| \tab |5 \tab 8 \tab 7 \tab 6|
      }
    }
    \item{Note that the permutations above are equivalent partitions of 2 groups of size 4 as only the last four elements are permuted. If we look at at the \eqn{25^{th}}{25th} lexicographical permutation, we observe our second distinct partition.
      \tabular{ccccccccc}{
      \tab \tab \tab Grp1 \tab \tab \tab Grp2 \tab \tab \cr
      \tab C1\tab C2\tab C3\tab C4\tab C5\tab C6\tab C7\tab C8\cr
      R24 \tab  |1 \tab 2 \tab 3 \tab 4| \tab |8 \tab 7 \tab 6 \tab 5|\cr
      \strong{R25} \tab  \strong{|1} \tab \strong{2} \tab \strong{3} \tab \strong{5|} \tab \strong{|4} \tab \strong{6} \tab \strong{7} \tab \strong{8|}\cr
      R26 \tab  |1 \tab 2 \tab 3 \tab 5| \tab |4 \tab 6 \tab 8 \tab 7|\cr
      R27 \tab  |1 \tab 2 \tab 3 \tab 5| \tab |4 \tab 7 \tab 6 \tab 8|\cr
      R28 \tab  |1 \tab 2 \tab 3 \tab 5| \tab |4 \tab 7 \tab 8 \tab 6|
      }
    }
    \item{Continuing on, we will reach the \eqn{3,457^{th}}{3,457th} lexicographical permutation, which represents the last result:
      \tabular{ccccccccc}{
      \tab \tab \tab Grp1 \tab \tab \tab Grp2 \tab \tab \cr
      \tab C1\tab C2\tab C3\tab C4\tab C5\tab C6\tab C7\tab C8\cr
      R3454 \tab  |1 \tab 6 \tab 7 \tab 5| \tab |8 \tab 3 \tab 4 \tab 2|\cr
      R3455 \tab  |1 \tab 6 \tab 7 \tab 5| \tab |8 \tab 4 \tab 2 \tab 3|\cr
      R3456 \tab  |1 \tab 6 \tab 7 \tab 5| \tab |8 \tab 4 \tab 3 \tab 2|\cr
      \strong{R3457} \tab  \strong{|1} \tab \strong{6} \tab \strong{7} \tab \strong{8|} \tab \strong{|2} \tab \strong{3} \tab \strong{4} \tab \strong{5|}\cr
      R3458 \tab  |1 \tab 6 \tab 7 \tab 8| \tab |2 \tab 3 \tab 5 \tab 4|
      }
    }
    \item{
    For this small example, the method above will not be that computationally expensive. In fact, there are only 35 total partitions of \code{1:8} into 2 groups of size 4 out of a possible \code{factorial(8) = 40320} permutations. However, just doubling the size of the vector will make this approach infeasible as there are over 10 trillion permutations of \code{1:16}.
    }
    \item{The algorithm in \code{comboGroups} avoids these duplicate partitions of groups by utlizing an efficient algorithm analogous to the \href{https://en.cppreference.com/w/cpp/algorithm/next_permutation}{std::next_permutation} found in the standard algorithm library in C++.
    }
  }
}
\value{
By default, a matrix is returned with column names corresponding to the associated group. If \code{retType = "3Darray"}, a 3D array is returned.
}
\author{
Joseph Wood
}
\note{
The maximum number of partitions of groups that can be generated at one time is \eqn{2^{31} - 1}{2^31 - 1}. Utilizing \code{lower} and \code{upper} makes it possible to generate additional combinations/permutations.
}
\examples{
## return a matrix
comboGroups(8, 2)

## or a 3 dimensional array
temp = comboGroups(8, 2, "3Darray")

## view the first partition
temp[1, , ]
}
