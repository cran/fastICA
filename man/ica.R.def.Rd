\name{ica.R.def}
\alias{ica.R.def}
\title{R code for FastICA using a deflation scheme}
\description{R code for FastICA using a deflation scheme in which the
  components are estimated one by one. This function is called by the
  fastICA function.
}
\usage{
ica.R.def(X, n.comp, tol, fun, alpha, maxit, verbose, w.init)
}
\arguments{
  \item{X}{data matrix}
  \item{n.comp}{number of components to be extracted}
  \item{tol}{a positive scalar giving the tolerance at which the
    un-mixing matrix is considered to have converged.}
  \item{fun}{the functional form of the G function used in the
    approximation to negentropy (see details)}
  \item{alpha}{constant in range [1,2] used in approximation to
    negentropy when \code{fun == "logcosh"}}
  \item{maxit}{maximum number of iterations to perform}
  \item{verbose}{a logical value indicating the level of output as the
    algorithm runs.}
  \item{w.init}{Initial value of un-mixing matrix.}
}
\value{The estimated un-mixing matrix W (see details)}

\details{
  See help on \code{\link{fastICA}} for details.
}
\author{
  J L Marchini and C Heaton 
}
\seealso{\code{\link{fastICA}}, \code{\link{ica.R.par}}}

\keyword{utilities}
