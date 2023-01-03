\name{wnr}
\alias{wnr}
\title{Maximum likelihood estimators for the Weibull distribution with type I right censorship}
\description{
  Obtaint maximum likelihood estimators for the Weibull distribution.
}
\usage{
wnr(y, d, init, niter = 30, tol = 0.0001, hist = TRUE)
}
\arguments{
  \item{y}{Failure time.}
  \item{d}{Status. Usually 1 means complete observation and 0 right censoring.}
  \item{init}{Initial values for both shape and scale parameters.}
  \item{maxiter}{Maximum number of iterations for the Newton-Raphson algorithm, the default is 30.}
  \item{epsilon}{Stopping criteria, the routine finishes when the changes from one iteration to other is smaller a given treshold (default is 0.0001).}
  \item{hist}{Displays iteration history, otherwise only returns the result.}
}

\details{
See vignette for further details.
}

\value{
  A list with 9 components.
  \item{x0}{Initial values.}
  \item{ops}{Stoppin criteria (`niter` and `tol`)}
  \item{n}{Number of observations.}
  \item{ti}{Failure times.}
  \item{di}{Status.}
  \item{cod}{Convergence code: 1 Successful, otherwise there were some problems in optimization.}
  \item{xfin}{Final values.}
  \item{vcov}{Covariance matrix of estimators}
  \item{hist}{History code.}
}

\examples{

 library(wnr)
 data(device)

 #Obtainning the maximum likelihood estimators
 params <- c(1720, 1.3)
 wnr(device$t, device$d, params)

}

\keyword{methods}
