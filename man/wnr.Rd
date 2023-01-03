\name{wnr}
\alias{wnr}
\title{Maximum likelihood estimators for the Weibull distribution with type I right censorship}
\description{
  Obtaint maximum likelihood estimators for the Weibull distribution.
}
\usage{
wnr(y, d, init = c(mean(y), 1), niter = 30, tol = 0.0001, hist = TRUE)
}
\arguments{
  \item{y}{Failure time.}
  \item{d}{Status. Usually 1 means complete observation and 0 right censoring.}
  \item{init}{Initial values for both shape and scale parameters.}
  \item{niter}{Maximum number of iterations for the Newton-Raphson algorithm, the default is 30.}
  \item{tol}{Stopping criteria, the routine finishes when the changes from one iteration to other is smaller a given treshold (default is 0.0001).}
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
 time = c(283,361,515,638,854,1024,1030,1045,1767,1777,1856,1951,1964,2884,5000)
 status = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,0)
 # Default values
 wnr(time,status)
 # Fast convergence
 # Obtainning the maximum likelihood estimators
 params <- c(1720, 1.3)
 wnr(time,status,init = params)
}

\keyword{methods}
