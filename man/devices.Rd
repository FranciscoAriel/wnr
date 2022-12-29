\name{device}
\alias{device}
\title{Temperature-Accelerated Life Test Data for Device-A}
\description{
	The dataset contains information about temperature-accelerated life test data on a certain device.
  Device samples were tested for up to 5000 hours at four separate temperatures.

  This dataset only contains information where temperature is 80 degrees celsius.
}

\usage{
  data(device)
}
\format{
  The data.frame `device' has two columns: 1)t and 2)d.
  t: Failure time.
  d: Status (1) failed, (0) censored.
}

\references{
  Hoopers, J. H. and Amster, S. J. (1990), Analysis and presentation of reliability data, in Handbook of Statistical Methods for Engineers and Scientists, McGraw-Hill, New York. Harrison M. Wadsworth, Editor.

  Meeker, W. Q. and Escobar, L. A. (1998), Statistical Methods for Reliability Data, New York, NY; John Wiley & Sons.
}

\keyword{datasets}

