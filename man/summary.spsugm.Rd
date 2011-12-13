\name{summary.spsugm}
\alias{summary.spsugm}
\title{summary method for class spse}

\description{
Method for summarizing the results of objects of class \code{'spsugm'}
}

\usage{

\method{summary}{spsugm}(object,...)
}

\arguments{
\item{object}{an object of class \code{'spse'}}
\item{...}{additional arguments to be passed}
}


\details{
The summary function \code{summary.spsugm} 
returns an objects of class 'spsugm' 
organized in a coefficient matrix.

The summary 
is generated looping over the number
of equations in the system.
}

\value{
An object of class \code{'summary.spsugm'}

}

\seealso{
  \code{\link{spseml}}, \code{\link{spsegm}}, 
}
\author{Gianfranco Piras\email{gpiras@mac.com}}


\keyword{spatial}