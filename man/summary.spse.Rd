\name{summary.spse}
\alias{summary.spse}
\title{summary method for class spse}

\description{
Method for summarizing the results of objects of class \code{'spse'}
}

\usage{

\method{summary}{spse}(object,...)
}

\arguments{
\item{object}{an object of class \code{'spse'}}
\item{...}{additional arguments to be passed}
}


\details{
The summary function \code{summary.spse} 
returns an objects of class 'spse' 
organized in a coefficient matrix.

The summary 
is generated looping over the number
of equations in the system.
}

\value{
An object of class \code{'summary.spse'}

}

\seealso{
  \code{\link{spseml}}, \code{\link{spsegm}}, 
}
\author{Gianfranco Piras\email{gpiras@mac.com}}

\examples{
data(Produc, package = "Ecdat")
data(usaww)
Produc <- Produc[Produc$year<1973, ]
eq1 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
eq2 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
eq3 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
formula<-list(tp1 = eq1, tp2 = eq2, tp3=eq3)
w<-mat2listw(usaww)
se<-spsegm(formula, data=Produc, w=w, panel= TRUE,
 lags=list(c(TRUE,TRUE,TRUE),c(TRUE,TRUE,TRUE),c(TRUE,TRUE,TRUE)),
 errors=list(FALSE,TRUE,FALSE),
 endogenous=list(c(FALSE,TRUE,FALSE),c(TRUE,FALSE,FALSE),c(TRUE,FALSE,FALSE)))
summary(se)

}

\keyword{spatial}
