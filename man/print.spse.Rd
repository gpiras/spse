\name{print.spse}
\alias{print.spse}
\title{print method for class splm}

\description{
Method to print objects of class \code{summary.spse} and \code{spse}
}

\usage{
\method{print}{spse}(x, digits = max(3,getOption("digits") -3), ...)
}

\arguments{
\item{x}{an object of class \code{spse}}
\item{digits}{minimal number of significant digits, see \code{\link{print.default}}}
\item{...}{additional arguments to be passed}
}


\details{
The summary function \code{summary.spse} 
returns an objects of class 'spse' 
organized in a coefficient matrix.

The \code{summary} 
is generated looping over the number
of equations in the system.
}

\seealso{
\code{\link{spsegm}} 
}
\author{ Gianfranco Piras\email{gpiras@mac.com}}

\examples{
data(Produc, package = "Ecdat")
data(usaww)
Produc <- Produc[Produc$year<1973, ]
eq1 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
eq2 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
eq3 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
formula<-list(tp1 = eq1, tp2 = eq2, tp3=eq3)
w<-mat2listw(usaww)
se<-spsegm(formula, data=Produc, w=w, panel= TRUE,lags=list(c(TRUE,TRUE,TRUE),c(TRUE,TRUE,TRUE),c(TRUE,TRUE,TRUE)), errors=list(FALSE,TRUE,FALSE),endogenous=list(c(FALSE,TRUE,FALSE),c(TRUE,FALSE,FALSE),c(TRUE,FALSE,FALSE)))
summary(se)
}
\keyword{spatial}