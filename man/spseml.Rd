\name{spseml}
\alias{spseml}
\alias{similar.listw}
\alias{similar.listw_spam}
\alias{similar.listw_Matrix}
\alias{can.be.simmed}
\alias{llsurerror}
\alias{llsurlag}



\title{Spatial SUR - Lag and Error}

\description{Maximum likelihood estimation of symultaneous systems of 
spatial SUR model. 

}

\usage{
spseml(formula, data=list(),panel=TRUE,index=NULL,w,method="eigen",
 quiet=NULL, model = c("lag","error"), zero.policy=NULL, interval=NULL,
 tol.solve=1.0e-10, trs=NULL, control=list(), initval=NULL)}


\arguments{
  \item{formula}{a list of objects of class \code{formula} }
  \item{data}{an object of class \code{\link{data.frame}} or \code{pdata.frame}. 
  An optional data frame containing the variables in the model. 
  When the obect is a data.frame, the first two columns may contain the indexes. See index }
  \item{panel}{logical. When TRUE, the data frame is a panel data set with cross-sectional and time observations}
  \item{index}{if not NULL (default), a character vector to identify the indexes among the columns of the \code{\link{data.frame}}}
  \item{w}{an object of class \code{listw} created for example by \code{nb2listw} or a \code{matrix}}
  \item{method}{\code{spatialsim}}
  \item{quiet}{default \code{NULL}}
  \item{model}{\code{lag} or \code{error} If \code{lag}, a spatial lag of the dependent variable is added in each equation. 
  If \code{error} a spatially autocorrelated error term is assumed in each equation.}
  \item{zero.policy}{See \code{lagsarlm} for details}
  \item{interval}{search interval for spatial parameters}
  \item{tol.solve}{tolerance}
  \item{trs}{default NULL}
  \item{control}{A list of control parameters for the optimization}
  \item{initval}{the initial values of the spatial parameters}
}


\details{
The function can be specified with any number of equations. 
The number of equations is determined through the \code{formula} object.
The data can also be specified as a panel data frame. The logical argument
PANEL should then be set to TRUE.

}


\value{
  An object of class \code{"splm"}.
  \item{coefficients}{ML coefficients of the model parameters (for all equations)}
  \item{vcov}{the variance covariance matrix of the estimated coefficients}
  \item{type}{'spsegm'}
  \item{model.data}{the matrix of the data used (responses in each equation are reported first, then the 
  explanatory variables)}
  \item{model}{\code{lag} or \code{error}}
  \item{N}{the number of cross-sectional observations}  
  \item{Eq}{the number of equations in the system}  
  \item{K}{the total number of explanatory variables in the system}
  \item{call}{the call used to create the object}
  \item{terms}{the \code{terms} object generated from \code{formula} and some optional arguments}
  \item{Xnames}{the names of the variables in the matrix of explanatory variables}
  \item{Ynames}{the names of the responses}
  \item{spec}{the number of explanatory variables in each equation}
  \item{lags}{a logical list}
  \item{errors}{a logical list}
  \item{endogenous}{a logical list}
  \item{rho}{the estimated spatial autocorrelation coefficients}
  \item{method}{"GM"}
}

\references{
  Anselin, L. (1988) 
Spatial Econometrics: Methods and Models,
    \emph{Kluwer Academic Publisher}.
}

\author{Gianfranco Piras\email{gpiras@mac.com}}

\seealso{
\code{\link{spsegm}}
}
\examples{
data(Produc, package = "Ecdat")
data(usaww)
Produc <- Produc[Produc$year<1973, ]
eq1 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
eq2 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
eq3 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
formula<-list(tp1 = eq1, tp2 = eq2, tp3=eq3)
listw<-mat2listw(usaww)
sur.error<-spseml(formula, data = Produc, w = listw, model = "error",
 method = "eigen", quiet = TRUE)
summary(sur.error)
sur.lag<-spseml(formula, data = Produc, w = listw, model = "lag",
 method = "eigen", quiet = FALSE)
summary(sur.lag)
}

\keyword{spatial}
