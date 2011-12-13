\name{spsepgm}
\alias{spsepgm}
\title{Spatial Simultaneous Equations Panel}

\description{Feasible generalized three stages least square estimator (FGS3SLS)
of symultaneous systems of spatially interrelated cross sectional equations of the form: 

\deqn{Y = Y B + X C +  \bar{Y} L + U}

\deqn{ U = \bar{U} R + E}

with \eqn{Y = (y_1,\dots,y_m)}, 
\eqn{X = (x_1,\dots,x_m)}, \eqn{U = (u_1,\dots,u_m)}, 
\eqn{\bar{Y} = (\bar{y}_1,\dots,\bar{y}_m)} 
and \eqn{\bar{y}_j = W y_j} \eqn{j=1,\dots,m}, 
\eqn{E = (e_1,\dots,e_m)}, \eqn{\bar{U} = (\bar{u}_1,\dots,\bar{u}_m)} 
and \eqn{\bar{u}_j = W u_j} \eqn{j=1,\dots,m}. 
B, C, L and \eqn{R=diag_{j=1}^m (\rho_j)} 
are matrix of parameters


}

\usage{
spsepgm(formula, data = list(), index = NULL, listw, model = c("within","random"), lag = NULL, 
spatial.error = NULL, moments = c("initial", "weights", "fullweights"), endog = NULL, 
instruments = NULL, verbose = FALSE, method = c("w2sls", "b2sls", "g2sls", "ec2sls"), control = list())
}
\arguments{
  \item{formula}{a list of objects of class \code{formula} }
  \item{data}{an object of class \code{\link{data.frame}} or \code{pdata.frame}.An optional data frame containing the variables in the model. An optional data frame containing the variables in the model. When the obect is a data.frame, the first two columns may contain the indexes. See index}
  \item{index}{if not NULL (default), a character vector to identify the indexes among the columns of the \code{\link{data.frame}}}
  \item{listw}{an object of class \code{listw} created for example by \code{nb2listw} or a \code{matrix}}
  \item{model}{\code{"within"} or \code{"random"}}
  \item{lag}{A logical list of length equal to the number of equations. If TRUE the spatial lag of the variable is included in the equation}
  \item{spatial.error}{A logical vector. When TRUE a spatially autocorrelated error term is included in the corresponding equation}
  \item{moments}{\code{c("initial", "weights", "fullweights")}}
  \item{endog}{Additional endogenous variables}
  \item{instruments}{Additional instruments}
  \item{verbose}{report the output of the optimization results}
  \item{method}{one of c("w2sls", "b2sls", "g2sls", "ec2sls")}
  \item{control}{a list of controls}
}


\details{
The function can be specified with any number of equations. 
The number of equations is determined through the \code{formula} object.
The data can also be specified as a panel data frame. The logical argument
PANEL should then be set to TRUE.

The logical list \code{lags} controls which spatial lags should be included in the equations. 
The logical list \code{errors} determines which equations should include an autoregressive term. 
The logical list \code{endogenous} determines which equations should include an autoregressive term. 

}


\value{
  An object of class \code{"splm"}.
  \item{coefficients}{FG3SLS coefficients estimate of the model parameters (for all equations)}
  \item{vcov}{the variance covariance matrix of the estimated coefficients}
  \item{residuals}{FG3SLS residuals}
  \item{fitted.values}{the fitted values}  
  \item{sigma2}{error variance}  
  \item{type}{"Seemingly unrelated regressions with spatial error components"}
  \item{rho}{a list with the spatial error coefficient and the variance component}
  \item{call}{the call used to create the object}
  \item{model.data}{the data from teh third step}
  \item{expv}{the number of explanatory variables in each equation}
  \item{namesx}{the names of the responses in all equation}
  \item{NT}{the number of observation}
  \item{N}{the number of cross-sections}
   \item{N}{the number of time periods}
  \item{eq}{the number of equations}
}

\references{
  Kelejian, H.H. and Prucha, I.R. (2004) 
Estimation of Simultaneous systems of spatially interrelated cross sectional equations,
    \emph{Journal of Econometrics}, \bold{118}, pages 27--50.

  Kelejian, H.H. and Prucha, I.R. (1999) 
A Generalized Moments Estimator for the Autoregressive Parameter in a Spatial Model,
    \emph{International Economic Review}, \bold{40}, pages 509--533.
    
      Kelejian, H.H. and Prucha, I.R. (1998) 
A Generalized Spatial Two Stage Least Square Procedure for Estimating a Spatial Autoregressive
Model with Autoregressive Disturbances,
    \emph{Journal of Real Estate Finance and Economics}, \bold{17}, pages 99--121.

}

\author{Gianfranco Piras\email{gpiras@mac.com}}

\seealso{
\code{\link{spseml}}
}
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