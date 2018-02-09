spsepgm<-function(formula, data = list(), index = NULL, listw, model = c("within","random"), lag = NULL, spatial.error = NULL, moments = c("initial", "weights", "fullweights"), endog = NULL, instruments = NULL, verbose = FALSE, method = c("w2sls", "b2sls", "g2sls", "ec2sls"), control = list()) {

#source("spregm.R")

    effects<-match.arg(model)
    moments<-match.arg(moments)

iindex<-index

    if (!is.null(index)) {
        #require(plm)
        data <- plm::plm.data(data, index)
    }
    
    index <- data[, 1]
    tindex <- data[, 2]
    cl <- match.call()
    if (dim(data)[[1]] != length(index))  stop("Non conformable arguments")
    
    
    mt <- terms(formula[[1]], data = data)
    mf <- lm(formula[[1]], data, na.action = na.fail, method = "model.frame")

    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)


    names(index) <- row.names(data)
    ind <- index[which(names(index) %in% row.names(x))]
    tind <- tindex[which(names(index) %in% row.names(x))]
    N <- length(unique(ind))
    T <- max(tapply(x[, 1], ind, length))

    NT <- length(ind)

    balanced <- N * T == NT
    if (!balanced) 
        stop("Estimation method unavailable for unbalanced panels")


index<-iindex
  	
    if (!inherits(formula, "list")) 
        stop("formula should be a list in simultaneous equation models")

        eq <- length(formula)

    if (!inherits(listw, "list")) {

    if (!inherits(listw, c("listw", "matrix"))) 
        stop("listw has to be of class 'listw' or 'matrix' if it is not a list")

else	{
	w.list <- vector("list", length = eq )
for ( i in 1:eq)    		w.list[[i]]<-listw
	
	}
    	
    	}
    	else{
    		w.list<-listw
    		}


    n <- dim(data)[[1]]

			
        if(!is.null(lag) & length(lag) != eq) stop("Argument lag incorrectly specified")
        if(is.null(lag)) lag<-rep(FALSE, eq)
        if(!is.null(endog) && is.null(instruments)) stop("Instruments not specified")
        
        
        est.res <- vector("list", length = eq)
        bigk<-0
        expv<-vector("numeric", length=0)
        rhos<-vector("numeric", length=0)
        namesx<-vector("list", length=0)

        for (i in 1:eq) {

est.res[[i]]<-splm::spgm(formula=formula[[i]], data=data, index = index, listw = w.list[[i]], 
model = model, lag = lag[[i]], spatial.error = spatial.error[[i]],  moments = moments, endog = endog, instruments = instruments, 
verbose = verbose, method = method, control = list())
#print(summary(est.res[[i]]))

bigk<-bigk + est.res[[i]]$rhs
expv<-c(expv,length(est.res[[i]]$coefficients))
namesx[[i]]<-names(est.res[[i]]$coefficients)
rhos<-c(rhos,est.res[[i]]$rho[1])


}


Umat <- matrix(,NT,eq)
bigy<-matrix(0,NT*eq,1)
bigx<-matrix(0,NT*eq,bigk)
model.data<-data.frame(matrix(0,NT,0))
counter<-0


        for (i in 1:eq) {
Umat[,i]<-est.res[[i]]$residuals
upl<-NT*i
lowl<-upl-NT+1
counter<-counter+est.res[[i]]$rhs
bigy[lowl:upl,]<-est.res[[i]]$coy
bigx[lowl:upl,((counter-est.res[[i]]$rhs + 1):counter)]<-est.res[[i]]$cox
model.data<-data.frame(cbind(model.data,est.res[[i]]$model))
}
    	
#pos<-c(1,10,13,22)
#print(bigx[1:10,])


ind <- rep(seq(1, N),T)
avctp <- function(x) rep(unlist(tapply(x, ind, mean, simplify = TRUE)), T)
Q1Umat <- apply(Umat, 2, avctp)
Q2Umat <- Umat - Q1Umat


Omv<-crossprod(Umat,Q2Umat)/(N*(T-1))
Om1<-crossprod(Umat,Q1Umat)/N

#print(Omv)
#print(Om1)

Omvinv<-invfrac(Omv, -0.5)
Om1inv<-invfrac(Om1, -0.5)

#print(Omvinv)
#print(Om1inv)



JT<-matrix(1/T,T,T)
IN<-Diagonal(N)
Q1<-kronecker(JT, IN)
IT<-Diagonal(T)
FP<-matrix((1-(1/T)),T,T)
Q2<-kronecker(FP, IN)

Omuinv<-kronecker(Om1inv, Q1) + kronecker(Omvinv, Q2)

#print(Omuinv)
finy<-Omuinv %*% bigy
finx<-Omuinv %*% bigx
finx<-as.matrix(finx)
#print(finx[1:10,])

#print(finx)

#finxpfinx<-crossprod(as.matrix(finx))
#print(cor(as.matrix(finxpfinx)))

#finxpfinxinv<-solve(as.matrix(finxpfinx), tol = 2e-20)
#pos<-c(1,10,13,22)
#print(finx[,pos])
modfin<-lm(as.matrix(finy)~as.matrix(finx)-1)

#betaFGLS<-finxpfinxinv %*% crossprod(finx,finy)
betaFGLS<-coefficients(modfin)

#print(betaFGLS)

#    fv <- as.vector(finx %*% betaFGLS)
fv<-fitted.values(modfin)
    egls <- finy - fv
    SGLS <- sum(egls^2)/(N - 1)
#    xfpxfNT <- (1/NT) * crossprod(finx)
 #   PSI <- solve(xfpxfNT, tol=2e-20)
    covbeta <- vcov(modfin)

type<- "Seemingly unrelated regressions with spatial error components"


#print(covbeta)

spmod<- list(coefficients= betaFGLS, vcov = covbeta, residuals = as.vector(egls), fitted.values = fv, 
sigma2 = SGLS, type = type, rho= rhos, model.data = model.data, call = cl, expv = expv, namesx = namesx, NT = NT, N = N, T = T, eq = eq)

class(spmod)<-"spsugm"

return(spmod)

    	}

summary.spsugm<-function(object,...){

	    coeff<-object$coefficients	    
	    eq <- object$eq
       expv <- object$expv	
       var<-object$vcov

       residuals<-object$residuals
       serr<-sqrt(diag(as.matrix(var)))
       tval<-coeff/serr
       pval<-pnorm(abs(as.matrix(tval)), lower.tail = FALSE) * 2
       
		 Xnam<-object$namesx
		 NT<-object$NT
		 N<-object$N
		 T<-object$T
		 NT<-object$NT
		 rho<-object$rho
		 tmp<-1:eq
		 tmp2<-rep(tmp,expv)

		 beta<-split(coeff,tmp2)
#print(coeff)		 
#print(beta)
		 se<-split(serr,tmp2)
		 tv<-split(tval,tmp2)
 		 pv<-split(pval,tmp2)


 		 namesx<-split(unlist(object$namesx),tmp2)

 		 tmp3<-rep(tmp,each=NT)
 		 resid<-split(residuals, tmp3) 		 
 		 object$beta<-beta
 		 object$se<-se
 		 object$tv<-tv 		 
 		 object$pv<-pv
 		 object$eq<-eq

		 object$namesx<-namesx
 		 object$Xnam<-Xnam
 		 object$resid<-resid
class(object)<- c("summary.spsugm","spsugm")

object
	}


print.summary.spsugm<-function(x, digits = max(3, getOption("digits") - 2), width = getOption("width"), 
    ...){
        cat("\nSeemingly unrelated regressions with spatial error components\n")
        cat("\nCall:\n")

        print(x$call)
        eq <- x$eq
        save.digits <- unlist(options(digits = digits))
        on.exit(options(digits = save.digits))
        namesx <- x$namesx


        for (i in 1:eq) {
            cat(" \n")
            cat(paste("Equation", i, sep = " ", collapse = ""), "\n", sep = "")
            
        cat("\nResiduals:\n")
        save.digits <- unlist(options(digits=digits))
        on.exit(options(digits=save.digits))

  sr <- summary(x$resid[[i]])
  srm <- mean(x$resid[[i]])
  if (abs(srm) < 1e-10) sr <- sr[c(1:3,5:6)]

        print(sr)

tables <- cbind(as.matrix(x$beta[[i]]), as.matrix(x$se[[i]]), as.matrix(x$tv[[i]]), as.matrix(x$pv[[i]]))


            dimnames(tables) <- list(namesx[[i]], c("Estimate", "Std.Error", "t value", "Pr(>|t|)"))
#                rownames(tables)


            if (i == eq) {
                legend = TRUE
            }
            else {
                legend = FALSE
            }
            
            printCoefmat(tables, digits = digits, signif.stars = TRUE, 
                signif.legend = legend)
                
            cat(" \n")
            if (!is.null(x$rho[i])) {
                    cat(paste("Spatial autoregressive parameter:", 
                    
                      round(x$rho[i], 4), sep = " ", collapse = ""), 
                      "\n", sep = "")
                }
            cat("\n _______________________________________________________ \n")
        }


    invisible(x)
}




invfrac<-function(Rmat, p){
	#if(p<0) Rmat<-solve(Rmat)
	#print(p)
	#print(solve(Rmat))
ee<-eigen(Rmat, symmetric = TRUE)
Vmat<-ee$vectors
Dmat<-diag(ee$values)
Dmatp<-Dmat^abs(p)
Rmatoh<-Vmat %*% Dmatp %*% solve(Vmat)
if(p<0) Rmatoh<-solve(Rmatoh)
return(Rmatoh)	
	}


`print.spsugm` <-
function(x, digits = max(3, getOption("digits") - 2), ...) {
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2,
                      quote = FALSE)
    } else {
        cat("No coefficients\n")
    }

    ## add printing of error variance parameters
    cat("\n")
    ec <- x$errcomp
    if (length(ec)) {
        cat("Error covariance parameters:\n")
        print.default(format(ec, digits = digits), print.gap = 2,
                      quote = FALSE)
    }

    else cat("No error covariance parameters\n")
    cat("\n")

    invisible(x)
}

