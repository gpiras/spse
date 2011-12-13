##still needs to check the SUR- lag model..... In particular check the expression for the likelihood function

spseml<-function(formula, data=list(),panel=TRUE,index=NULL,w,method="eigen", quiet=NULL, model = c("lag","error"), zero.policy=NULL, interval=NULL, tol.solve=1.0e-10, trs=NULL, control=list(), initval=NULL ){

#####################################################################
        timings <- list()
        .ptime_start <- proc.time()
        con <- list(tol.opt=.Machine$double.eps^0.5,
            fdHess=NULL, optimHess=FALSE, compiled_sse=FALSE, Imult=2,
            cheb_q=5, MC_p=16, MC_m=30, super=FALSE)
        nmsC <- names(con)
        con[(namc <- names(control))] <- control
    
    if (length(noNms <- namc[!namc %in% nmsC])) 
            warning("unknown names in control: ", paste(noNms, collapse = ", "))

    if (is.null(quiet)) 
	quiet <- !get("verbose", env = spdep:::.spdepOptions)
    stopifnot(is.logical(quiet))

	if (is.null(zero.policy))
            zero.policy <- get.ZeroPolicyOption()
        stopifnot(is.logical(zero.policy))

if (!inherits(formula, "list")) stop("formula should be a list in simultaneous equation models")

if(panel){
	
	if (!is.null(index)) {
        require(plm)
        data <- plm.data(data, index)
    }
    
    cl <- match.call()
    require(nlme)
    if (is.matrix(w)) {
        if ("listw" %in% class(w)) {
            require(spdep)
            w <- listw2mat(w)
        }
        else {
            stop("w has to be either a 'matrix' or a 'listw' object")
        }
    }


listw<-w


	if (!inherits(listw, "listw")) stop("No neighbourhood list")
        if (is.null(con$fdHess)) con$fdHess <- method != "eigen"
        stopifnot(is.logical(con$fdHess))
	can.sim <- FALSE
	if (listw$style %in% c("W", "S")) 
		can.sim <- spdep:::can.be.simmed(listw)

    
    
    index <- data[, 1]
    tindex <- data[, 2]
    shortt<-unique(tindex)
    shorti<-unique(index)
    n<-length(shorti)
    eq<-t<-length(shortt)

if (n != length(listw$neighbours)) stop("listw objects and variables have different dimension")
if(eq != length(formula)) stop("Number of equations and time periods in the data should correspond")    
    balanced <- (n*t) == nrow(data)
    if (!balanced) 
        stop("Estimation method unavailable for unbalanced data")

#print(listw)
xlist<-vector("list",length=eq)
ylist<-vector("list",length=eq)
Wxlist<-vector("list",length=eq)
Wylist<-vector("list",length=eq)
K<-vector("numeric",length=eq)

for (i in 1:eq){
    xlist[[i]] <- model.matrix(formula[[i]], data = data[tindex == shortt[i], ])
   colnames(xlist[[i]]) <- attributes(model.matrix(formula[[i]], data = data[tindex == shortt[i], ]) )$dimnames[[2]]
    ylist[[i]] <- model.response(model.frame(formula[[i]], data = data[tindex == shortt[i], ]))    
    names(ylist[[i]])[1]<-attributes(model.frame(formula[[i]], data = data[tindex == shortt[i], ]))$names[1] 

    Wylist[[i]] <- lag.listw(listw, ylist[[i]])
    K[i] <- dim(xlist[[i]])[[2]]
    Wxlist[[i]] <- lag.listw(listw, xlist[[i]])
}
	}
	
else{

n<-dim(data)[[1]]
if (n != length(listw$neighbours)) stop("listw objects and variables have different dimension")
eq<-length(formula)

xlist<-vector("list",length=eq)
ylist<-vector("list",length=eq)
Wxlist<-vector("list",length=eq)
Wylist<-vector("list",length=eq)
K<-vector("numeric",length=eq)

for (i in 1:eq){
    xlist[[i]] <- model.matrix(formula[[i]], data = data)
    Wxlist[[i]] <- lag.listw(listw, xlist[[i]])
    ylist[[i]] <- as.matrix(model.frame(formula[[i]], data = data)[1])
    Wylist[[i]] <- lag.listw(listw, ylist[[i]])
    K[i] <- dim(xlist[[i]])[[2]]
}
}

allnames<-NULL
Xnames<-vector("list",length=eq)
Ynames<-NULL

for (i in 1:eq) {
hold<-colnames(xlist[[i]])
Xnames[[i]]<-hold
if(panel) holdy<-names(ylist[[i]])[1]
else holdy<-colnames(ylist[[i]])
Ynames<-c(Ynames, holdy)
	allnames<-c(allnames,hold)
	}
xall<-matrix(unlist(xlist),n,sum(K))
colnames(xall)	<-allnames

	switch(model, lag = if (!quiet) cat("\nSimultaneous system of equations with spatial lags\n"),
	    error = if (!quiet) cat("\nSimultaneous system of equations with error correlation\n"),
	    stop("\nUnknown model type\n"))


yvec<-matrix(unlist(ylist), n*eq,1)
Wyvec<-matrix(unlist(Wylist), n*eq,1)
ncolvecx<-sum(K)
pos<-c(0,cumsum(K))
xvec<-Matrix(0,n*eq,ncolvecx)
Wxvec<-Matrix(0,n*eq,ncolvecx)
for (i in 1:eq){
	lower<- n*(i-1)+1
	upper<- n*i
 	xvec[lower:upper,((pos[i]+1):pos[i+1])]<-xlist[[i]]
 	Wxvec[lower:upper,((pos[i]+1):pos[i+1])]<-Wxlist[[i]]
 	}

##ols estimation eq by eq

bols<-vector("list",length=eq)
eols<-matrix(,nrow=n,ncol=eq)
for (i in 1:eq){
	bols[[i]]<-solve(crossprod(xlist[[i]]),crossprod(xlist[[i]], ylist[[i]]) )
	eols[,i]<-ylist[[i]] - xlist[[i]] %*% bols[[i]]
}

sigma<-crossprod(eols)
sigmainv<-solve(sigma)
sigmasys<-kronecker(sigmainv,Diagonal(n))


sigmasysX<-sigmasys %*% xvec  
sigmasysY<-sigmasys %*% yvec
bgls<- as.matrix(solve(crossprod(xvec,sigmasysX), crossprod(xvec,sigmasysY)))
yf<-as.matrix(xvec%*%bgls)
egls<-yvec-yf
eglsmat<-matrix(as.numeric(egls),n,eq)

if (!is.null(initval)) lambda<-initval
else lambda<-rep(0.2,eq) ##initial values for lambda improve

env <- new.env(parent=globalenv())
assign("yvec",yvec, envir=env)
assign("xvec",xvec, envir=env)
assign("eglsmat",eglsmat, envir=env)
assign("ylist",ylist, envir=env)
assign("xlist",xlist, envir=env)
assign("Wylist",Wylist, envir=env)
assign("Wxlist",Wxlist, envir=env)
assign("eq",eq, envir=env)
assign("n",n, envir=env)
assign("K",K, envir=env)
assign("family", "SAR", envir=env)
assign("verbose", !quiet, envir=env)
assign("can.sim", can.sim, envir=env)
assign("listw", listw, envir=env)
assign("compiled_sse", con$compiled_sse, envir=env)
assign("similar", FALSE, envir=env)
timings[["set_up"]] <- proc.time() - .ptime_start
.ptime_start <- proc.time()



	if (!quiet) cat("Jacobian calculated using ")
	switch(method, 
		eigen = {
                    if (!quiet) cat("neighbourhood matrix eigenvalues\n")
                    eigen_setup(env)
                    er <- get("eig.range", envir=env)
                    if (is.null(interval)) 
                        interval <- c(er[1]+.Machine$double.eps, 
                            er[2]-.Machine$double.eps)
                },
	        Matrix = {
		    if (listw$style %in% c("W", "S") && !can.sim)
		    stop("Matrix method requires symmetric weights")
		    if (listw$style %in% c("B", "C") && 
			!(is.symmetric.glist(listw$neighbours, listw$weights)))
		    stop("Matrix method requires symmetric weights")
                    if (listw$style == "U") stop("U style not permitted, use C")
		    if (!quiet) cat("sparse matrix Cholesky decomposition\n")
	            Imult <- con$Imult
                    if (is.null(interval)) {
	                if (listw$style == "B") {
                            Imult <- ceiling((2/3)*max(sapply(listw$weights,
                                sum)))
	                    interval <- c(-0.5, +0.25)
	                } else interval <- c(-1, 0.999)
                    }
                    Matrix_setup(env, Imult, con$super)
                    W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
        	    I <- as_dsCMatrix_I(n)
		},
	        spam = {
                    if (!require(spam)) stop("spam not available")
		    if (listw$style %in% c("W", "S") && !can.sim)
		    stop("spam method requires symmetric weights")
		    if (listw$style %in% c("B", "C", "U") && 
			!(is.symmetric.glist(listw$neighbours, listw$weights)))
		    stop("spam method requires symmetric weights")
		    if (!quiet) cat("sparse matrix Cholesky decomposition\n")
                    spam_setup(env)
                    W <- as.spam.listw(get("listw", envir=env))
                    if (is.null(interval)) interval <- c(-1,0.999)
		},
                Chebyshev = {
		    if (listw$style %in% c("W", "S") && !can.sim)
		        stop("Chebyshev method requires symmetric weights")
		    if (listw$style %in% c("B", "C", "U") && 
			!(is.symmetric.glist(listw$neighbours, listw$weights)))
		        stop("Chebyshev method requires symmetric weights")
		    if (!quiet) cat("sparse matrix Chebyshev approximation\n")
                    cheb_setup(env, q=con$cheb_q)
                    W <- get("W", envir=env)
        	    I <- as_dsCMatrix_I(n)
                    if (is.null(interval)) interval <- c(-1,0.999)
                },
                MC = {
		    if (!listw$style %in% c("W"))
		        stop("MC method requires row-standardised weights")
		    if (!quiet) cat("sparse matrix Monte Carlo approximation\n")
                    mcdet_setup(env, p=con$MC_p, m=con$MC_m)
                    W <- get("W", envir=env)
        	    I <- as_dsCMatrix_I(n)
                    if (is.null(interval)) interval <- c(-1,0.999)
                },
                LU = {
		    if (!quiet) cat("sparse matrix LU decomposition\n")
                    LU_setup(env)
                    W <- get("W", envir=env)
                    I <- get("I", envir=env)
                    if (is.null(interval)) interval <- c(-1,0.999)
                },
		stop("...\nUnknown method\n"))

nm <- paste(method, "set_up", sep="_")
timings[[nm]] <- proc.time() - .ptime_start
.ptime_start <- proc.time()


##likelihood optimization
if(model=="error"){ 
spatpar <- nlminb(lambda, llsurerror, env=env)
bgls<-get("bgls", envir=env)
egls<-get("egls", envir=env)
yfgls<-get("yfgls", envir=env)
lambdas<-spatpar$par
llik<-spatpar$objective
optres<-spatpar
S2<-crossprod(egls)/n
xVxinv<-get("xVxinv", envir=env)
VC<- xVxinv
covPRL <- solve(-fdHess(lambdas, function(x) -llsurerror(x,env=env) )$Hessian )
}
else{
spatpar <- nlminb(lambda, llsurlag, env=env)
bgls<-get("bols", envir=env)
egls<-get("res", envir=env)
yfgls<-get("yfols", envir=env)
lambdas<-spatpar$par
#print(lambdas)
llik<-spatpar$objective
optres<-spatpar
S2<-crossprod(egls)/n
vcgls<-get("vcgls", envir=env)

#print(-fdHess(lambdas, function(x) -llsurlag(x,env=env) )$Hessian )
covPRL <- solve(-fdHess(lambdas, function(x) -llsurlag(x,env=env) )$Hessian )

####
ytlist<-ylist
for (i in 1:eq) ytlist[[i]]<-ylist[[i]] - lambdas[i]*Wylist[[i]]

ytvec<-matrix(unlist(ytlist), n*eq,1)
ncolxvec<-sum(K)
pos<-c(0,cumsum(K))
xvec<-Matrix(0,n*eq,ncolxvec)
for (i in 1:eq){
	lower<- n*(i-1)+1
	upper<- n*i
	xvec[lower:upper,((pos[i]+1):pos[i+1])]<-xlist[[i]]
 	}
mod.filt<-lm(as.matrix(ytvec)~as.matrix(xvec)-1)
res.filt<-residuals(mod.filt)
s2.filt<-crossprod(res.filt)/(n-sum(K))
xVx<-crossprod(xvec)
xVxinv<-as.numeric(s2.filt)*solve(xVx)
####
VC<- xVxinv
	}



type<-"spseml"
model.data <- list(ylist,xlist)


spmod <- list(method=method, coefficients=bgls, errcomp=NULL, vcov=VC, 
			  vcov.errcomp= covPRL, residuals=NULL, fitted.values=NULL,
			  sigma2=NULL, type=type, model= model, model.data=model.data,  N=n,
			  EQ=eq,K=sum(K), call=cl,terms=NULL, Xnames=Xnames,Ynames=Ynames, 
			  spec=K,lags=NULL, errors=NULL, endogenous=NULL, rho=lambdas )

class(spmod)<- "spse"
return(spmod)

#
#
#



#res<-list(bgls,egls,spatpar)
	}

#
#
llsurlag<-function(lambda,env){
	eglsmat <- get("eglsmat", envir=env)
	yvec <- get("yvec", envir=env)
	xvec <- get("xvec", envir=env)
	ylist <- get("ylist", envir=env)
	xlist <- get("xlist", envir=env)
	Wylist <- get("Wylist", envir=env)
	Wxlist <- get("Wxlist", envir=env)
	eq <- get("eq", envir=env)
	n <- get("n", envir=env)
	K <- get("K", envir=env)
	

ytlist<-ylist
for (i in 1:eq) ytlist[[i]]<-ylist[[i]] - lambda[i]*Wylist[[i]]

bols<-vector("list",length=eq)
for (i in 1:eq) bols[[i]]<- coefficients(lm(ytlist[[i]]~ xlist[[i]]-1))
assign("bols",matrix(unlist(bols)), envir=env)

Xbeta<-vector("list", length=eq)
for (i in 1:eq) Xbeta[[i]]<- xlist[[i]] %*% bols[[i]]

yfols<-matrix(unlist(Xbeta),eq*n,1)
assign("yfols",yfols,envir=env)

res<- matrix(unlist(ytlist),eq*n,1) - yfols
assign("res",res,envir=env)
resmat<-matrix(res,n,eq)

ZpZ<-crossprod(resmat)
Sigma <-ZpZ/n
Sigmai<-solve(Sigma)
vcgls<-kronecker(Sigmai,Diagonal(n))
#print(dim(vcgls))
assign("vcgls",vcgls, envir=env)

deti<-c(0,eq)
for(i in 1:eq) 	deti[i] <- do_ldet(lambda[i], env)
#print(deti)
deter<-sum(deti) 
#print(class(deter))
#print(class(ZpZ))	
const <- - (n/2) - (n/2) * log(2 * pi)
obj<-  const - deter + (n/2)* log(det(Sigma)) + (1/2)*crossprod(res)

#print(as.numeric(obj))
as.numeric(obj)
	}
	




llsurerror<-function(lambda,env){
	eglsmat <- get("eglsmat", envir=env)
	yvec <- get("yvec", envir=env)
	xvec <- get("xvec", envir=env)
	ylist <- get("ylist", envir=env)
	xlist <- get("xlist", envir=env)
	Wylist <- get("Wylist", envir=env)
	Wxlist <- get("Wxlist", envir=env)
	eq <- get("eq", envir=env)
	n <- get("n", envir=env)
	K <- get("K", envir=env)
	listw<-get("listw", envir=env)
	
Z<-matrix(,nrow(eglsmat),ncol(eglsmat))
for(i in 1:eq) Z[,i]<- eglsmat[,i] - lambda[i]* lag.listw(listw,eglsmat[,i])
ZpZ<-crossprod(Z)
Sigma <-ZpZ/n
Sigmai<-solve(Sigma)

deti<-c(0,eq)
for(i in 1:eq) 	deti[i] <- do_ldet(lambda[i], env)
#print(deti)
deter<-sum(deti) 
	
vcgls<-kronecker(Sigmai,Diagonal(n))


ytlist<-ylist
xtlist<-xlist
for (i in 1:eq) ytlist[[i]]<-ylist[[i]] - lambda[i]*Wylist[[i]]
for (i in 1:eq) xtlist[[i]]<-xlist[[i]] - lambda[i]*Wxlist[[i]]

ytvec<-matrix(unlist(ytlist), n*eq,1)
ncolxvec<-sum(K)
pos<-c(0,cumsum(K))
xtvec<-Matrix(0,n*eq,ncolxvec)
for (i in 1:eq){
	lower<- n*(i-1)+1
	upper<- n*i
	xtvec[lower:upper,((pos[i]+1):pos[i+1])]<-xtlist[[i]]
 	}
VX<-vcgls %*% xtvec 
VY<- vcgls %*% ytvec

#xt<-Matrix(0,n*eq,2*k)
#for (i in 1:eq) xt[(i*n-n+1):(i*n),(i*k-k+1):(i*k)]<-x[(i*n-n+1):(i*n),] 
#
xVx<-crossprod(xtvec,as.matrix(VX))
xVy<-crossprod(xtvec,as.matrix(VY))
xVxinv<-solve(as.matrix(xVx))
assign("xVxinv",xVxinv, envir=env)
bgls<- xVxinv %*% xVy
assign("bgls",bgls, envir=env)
#print(bGLS)
yfgls<-xvec %*% bgls
assign("yfgls",yfgls, envir=env)
egls<-as.numeric(yvec - yfgls)
assign("egls",egls, envir=env)
for(i in 1:eq) Z[,i]<- eglsmat[,i] - lambda[i]* lag.listw(listw,eglsmat[,i])
vcglsee<-vcgls %*% matrix(Z,n*eq,1)
ter<-crossprod(matrix(Z,n*eq,1),vcglsee)
const <- - (n)/2 * log(2 * pi)
obj<-  const - deter +  (1/2) *ter + (n/2)* log(det(Sigma)) 
#print(as.numeric(obj))
as.numeric(obj)
	}
	
	
	
	
