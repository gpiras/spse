`spsegm` <-
function(formula,data=list(), panel=FALSE,index=NULL,  w, method='spatialsim', lags=NULL, errors=NULL, endogenous=NULL, zero.policy = FALSE){

#source("tslssp.R")
#source("Ggsararsp.R")
#source("searg.R")
#source("summary.splm.R")
#source("print.summary.splm.R")

type<-'spsegm'

if(!inherits(formula, "list")) stop("formula should be a list in simultaneous equation models")

if(!inherits(w, c("listw","matrix"))) stop("w has to be an object of class 'listw' of 'matrix'")

    if(is.matrix(w)) {
            require(spdep)
            w <- mat2listw(w)
        }

listw<-w        
cl <- match.call()

if(panel){
	
	if (!is.null(index)) {
        require(plm)
        data <- plm.data(data, index)
    }
    
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
#Wxlist<-vector("list",length=eq)
#WWxlist<-vector("list",length=eq)
Wylist<-vector("list",length=eq)
K<-vector("numeric",length=eq)

#print(attributes(model.frame(formula[[i]], data = data[tindex == shortt[i], ]))$names[1] )
for (i in 1:eq){
    xlist[[i]] <- model.matrix(formula[[i]], data = data[tindex == shortt[i], ])
   colnames(xlist[[i]]) <- attributes(model.matrix(formula[[i]], data = data[tindex == shortt[i], ]) )$dimnames[[2]]
    ylist[[i]] <- model.response(model.frame(formula[[i]], data = data[tindex == shortt[i], ]))    
    names(ylist[[i]])[1]<-attributes(model.frame(formula[[i]], data = data[tindex == shortt[i], ]))$names[1] 

    Wylist[[i]] <- lag.listw(listw, ylist[[i]])
    K[i] <- dim(xlist[[i]])[[2]]
#int <- ifelse(colnames(xlist[[i]])[1] == "(Intercept)", 2, 1)
#	wx <- matrix(nrow = n, ncol = (K[i] - (int - 1)))
 #       for (j in int : K[i]) {
  #          wx[,j- (int - 1)] <- lag.listw(listw, xlist[[i]][,j])
   #     }    
#Wxlist[[i]]<- wx
#WWxlist[[i]]<- lag.listw(listw, Wxlist[[i]])
}

	}
	
else{

n<-dim(data)[[1]]
if (n != length(listw$neighbours)) stop("listw objects and variables have different dimension")
eq<-length(formula)

xlist<-vector("list",length=eq)
ylist<-vector("list",length=eq)
#Wxlist<-vector("list",length=eq)
#WWxlist<-vector("list",length=eq)
Wylist<-vector("list",length=eq)
K<-vector("numeric",length=eq)

for (i in 1:eq){
    xlist[[i]] <- model.matrix(formula[[i]], data = data)
#    Wxlist[[i]] <- lag.listw(listw, xlist[[i]])
    ylist[[i]] <- as.matrix(model.frame(formula[[i]], data = data)[1])
    Wylist[[i]] <- lag.listw(listw, ylist[[i]])
    K[i] <- dim(xlist[[i]])[[2]]
#int <- ifelse(colnames(xlist[[i]])[1] == "(Intercept)", 2, 1)
#	wx <- matrix(nrow = n, ncol = (K[i] - (int - 1)))
 #       for (j in int : K[i]) {
  #          wx[,j- (int - 1)] <- lag.listw(listw, xlist[[i]][,j])
   #     }    
#Wxlist[[i]]<- wx
# WWxlist[[i]]<- lag.listw(listw, Wxlist[[i]])
}
}

if(!is.null(lags)){
if(length(lags)!= eq) stop("The length of lags should equal the number of equations")
	}

if(!is.null(errors)){
if(length(errors)!= eq) stop("The length of errors should equal the number of equations")
	}

if(!is.null(endogenous)){
#	print(endogenous)
if(length(endogenous)!= eq) stop("The length of errors should equal the number of equations")
for(i in 1:eq){
	if(endogenous[[i]][i]==TRUE){
		endogenous[[i]][i]<-FALSE 		
		warning(paste("Dependent variable cannot be on the right hand side in equation",i))
#		print(endogenous[[i]][i])
		} 
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

allnames<-unique(allnames)
sel<-vector("numeric",length=length(allnames))
for (j in 1:length(allnames)) sel[j]<-unique(which(colnames(xall)==allnames[j]))[1]
xall<-xall[,sel]
if(!is.null(which(colnames(xall)=="(Intercept)"))) {
	xallni<-xall[,-which(colnames(xall)=="(Intercept)")]
	}
Wxall<-lag.listw(listw,xallni)
WWxall<-lag.listw(listw,Wxall)

inst<-cbind(Wxall,WWxall)

##########################################################
#if(is.null(endogenous) && is.null(lags)){
#2sls<-FALSE	
#	}




if(is.null(endogenous)){
endogenous<-vector("list",eq)	
	for (i in 1:eq){
		endogenous[[i]]<-rep(TRUE,(eq))
		endogenous[[i]][i]<-FALSE
		} 
	}

if(is.null(lags)){
lags<-vector("list",eq)	
	for (i in 1:eq) lags[[i]]<-rep(TRUE,(eq))
	}

ymat<-matrix(unlist(ylist), n,eq)
Wymat<-matrix(unlist(Wylist),n,eq)

####FIRST STEP ESTIMATION
r<-matrix(,n,eq)
b<-vector(mode="list",eq)

		for (i in 1:eq){
Wy<-Wymat[,which(lags[[i]]==TRUE)]
#yhold<-matrix(ymat[,-i],nrow=nrow(ymat),ncol=ncol(ymat)-1)
end<-ymat[,which(endogenous[[i]] ==TRUE)]
yend<-cbind(Wy,end)
if(!is.null(yend)){
est<-tslssp(ylist[[i]],yend,xlist[[i]],inst)
r[,i]<-est$residuals
b[[i]]<-est$coefficients
#print(b[[i]])
}
else{
est<-lm(ylist[[i]]~xlist[[i]])	
r[,i]<-est$residuals
b[[i]]<-est$coefficients
	}
	} 
#print(b)
#once the first step coefficients are obtained, one can proceed to perform the GM estimator

if(is.null(errors)) errors<-rep(TRUE,eq)

et<-length(which(errors==TRUE))
rho<-vector("numeric",length=et)
sigma<-vector("numeric",length=eq)
i=1
while (i <=eq){
if(errors[[i]]){	
	mom<-Ggsararsp(u=r[,i],W=listw)
	pars<-c(0,0)
	estim <- nlminb(pars, searg, v = mom, verbose =FALSE)
	rho[i]<-estim$par[1]
	sigma[i]<-estim$par[2]
	}
	i=i+1
}

#print(rho)




H<-cbind(xall,Wxall,WWxall)
HH<-crossprod(H,H)
HHinv<-solve(HH)
Hp<-t(H)
P<-H%*%HHinv%*%Hp


##transform the y's and Wy's
xtlist<-xlist
ytlist<-ylist
Wytlist<-Wylist

r2<-matrix(,n,eq)
b2<-vector(mode="list",eq)

nlags<-length(which(unlist(lags)==TRUE))
nend<-length(which(unlist(endogenous)==TRUE))
Z<-Matrix(0,eq*n,sum(K) +nlags+ nend)
PZ<-Matrix(0,eq*n,sum(K) +nlags+ nend)
Yt<-matrix(,n*eq,1)
pos1<-0



for (i in 1:eq){
if(errors[[i]]){
	xtlist[[i]] <- xlist[[i]] - as.numeric(rho[i]) * lag.listw(listw,xlist[[i]]) 
for (t in 1:eq){	
	ytlist[[t]]<- ylist[[t]] - as.numeric(rho[i]) * lag.listw(listw,ylist[[t]])
	Wytlist[[t]] <- lag.listw(listw,ytlist[[t]])
	}
}
Wytmat<-matrix(unlist(Wytlist), n,eq)	
ytmat<-matrix(unlist(ytlist), n,eq)	
Wyt<-Wytmat[,which(lags[[i]]==TRUE)]
endt<-ytmat[,which(endogenous[[i]] ==TRUE)]
ytend<-cbind(Wyt,endt)
if(!is.null(ytend)){
	est<-tslssp(ytmat[,i],ytend,xtlist[[i]],inst)
	r2[,i]<-est$residuals
	b2[[i]]<-est$coefficients

}
else{
est<-lm(ytlist[[i]]~xtlist[[i]])	
r2[,i]<-est$residuals
b2[[i]]<-est$coefficients
	}
	lower<- n*(i-1)+1
	upper<- n*i
	
	hold<-xtlist[[i]]
	ytend<-cbind(Wyt,endt)
if(!is.null(ytend))	hold<-cbind(endt,hold,Wyt)
 pos1<-pos1+ ncol(hold)  
 	Z[lower:upper,((pos1-ncol(hold)+1):pos1)]<-hold
 	PZ[lower:upper,((pos1-ncol(hold)+1):pos1)]<-P%*%hold
	Yt[lower:upper,]<-ytmat[,i]
	
	} 
#print(b2)
SIGMA<-crossprod(r2)
SIGMAinv<-solve(SIGMA)


AZ <- Matrix(0,n*eq,ncol(Z))
for(i in 1:ncol(Z)){
	tmp <- matrix(Z[,i],n,eq)
tmpt <- t(tmp)
tmp2 <- SIGMAinv%*%tmpt
tmp3 <- matrix(t(tmp2),nrow=n*eq,ncol=1)
AZ[,i]<-tmp3
}

AZH <- matrix(,n*eq,ncol(PZ))
for(i in 1:ncol(PZ)){
	tmp <- matrix(PZ[,i],n,eq)
tmpt <- t(tmp)
tmp2 <- SIGMAinv%*%tmpt
tmp3 <- matrix(t(tmp2),nrow=n*eq,ncol=1)
AZH[,i]<-tmp3
}

#print(AZH)


Ytm <- matrix(Yt,n,eq)
Ytmt <- t(Ytm)
Ay1 <- SIGMAinv%*%Ytmt
Ay <- matrix(t(Ay1),nrow=n*eq,ncol=1)

#print(Ytm)
ZHAZ<-crossprod(PZ,AZ)
ZHAZinv<-solve(ZHAZ)

#print(ZHAZinv)

ZAy<-crossprod(PZ,Ay)
#print(ZAy)

delta<-ZHAZinv%*%ZAy
ZHAZH<-crossprod(PZ,AZH)
VC<-solve(ZHAZH)
model.data <- list(ylist,xlist)

spmod <- list(method=method, coefficients=delta, errcomp=NULL, vcov=VC, 
			  vcov.errcomp= NULL, residuals=NULL, fitted.values=NULL,
			  sigma2=NULL, type=type, model= model.data,  N=n,
			  EQ=eq,K=sum(K), call=cl,terms=NULL, Xnames=Xnames,Ynames=Ynames, 
			  spec=K,lags=lags, errors=errors, endogenous=endogenous, rho=rho )

class(spmod)<- "spse"
return(spmod)
}

