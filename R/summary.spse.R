`summary.spse` <-
function(object,...){


    	coeff<-object$coefficients
		eq<-object$EQ
		var<-as.matrix(object$vcov)
		ser<-sqrt(diag(as.matrix(var)))
		tr<-coeff/ser
		pr<-pnorm(abs(as.matrix(tr)), lower.tail=FALSE)*2
		Xnam<-object$Xnames
		Ynam<-object$Ynames
		reg<-object$K
		sp<-object$spec
		lags<-object$lags
		errors<-object$errors
		endogenous<-object$endogenous
		rho<-object$rho

if(object$type=='spseml') {
	std.errl <- sqrt(diag(object$vcov.errcomp))
	pvlam<-2*pnorm(abs(rho/std.errl),lower.tail=FALSE)
	}

if(object$type=='spsegm'){		
		nlags<-length(which(unlist(lags)==TRUE))
      nend<-length(which(unlist(endogenous)==TRUE))
		numx<-vector("numeric",eq)
for (i in 1:eq) numx[i]<- sp[i] + length(which(lags[[i]]==TRUE)) + length(which(endogenous[[i]]==TRUE))
}
else numx<-sp
		tmp<-seq(1,eq)
      tmp2<-rep(tmp,numx)
		b<-split(coeff,tmp2)
		se<-split(as.matrix(ser), tmp2)
		t<-split(tr,tmp2)
		p<-split(pr,tmp2)
		object$b <- b
		object$se<-se
		object$t<-t
		object$p<-p
		object$eq<-eq
		object$xnam<-Xnam
		object$ynam<-Ynam
		object$sp<-object$spec
		object$lags<-object$lags
		object$errors<-object$errors
		object$endogenous<-object$endogenous
if(object$type=='spseml') object$pvlam<-pvlam
		object$type<-object$type
		object$model<-object$model
		class(object)<- c("summary.spse","spse")
		object
  	} 
