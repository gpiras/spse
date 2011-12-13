print.summary.spse <- function(x,digits= max(3, getOption("digits") - 2),width=getOption("width"),...) {

        cat("\nSimultaneous Equations Model:\n")
        cat("\nCall:\n")
        print(x$call)
        eq<-x$eq
        save.digits <- unlist(options(digits=digits))
        on.exit(options(digits=save.digits))
		  numx<-vector("numeric",eq)
		  
        for (i in 1:eq) {
            cat(" \n" )
            cat(paste('Equation', i, sep = " ", collapse = ""),"\n", sep = "")
            #cat("\n _______________________________________________________ \n")
            tables<- cbind( x$b[[i]], x$se[[i]], x$t[[i]], x$p[[i]] )
if(x$type=='spsegm'){
if(length(which(x$lags[[i]]==TRUE)) !=0 ){
	Wynames<-paste("W",x$ynam[which(x$lags[[i]]==TRUE)])
   rn<-c( x$ynam[which(x$endogenous[[i]]==TRUE)], x$xnam[[i]], Wynames )
	} 
	else rn<-c(x$ynam[which(x$endogenous[[i]]==TRUE)], x$xnam[[i]] )
        }
        
if(x$type=='spseml'){
	if(!is.null(x$rho[i])) rn<- x$xnam[[i]]
	else rn<- c(paste("W",x$ynam[[i]]), x$xnam[[i]] )
	
	}
            dimnames(tables)<-list(rn ,c("Estimate", "Std.Error", "t value", "Pr(>|t|)"))


            if(i==eq) {
                legend=TRUE
            } else {
                legend=FALSE
            }

            printCoefmat(tables,digits=digits, signif.stars=TRUE,signif.legend=legend)
       cat(" \n" )     
     if(!is.null(x$rho[i])){
     	if(x$type=='spseml'){
     		if(x$model=="lag")cat(paste('Spatial autoregressive parameter:', round(x$rho[i],4),"Pr(>|t|)", round(x$pvlam[i],4), sep = " ", collapse = ""),"\n", sep = "")
     		else cat(paste('Spatial autocorrelation coefficient:', round(x$rho[i],4),"Pr(>|t|)", round(x$pvlam[i],4), sep = " ", collapse = ""),"\n", sep = "")
     		} 
     	  else{
            if(x$errors[i] != FALSE) cat(paste('Spatial autoregressive parameter:', round(x$rho[i],4), sep = " ", collapse = ""),"\n", sep = "")
            }
}
            cat("\n _______________________________________________________ \n")
            }
}