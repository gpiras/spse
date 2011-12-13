`searg` <-function (rhopar, v, verbose = FALSE) {
    vv <- v$bigG %*% c(rhopar[1], rhopar[1]^2, rhopar[2]) - v$litg
    value <- sum(vv^2)
    if (verbose) 
        cat("function:", value, "rho:", rhopar[1], "sig2:", 
            rhopar[2], "\n")
    value
}


`tslssp` <-function(y,yend,X,Zinst,robust=FALSE) {
	Q <- cbind(X,Zinst)
	Z <- cbind(yend,X)
	df <- nrow(Z) - ncol(Z)
	QQ <- crossprod(Q,Q)
	Qye <- crossprod(Q,yend)
	bz <- solve(QQ,Qye)
	yendp <- Q %*% bz
	Zp <- cbind(yendp,X)
	ZpZp <- crossprod(Zp,Zp)
	ZpZpi <- solve(ZpZp)
	Zpy <- crossprod(Zp,y)
	biv <- crossprod(ZpZpi,Zpy)
	#print(biv)
	rownames(biv)<-c(colnames(yend), colnames(X))
	yp <- Z %*% biv
	e <- y - yp
    	s2 <- crossprod(e,e) / df
	    varb <- ZpZpi * s2[1,1]
	    sebiv <- sqrt(diag(varb))
	    tbiv <- biv / sebiv
	    pbiv <- pnorm(abs(tbiv),lower.tail=FALSE) * 2
	    result <- list(coefficients=biv,se=sebiv,t=tbiv,
	          p=pbiv,var=varb,s2=s2,
	          residuals=e,yhat=yp)
	result
}


`modtslssp` <-function(y,yend,X,Zinst) {
	Q <- Zinst
	Z <- cbind(yend,X)
	df <- nrow(Z) - ncol(Z)
	QQ <- crossprod(Q,Q)
	Qye <- crossprod(Q,yend)
	bz <- solve(QQ,Qye)
	yendp <- Q %*% bz
	Zp <- cbind(yendp,X)
	ZpZp <- crossprod(Zp,Zp)
	ZpZpi <- solve(ZpZp)
	Zpy <- crossprod(Zp,y)
	biv <- crossprod(ZpZpi,Zpy)
	#rownames(biv)<-c(colnames(yend), colnames(X))
	yp <- Z %*% biv
	e <- y - yp
    	s2 <- crossprod(e,e) / df
	    varb <- ZpZpi * s2[1,1]
	    sebiv <- sqrt(diag(varb))
	    tbiv <- biv / sebiv
	    pbiv <- pnorm(abs(tbiv),lower.tail=FALSE) * 2
	    result <- list(coefficients=biv,se=sebiv,t=tbiv,
	          p=pbiv,var=varb,s2=s2,
	          residuals=e,yhat=yp)
	    
	result
}


`Ggsararsp` <-
function (W, u, zero.policy = FALSE) 
{
      n <- length(u)
      tt<-matrix(0,n,1)
      tr<-sum(unlist(W$weights)^2)
      wu<-lag.listw(W,u)
      wwu<-lag.listw(W,wu)
    	uu <- crossprod(u, u)
    	uwu <- crossprod(u, wu)
 	uwpuw <- crossprod(wu, wu)
    	uwwu <- crossprod(u, wwu)
    	wwupwu <- crossprod(wwu, wu)
    	wwupwwu <- crossprod(wwu, wwu)
    	bigG <- matrix(0, 3, 3)
    	bigG[, 1] <- c(2 * uwu, 2 * wwupwu, (uwwu + uwpuw))/n
    	bigG[, 2] <-  -c(uwpuw, wwupwwu, wwupwu)/n
    	bigG[, 3] <- c(1, tr/n, 0)
    	litg <- c(uu, uwpuw, uwu)/n
    	list(bigG = bigG, litg = litg)
}
