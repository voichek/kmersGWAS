"polylik" <- function(coeffs,y,desmat,relmat,eigenRes,fixh2,trait.type,
		startlik=.Machine$double.xmax-2,scaleh2=1) {
	#print(coeffs)
	if (!missing(fixh2)) {
		h2 <- fixh2
		fixeff <- coeffs[1:(length(coeffs)-1)]
	} else {
		h2 <- coeffs[(length(coeffs)-1)] # coeffs[(length(coeffs))] 
		fixeff <- coeffs[1:(length(coeffs)-2)]
	}
	h2 <- h2*scaleh2
	tvar <- coeffs[(length(coeffs))] # var(y) 
	if (h2<1e-8 || h2>(1-1e-8) || tvar<1e-6) return(startlik+1)
	nids <- length(y)
	sigma <- h2*tvar*relmat + (1-h2)*tvar*diag(x=1,ncol=nids,nrow=nids)
	if (trait.type=="binomial") {
		expb <- exp(desmat %*% fixeff)
		qt <- y - expb/(1.+expb)
	} else {
		qt <- y - desmat %*% fixeff
	}
	es <- 1./diag(t(eigenRes$vectors) %*% (sigma) %*% eigenRes$vectors)
	ginvsig <- eigenRes$vectors %*% diag(es,ncol=length(qt)) %*% t(eigenRes$vectors)
	a <- determinant(sigma,logarithm=T)
	a <- a$modulus * a$sign
	b <- (t(qt) %*% ginvsig %*% qt)
	loglik <- a+b
	loglik[1,1]
}
"polylik_eigen" <- function(coeffs,y,desmat,relmat,eigenRes,fixh2,trait.type,
		startlik=.Machine$double.xmax-2,scaleh2=1) {
	if (!missing(fixh2)) {
		h2 <- fixh2
		fixeff <- coeffs[1:(length(coeffs)-1)]
	} else {
		h2 <- coeffs[(length(coeffs)-1)] # coeffs[(length(coeffs))] 
		fixeff <- coeffs[1:(length(coeffs)-2)]
	}
	h2 <- h2*scaleh2
	tvar <- coeffs[(length(coeffs))]
	if (h2<1e-8 || h2>(1-1e-8) || tvar<1e-6) return(startlik+1)
	if (trait.type=="binomial") {
		expb <- exp(desmat %*% fixeff)
		qt <- y - expb/(1.+expb)
	} else {
		qt <- y - desmat %*% fixeff
	}
	y_0 <- as.vector(crossprod(eigenRes$vectors, qt))
	eigenvalueCOV <- tvar*(eigenRes$values*h2+1.-h2)
	loglik <- sum(log(abs(eigenvalueCOV))) + sum(y_0^2/eigenvalueCOV)
	loglik
}
