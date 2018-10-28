#' Estimation of polygenic model
#' 
#' This function maximises the likelihood of the data under polygenic 
#' model with covariates an reports twice negative maximum likelihood estimates 
#' and the inverse of the variance-covariance matrix at the point of ML. 
#' 
#' One of the major uses of this function is to estimate residuals of the 
#' trait and the inverse of the variance-covariance matrix for 
#' further use in analysis with \code{\link{mmscore}} and 
#' \code{\link{grammar}}.
#' 
#' Also, it can be used for a variant of GRAMMAR analysis, which 
#' allows for permutations for GW significance by use of 
#' environmental residuals as an analysis trait with \code{\link{qtscore}}.
#' 
#' "Environmental residuals" (not to be mistaken with just "residuals") are 
#' the residual where both the effect of covariates AND the estimated 
#' polygenic effect (breeding values) are factored out. This thus 
#' provides an estimate of the trait value contributed by environment
#' (or, turning this other way around, the part of the trait not explained 
#' by covariates and by the polygene). Polygenic residuals are estimated 
#' as
#' 
#' \deqn{
#' \sigma^2 V^{-1} (Y - (\hat{\mu} + \hat{\beta} C_1 + ...))
#' }
#' 
#' where \eqn{sigma^2} is the residual variance, \eqn{V^{-1}} is the 
#' InvSigma (inverse of the var-cov matrix at the maximum of 
#' polygenic model) and 
#' \eqn{(Y - (\hat{\mu} + \hat{\beta} C_1 + ...))} is the trait 
#' values adjusted for covariates (also at at the maximum of 
#' polygenic model likelihood). 
#' 
#' It can also be used for heritability analysis.
#' If you want to test significance of heritability, 
#' estimate the model and write down 
#' the function minimum reported at the "h2an" element of the output 
#' (this is twice the negative MaxLikelihood). Then do a next round of 
#' estimation, but set fixh2=0. The difference between your function minima 
#' gives a test distributed as chi-squared with 1 d.f.
#' 
#' The way to compute the likelihood is partly based on 
#' the paper of Thompson (see refs), namely instead of 
#' taking the inverse of the var-cov matrix every time, 
#' eigenvectors of the inverse of G (taken only once) 
#' are used.
#' 
#' 
#' @param formula Formula describing fixed effects to be used in the analysis, e.g. 
#' y ~ a + b means that outcome (y) depends on two covariates, a and b. 
#' If no covariates used in the analysis, skip the right-hand side of the 
#' equation.
#' @param kinship.matrix Kinship matrix, as provided by e.g. ibs(,weight="freq"), 
#' or estimated outside of GenABEL from pedigree data.
#' @param data An (optional) object of \code{\link{gwaa.data-class}} or a data frame with 
#' outcome and covariates
#' @param fixh2 Optional value of heritability to be used, instead of maximisation. 
#' The uses of this option are two-fold: (a) testing significance of 
#' heritability and (b) using a priori known heritability to derive the 
#' rest of MLEs and var.-cov. matrix.
#' @param starth2 Starting value for h2 estimate
#' @param trait.type "gaussian" or "binomial"
#' @param opt.method "nlm" or "optim". These two use different optimisation functions. 
#' We suggest using the default \code{\link{nlm}}, although 
#' \code{\link{optim}} may give better results in some situations
#' @param scaleh2 Only relevant when "nlm" optimisation function is used. 
#' "scaleh2" is the heritability 
#' scaling parameter, regulating how "big" are parameter changes in h2 with
#' respect to changes in other parameters. As other parameters are estimated 
#' from previous regression, these are expected to change little from the 
#' initial estimate. The default value of 1000 proved to work rather well under a 
#' range of conditions.
#' @param quiet If FALSE (default), details of optimisation process are reported
#' @param steptol steptal parameter of "nlm"
#' @param gradtol gradtol parameter of "nlm" 
#' @param optimbou fixed effects boundary scale parameter for 'optim'
#' @param fglschecks additional check for convergence on/off (convergence 
#' between estimates obtained and that from FGLS)
#' @param maxnfgls number of fgls checks to perform
#' @param maxdiffgls max difference allowed in fgls checks 
#' @param patchBasedOnFGLS if FGLS checks not passed, 'patch' fixed 
#' effect estimates based on FGLS expectation
#' @param llfun function to compute likelihood (default 'polylik_eigen', also 
#' available -- but not recommended -- 'polylik')
#' @param eigenOfRel results of eigen(relationship matrix = 2*kinship.matrix).  
#' Passing this can decrease computational time substantially if multiple traits 
#' are analysed using the same kinship matrix. This option will not work if any 
#' NA's are found in the trait and/or covariates and if the dimensions of the 
#' 'eigen'-object, trait, covariates, kinship do not match. 
#' @param ... Optional arguments to be passed to \code{\link{nlm}} or (\code{\link{optim}}) 
#' minimisation function
#' 
#' @return 
#' A list with values 
#' \item{h2an}{A list supplied by the \code{\link{nlm}} minimisation routine. 
#' Of particular interest are elements "estimate" containing parameter 
#' maximal likelihood estimates (MLEs) (order: mean, betas for covariates, 
#' heritability, (polygenic + residual variance)). The value of 
#' twice negative maximum log-likelihood
#' is returned as h2an\$minimum.}
#' \item{esth2}{Estimate (or fixed value) of heritability}
#' \item{residualY}{Residuals from analysis, based on covariate effects only; 
#' NOTE: these are NOT grammar "environmental residuals"!}
#' \item{pgresidualY}{Environmental residuals from analysis, based on covariate effects 
#' and predicted breeding value.
#' }
#' \item{grresidualY}{GRAMMAR+ transformed trait residuals}
#' \item{grammarGamma}{list with GRAMMAR-gamma correction factors}
#' \item{InvSigma}{Inverse of the variance-covariance matrix, computed at the 
#' MLEs -- these are used in \code{\link{mmscore}} and \code{\link{grammar}}
#' functions.}
#' \item{call}{The details of call}
#' \item{measuredIDs}{Logical values for IDs who were used in analysis 
#' (traits and all covariates measured) == TRUE}
#' \item{convFGLS}{was convergence achieved according to FGLS criterionE}
#' 
#' 
#' @references 
#' Thompson EA, Shaw RG (1990) Pedigree analysis for quantitative 
#' traits: variance components without matrix inversion. Biometrics 
#' 46, 399-413.
#' 
#' for original GRAMMAR
#' 
#' Aulchenko YS, de Koning DJ, Haley C. Genomewide rapid association using mixed model 
#' and regression: a fast and simple method for genome-wide pedigree-based quantitative 
#' trait loci association analysis. Genetics. 2007 177(1):577-85.
#' 
#' for GRAMMAR-GC
#' 
#' Amin N, van Duijn CM, Aulchenko YS. A genomic background based method for 
#' association analysis in related individuals. PLoS ONE. 2007 Dec 5;2(12):e1274.
#' 
#' for GRAMMAR-Gamma
#' 
#' Svischeva G, Axenovich TI, Belonogova NM, van Duijn CM, Aulchenko YS. Rapid 
#' variance components-based method for whole-genome association analysis. 
#' Nature Genetics. 2012 44:1166-1170. doi:10.1038/ng.2410 
#' 
#' for GRAMMAR+ transformation
#' 
#' Belonogova NM, Svishcheva GR, van Duijn CM, Aulchenko YS, Axenovich TI (2013) 
#' Region-Based Association Analysis of Human Quantitative Traits in Related Individuals. 
#' PLoS ONE 8(6): e65395. doi:10.1371/journal.pone.0065395 
#' 
#' @author Yurii Aulchenko, Gulnara Svischeva
#' 
#' @note 
#' Presence of twins may complicate your analysis. Check the kinship matrix for 
#' singularities, or rather use \code{\link{check.marker}} for identification 
#' of twin samples. Take special care in interpretation.
#' 
#' If a trait (no covariates) is used, make sure that the order of IDs in the
#' kinship.matrix is exactly the same as in the outcome
#' 
#' Please note that there is alternative to 'polygenic', 
#' \code{\link{polygenic_hglm}}, which is faster than 
#' polygenic() with the llfun='polylik' option, but slightly slower than the
#' default polygenic().
#' 
#' @seealso 
#' \code{\link{polygenic_hglm}},
#' \code{\link{mmscore}},
#' \code{\link{grammar}}
#' 
#' @examples 
#' # note that procedure runs on CLEAN data
#' data(ge03d2ex.clean)
#' gkin <- ibs(ge03d2ex.clean,w="freq")
#' h2ht <- polygenic(height ~ sex + age, kin=gkin, ge03d2ex.clean)
#' # estimate of heritability
#' h2ht$esth2
#' # other parameters
#' h2ht$h2an
#' # the minimum twice negative log-likelihood
#' h2ht$h2an$minimum
#' # twice maximum log-likelihood
#' -h2ht$h2an$minimum
#' 
#' # for binary trait (experimental)
#' h2dm <- polygenic(dm2 ~ sex + age, kin=gkin, ge03d2ex.clean, trait="binomial")
#' # estimated parameters
#' h2dm$h2an
#' 
#' @keywords htest
#' 
#' 
"polygenic" <-
	function(formula,kinship.matrix,data,fixh2,starth2=0.3,trait.type="gaussian",
		 opt.method="nlm",scaleh2=1,quiet=FALSE,  
		 steptol=1e-8, gradtol = 1e-8, optimbou = 8, 
		 fglschecks=TRUE,maxnfgls=8,maxdiffgls=1e-4, patchBasedOnFGLS = TRUE, 
		 llfun = "polylik_eigen", eigenOfRel, ...) {
		if (!missing(data)) if (is(data,"gwaa.data")) 
		{
			checkphengen(data)
			data <- phdata(data)
		}
		if (!missing(data)) if (!is(data,"data.frame")) stop("data should be of gwaa.data or data.frame class")
		if (starth2<0 || starth2>1) stop("Starting value of h2 (starth2) must be between 0 and 1")
		ttargs <- c("gaussian","binomial")
		if (!(match(trait.type,ttargs,nomatch=0)>0)) {
			out <- paste("trait.type argument should be one of",ttargs,"\n")
			stop(out)
		}
		if (trait.type=="gaussian") fam <- gaussian()
		if (trait.type=="binomial") fam <- binomial()
		optargs <- c("nlm","optim")
		if (!(match(opt.method,optargs,nomatch=0)>0)) {
			out <- paste("opt.method argument should be one of",optargs,"\n")
			stop(out)
		}
		if (llfun == "polylik_eigen") llFUN <- polylik_eigen
		else if (llfun == "polylik") llFUN <- polylik
		else stop("llfun should be 'polylik' or 'polylik_eigen'")

		if ( is(try(formula,silent=TRUE),"try-error") ) { 
			formula <- data[[as(match.call()[["formula"]],"character")]] 
		}
		if (is(formula, "formula")){
			mf <- model.frame(formula, data, na.action = na.omit, 
					  drop.unused.levels = TRUE)
			ciccio=function(x)nlevels(as.factor(x))    
			livelli=apply(mf,2,ciccio)
			ch=length(livelli)
			bad=which(livelli<2)
			if(length(bad)>0) warning(paste("Variable",names(mf)[bad],"has only one value: Dropped"))
			livelli=livelli[livelli>1]
			mf=mf[names(livelli)]
			if(length(livelli)>1){
				if(length(livelli)==2){ formula=as.formula(paste(names(mf)[1],"~",names(mf)[2]))

				}else{
					formula=(paste(names(mf)[1],"~",names(mf)[2]))
					for(i in 3:length(livelli)){

						formula=paste(formula,"+",names(mf)[i])

					}
					formula=as.formula(formula)	
				}
			}else{
				stop("All covariates have only one value")	
			}
		}
		# end patch bug #1322

		if (is(formula,"formula")) {
			clafo <- "formula"
			mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
			y <- model.response(mf)
			sdy <- sd(y)
			meany <- mean(y)
			if (trait.type=="gaussian") y <- (y-meany)/sdy
			desmat <- model.matrix(formula,mf)
			rglm <- glm(y~0+desmat,family=fam)
			sglm <- summary(rglm)
			iniest <- sglm$coef[,1]
			inierr <- sglm$coef[,2]
			allids <- rownames(data)
			phids <- allids[rownames(data) %in% rownames(mf)]
			#print("bb")
			#print(phids)
			relmat <- kinship.matrix[phids,phids]*2.0
			#print("bbb")
			mids <- (rownames(data) %in% rownames(mf))
			if (trait.type=="binomial") {
				tvar <- var(rglm$resid)
			} else {
				tvar <- var(rglm$resid)
			}
		} else {
			clafo <- "NOT"
			y <- formula
			if (length(y) != dim(kinship.matrix)[1]) stop("dimension of outcome and kinship.matrix do not match")
			mids <- (!is.na(y))
			allids <- data$id
			phids <- allids[mids]
			y <- y[mids]
			relmat <- kinship.matrix[mids,mids]*2.0
			sdy <- sd(y)
			meany <- mean(y)
			if (trait.type=="gaussian") y <- (y-meany)/sdy
			desmat <- matrix(1,nrow=length(y))
			if (trait.type=="binomial") {
				tvar <- var(y)
				tmp <- mean(y)
				iniest <- c(log(tmp/(1.-tmp)))
				inierr <- sqrt(tmp*(1-tmp)/length(y))
			} else {
				tvar <- var(y)
				iniest <- c(mean(y))
				inierr <- sqrt(var(y)/length(y))
			}
		}

		#	if (!missing(data)) detach(data)
		tmp <- t(relmat)
		relmat[upper.tri(relmat)] <- tmp[upper.tri(tmp)]
		rm(tmp);gc()
		if (llfun=="polylik") eigres <- eigen(ginv(relmat),symmetric=TRUE)
		else if (llfun=="polylik_eigen") {
			if (!missing(eigenOfRel))
				eigres <- eigenOfRel
			else
				eigres <- eigen(relmat,symmetric=TRUE)
		if (any(eigres$values<1e-8)) {
			eigres$values[eigres$values<1e-8] <- 1e-8
			msg <- paste("some eigenvalues close/less than 1e-8, setting them to 1e-8\nyou can also try option llfun='polylik' instead")
			warning(msg)
		}
		} else stop("cannot be here...")
		if (!quiet) {
			cat("LM estimates of fixed parameters:\n")
			print(iniest);
			flush.console()
		}
		iniest <- iniest# + 0.001*iniest

		# check dimensions and no NA if eigenOfRel is specified
		if (!missing(eigenOfRel)) {
			if (!(is(eigenOfRel,"list") & all(names(eigenOfRel) == c("values","vectors")) )) 
				stop('eigenOfRel does not look like and eigen-generated object!')
			nIds <- length(y)
			if (dim(desmat)[1] != nIds) stop('dimensions of Y and X do not match')
			if (dim(relmat)[1] != nIds) stop('dimensions of Y and kinship matrix do not match')
			if (dim(eigres$vectors)[1] != nIds) stop('dimension 1 of eigen-vectors do not match other data')
			if (dim(eigres$vectors)[2] != nIds) stop('dimension 2 of eigen-vectors do not match other data')
			if (length(eigres$values) != nIds) stop('number of eigen-values does not match other data')
		}
		# END check dimensions and no NA if eigenOfRel is specified



		convFGLS <- NULL;

		if (!missing(fixh2)) {
			startlik <- llFUN(c(iniest,tvar),y=y,desmat=desmat,relmat=relmat,
					  eigenRes=eigres,fixh2=(fixh2/scaleh2),trait.type=trait.type,
					  scaleh2=scaleh2)
			if (opt.method=="nlm") {
				prnlev <- 0; if (!quiet) prnlev <- 2;
				h2an <- nlm(llFUN,p=c(iniest,tvar),y=y,desmat=desmat,relmat=relmat,eigenRes=eigres,
					    fixh2=(fixh2/scaleh2),trait.type=trait.type,startlik=startlik,scaleh2=scaleh2,
					    print.level=prnlev,steptol=steptol,gradtol=gradtol,...)
			} else {
				lower <- c(iniest-inierr*optimbou,1.e-4)
				upper <- c(iniest+inierr*optimbou,1)
				cntrl <- list(); if (!quiet) cntrl <- list(trace=6,REPORT=1)
				h2an <- optim(fn=llFUN,par=c(iniest,tvar),method="L-BFGS-B",lower=lower,upper=upper,
					      y=y,desmat=desmat,relmat=relmat,eigenRes=eigres,fixh2=(fixh2),trait.type=trait.type,
					      control=cntrl,scaleh2=1,...)
			}
		} else {
			nfgls <- 0
			diffgls <- maxdiffgls + 1
			parsave <- NULL	
			while (nfgls < maxnfgls & diffgls > maxdiffgls)
			{
				if (is.null(parsave)) {
					if (opt.method=="nlm") parsave <- c(iniest,starth2/scaleh2,tvar)
					else parsave <- c(iniest,starth2,tvar)
				}
				startlik<-llFUN(parsave,y=y,desmat=desmat,relmat=relmat,eigenRes=eigres,
						trait.type=trait.type,scaleh2=scaleh2)
				if (opt.method=="nlm") {
					#print(parsave)
					prnlev <- 0; if (!quiet) prnlev <- 2;
					h2an <- nlm(llFUN,p=parsave,y=y,desmat=desmat,relmat=relmat,eigenRes=eigres,
						    trait.type=trait.type,startlik=startlik,scaleh2=scaleh2,
						    print.level=prnlev,steptol=steptol,gradtol=gradtol,...)
					parsave <- h2an$estimate
				} else {
					#print(parsave)
					lower <- c(iniest-inierr*optimbou,1.e-4,1.e-4)
					upper <- c(iniest+inierr*optimbou,0.98,1)
					cntrl <- list(); if (!quiet) cntrl <- list(trace=6,REPORT=1)
					h2an <- optim(fn=llFUN,par=parsave,method="L-BFGS-B",lower=lower,upper=upper,
						      y=y,desmat=desmat,relmat=relmat,eigenRes=eigres,trait.type=trait.type,
						      control=cntrl,scaleh2=1,...)
					parsave <- h2an$par
				}


				if (fglschecks && missing(fixh2)) {
					npar <- length(parsave)
					h2 <- parsave[npar-1]*scaleh2
					if (llfun=="polylik") {
						if (!missing(eigenOfRel))
							eigres <- eigenOfRel
						else # ensure eigres contain eigen of RelMat (not Inv(RelMat))
							eigres <- eigen(relmat,symmetric=TRUE) 
					}
					es <- (eigres$value*h2+1.-h2)*parsave[npar]*sdy*sdy
					iSigma <- (eigres$vec) %*% diag(1./es,ncol=length(es)) %*% t(eigres$vec)
					# END new
					LHest <- parsave[1:(npar-2)]
					betaFGLS <- as.vector(ginv(t(desmat) %*% iSigma %*% desmat) %*% 
							      (t(desmat) %*% iSigma %*% y))
					difFGLS <- abs(betaFGLS-LHest)
					if (!quiet) {
						cat("difFGLS:\n")
						print(difFGLS)
					}
					diffgls <- max(difFGLS)
					steptol <- steptol/10;
					gradtol <- gradtol/10;
					if ((h2<1e-4 || h2>(1-1e-4)) && nfgls > floor(maxnfgls/2)) {
						parsave[npar-1] <- runif(1,min=0.05,max=0.95)/scaleh2
					}
					if (patchBasedOnFGLS && diffgls > maxdiffgls) {
						if (!quiet) {cat("fixed effect betas changed to FGLS-betas for re-estimation\n")}
						parsave[1:(npar-2)] <- betaFGLS;
					}
				} else {
					nfgls <- maxnfgls+1	
				}
				nfgls <- nfgls + 1
				#print(c("NFGLS=",nfgls))

			}

			if (diffgls<maxdiffgls) {
				convFGLS <- TRUE; 
				if (!quiet) {
					cat("\n******************************************\n");
					cat("*** GOOD convergence indicated by FGLS ***\n");
					cat("******************************************\n");
				}
			} else {
				convFGLS <- FALSE; 
				if (!quiet) {
					cat("\n***********************************************\n");
					cat("*** !!!BAD!!! convergence indicated by FGLS ***\n");
					cat("***********************************************\n");
				}
			}


		}

		if (opt.method=="optim") {
			scaleh2 <- 1
			h2an$estimate <- h2an$par; h2an$par
			h2an$minimum <- h2an$value; h2an$value
		}
		out <- list();
		npar <- length(h2an$estimate)
		if (trait.type=="gaussian") {
			if (!missing(fixh2)) {
				h2an$estimate[1:(npar-1)] <- h2an$estimate[1:(npar-1)]*sdy
				h2an$estimate[npar] <- h2an$estimate[npar]*sdy*sdy
				h2an$estimate[1] <- h2an$estimate[1]+meany
			} else {
				h2an$estimate[1:(npar-2)] <- h2an$estimate[1:(npar-2)]*sdy
				h2an$estimate[(npar-1)] <- h2an$estimate[(npar-1)]*scaleh2
				h2an$estimate[npar] <- h2an$estimate[npar]*sdy*sdy
				h2an$estimate[1] <- h2an$estimate[1]+meany
			}
		} else {
			if (!missing(fixh2)) {
			} else {
				h2an$estimate[(npar-1)] <- h2an$estimate[(npar-1)]*scaleh2
			}
		}
		out$h2an <- h2an
		if (trait.type=="gaussian") {scay <- y*sdy+meany} else {scay<-y}
		if (clafo == "formula") {
			if (!missing(fixh2)) {fxeff <- h2an$est[1:(npar-1)]} else {fxeff <- h2an$est[1:(npar-2)]}
			if (trait.type=="gaussian") {eY <- desmat %*% fxeff} else {ee <- exp(desmat %*% fxeff); eY <- ee/(1.+ee);}
			resY <- scay - eY
		} else {
			if (trait.type=="gaussian") {eY <- h2an$estimate[1]} else {ee <- exp(h2an$estimate[1]); eY <- ee/(1.+ee);}
			resY <- scay - eY
		}
		if (!missing(fixh2)) {
			h2 <- fixh2
		} else {
			h2 <- h2an$estimate[npar-1]
		}
		out$esth2 <- h2
		tvar <- h2an$estimate[npar]
		# new implementation of InvSigma
		if (fglschecks && missing(fixh2)) { 
			ginvsig <- iSigma # already computed it in FGLS checks
		} else {
			if (llfun=="polylik") {
				if (!missing(eigenOfRel))
					eigres <- eigenOfRel
				else
					eigres <- eigen(relmat,symmetric=TRUE) # ensure eigres contain eigen of RelMat (not Inv(RelMat))
			}
			es <- tvar*(eigres$value*h2+1.-h2)
			#print(es[1:5])
			#print(eigres$vec[1:5,1:5])
			#print(((diag(1./es,ncol=length(es))))[1:5,1:5])
			ginvsig <- (eigres$vec) %*% diag(1./es,ncol=length(es)) %*% t(eigres$vec)
			#print(ginvsig[1:5,1:5])
		}
		out$InvSigma <- ginvsig #ginvsig
		# END new implementation of InvSigma
		rownames(out$InvSigma) <- phids
		colnames(out$InvSigma) <- phids
		out$measuredIDs <- mids
		# compute 'environmental' residuals for GRAMMAR-2007	
		InvSigma_x_residualY <- (out$InvSigma %*% resY)
		pgresY <- as.vector((1.-h2) * tvar * InvSigma_x_residualY)
		out$pgresidualY <- rep(NA,length(mids))
		out$pgresidualY[mids] <- pgresY
		names(out$pgresidualY) <- allids #phids
		# compute residuals 	
		out$residualY <- rep(NA,length(mids))
		out$residualY[mids] <- resY
		names(out$residualY) <- allids #phids

		# compute GRAMMAR+ (G. SvISCHEVA) residuals
		eigValInvSig<- as.vector(1./es)
		zu <- mean(eigValInvSig)*tvar
		fi <- (1-(1-out$esth2)*zu)/(out$esth2 * tvar)
		# this is most expensive operation
		Bu <- eigres$vec %*% diag(sqrt(eigValInvSig)) %*% t(eigres$vec)
		# GRAMMAR+ transformed outcome
		grresY <- as.vector((1/sqrt(fi)) * (Bu %*% resY))
		out$grresidualY <- rep(NA,length(mids))
		out$grresidualY[mids] <- grresY
		names(out$grresidualY) <- allids #phids
		# GRAMMAR+ correction coefficients
		VarYG1 <- mean(pgresY^2)
		z <- 1-(zu-1)*(1-out$esth2)/out$esth2
		out$grammarGamma <- list()
		out$grammarGamma$Beta <- z*(1-out$esth2)
		out$grammarGamma$Test <- z*((1-out$esth2)^2)*tvar/VarYG1
		# END GRAMMAR+ computations

		out$call <- match.call()
		out$convFGLS <- convFGLS

		# old variant
		#	time0 <- proc.time()
		#	a <- determinant(sigma,logarithm=T)
		#	a <- a$modulus * a$sign
		#	print(proc.time()-time0)
		# end old variant
		# new variant
		a <- sum(log(es)) # logarithm of determinant of sigma as sum of logarithm of eigenvalues
		# end new variant
		b <- (t(resY) %*% InvSigma_x_residualY)
		# this is -2*lnLikelihood
		loglik <- a+b
		out$h2an$minimum <- as.vector(loglik)


		class(out) <- "polygenic"
		out
	}

