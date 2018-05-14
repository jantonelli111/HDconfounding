#' Estimate causal effect allowing for heterogeneous treatment effects
#' 
#' This function will take in the observed data and estimate a treatment effect allowing for 
#' heterogeneous treatment effects. y,x, and z must
#' all be supplied, though all other parameters have pre-set values the user can proceed with unless
#' they wish to change the prior specification. We recommend using the EB option to estimate lambda0
#' as results can be very sensitive to this parameter and choosing it a-priori is a difficult task.
#' Note, however, that the EB option will take longer as it requires running multiple MCMCs instead
#' of just one
#'
#' @param y              The outcome to be analyzed
#' @param z              The treatment whose causal effect is to be estimated
#' @param x              An n by p matrix of potential confounders
#' @param lambda0        Either a numeric value to be used for the value of lambda0
#'                       or "EB" is specified to indicate that it will be estimated
#'                       via empirical Bayes
#' @param lambda1        A numeric value for lambda1
#' @param nScans         The number of MCMC scans to run
#' @param burn           The number of MCMC scans that will be dropped as a burn-in
#' @param thin           This number represents how many iterations between each scan 
#'                       that is kept
#' @param thetaA         The first parameter of the beta prior on the overall sparsity level
#' @param thetaB         The second parameter of the beta prior on the overall sparsity level
#' @param weight         The weight used to prioritize variables associated with treatment. This
#'                       parameter has a default of NULL, and should only be used when lambdo0
#'                       is provided instead of estimated using empirical Bayes
#' @param kMax           The maximum number of covariates to be prioritized due to association with treatment
#' @param EBiterMax      The maximum number of iterations to update EB algorithm. The algorithm is updated and
#'                       checked for convergence every 50th MCMC scan. We recommend a high value such as 300 or
#'                       500 for this parameter to ensure convergence, though the program will stop well short
#'                       of 300 in most applications if it has converged.
#'                                          
#'
#' @return A list of values that contain the treatment effect, confidence interval for the 
#'         treatment effect, and full posterior draws for the treatment effect. It also
#'         contains posterior means, confidence intervals for the regression coefficients for confounders, and
#'         the posterior mean of the binary variables indicating whether a covariate is 
#'         important for both the treated and control group models.
#'
#' @export
#' @examples
#'
#' ## p can be larger than n, but we keep the number of covariates small here
#' ## just for illustration so that the code will finish faster
#' ## n=200 and p=200 takes ~5 minutes to run for example
#' n = 200
#' p = 20
#' x = matrix(rnorm(n*p), n, p)
#' z = rbinom(n, 1, p=pnorm(0.5 + 0.7*x[,1] + 0.3*x[,2]))
#' y = rnorm(n, mean=z + 0.3*x[,1] + 0.6*x[,2] + 0.5*x[,3] + 0.5*z*x[,1], sd=1)
#' 
#' ssl = SSLhetero(y=y, z=z, x=x, nScans=3000, burn=1000, thin=2)
#' ## Output treatment effect and credible interval
#' print(ssl$TreatEffect)
#' print(ssl$TreatEffectCI)


SSLhetero = function(nScans = 20000, burn = 10000, thin = 10,
               y, x, z, lambda1 = 0.1, thetaA = 1, thetaB = 0.2*dim(x)[2],
               lambda0 = "EB", weight=NULL, kMax=20, EBiterMax = 300) {
  
  n = dim(x)[1]
  p = dim(x)[2]
  
  x = scale(x)
  
  if (length(unique(z)) == 2) {
    fit.lasso.treat <- glmnet::cv.glmnet(x=as.matrix(x), y=z, intercept=TRUE, family="binomial")
    activeX <- which(coef(fit.lasso.treat, s='lambda.1se')[-1] != 0)
    
    if (length(activeX) > kMax) {
      activeX <- order(abs(coef(fit.lasso.treat, s='lambda.1se')[-1]), 
                       decreasing=TRUE)[1 : kMax]
    } 
  }  else {
    stop("z must be binary for heterogeneous treatment effects")
  }
  
  if (lambda0 == "EB") {
    w = rep(1,p)
    
    print("Running initial empirical Bayes estimates to calculate weights for the treated group")
    
    EBresults1 = BayesSSLemHetero(p = ncol(x), y = y[z==1],
                                  x = x[z==1,], lambda1 = lambda1, lambda0start = 20,
                                  numBlocks = 10, w=w, EBiterMax = EBiterMax)
    
    
    thetaEst1 = EBresults1$thetaEst
    lambda0est1 = EBresults1$lambda0est
    wj_vec = seq(0, 1, length=2000)
    
    wj_final1 = wj_vec[which(pStar(beta=0, thetaEst1^wj_vec,
                                   lambda1=lambda1, lambda0=lambda0est1) < 0.1)[1]]
    
    
    
    if (length(activeX) == 0) {
      w= rep(1,p)
    } else {
      w = rep(1, p)
      w[activeX] = wj_final1
    }
    
    ## Now estimate lambda0 conditional on the weights
    
    print("Now estimating empirical Bayes estimates of Lambda0 conditional on weights for treated group")
    
    EBresults1.2 = BayesSSLemHetero(p = ncol(x), y = y[z==1],
                                    x = x[z==1,], lambda1 = lambda1, lambda0start = 20,
                                    numBlocks = 10, w=w, EBiterMax = EBiterMax)
    
    lambda0est1.2 = EBresults1.2$lambda0est
    
    ## Now do final analysis
    
    print("Running main analysis for treated group now")
    
    MainAnalysisBayes1 = BayesSSLHetero(nScans=nScans, burn=burn, thin=thin, p = ncol(x), y = y[z==1],
                                        x = x[z==1,], lambda1 = lambda1, lambda0 = lambda0est1.2,
                                        numBlocks = 10, w=w)
    
    
    
    ####################################################################################
    ################################ Now do it for the controls ########################
    ####################################################################################
    
    
    ## first estimate lambda0 and theta when weights are all 1
    
    print("Running initial empirical Bayes estimates to calculate weights for the control group")
    
    EBresults0 = BayesSSLemHetero(p = ncol(x), y = y[z==0],
                                  x = x[z==0,], lambda1 = lambda1, lambda0start = 20,
                                  numBlocks = 10, w=w, EBiterMax = EBiterMax)
    
    
    thetaEst0 = EBresults0$thetaEst
    lambda0est0 = EBresults0$lambda0est
    wj_vec = seq(0, 1, length=2000)
    
    wj_final0 = wj_vec[which(pStar(beta=0, thetaEst0^wj_vec,
                                   lambda1=lambda1, lambda0=lambda0est0) < 0.1)[1]]
    
    
    
    if (length(activeX) == 0) {
      w= rep(1,p)
    } else {
      w = rep(1, p)
      w[activeX] = wj_final0
    }
    
    ## Now estimate lambda0 conditional on the weights
    
    print("Now estimating empirical Bayes estimates of Lambda0 conditional on weights for control group")
    
    EBresults0.2 = BayesSSLemHetero(p = ncol(x), y = y[z==0],
                                    x = x[z==0,], lambda1 = lambda1, lambda0start = 20,
                                    numBlocks = 10, w=w, EBiterMax = EBiterMax)
    
    lambda0est0.2 = EBresults0.2$lambda0est
    
    ## Now do final analysis
    
    print("Running main analysis for control group now")
    
    MainAnalysisBayes0 = BayesSSLHetero(nScans=nScans, burn=burn, thin=thin, p = ncol(x), y = y[z==0],
                                        x = x[z==0,], lambda1 = lambda1, lambda0 = lambda0est0.2,
                                        numBlocks = 10, w=w)
    
    
    
    ####################################################################################
    ###################### Now combine them to estimate ATE ############################
    ####################################################################################
    

    Design = cbind(rep(1,n), x)
    
    atePost = rep(NA, dim(MainAnalysisBayes1$beta)[1])
    for (ni in 1 : dim(MainAnalysisBayes1$beta)[1]) {
      atePost[ni] = mean(cbind(rep(1,n), x) %*% MainAnalysisBayes1$beta[ni,] -
                           cbind(rep(1,n), x) %*% MainAnalysisBayes0$beta[ni,])
    }
    
    nBoot = 500
    ateBoot = matrix(NA, nBoot, dim(MainAnalysisBayes1$beta)[1])
    for (bi in 1 : nBoot) {
      uvec = c(0,sort(runif(n-1)),1)
      gvec = diff(uvec)
      for (ni in 1 : dim(MainAnalysisBayes1$beta)[1]) {
        ateBoot[bi,ni] = weighted.mean(cbind(rep(1,n), x) %*% MainAnalysisBayes1$beta[ni,] -
                                cbind(rep(1,n), x) %*% MainAnalysisBayes0$beta[ni,], w=gvec)
      }
    }
        

  } else {
    
    if (is.null(weight)) {
      stop("Weight must be provided if empirical Bayes is not used")
    } else if (weight > 1 | weight <= 0) {
      stop("Weight must be between 0 and 1")
    } else {
      w = rep(1, p)
      
      if (length(activeX) == 0) {
        w= rep(1,p)
      } else {
        w = rep(1, p)
        w[activeX] = weight
      }
      
      print("Running main analysis for treated group now")
      
      MainAnalysisBayes1 = BayesSSLHetero(nScans=nScans, burn=burn, thin=thin, p = ncol(x), y = y[z==1],
                                          x = x[z==1,], lambda1 = lambda1, lambda0 = lambda0,
                                          numBlocks = 10, w=w)
      
      
      
      ####################################################################################
      ################################ Now do it for the controls ########################
      ####################################################################################
      
      print("Running main analysis for control group now")
      
      MainAnalysisBayes0 = BayesSSLHetero(nScans=nScans, burn=burn, thin=thin, p = ncol(x), y = y[z==0],
                                          x = x[z==0,], lambda1 = lambda1, lambda0 = lambda0,
                                          numBlocks = 10, w=w)
      
      
      
      ####################################################################################
      ###################### Now combine them to estimate ATE ############################
      ####################################################################################
      
      
      Design = cbind(rep(1,n), x)
      
      atePost = rep(NA, dim(MainAnalysisBayes1$beta)[1])
      for (ni in 1 : dim(MainAnalysisBayes1$beta)[1]) {
        atePost[ni] = mean(cbind(rep(1,n), x) %*% MainAnalysisBayes1$beta[ni,] -
                             cbind(rep(1,n), x) %*% MainAnalysisBayes0$beta[ni,])
      }
      
      nBoot = 500
      ateBoot = matrix(NA, nBoot, dim(MainAnalysisBayes1$beta)[1])
      for (bi in 1 : nBoot) {
        uvec = c(0,sort(runif(n-1)),1)
        gvec = diff(uvec)
        for (ni in 1 : dim(MainAnalysisBayes1$beta)[1]) {
          ateBoot[bi,ni] = weighted.mean(cbind(rep(1,n), x) %*% MainAnalysisBayes1$beta[ni,] -
                                  cbind(rep(1,n), x) %*% MainAnalysisBayes0$beta[ni,], w=gvec)
        }
      }
      
      
    }
  }
  
  l = list(TreatEffect = mean(atePost),
           TreatEffectCI = quantile(ateBoot, c(.025, .975)),
           betaPostMean1 = apply(MainAnalysisBayes1$beta[,2:(p+1)], 2, mean),
           betaPostCI1 = apply(MainAnalysisBayes1$beta[,2:(p+1)], 2, quantile, c(.025, .975)),
           gammaPostMean1 = apply(MainAnalysisBayes1$gamma, 2, mean),
           betaPostMean0 = apply(MainAnalysisBayes0$beta[,2:(p+1)], 2, mean),
           betaPostCI0 = apply(MainAnalysisBayes0$beta[,2:(p+1)], 2, quantile, c(.025, .975)),
           gammaPostMean0 = apply(MainAnalysisBayes0$gamma, 2, mean))
  
  return(l)
}




#' Estimate causal effect assuming homogeneous treatment effect
#'
#' This function will take in the observed data and estimate a treatment effect assuming a 
#' homogeneous treatment effect. y,x, and z must
#' all be supplied, though all other parameters have pre-set values the user can proceed with unless
#' they wish to change the prior specification. We recommend using the EB option to estimate lambda0
#' as results can be very sensitive to this parameter and choosing it a-priori is a difficult task.
#' Note, however, that the EB option will take longer as it requires running multiple MCMCs instead
#' of just one
#'
#' @param y              The outcome to be analyzed
#' @param z              The treatment whose causal effect is to be estimated
#' @param x              An n by p matrix of potential confounders
#' @param z_type         String indicating whether treatment is binary or continuous. The default
#'                       value is z_type="binary", though it should be set to z_type = "continuous"
#'                       if the treatment is continuous
#' @param lambda0        Either a numeric value to be used for the value of lambda0
#'                       or "EB" is specified to indicate that it will be estimated
#'                       via empirical Bayes
#' @param lambda1        A numeric value for lambda1
#' @param nScans         The number of MCMC scans to run
#' @param burn           The number of MCMC scans that will be dropped as a burn-in
#' @param thin           This number represents how many iterations between each scan 
#'                       that is kept
#' @param thetaA         The first parameter of the beta prior on the overall sparsity level
#' @param thetaB         The second parameter of the beta prior on the overall sparsity level
#' @param weight         The weight used to prioritize variables associated with treatment. This
#'                       parameter has a default of NULL, and should only be used when lambdo0
#'                       is provided instead of estimated using empirical Bayes
#' @param kMax           The maximum number of covariates to be prioritized due to association with treatment
#' @param EBiterMax      The maximum number of iterations to update EB algorithm. The algorithm is updated and
#'                       checked for convergence every 50th MCMC scan. We recommend a high value such as 300 or
#'                       500 for this parameter to ensure convergence, though the program will stop well short
#'                       of 300 in most applications if it has converged.
#'                    
#'
#' @return A list of values that contain the treatment effect, confidence interval for the 
#'         treatment effect, full posterior draws for the treatment effect, posterior means
#'         and confidence intervals for the regression coefficients for confounders, and
#'         the posterior mean of the binary variables indicating whether a covariate is 
#'         important.
#'
#' @export
#' @examples
#'
#' ## p can be larger than n, but we keep the number of covariates small here
#' ## just for illustration so that the code will finish faster
#' ## n=200 and p=200 takes ~5 minutes to run for example
#' n = 200
#' p = 20
#' x = matrix(rnorm(n*p), n, p)
#' z = rbinom(n, 1, p=pnorm(0.7*x[,1] + 0.3*x[,2]))
#' y = rnorm(n, mean=z + 0.3*x[,1] + 0.6*x[,2] + 0.5*x[,3], sd=1)
#' 
#' ssl = SSL(y=y, z=z, x=x, nScans=3000, burn=1000, thin=2)
#' ## Output treatment effect and credible interval
#' print(ssl$TreatEffect)
#' print(ssl$TreatEffectCI)
#' 
#' ## Print the posterior inclusion probabilities for confounders
#' print(ssl$gammaPostMean)

SSL = function(nScans = 20000, burn = 10000, thin = 10,
               y, x, z, z_type = "binary", lambda1 = 0.1, thetaA = 1, thetaB = 0.2*dim(x)[2],
               lambda0 = "EB", weight=NULL, kMax=20, EBiterMax=300) {
  
  n = dim(x)[1]
  p = dim(x)[2]
  
  x = scale(x)
  
  if (z_type == "binary") {
    fit.lasso.treat <- glmnet::cv.glmnet(x=as.matrix(x), y=z, intercept=TRUE, family="binomial")
    activeX <- which(coef(fit.lasso.treat, s='lambda.1se')[-1] != 0)
    
    if (length(activeX) > kMax) {
      activeX <- order(abs(coef(fit.lasso.treat, s='lambda.1se')[-1]), 
                       decreasing=TRUE)[1 : kMax]
    } 
  } else if (z_type == "continuous"){
    fit.lasso.treat <- glmnet::cv.glmnet(x=as.matrix(x), y=z, intercept=TRUE)
    activeX <- which(coef(fit.lasso.treat, s='lambda.1se')[-1] != 0)
    
    if (length(activeX) > kMax) {
      activeX <- order(abs(coef(fit.lasso.treat, s='lambda.1se')[-1]), 
                       decreasing=TRUE)[1 : kMax]
    }
  } else {
    stop("z_type must be set to binary or continuous")
  }
  
  if (lambda0 == "EB") {
    w = rep(1, p)
    
    print("Running initial empirical Bayes estimates to calculate weights")
    
    EBresults = BayesSSLem(n=n, p = ncol(x), y = y,
                           x = x, z=z, lambda1 = lambda1, lambda0start = 10,
                           numBlocks = 10, w=w, EBiterMax = EBiterMax)
    
    thetaEst = EBresults$thetaEst
    lambda0est = EBresults$lambda0est
    wj_vec = seq(0, 1, length=2000)
    
    print(lambda0est)
    
    wj_final = wj_vec[which(pStar(beta=0, thetaEst^wj_vec, 
                                  lambda1=lambda1, lambda0=lambda0est) < 0.1)[1]]
    
    
    if (length(activeX) == 0) {
      w= rep(1,p)
    } else {
      w = rep(1, p)
      w[activeX] = wj_final
    }
    
    ## Now estimate lambda0 conditional on the weights
    
    print("Now estimating empirical Bayes estimates of Lambda0 conditional on weights")
    
    EBresults2 = BayesSSLem(n=n, p = ncol(x), y = y,
                            x = x, z=z, lambda1 = lambda1, lambda0start = 10,
                            numBlocks = 10, w=w, EBiterMax = EBiterMax)
    
    lambda0est = EBresults2$lambda0est
    print(lambda0est)
    
    ## Now do final analysis
    
    print("Running final analysis now")
    
    MainAnalysisBayes = BayesSSL(nScans=nScans, burn=burn, thin=thin, n=n, p = ncol(x), y = y, 
                                 x = x, z=z, lambda1 = lambda1, lambda0 = lambda0est,
                                 numBlocks = 10, w=w)
  } else {
    
    if (is.null(weight)) {
      stop("Weight must be provided if empirical Bayes is not used")
    } else if (weight > 1 | weight <= 0) {
      stop("Weight must be between 0 and 1")
    } else {
      w = rep(1, p)
      
      if (length(activeX) == 0) {
        w= rep(1,p)
      } else {
        w = rep(1, p)
        w[activeX] = weight
      }
      
      ## Now do final analysis
      
      print("Running final analysis now")
      
      MainAnalysisBayes = BayesSSL(nScans=nScans, burn=burn, thin=thin, n=n, p = ncol(x), y = y, 
                                   x = x, z=z, lambda1 = lambda1, lambda0 = lambda0,
                                   numBlocks = 10, w=w)
    }
  }
  
  l = list(TreatEffect = mean(MainAnalysisBayes$beta[,2]),
           TreatEffectCI = quantile(MainAnalysisBayes$beta[,2], c(.025, .975)),
           TreatEffectPost = MainAnalysisBayes$beta[,2],
           betaPostMean = apply(MainAnalysisBayes$beta[,3:(p+2)], 2, mean),
           betaPostCI = apply(MainAnalysisBayes$beta[,3:(p+2)], 2, quantile, c(.025, .975)),
           gammaPostMean = apply(MainAnalysisBayes$gamma, 2, mean))
  
  return(l)
}
