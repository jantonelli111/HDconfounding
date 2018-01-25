#' Main function for analyzing the effect of an exposure on a continuous outcome
#'
#' This function will take in the observed data and estimate a treatment effect. y,x, and z must
#' all be supplied, though all other parameters have pre-set values the user can proceed with unless
#' they wish to change the prior specification. We recommend using the EM option to estimate lambda0
#' as results can be very sensitive to this parameter and choosing it a-priori is a difficult task.
#'
#' @param y              The outcome to be analyzed
#' @param z              The treatment whose causal effect is to be estimated
#' @param x              An n by p matrix of potential confounders
#' @param lambda0        Either a numeric value to be used for the value of lambda0
#'                       or "EM" is specified to indicate that it will be estimated
#'                       via empirical Bayes
#' @param lambda1        A numeric value for lambda1
#' @param nScans         The number of MCMC scans to run
#' @param burn           The number of MCMC scans that will be dropped as a burn-in
#' @param thin           This number represents how many iterations between each scan 
#'                       that is kept
#' @param thetaA         The first parameter of the beta prior on the overall sparsity level
#' @param thetaB         The second parameter of the beta prior on the overall sparsity level
#'                    
#'
#' @return An list of values that contain the treatment effect, confidence interval for the 
#'         treatment effect, as well as the corresponding quantities for the regression parameters
#'         for the confounders.
#'
#' @export
#' @examples
#'
#' put example run here

SSL = function(nScans = 20000, burn = 10000, thin = 10,
               y, x, z, lambda1 = 0.1, thetaA = 1, thetaB = 0.2*dim(x)[2],
               lambda0 = "EM") {

  n = dim(x)[1]
  p = dim(x)[2]
  
  fit.lasso.treat <- glmnet::cv.glmnet(x=as.matrix(x), y=z, intercept=TRUE, family="binomial")
  activeX <- which(coef(fit.lasso.treat, s='lambda.1se')[-1] != 0)

  if (lambda0 == "EM") {
    w = rep(1, p)
    
    print("Running initial empirical Bayes estimates to calculate weights")
    
    EMresults = BayesSSLem(nScans=nScans, burn=burn, thin=thin, n=n, p = ncol(x), y = y,
                           x = x, z=z, lambda1 = 0.1, lambda0start = 8,
                           numBlocks = 10, w=w)
    
    thetaEst = EMresults$thetaEst
    lambda0est = EMresults$lambda0est
    wj_vec = seq(0, 1, length=2000)
    
    wj_final = wj_vec[which(pStar(beta=0, thetaEst^wj_vec, 
                                  lambda1=0.1, lambda0=lambda0est) < 0.1)[1]]
    
    
    if (length(activeX) == 0) {
      w= rep(1,p)
    } else {
      w = rep(1, p)
      w[activeX] = wj_final
    }
    
    ## Now estimate lambda0 conditional on the weights
    
    print("Now estimating empirical Bayes estimates of Lambda0 conditional on weights")
    
    EMresults2 = BayesSSLem(nScans=nScans, burn=burn, thin=thin, n=n, p = ncol(x), y = y,
                            x = x, z=z, lambda1 = 0.1, lambda0start = 8,
                            numBlocks = 10, w=w)
    
    lambda0est = EMresults2$lambda0est
    
    ## Now do final analysis
    
    print("Running final analysis now")
    
    MainAnalysisBayes = BayesSSL(nScans=nScans, burn=burn, thin=thin, n=n, p = ncol(x), y = y, 
                                 x = x, z=z, lambda1 = 0.1, lambda0 = lambda0est,
                                 numBlocks = 10, w=w)
  } else {
    w = rep(1, p)
    lambda0est = lambda0
    wj_vec = seq(0, 1, length=2000)
    
    wj_final = wj_vec[which(pStar(beta=0, thetaEst^wj_vec, 
                                  lambda1=0.1, lambda0=lambda0est) < 0.1)[1]]
    
    wjSave[r] = wj_final
    
    
    if (length(activeX) == 0) {
      w= rep(1,p)
    } else {
      w = rep(1, p)
      w[activeX] = wj_final
    }
    
    ## Now do final analysis
    
    print("Running final analysis now")
    
    MainAnalysisBayes = BayesSSL(nScans=nScans, burn=burn, thin=thin, n=n, p = ncol(x), y = y, 
                                 x = x, z=z, lambda1 = 0.1, lambda0 = lambda0est,
                                 numBlocks = 10, w=w)
  }

}

BayesSSLem = function(nScans = 20000, burn = 10000, thin = 10,
                      n, p, y, x, z, lambda1 = 0.1,
                      lambda0start = 20, numBlocks = 10, w,
                      thetaA = 1, thetaB = .2*p) {

  betaPost = matrix(NA, nScans, p+2)
  gammaPost = matrix(NA, nScans, p)
  sigma2Post = rep(NA, nScans)
  thetaPost = rep(NA, nScans)
  tauPost = matrix(NA, nScans, p)
  lambda0Post = rep(NA, nScans)

  betaPost[1,] = rnorm(p+2, sd=0.2)
  gammaPost[1,] = rep(0,p)
  sigma2Post[1] = 1
  thetaPost[1] = 0.02
  tauPost[1,] = rgamma(p, 2)

  sigmaA = 0.001
  sigmaB = 0.001
  lambda0 = lambda0start
  lambda0Post[1] = lambda0

  K = 10000

  Design = as.matrix(cbind(rep(1, n), z, x))

  accTheta = 0

  for (i in 2 : nScans) {

    if (i %% 100 == 0) print(i)

    D = diag(c(K,K,tauPost[i-1,]))
    Dinv = diag(1/diag(D))

    ## SIGMA^2
    SS = sum((y - Design %*% betaPost[i-1,])^2)
    sigma2Post[i] = 1/rgamma(1, (n-1)/2 + p/2,
                             SS/2 + t(betaPost[i-1,]) %*% Dinv %*% betaPost[i-1,]/2)

    ## THETA using MH update
    BoundaryLow2 = max(0, thetaPost[i-1] - 0.02)
    BoundaryAbove2 = min(1, thetaPost[i-1] + 0.02)
    thetaNew = runif(1, BoundaryLow2, BoundaryAbove2)

    BoundaryLow1 = max(0, thetaNew - 0.02)
    BoundaryAbove1 = min(1, thetaNew + 0.02)

    logAR = LogTheta(thetaNew, a=thetaA, b=thetaB, w=w, gamma=gammaPost[i-1,]) +
      dunif(thetaPost[i-1], BoundaryLow1, BoundaryAbove1, log=TRUE) -
      LogTheta(thetaPost[i-1], a=thetaA, b=thetaB, w=w, gamma=gammaPost[i-1,]) -
      dunif(thetaNew, BoundaryLow2, BoundaryAbove2, log=TRUE)

    if (logAR > log(runif(1))) {
      thetaPost[i] = thetaNew
      accTheta = accTheta + 1
    } else {
      thetaPost[i] = thetaPost[i-1]
    }

    ## BETA
    jumps = floor((p+2) / numBlocks)
    betaPost[i,] = betaPost[i-1,]
    for (jj in 1 : numBlocks) {
      wp = ((jj-1)*jumps + 1) : (jj*jumps)
      if (jj == numBlocks) {
        wp = ((jj-1)*jumps + 1) : (p+2)
      }

      tempY = y - (Design[,-wp] %*% betaPost[i,-wp])
      tempV = sigma2Post[i]*solve((t(Design[,wp]) %*% Design[,wp]) + Dinv[wp,wp])
      tempMU = (1/sigma2Post[i])*tempV %*% t(Design[,wp]) %*% tempY
      betaPost[i,wp] = mvtnorm::rmvnorm(1, mean=tempMU, sigma=tempV)
    }

    ## GAMMA
    RandOrder = sample(1:p, p, replace=FALSE)
    for (j in RandOrder) {
      tempProb = ((lambda1/sqrt(sigma2Post[i]))*(thetaPost[i]^w[j])*
                    exp(-(lambda1/sqrt(sigma2Post[i]))*abs(betaPost[i,j+2]))) /
        (((lambda1/sqrt(sigma2Post[i]))*(thetaPost[i]^w[j])*
            exp(-(lambda1/sqrt(sigma2Post[i]))*abs(betaPost[i,j+2]))) +
           ((lambda0Post[i-1]/sqrt(sigma2Post[i]))*(1 - thetaPost[i]^w[j])*
              exp(-(lambda0Post[i-1]/sqrt(sigma2Post[i]))*abs(betaPost[i,j+2]))))
      gammaPost[i,j] = rbinom(1, 1, tempProb)
    }

    ## TAU
    RandOrder = sample(1:p, p, replace=FALSE)
    for (j in RandOrder) {
      tempLambda = gammaPost[i,j]*lambda1 + (1 - gammaPost[i,j])*lambda0Post[i-1]
      lambdaPrime = tempLambda^2
      muPrime = sqrt(tempLambda^2 * sigma2Post[i] / betaPost[i,j+2]^2)
      tauPost[i,j] = 1 / statmod::rinvgauss(1, muPrime, lambdaPrime)
    }

    ## LAMBDA0
    lambda0Post[i] =  lambda0
    if (i %% 50 == 0 & i > 500) {
      wut1 = apply(gammaPost[(i-49):i, ], 1, sum)
      wut2 = sum(apply(tauPost[(i-49):i,] * (gammaPost[(i-49):i, ] == 0), 2, mean))

      lambda0 = sqrt(2*(p - mean(wut1)) / mean(wut2))
      diff = lambda0 - lambda0Post[i]
      print(c(lambda0, diff))
      lambda0Post[i] = lambda0

    }
  }

  keep = seq((burn + 1), nScans, by=thin)
  return(list(lambda0est = mean(lambda0Post[keep], na.rm=TRUE),
              thetaEst = mean(thetaPost[keep])))
}



expit = function(x) {
  return(exp(x) / (1 + exp(x)))
}

pStar = function(beta, theta, lambda1, lambda0) {
  part1 = theta*lambda1*exp(-lambda1*abs(beta))
  part0 = (1-theta)*lambda0*exp(-lambda0*abs(beta))
  return(part1 / (part1 + part0))
}

LogTheta = function(theta, a, b, w, gamma) {
  part1 = (a-1)*log(theta)
  part2 = (b-1)*log(1 - theta)
  part3 = sum(w*gamma)*log(theta)
  part4 = sum(log((1 - theta^w)^(1-gamma)))
  return(part1 + part2 + part3 + part4)
}
