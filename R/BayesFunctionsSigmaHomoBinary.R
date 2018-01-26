BayesSSLemBinary = function(n, p, y, x, z, lambda1 = 0.1, 
                            lambda0start = 20, numBlocks = 10, w,
                            thetaA = 1, thetaB = .2*p, EMiterMax=300) {
  
  EMiterMax = max(20, EMiterMax)
  
  nScans = EMiterMax*60
  
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
  
  diffCounter = c()
  
  K = 10000
  
  Design = as.matrix(cbind(rep(1, n), z, x))
  
  accTheta = 0
  
  i = 2
  EMconverge = FALSE
  counter = 1
  
  while (counter < EMiterMax & EMconverge == FALSE) {
    
    if (i %% 1000 == 0) print(paste(i, "MCMC scans have finished"))
    
    ## latent continuous gibbs step
    Zy = rep(NA, n)
    
    meanZy = Design %*% betaPost[i-1,]
    
    Zy[y==1] = truncnorm::rtruncnorm(sum(y==1), a=0, mean = meanZy[y==1], sd=1)
    Zy[y==0] = truncnorm::rtruncnorm(sum(y==0), b=0, mean = meanZy[y==0], sd=1)
    
    D = diag(c(K,K,tauPost[i-1,]))
    Dinv = diag(1/diag(D))
    
    ## SIGMA^2
    SS = sum((Zy - Design %*% betaPost[i-1,])^2)
    sigma2Post[i] = 1
    
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
      
      tempY = Zy - (Design[,-wp] %*% betaPost[i,-wp])
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
    if (i %% 50 == 0 & i > 200) {
      wut1 = apply(gammaPost[(i-49):i, ], 1, sum)
      wut2 = sum(apply(tauPost[(i-49):i,] * (gammaPost[(i-49):i, ] == 0), 2, mean))
      
      lambda0 = sqrt(2*(p - mean(wut1)) / mean(wut2))
      diff = lambda0 - lambda0Post[i]
      lambda0Post[i] = lambda0
      diffCounter = c(diffCounter, lambda0)
      counter = counter + 1
      
      ## test if it has converged yet. Only test after 10 iterations
      if (counter > 10) {
        lD = length(diffCounter)
        mainSign = sign(diffCounter[lD] - diffCounter[1])
        if (sign(diffCounter[lD] - diffCounter[lD-5]) != mainSign) EMconverge = TRUE 
      }
    }
    i = i + 1
  }
  
  return(list(lambda0est = lambda0,
              thetaEst = thetaPost[i-1]))
}


BayesSSLBinary = function(nScans, burn, thin, 
                    n, p, y, x, z, lambda1 = 0.1, lambda0,
                    numBlocks = 10, w, thetaA = 1, thetaB = .2*p) {
  
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
  lambda0Post[1] = lambda0
  
  K = 10000
  
  Design = as.matrix(cbind(rep(1, n), z, x))
  
  accTheta = 0
  
  for (i in 2 : nScans) {
    
    if (i %% 1000 == 0) print(paste(i, "MCMC scans have finished"))
    
    ## latent continuous gibbs step
    Zy = rep(NA, n)
    
    meanZy = Design %*% betaPost[i-1,]
    
    Zy[y==1] = truncnorm::rtruncnorm(sum(y==1), a=0, mean = meanZy[y==1], sd=1)
    Zy[y==0] = truncnorm::rtruncnorm(sum(y==0), b=0, mean = meanZy[y==0], sd=1)
    
    D = diag(c(K,K,tauPost[i-1,]))
    Dinv = diag(1/diag(D))
    
    ## SIGMA^2
    SS = sum((Zy - Design %*% betaPost[i-1,])^2)
    sigma2Post[i] = 1
    
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
      
      tempY = Zy - (Design[,-wp] %*% betaPost[i,-wp])
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
           ((lambda0/sqrt(sigma2Post[i]))*(1 - thetaPost[i]^w[j])*
              exp(-(lambda0/sqrt(sigma2Post[i]))*abs(betaPost[i,j+2]))))
      gammaPost[i,j] = rbinom(1, 1, tempProb)
    }
    
    ## TAU
    RandOrder = sample(1:p, p, replace=FALSE)
    for (j in RandOrder) {
      tempLambda = gammaPost[i,j]*lambda1 + (1 - gammaPost[i,j])*lambda0
      lambdaPrime = tempLambda^2
      muPrime = sqrt(tempLambda^2 * sigma2Post[i] / betaPost[i,j+2]^2)
      tauPost[i,j] = 1 / statmod::rinvgauss(1, muPrime, lambdaPrime)
    }
  }
  
  keep = seq((burn + 1), nScans, by=thin)
  return(list(beta = betaPost[keep,],
              gamma = gammaPost[keep,]))
}




