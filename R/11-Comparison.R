package.install = function(pack) {
  local({r <- getOption("repos");r["CRAN"] <- "http://cran.r-project.org"; options(repos=r)})
  
  # name of package to install / load
  pack = pack
  
  if (pack %in% rownames(installed.packages())) {
    library(pack, character.only=T)
  } else {
    if (pack %in% rownames(installed.packages(lib.loc='~/Rpackages'))) {
      library(pack, lib.loc='~/Rpackages', character.only=T)
    } else {
      install.packages(pack, lib='~/Rpackages')
      library(pack, lib.loc='~/Rpackages', character.only=T)
    }
  }
}

expit = function(x) {
  return(exp(x) / (1 + exp(x)))
}

pStar = function(beta, theta, lambda1, lambda0) {
  part1 = theta*lambda1*exp(-lambda1*abs(beta))
  part0 = (1-theta)*lambda0*exp(-lambda0*abs(beta))
  return(part1 / (part1 + part0))
}

package.install("mvtnorm")
package.install("statmod")
package.install("glmnet")
package.install("balanceHD")

index <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(index)
R = 5

source("~/SSLconfounding/Code5/BayesHomo/BayesFunctionsSigmaHomo.R")
source("~/SSLconfounding/Code5/BayesHetero/BayesFunctionsSigmaHetero.R")

n <- 200
p <- 500

StrongConfoundersX = StrongConfoundersY = c(1, -1)

WeakConfoundersX = StrongConfoundersX
WeakConfoundersY = c(0.3, -0.3)

InstrumentsX = StrongConfoundersX
InstrumentsY = rep(0,2)

PredictorsX = rep(0,2)
PredictorsY = StrongConfoundersY

beta.c <- c(StrongConfoundersY, WeakConfoundersY,
            InstrumentsY, PredictorsY, rnorm(p-8, sd=0.1))

gamma <- c(StrongConfoundersX, WeakConfoundersX,
           InstrumentsX, PredictorsX,rep(0, p-8))

beta <- 1

sigma = matrix(0.6, p,p)
diag(sigma) = 1

BetaSave = array(NA, dim=c(R, 13))
IntervalSave = array(NA, dim=c(R, 13))
IntervalWidthSave = array(NA, dim=c(R, 13))
wjSave = rep(NA, R)
StandardErrorSave = array(NA, dim=c(R,13))

for (r in 1 : R) {
  
  x = rmvnorm(n, sigma=sigma)
  z <- 1*(0 + x %*% gamma + rlogis(n) > 0)
  y <- 200 + z + x %*% beta.c + rnorm(n, sd=1)
  
  fit.lasso.treat <- cv.glmnet(x=as.matrix(x), y=z, intercept=TRUE, family="binomial")
  activeX <- which(coef(fit.lasso.treat, s='lambda.1se')[-1] != 0)
  
  # Need a new design
  xtilde = cbind(x, x*z[1:n])
  
  ## What does the usual LASSO give you
  fit.lasso <- cv.glmnet(cbind(z, x), y, intercept = TRUE, penalty.factor = rep(c(0,1), c(1, p)))
  active <- which(coef(fit.lasso)[-c(1,2)] != 0)
  BetaSave[r,1] = coef(fit.lasso)[2]
  
  ## post selection lasso
  ModelLasso = lm(y ~ z + x[,active])
  BetaSave[r,2] = ModelLasso$coef[2]
  LassoInterval = c(BetaSave[r,2] - 1.96*sqrt(vcov(ModelLasso)[2,2]),
                    BetaSave[r,2] + 1.96*sqrt(vcov(ModelLasso)[2,2]))
  IntervalSave[r,2] = as.numeric(beta < LassoInterval[2] & beta > LassoInterval[1])
  IntervalWidthSave[r,2] = LassoInterval[2] - LassoInterval[1]
  StandardErrorSave[r,2] = sqrt(vcov(ModelLasso)[2,2])
  
  ## What does the usual LASSO give you heterogeneous
  fit.lasso <- cv.glmnet(cbind(z, xtilde), y, intercept = TRUE, penalty.factor = rep(c(0,1), c(1, 2*p)))
  active <- which(coef(fit.lasso) != 0)
  BetaSave[r,3] = mean(cbind(rep(1,n), rep(1,n), x, x)[,active] %*% coef(fit.lasso)[active] - 
                         cbind(rep(1,n), rep(0,n), x, matrix(0,n,p))[,active] %*% coef(fit.lasso)[active])
  
  
  ## post selection lasso heterogeneous
  active <- which(coef(fit.lasso)[-c(1,2)] != 0)
  ModelLasso = lm(y ~ z + xtilde[,active])
  BetaSave[r,4] = mean(cbind(rep(1,n), rep(1,n), cbind(x,x)[,active]) %*% coef(ModelLasso) - 
                         cbind(rep(1,n), rep(0,n), 
                               cbind(x,matrix(0,n,p))[,active]) %*% coef(ModelLasso))
  
  ## Double selection method with 1se CV
  
  fit.lasso <- cv.glmnet(cbind(z, x), y, intercept = TRUE, penalty.factor = rep(c(0,1), c(1, p)))
  fit.lasso.treat <- cv.glmnet(x=as.matrix(x), y=z, intercept=TRUE, family="binomial")
  
  activeX <- which(coef(fit.lasso.treat, s='lambda.1se')[-1] != 0)
  activeY <- which(coef(fit.lasso, s='lambda.1se')[-c(1,2)] != 0)
  DoubleActive = union(activeX, activeY)
  ModelDouble = lm(y ~ z + x[,DoubleActive])
  BetaSave[r,5] = ModelDouble$coef[2]
  DoubleInterval = c(BetaSave[r,5] - 1.96*sqrt(vcov(ModelDouble)[2,2]),
                     BetaSave[r,5] + 1.96*sqrt(vcov(ModelDouble)[2,2]))
  IntervalSave[r,5] = as.numeric(beta < DoubleInterval[2] & beta > DoubleInterval[1])
  IntervalWidthSave[r,5] = DoubleInterval[2] - DoubleInterval[1]
  StandardErrorSave[r,5] = sqrt(vcov(ModelDouble)[2,2])
  
  ## Double selection method with 1se CV heterogeneous
  
  fit.lasso <- cv.glmnet(cbind(z, xtilde), y, intercept = TRUE, penalty.factor = rep(c(0,1), c(1, 2*p)))
  
  activeX <- which(coef(fit.lasso.treat, s='lambda.1se')[-1] != 0)
  activeY <- which(coef(fit.lasso, s='lambda.1se')[-c(1,2)] != 0)
  DoubleActive = union(activeX, activeY)
  ModelDoubleY = lm(y ~ z + xtilde[,DoubleActive])
  ModelDoubleX = glm(z ~ x[,activeX], family=binomial)
  
  h1 = cbind(rep(1,n), rep(0,n), cbind(x,matrix(0,n,p))[,DoubleActive]) %*% ModelDoubleY$coefficients
  h2 = cbind(rep(1,n), rep(1,n), cbind(x,x)[,DoubleActive]) %*% ModelDoubleY$coefficients
  h3 = expit(cbind(rep(1,n), x[,activeX]) %*% ModelDoubleX$coefficients)
  
  BetaSave[r,6] = mean(h2 - h1 + (z*(y - h2)/h3) - ((1-z)*(y - h1)/(1-h3)))
  varDouble = mean((BetaSave[r,6] - (h2 - h1 + (z*(y - h2)/h3) - ((1-z)*(y - h1)/(1-h3))))^2)/n
  seDouble = sqrt(varDouble)
  
  DoubleInterval = c(BetaSave[r,6] - 1.96*seDouble,
                     BetaSave[r,6] + 1.96*seDouble)
  IntervalSave[r,6] = as.numeric(beta < DoubleInterval[2] & beta > DoubleInterval[1])
  IntervalWidthSave[r,6] = DoubleInterval[2] - DoubleInterval[1]
  StandardErrorSave[r,6] = seDouble
  
  ## Approximate residual debiasing
  
  fit.lasso <- cv.glmnet(cbind(z, x), y, intercept = TRUE, penalty.factor = rep(c(0,1), c(1, p)))
  
  
  tau.hat = residualBalance.ate(x, y, z, estimate.se = TRUE)
  BetaSave[r,7] = tau.hat[1]
  IntervalSave[r,7] = (tau.hat[1] + 1.96*tau.hat[2] > beta &
                   tau.hat[1] - 1.96*tau.hat[2] < beta)
  IntervalWidthSave[r,7] = 2*1.96*tau.hat[2]
  StandardErrorSave[r,7] = tau.hat[2]
  
  ## Farrell (2015) estimator
  
  fit.x <- cv.glmnet(x, z, intercept = TRUE, family = 'binomial')
  fit.y2 <- cv.glmnet(cbind(z,x), y, intercept = TRUE, penalty.factor = rep(c(0,1), c(1, p)))
  
  active.x <- which(coef(fit.x)[-1] != 0)
  active.y <- which(coef(fit.y2)[-c(1,2)] != 0)
  
  if (length(active.x) == 0) {
    ps.mod.farrell = glm(z ~ 1, family='binomial')
  } else {
    ps.mod.farrell = glm(z ~ x[,active.x], family='binomial')
  }
  
  if (length(active.y) == 0) {
    out.mod.farrell = lm(y ~ z)
  } else {
    out.mod.farrell = lm(y ~ z + x[,active.y]) 
  }
  
  ps1.farrell = ps.mod.farrell$fitted.values
  ps0.farrell = 1 - ps1.farrell
  
  out1.farrell = cbind(rep(1,n), rep(1,n), x[,active.y]) %*% out.mod.farrell$coefficients
  out0.farrell = cbind(rep(1,n), rep(0,n), x[,active.y]) %*% out.mod.farrell$coefficients
  
  BetaSave[r,8] = sum(z*(y - out1.farrell)/ps1.farrell + out1.farrell)/n -
    sum((1-z)*(y - out0.farrell)/ps0.farrell + out0.farrell)/n
  
  part1 = mean((z*(y - out1.farrell)^2)/(ps1.farrell^2))
  part2 = mean((out1.farrell - mean(y[z==1]))^2)
  var1 = part1 + part2
  
  part2 = mean((out1.farrell - mean(y[z==1])) * (out0.farrell - mean(y[z==0])))
  cov01 = part2
  
  part1 = mean(((1-z)*(y - out0.farrell)^2)/(ps0.farrell^2))
  part2 = mean((out0.farrell - mean(y[z==0]))^2)
  var0 = part1 + part2
  
  CovMat = matrix(c(var0, cov01, cov01, var1), nrow=2)/n
  
  totalVar = c(-1,1) %*% CovMat %*% c(-1,1)
  totalSE = sqrt(totalVar)
  
  IntervalSave[r,8] = (BetaSave[r,8] + 1.96*totalSE > beta &
                         BetaSave[r,8] - 1.96*totalSE < beta)
  IntervalWidthSave[r,8] = 2*1.96*totalSE
  StandardErrorSave[r,8] = totalSE
  
  
  ## Farrell (2015) estimator heterogeneous
  
  fit.x <- cv.glmnet(x, z, intercept = TRUE, family = 'binomial')
  fit.y0 <- cv.glmnet(x[z==0,], y[z==0], intercept = TRUE, penalty.factor = rep(1, p))
  fit.y1 <- cv.glmnet(x[z==1,], y[z==1], intercept = TRUE, penalty.factor = rep(1, p))
  
  active.x <- which(coef(fit.x)[-1] != 0)
  active.y0 <- which(coef(fit.y0)[-c(1)] != 0)
  active.y1 <- which(coef(fit.y1)[-c(1)] != 0)
  
  if (length(active.x) == 0) {
    ps.mod.farrell = glm(z ~ 1, family='binomial')
  } else {
    ps.mod.farrell = glm(z ~ x[,active.x], family='binomial')
  }
  
  if (length(active.y0) == 0) {
    out0.mod.farrell = lm(y[z==0] ~ 1)
  } else {
    out0.mod.farrell = lm(y[z==0] ~ x[z==0,active.y0]) 
  }
  
  if (length(active.y1) == 0) {
    out1.mod.farrell = lm(y[z==1] ~ 1)
  } else {
    out1.mod.farrell = lm(y[z==1] ~ x[z==1,active.y1]) 
  }
  
  ps1.farrell = ps.mod.farrell$fitted.values
  ps0.farrell = 1 - ps1.farrell
  
  out1.farrell = cbind(rep(1,n), x[,active.y1]) %*% out1.mod.farrell$coefficients
  out0.farrell = cbind(rep(1,n), x[,active.y0]) %*% out0.mod.farrell$coefficients
  
  BetaSave[r,9] = sum(z*(y - out1.farrell)/ps1.farrell + out1.farrell)/n -
    sum((1-z)*(y - out0.farrell)/ps0.farrell + out0.farrell)/n
  
  part1 = mean((z*(y - out1.farrell)^2)/(ps1.farrell^2))
  part2 = mean((out1.farrell - mean(y[z==1]))^2)
  var1 = part1 + part2
  
  part2 = mean((out1.farrell - mean(y[z==1])) * (out0.farrell - mean(y[z==0])))
  cov01 = part2
  
  part1 = mean(((1-z)*(y - out0.farrell)^2)/(ps0.farrell^2))
  part2 = mean((out0.farrell - mean(y[z==0]))^2)
  var0 = part1 + part2
  
  CovMat = matrix(c(var0, cov01, cov01, var1), nrow=2)/n
  
  totalVar = c(-1,1) %*% CovMat %*% c(-1,1)
  totalSE = sqrt(totalVar)
  
  IntervalSave[r,9] = (BetaSave[r,9] + 1.96*totalSE > beta &
                         BetaSave[r,9] - 1.96*totalSE < beta)
  IntervalWidthSave[r,9] = 2*1.96*totalSE
  StandardErrorSave[r,9] = totalSE
  
  
  
  ## EM-SSL approach
  
  w = rep(1, p)
  
  EMresults = BayesSSLem(nScans=30000, burn=25000, thin=5, n=n, p = ncol(x), y = y,
                         x = x, z=z, lambda1 = 0.1, lambda0start = 8,
                         numBlocks = 10, w=w)
  
  
  thetaEst = EMresults$thetaEst
  lambda0est = EMresults$lambda0est
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
  
  ## Now estimate lambda0 conditional on the weights
  
  EMresults2 = BayesSSLem(nScans=30000, burn=25000, thin=5, n=n, p = ncol(x), y = y,
                          x = x, z=z, lambda1 = 0.1, lambda0start = 8,
                          numBlocks = 10, w=w)
  
  lambda0est = EMresults2$lambda0est
  
  ## Now do final analysis
  
  MainAnalysisBayes = BayesSSL(nScans=10000, burn=2000, thin=8, n=n, p = ncol(x), y = y, 
                               x = x, z=z, lambda1 = 0.1, lambda0 = lambda0est,
                               numBlocks = 10, w=w)
  
  BetaSave[r,10] = mean(MainAnalysisBayes[,2], na.rm=TRUE)
  intervalBETA = quantile(MainAnalysisBayes[,2], c(.025, .975), na.rm=TRUE)
  IntervalSave[r,10] = as.numeric(beta < intervalBETA[2] & beta > intervalBETA[1])
  IntervalWidthSave[r,10] = intervalBETA[2] - intervalBETA[1]
  StandardErrorSave[r,10] = sd(MainAnalysisBayes[,2], na.rm=TRUE)

  ## EM-SSL approach

  ## first estimate lambda0 and theta when weights are all 1

  w = rep(1,p)
  EMresults1 = BayesSSLemHetero(nScans=30000, burn=25000, thin=5, p = ncol(x), y = y[z==1],
                          x = x[z==1,], lambda1 = 0.1, lambda0start = 8,
                          numBlocks = 10, w=w)


  thetaEst1 = EMresults1$thetaEst
  lambda0est1 = EMresults1$lambda0est
  wj_vec = seq(0, 1, length=2000)

  wj_final1 = wj_vec[which(pStar(beta=0, thetaEst1^wj_vec,
                                 lambda1=0.1, lambda0=lambda0est1) < 0.1)[1]]



  if (length(activeX) == 0) {
    w= rep(1,p)
  } else {
    w = rep(1, p)
    w[activeX] = wj_final1
  }

  ## Now estimate lambda0 conditional on the weights

  EMresults1.2 = BayesSSLemHetero(nScans=30000, burn=25000, thin=5, p = ncol(x), y = y[z==1],
                            x = x[z==1,], lambda1 = 0.1, lambda0start = 8,
                            numBlocks = 10, w=w)

  lambda0est1.2 = EMresults1.2$lambda0est

  ## Now do final analysis

  MainAnalysisBayes1 = BayesSSLHetero(nScans=10000, burn=2000, thin=8, p = ncol(x), y = y[z==1],
                                x = x[z==1,], lambda1 = 0.1, lambda0 = lambda0est1.2,
                                numBlocks = 10, w=w)



  ####################################################################################
  ################################ Now do it for the controls ########################
  ####################################################################################


  ## first estimate lambda0 and theta when weights are all 1

  EMresults0 = BayesSSLemHetero(nScans=30000, burn=25000, thin=5, p = ncol(x), y = y[z==0],
                          x = x[z==0,], lambda1 = 0.1, lambda0start = 8,
                          numBlocks = 10, w=w)


  thetaEst0 = EMresults0$thetaEst
  lambda0est0 = EMresults0$lambda0est
  wj_vec = seq(0, 1, length=2000)

  wj_final0 = wj_vec[which(pStar(beta=0, thetaEst0^wj_vec,
                                 lambda1=0.1, lambda0=lambda0est0) < 0.1)[1]]



  if (length(activeX) == 0) {
    w= rep(1,p)
  } else {
    w = rep(1, p)
    w[activeX] = wj_final0
  }

  ## Now estimate lambda0 conditional on the weights

  EMresults0.2 = BayesSSLemHetero(nScans=30000, burn=25000, thin=5, p = ncol(x), y = y[z==0],
                            x = x[z==0,], lambda1 = 0.1, lambda0start = 8,
                            numBlocks = 10, w=w)

  lambda0est0.2 = EMresults0.2$lambda0est

  ## Now do final analysis

  MainAnalysisBayes0 = BayesSSLHetero(nScans=10000, burn=2000, thin=8, p = ncol(x), y = y[z==0],
                                x = x[z==0,], lambda1 = 0.1, lambda0 = lambda0est0.2,
                                numBlocks = 10, w=w)



  ####################################################################################
  ###################### Now combine them to estimate ATE ############################
  ####################################################################################


  Design = cbind(rep(1,n), x)

  atePost = rep(NA, dim(MainAnalysisBayes1)[1])

  for (i in 1 : dim(MainAnalysisBayes1)[1]) {
    atePost[i] = mean(Design %*% MainAnalysisBayes1[i,] -
                        Design %*% MainAnalysisBayes0[i,])
  }

  BetaSave[r,11] = mean(atePost, na.rm=TRUE)
  intervalBETA = quantile(atePost, c(.025, .975), na.rm=TRUE)
  IntervalSave[r,11] = as.numeric(beta < intervalBETA[2] & beta > intervalBETA[1])
  IntervalWidthSave[r,11] = intervalBETA[2] - intervalBETA[1]
  StandardErrorSave[r,11] = sd(atePost, na.rm=TRUE) 
  


  ####################################################################################
  ###################### Now only allow up to 10 in activeX ##########################
  #################################################################################### 

  fit.lasso.treat <- cv.glmnet(x=as.matrix(x), y=z, intercept=TRUE, family="binomial")
  nNonZero = length(which(coef(fit.lasso.treat, s='lambda.1se')[-1] != 0))
  if (nNonZero > 0) {
    activeX <- order(abs(coef(fit.lasso.treat, s='lambda.1se')[-1]), 
                     decreasing=TRUE)[1 : min(10, nNonZero)]
  } else {
    activeX = c()
  }

  ## EM-SSL approach

  w = rep(1, p)

  EMresults = BayesSSLem(nScans=30000, burn=25000, thin=5, n=n, p = ncol(x), y = y,
                         x = x, z=z, lambda1 = 0.1, lambda0start = 8,
                         numBlocks = 10, w=w)


  thetaEst = EMresults$thetaEst
  lambda0est = EMresults$lambda0est
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

  ## Now estimate lambda0 conditional on the weights

  EMresults2 = BayesSSLem(nScans=30000, burn=25000, thin=5, n=n, p = ncol(x), y = y,
                          x = x, z=z, lambda1 = 0.1, lambda0start = 8,
                          numBlocks = 10, w=w)

  lambda0est = EMresults2$lambda0est

  ## Now do final analysis

  MainAnalysisBayes = BayesSSL(nScans=10000, burn=2000, thin=8, n=n, p = ncol(x), y = y,
                               x = x, z=z, lambda1 = 0.1, lambda0 = lambda0est,
                               numBlocks = 10, w=w)

  BetaSave[r,12] = mean(MainAnalysisBayes[,2], na.rm=TRUE)
  intervalBETA = quantile(MainAnalysisBayes[,2], c(.025, .975), na.rm=TRUE)
  IntervalSave[r,12] = as.numeric(beta < intervalBETA[2] & beta > intervalBETA[1])
  IntervalWidthSave[r,12] = intervalBETA[2] - intervalBETA[1]
  StandardErrorSave[r,12] = sd(MainAnalysisBayes[,2], na.rm=TRUE)

  ## EM-SSL approach

  ## first estimate lambda0 and theta when weights are all 1

  w = rep(1,p)
  EMresults1 = BayesSSLemHetero(nScans=30000, burn=25000, thin=5, p = ncol(x), y = y[z==1],
                          x = x[z==1,], lambda1 = 0.1, lambda0start = 8,
                          numBlocks = 10, w=w)


  thetaEst1 = EMresults1$thetaEst
  lambda0est1 = EMresults1$lambda0est
  wj_vec = seq(0, 1, length=2000)

  wj_final1 = wj_vec[which(pStar(beta=0, thetaEst1^wj_vec,
                                 lambda1=0.1, lambda0=lambda0est1) < 0.1)[1]]



  if (length(activeX) == 0) {
    w= rep(1,p)
  } else {
    w = rep(1, p)
    w[activeX] = wj_final1
  }

  ## Now estimate lambda0 conditional on the weights

  EMresults1.2 = BayesSSLemHetero(nScans=30000, burn=25000, thin=5, p = ncol(x), y = y[z==1],
                            x = x[z==1,], lambda1 = 0.1, lambda0start = 8,
                            numBlocks = 10, w=w)

  lambda0est1.2 = EMresults1.2$lambda0est

  ## Now do final analysis

  MainAnalysisBayes1 = BayesSSLHetero(nScans=10000, burn=2000, thin=8, p = ncol(x), y = y[z==1],
                                x = x[z==1,], lambda1 = 0.1, lambda0 = lambda0est1.2,
                                numBlocks = 10, w=w)



  ####################################################################################
  ################################ Now do it for the controls ########################
  ####################################################################################


  ## first estimate lambda0 and theta when weights are all 1

  EMresults0 = BayesSSLemHetero(nScans=30000, burn=25000, thin=5, p = ncol(x), y = y[z==0],
                          x = x[z==0,], lambda1 = 0.1, lambda0start = 8,
                          numBlocks = 10, w=w)


  thetaEst0 = EMresults0$thetaEst
  lambda0est0 = EMresults0$lambda0est
  wj_vec = seq(0, 1, length=2000)

  wj_final0 = wj_vec[which(pStar(beta=0, thetaEst0^wj_vec,
                                 lambda1=0.1, lambda0=lambda0est0) < 0.1)[1]]



  if (length(activeX) == 0) {
    w= rep(1,p)
  } else {
    w = rep(1, p)
    w[activeX] = wj_final0
  }

  ## Now estimate lambda0 conditional on the weights

  EMresults0.2 = BayesSSLemHetero(nScans=30000, burn=25000, thin=5, p = ncol(x), y = y[z==0],
                            x = x[z==0,], lambda1 = 0.1, lambda0start = 8,
                            numBlocks = 10, w=w)

  lambda0est0.2 = EMresults0.2$lambda0est

  ## Now do final analysis

  MainAnalysisBayes0 = BayesSSLHetero(nScans=10000, burn=2000, thin=8, p = ncol(x), y = y[z==0],
                                x = x[z==0,], lambda1 = 0.1, lambda0 = lambda0est0.2,
                                numBlocks = 10, w=w)



  ####################################################################################
  ###################### Now combine them to estimate ATE ############################
  ####################################################################################


  Design = cbind(rep(1,n), x)

  atePost = rep(NA, dim(MainAnalysisBayes1)[1])

  for (i in 1 : dim(MainAnalysisBayes1)[1]) {
    atePost[i] = mean(Design %*% MainAnalysisBayes1[i,] -
                        Design %*% MainAnalysisBayes0[i,])
  }

  BetaSave[r,13] = mean(atePost, na.rm=TRUE)
  intervalBETA = quantile(atePost, c(.025, .975), na.rm=TRUE)
  IntervalSave[r,13] = as.numeric(beta < intervalBETA[2] & beta > intervalBETA[1])
  IntervalWidthSave[r,13] = intervalBETA[2] - intervalBETA[1]
  StandardErrorSave[r,13] = sd(atePost, na.rm=TRUE)
  
  save(BetaSave, file=paste("~/SSLconfounding/Code5/Output/11-beta-", 
                            index, ".dat", sep=''))
  save(wjSave, file=paste("~/SSLconfounding/Code5/Output/11-wj-", 
                          index, ".dat", sep=''))
  save(IntervalSave, file=paste("~/SSLconfounding/Code5/Output/11-interval-",
                                index, ".dat", sep=''))
  save(IntervalWidthSave, file=paste("~/SSLconfounding/Code5/Output/11-interval-width-",
                                     index, ".dat", sep=''))
  save(StandardErrorSave, file=paste("~/SSLconfounding/Code5/Output/11-standard-error-",
                                     index, ".dat", sep=''))
  
}






