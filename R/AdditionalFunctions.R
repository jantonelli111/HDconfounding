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


WAIChomo = function(y, x, t, Post,
                       type="continuous",
                       totalScans) {
  
  n = length(y)
  
  LHood = array(NA, dim=c(totalScans, n))
  
  designY = cbind(rep(1,n), t, x)
  
  if (type == "continuous") {
    for (ni in 1 : totalScans) {
      linPred = designY %*% Post$beta[ni,]
      LHood[ni,] = (1/(sqrt(2*pi*Post$sigma[ni]))) *
        exp((-(y - linPred)^2)/(2*Post$sigma[ni]))
    }
  } else {
    for (ni in 1 : totalScans) {
      linPred = pnorm(designY %*% Post$beta[ni,])
      LHood[ni,] = (linPred^y) * (1 - linPred)^(1-y)
    }
  }
  
  return(-2*(sum(log(apply(LHood, 2, mean))) - sum(apply(log(LHood), 2, sd)^2)))
}


WAIChetero = function(y, x, t, Post1, Post0,
                    type="continuous",
                    totalScans) {
  
  n = length(y)
  n1 = sum(t)
  n0 = n - n1
  
  y1 = y[t==1]
  y0 = y[t==0]
  
  LHood = array(NA, dim=c(totalScans, n))
  
  designY1 = cbind(rep(1,n1), x[t==1,])
  designY0 = cbind(rep(1,n0), x[t==0,])
  
  if (type == "continuous") {
    for (ni in 1 : totalScans) {
      linPred1 = designY1 %*% Post1$beta[ni,]
      linPred0 = designY0 %*% Post0$beta[ni,]
      LHood[ni,t==1] = (1/(sqrt(2*pi*Post1$sigma[ni]))) *
        exp((-(y1 - linPred1)^2)/(2*Post1$sigma[ni]))
      LHood[ni,t==0] = (1/(sqrt(2*pi*Post0$sigma[ni]))) *
        exp((-(y0 - linPred0)^2)/(2*Post0$sigma[ni]))
    }
  } else {
    for (ni in 1 : totalScans) {
      linPred1 = pnorm(designY1 %*% Post1$beta[ni,])
      linPred0 = pnorm(designY0 %*% Post0$beta[ni,])
      LHood[ni,t==1] = (linPred1^y1) * (1 - linPred1)^(1-y1)
      LHood[ni,t==0] = (linPred0^y0) * (1 - linPred0)^(1-y0)
    }
  }
  
  return(-2*(sum(log(apply(LHood, 2, mean))) - sum(apply(log(LHood), 2, sd)^2)))
}

