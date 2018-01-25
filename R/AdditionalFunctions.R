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

