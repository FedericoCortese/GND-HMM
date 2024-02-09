dGND=function(theta,eta,y1,y2,logd=T){
  
  # dGND computes density for a bivariate generalized normal distributed vector
  # theta is the vector of parameters
  # eta is the normalizing constant
  # y1 and y2 are vectors of time-series
  # logd=T returns log of the density
  
  theta1=theta[1]
  theta2=theta[2]
  theta3=theta[3]
  theta4=theta[4]
  theta5=theta[5]
  theta6=theta[6]
  
  f=exp((-1/2)*theta1*y1^2-(1/2)*theta2*y2^2+theta3
        *y1*y2+theta4*y1*y2^2+theta5
        *y1^2*y2-theta6*y1^2*y2^2-eta)
  
  if(logd){
    f=log(f)
  }
  
  return(f)
}

compute_theta=function(rho,theta4=0,theta5=0,theta6=0){
  
  # compute_theta computes theta1, theta2, and theta3 starting from the linear correlation coefficient rho
  
  theta=rep(0,6)
  theta[1]=1/(1-rho^2)
  theta[2]=theta[1]
  theta[3]=rho/(1-rho^2)
  theta[4]=theta4
  theta[5]=theta5
  theta[6]=theta6
  
  return(theta)
}

cond_mom=function(yA,yB,theta){
  
  # cond_mom computes conditional first and second order moments
  # yA is the vector of observations
  # yB is the vector of observations for the conditioning variable
  # theta is the vector of the GND parameters
  
  check_cond1=theta[5]^2-2*theta[6]*theta[1]<0
  check_cond2=theta[4]^2-2*theta[6]*theta[2]<0
  
  if(check_cond1&check_cond2){
    sig2AB=1/(theta[1]+2*theta[6]*yB^2)
    muAB=(theta[3]*yB+theta[4]*yB^2)*sig2AB
    return(list(sig2AB=sig2AB,
                muAB=muAB))
  }
  else{
    stop("Inequality conditions not satisfied")
  }
}

gfun_GND=function(y,theta,sig2AB){
  # This function returns f(y)*exp(eta)
  g=(2*pi)/sqrt(theta[2])*exp((1/2)*(theta[3]*y+theta[4]*y^2)^2*sig2AB)*sqrt(sig2AB*dnorm(y,sd=1/theta[2]))
  return(g)
}

get_eta=function(y1,y2,theta){
  
  # get_eta computes the normalizing constant of the GND density
  
  sig2AB=cond_mom(y1,y2,theta)
  sig2AB=sig2AB$sig2AB
  
  # Set integral boundaries to -5 and 5 for convergence 
  int=integrate(f=gfun_GND,lower = -5,upper = 5,theta=theta,sig2AB=sig2AB)
  eta=log(int$value)
  return(eta)
}

# EM ---------------------------------------------------------------------

EstGNDHMM=function (y1,y2, theta, Q){
  n = dim(y)[1]
  r = dim(Q)[2]
  eta_bar = matrix(0, n, r)
  eta = matrix(0, n, r)
  lambda = matrix(0, n, r)
  c = matrix(0, n, r)
  Lambda = array(0, c(r, r, n))
  M = matrix(0, r, r)
  
  for(j in 1:reg){
    c[, j] = dGND(theta[j,],et[j],y1,y2,logd=T)
  }

  ###
  eta_bar[n, ] = 1/r
  for (k in 1:(n - 1)) {
    i = n - k
    j = i + 1
    v = (eta_bar[j, ] * c[j, ]) %*% t(Q)
    eta_bar[i, ] = v/sum(v)
  }
  eta0 = rep(1, r)/r
  v = (eta0 %*% Q) * c[1, ]
  eta[1, ] = v/sum(v)
  for (i in 2:n) {
    v = (eta[i - 1, ] %*% Q) * c[i, ]
    eta[i, ] = v/sum(v)
  }
  v = eta * eta_bar
  sv0 = rowSums(v)
  for (j in 1:r) {
    lambda[, j] = v[, j]/sv0
  }
  gc = eta_bar * c
  M = Q * (as.matrix(eta0) %*% gc[1, ])
  MM = sum(M)
  Lambda[, , 1] = M/MM
  for (i in 2:n) {
    M = Q * (as.matrix(eta[i - 1, ]) %*% gc[i, ])
    MM = sum(M)
    Lambda[, , i] = M/MM
  }
  nu = colMeans(lambda)
  Qnew = Q
  if (r >= 2) {
    for (j in 1:r) {
      sv = rowSums(Lambda[j, , ], dims = 1)
      ssv = sum(sv)
      Qnew[j, ] = sv/ssv
    }
  }
  
  fun = function(thetaa) {
    if (r < 2) {
      log_likelihood = -sum(lambda[, 1] * dGND(theta[1,],eta[1],y1,y2,logd=T)
      )
    } else if (r > 1) {
      log_likelihood = -sum(lambda[, 1] * dGND(theta[1,],eta[1],y1,y2,logd=T))
      
      for (i in 2:r) {
        log_likelihood = log_likelihood + (-sum(lambda[, 
                                                       i] * dGND(theta[i,],eta[i],y1,y2,logd=T)
        )
        )
      }
      return(log_likelihood)
    }
  }
  if (r >= 2) {
    theta_new = stats::optim(par = theta, fun, method = "Nelder-Mead")$par
  } else if (r == 1) {
    theta_new = stats::optim(par = theta, fun, method = "L-BFGS-B")$par
  }
  
  out = list(nu = nu, theta_new = theta_new, Qnew = Qnew, eta = eta, 
             eta_bar = eta_bar, lambda = lambda, Lambda = Lambda)
  return(out)
}


