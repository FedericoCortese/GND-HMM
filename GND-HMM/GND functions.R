
# GND functions ---------------------------------------------------------------------

dGND=function(theta,eta,y1,y2,logd=T){
  
  # dGND computes density for a bivariate generalized normal distributed vector
  # theta is the vector of parameters, first element is the correlation coefficient, the second, third, and fourth are the parameters corresponding to theta4, theta5, and theta6 of the GND
  # eta is the normalizing constant
  # y1 and y2 are vectors of time-series
  # logd=T returns log of the density
  
  rho=theta[1]
  theta1=1/(1-rho^2)
  theta2=theta1
  theta3=rho/(1-rho^2)
  # theta1=theta[1]
  # theta2=theta[2]
  # theta3=theta[3]
  theta4=theta[2]
  theta5=theta[3]
  theta6=theta[4]
  
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
  # It returns a vector of parameters of dimension 6
  
  theta=rep(0,6)
  theta[1]=1/(1-rho^2)
  theta[2]=theta[1]
  theta[3]=rho/(1-rho^2)
  theta[4]=theta4
  theta[5]=theta5
  theta[6]=theta6
  
  return(theta)
}

cond_mom=function(y1,y2,theta){
  
  # cond_mom computes conditional first and second order moments
  # y1 is the vector of observations
  # y2 is the vector of observations for the conditioning variable
  # theta is the vector of the GND parameters
  
  check_cond1=theta[5]^2-2*theta[6]*theta[1]<0
  check_cond2=theta[4]^2-2*theta[6]*theta[2]<0
  
  if(check_cond1&check_cond2){
    sig2AB=1/(theta[1]+2*theta[6]*y2^2-2*theta[5]*y)
    muAB=(theta[3]*y2+theta[4]*y2^2)*sig2AB
    return(list(sig2AB=sig2AB,
                muAB=muAB))
  }
  else{
    stop("Inequality conditions not satisfied")
  }
}

gfun_GND=function(y,theta
                  #,sig2AB
){
  # This function returns f(y)*exp(eta)
  # y is the vector of observations 
  # theta is the vector of parameters with 6 components, the first is the correlation coefficient, the second, third, and fourth are the parameters corresponding to theta4, theta5, and theta6 of the GND
  
  rho=theta[1]
  theta1=1/(1-rho^2)
  theta2=theta1
  theta3=rho/(1-rho^2)
  theta4=theta[2]
  theta5=theta[3]
  theta6=theta[4]
  
  check_cond1=theta5^2-2*theta6*theta1<0
  check_cond2=theta4^2-2*theta6*theta2<0
  
  
  if(check_cond1&check_cond2){
    sig2AB=1/(theta1+2*theta6*y^2-2*theta5*y)
    g=(2*pi)/sqrt(theta2)*exp((1/2)*(theta3*y+theta4*y^2)^2*sig2AB)*
      sqrt(sig2AB*dnorm(y,sd=1/theta2))
  }
  else{
    stop("Inequality conditions not satisfied")
  }
  return(g)
}

get_eta=function(
    #y1,y2,
  theta){
  
  # get_eta computes the normalizing constant of the GND density
  # y1 and y2 are vectors of observations
  # theta is the vector of parameters with 6 components
  
  # sig2AB=cond_mom(y1,y2,theta)
  # sig2AB=sig2AB$sig2AB
  
  # Set integral boundaries to -5 and 5 for convergence 
  int=integrate(f=gfun_GND,lower = -5,upper = 5,theta=theta
                #,sig2AB=sig2AB
                )
  eta=log(int$value)
  return(eta)
}


# Simulations ---------------------------------------------------------------------

GND_sim<-function(theta,T){
  
  #MC Simulations for the GND model
  #theta is the vector of parameters with the following components:
  #theta[1] is the correlation coefficient
  #theta[2:4] are theta_4, theta_5, and theta_6
  #T is the number of observations
  
  # theta1=theta[1]
  # theta2=theta[2]
  # theta3=theta[3]
  # 
  theta4=theta[2]
  theta5=theta[3]
  theta6=theta[4]
  rho=theta[1]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  stheta1<-stheta2<-1/sqrt(theta1)
  
  #Griddy-Gibbs for Y2:
  ngrid=10000 
  grid=seq(from=-8,to=8,length.out=ngrid)
  coef1=theta1+2*theta6*grid^2-2*theta5*grid
  sigma2_12=1/coef1
  ee1=exp(.5*(theta3*grid+theta4*grid^2)^2*sigma2_12 )
  pp=dnorm(grid,mean=0,sd=stheta1)*ee1*sqrt(sigma2_12)
  pp=pp/sum(pp)
  i = sample(ngrid, T, replace = TRUE, prob = pp)
  dx = diff(grid[1:2])
  Y2 = grid[i] + runif(T) * dx - dx/2
  coef1=theta1+2*theta6*Y2^2-2*theta5*Y2
  sigma2_12=1/coef1
  mu_12=sigma2_12*(theta3*Y2+theta4*Y2^2)
  Y1=rnorm(T,mean=mu_12,sd=sqrt(sigma2_12))
  return(cbind(Y1,Y2))
}



GNDHMM_sim=function(n,THETA, Q, init,seed){
  # This function simulates a GND-HMM model
  # n is the number of observations
  # THETA is a matrix of dimension 4 x k where k is the number of regimes, so that each column is the state specific vector of parameters 
  # Q is the transition matrix
  # init is the initial distribution
  # seed is the seed for the random number generator
  
  #markov chain simulation
  reg = dim(Q)[1]
  x <- numeric(n)
  set.seed(seed)
  x[1] <- sample(1:reg, 1, prob = init)
  for(i in 2:n){
    x[i] <- sample(1:reg, 1, prob = Q[x[i - 1], ])
  }
  
  d=2
  Sim = matrix(0, n, d * reg)
  SimData = matrix(0, n, d)
  
  #pseudo-observations simulation
  for (k in 1:reg) {
    u=GND_sim(THETA[,k],n)
    #u = rCopula(n, copula::tCopula(param=P2p(R[,,k]), dim = d,df=nu[k],dispstr = "un"))
    Sim[, (d * k - d + 1):(d * k)] = u
  }
  
  for (i in 1:n) {
    k = x[i]
    SimData[i, ] = Sim[i, (d * k - d + 1):(d * k)]
  }
  return(list(SimData=SimData,states=x))
  
}


# EM ---------------------------------------------------------------------

EstGNDHMM=function (y1,y2, THETA, Q){
  n = dim(y)[1]
  r = dim(Q)[2]
  eta_bar = matrix(0, n, r)
  eta = matrix(0, n, r)
  lambda = matrix(0, n, r)
  c = matrix(0, n, r)
  Lambda = array(0, c(r, r, n))
  M = matrix(0, r, r)
  
  for(j in 1:reg){
    c[, j] = dGND(THETA[j,],et[j],y1,y2,logd=T)
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
      log_likelihood = -sum(lambda[, 1] * dGND(theta[1,],et[1],y1,y2,logd=T)
      )
    } else if (r > 1) {
      log_likelihood = -sum(lambda[, 1] * dGND(theta[1,],et[1],y1,y2,logd=T))
      
      for (i in 2:r) {
        log_likelihood = log_likelihood + (-sum(lambda[, 
                                                       i] * dGND(theta[i,],et[i],y1,y2,logd=T)
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

