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
  theta4=theta[2]
  #theta5=theta[3]
  theta6=theta[3]
  
  f=exp((-1/2)*theta1*y1^2-(1/2)*theta2*y2^2+theta3
        *y1*y2+theta4*y1*y2^2
        #+theta5*y1^2*y2
        -theta6*y1^2*y2^2-eta)
  
  if(logd){
    f=log(f)
  }
  
  return(f)
}

# compute_theta=function(rho,theta4=0,theta5=0,theta6=0){
#   
#   # compute_theta computes theta1, theta2, and theta3 starting from the linear correlation coefficient rho
#   # It returns a vector of parameters of dimension 6
#   
#   theta=rep(0,6)
#   theta[1]=1/(1-rho^2)
#   theta[2]=theta[1]
#   theta[3]=rho/(1-rho^2)
#   theta[4]=theta4
#   theta[5]=theta5
#   theta[6]=theta6
#   
#   return(theta)
# }

cond_mom=function(y1,y2,theta){
  
  # cond_mom computes conditional first and second order moments
  # y1 is the vector of observations
  # y2 is the vector of observations for the conditioning variable
  # theta is the 4-dim vector of the GND parameters
  
  rho=theta[1]
  theta1=1/(1-rho^2)
  theta2=theta1
  theta3=rho/(1-rho^2)
  theta4=theta[2]
  #theta5=theta[3]
  theta6=theta[3]
  
  #check_cond1=theta5^2-2*theta6*theta1<0
  check_cond2=theta4^2-2*theta6*theta2<0
  
  if(
    #check_cond1&
     check_cond2){
    sig2AB=1/(theta1+2*theta6*y2^2
              #-2*theta5*y
              )
    muAB=(theta[3]*y2+theta4*y2^2)*sig2AB
    return(list(sig2AB=sig2AB,
                muAB=muAB))
  }
  else{
    stop("Inequality conditions not satisfied")
  }
}

gfun_GND=function(y,theta){
  # This function returns f(y)*exp(eta)
  # y is the vector of observations 
  # theta is the vector of parameters with 4 components, the first is the correlation coefficient, the second, third, and fourth are the parameters corresponding to theta4, theta5, and theta6 of the GND
  
  rho=theta[1]
  theta1=1/(1-rho^2)
  theta2=theta1
  theta3=rho/(1-rho^2)
  theta4=theta[2]
  #theta5=theta[3]
  theta6=theta[3]
  
  #check_cond1=theta5^2-2*theta6*theta1<0
  check_cond2=theta4^2-2*theta6*theta2<0
  
  
  if(
    #check_cond1&
     check_cond2){
    sig2AB=1/(theta1+2*theta6*y^2
             # -2*theta5*y
              )
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
  int=integrate(f=gfun_GND,lower = -5,upper = 5,theta=theta)
  
  eta=log(int$value)
  return(eta)
}


# Simulations ---------------------------------------------------------------------

GND_sim<-function(theta,T,seed=1234,ngrid=10000){
  
  #MC Simulations for the GND model
  #theta is the vector of parameters with the following components:
  #theta[1] is the correlation coefficient
  #theta[2:4] are theta_4, theta_5, and theta_6
  #T is the number of observations
  #seed is the seed for the random number generator
  #ngrid is the number of grid points for the Griddy-Gibbs algorithm
  
  set.seed(seed)
  # theta1=theta[1]
  # theta2=theta[2]
  # theta3=theta[3]
  # 
  theta4=theta[2]
  #theta5=theta[3]
  theta6=theta[3]
  rho=theta[1]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  stheta1<-stheta2<-1/sqrt(theta1)
  
  #Griddy-Gibbs for Y2:
  #ngrid=10000 
  grid=seq(from=-8,to=8,length.out=ngrid)
  coef1=theta1+2*theta6*grid^2
  #-2*theta5*grid
  sigma2_12=1/coef1
  ee1=exp(.5*(theta3*grid+theta4*grid^2)^2*sigma2_12 )
  pp=dnorm(grid,mean=0,sd=stheta1)*ee1*sqrt(sigma2_12)
  pp=pp/sum(pp)
  i = sample(ngrid, T, replace = TRUE, prob = pp)
  dx = diff(grid[1:2])
  Y2 = grid[i] + runif(T) * dx - dx/2
  coef1=theta1+2*theta6*Y2^2
  #-2*theta5*Y2
  sigma2_12=1/coef1
  mu_12=sigma2_12*(theta3*Y2+theta4*Y2^2)
  Y1=rnorm(T,mean=mu_12,sd=sqrt(sigma2_12))
  return(cbind(Y1,Y2))
}


# run ---------------------------------------------------------------------

source("GND functions.R")

thetaA=c(0.5,1,2)

N=5000
Y=GND_sim(thetaA,N,seed=1,ngrid = 10^6)
plot(Y)

fun = function(theta,Y) {
  
  theta[1]=(exp(theta[1])-1)/(exp(theta[1])+1)
  theta[3]=exp(theta[3])
  
  #theta[3]=0
  
  y1=Y[,1]
  y2=Y[,2]
  #theta=compute_theta(theta[1],theta[2],theta[3],theta[4])  
  eta=get_eta(theta)
  log_likelihood = -sum(dGND(theta,eta,y1,y2,logd=T))
  # for (i in 2:r) {
  #   eta=get_eta(THETA[,i])
  #   log_likelihood = log_likelihood + (-sum(dGND(THETA[,i],eta,y1,y2,logd=T)))
  # }
  
  return(log_likelihood)
}

fun(thetaA,Y)

ineq2=function(theta,Y){
  
  #theta[1]=(exp(theta[1])-1)/(exp(theta[1])+1)
  
  #rho=theta[1]
  rho=(exp(theta[1])-1)/(exp(theta[1])+1)
  theta1=1/(1-rho^2)
  theta2=theta1
  theta3=rho/(1-rho^2)
  # theta1=theta[1]
  # theta2=theta[2]
  # theta3=theta[3]
  theta4=theta[2]
  
  #theta[3]=0
  #theta5=theta[3]
  #theta6=theta[4]
  theta6=exp(theta[3])
  
  z11=#theta5
  -sqrt(2*theta6*theta1)
  z12=#theta5
  +sqrt(2*theta6*theta1)
  z21=theta4-sqrt(2*theta6*theta2)
  z22=theta4+sqrt(2*theta6*theta2)
  # z3=theta[1]
  # z4=theta[2]
  return(c(z11,z12,z21,z22))
}


# theta_new=Rsolnp::solnp(c(.5,0,0,5),fun,
#                         ineqfun = ineq2,ineqUB = c(0,Inf,0,Inf),ineqLB = c(-Inf,0,-Inf,0),
#                         LB=c(-.98,-2,-2,0),UB=c(.98,2,2,Inf),Y=Y)

theta_new=Rsolnp::solnp(c(0,1,0),fun,
                        ineqfun = ineq2,ineqUB = c(0,Inf,0,Inf),ineqLB = c(-Inf,0,-Inf,0),
                        LB=c(-10,-3,-10),UB=c(10,3,10),Y=Y)

theta_new$convergence
# 0 means convergence
theta_newpars=theta_new$par
theta_newpars[1]=(exp(theta_newpars[1])-1)/(exp(theta_newpars[1])+1)
theta_newpars[3]=exp(theta_newpars[3])
theta_newpars
thetaA
