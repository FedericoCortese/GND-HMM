source("GND functions.R")

thetaA=c(0.5,1,0,2)

N=5000
Y=GND_sim(thetaA,N,seed=1,ngrid = 10^6)
plot(Y)

fun = function(theta,Y) {
  
  theta[1]=(exp(theta[1])-1)/(exp(theta[1])+1)
  theta[4]=exp(theta[4])
  
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
  theta5=theta[3]
  #theta6=theta[4]
  theta6=exp(theta[4])
  
  z11=theta5-sqrt(2*theta6*theta1)
  z12=theta5+sqrt(2*theta6*theta1)
  z21=theta4-sqrt(2*theta6*theta2)
  z22=theta4+sqrt(2*theta6*theta2)
  # z3=theta[1]
  # z4=theta[2]
  return(c(z11,z12,z21,z22))
}


# theta_new=Rsolnp::solnp(c(.5,0,0,5),fun,
#                         ineqfun = ineq2,ineqUB = c(0,Inf,0,Inf),ineqLB = c(-Inf,0,-Inf,0),
#                         LB=c(-.98,-2,-2,0),UB=c(.98,2,2,Inf),Y=Y)

theta_new=Rsolnp::solnp(c(0,1,1,0),fun,
                        ineqfun = ineq2,ineqUB = c(0,Inf,0,Inf),ineqLB = c(-Inf,0,-Inf,0),
                        LB=c(-10,-3,-3,-10),UB=c(10,3,3,10),Y=Y)

theta_new$convergence
# 0 means convergence
theta_newpars=theta_new$par
theta_newpars[1]=(exp(theta_newpars[1])-1)/(exp(theta_newpars[1])+1)
theta_newpars[4]=exp(theta_newpars[4])
theta_newpars
thetaA
