source("GND functions.R")

# theta1=1/(1-thetaA[1]^2)
# theta2=theta1
# theta3=thetaA[1]/(1-thetaA[1]^2)
# theta4=thetaA[2]
# theta5=thetaA[3]
# theta6=thetaA[4]
# theta5^2-2*theta6*theta1<0
# theta4^2-2*theta6*theta2<0

fun = function(theta,Y,bounds=c(-5,5)) {
  
  # theta[1]=(exp(theta[1])-1)/(exp(theta[1])+1)
  # theta[4]=exp(theta[4])
  
  #theta[3]=0
  
  y1=Y[,1]
  y2=Y[,2]
  #theta=compute_theta(theta[1],theta[2],theta[3],theta[4])  
  eta=get_eta(theta,bounds)
  log_likelihood = -sum(dGND(theta,eta,y1,y2,logd=T))
  # for (i in 2:r) {
  #   eta=get_eta(THETA[,i])
  #   log_likelihood = log_likelihood + (-sum(dGND(THETA[,i],eta,y1,y2,logd=T)))
  # }
  
  return(log_likelihood)
}

# theta_inv=thetaA
# theta_inv[1]=log((1+thetaA[1])/(1-thetaA[1]))
# theta_inv[4]=log(thetaA[4])
# theta_inv
# thetaA
# #fun(theta_inv,Y)
fun(thetaA,Y)

ineq2=function(theta,Y,bounds){
  
  #theta[1]=(exp(theta[1])-1)/(exp(theta[1])+1)
  
  rho=theta[1]
  theta1=1/(1-rho^2)
  theta2=theta1
  theta3=rho/(1-rho^2)
  theta4=theta[2]
  
  theta5=theta[3]
  theta6=theta[4]
  #theta6=exp(theta[4])
  
  z11=theta5-sqrt(2*theta6*theta1)
  z12=theta5+sqrt(2*theta6*theta1)
  z21=theta4-sqrt(2*theta6*theta2)
  z22=theta4+sqrt(2*theta6*theta2)
  # z3=theta[1]
  # z4=theta[2]
  return(c(z11,z12,z21,z22))
}

thetaA=c(-0.5,-.1,0.5,2)
N=5000
Y=GND_sim(thetaA,N,seed=1,ngrid = 10^6)
plot(Y)
thetaA_new=Rsolnp::solnp(c(-0.5,0,0,1),
                        fun,
                        ineqfun = ineq2,
                        ineqUB = c(0,Inf,0,Inf),
                        ineqLB = c(-Inf,0,-Inf,0),
                        LB=c(-.98,-2.8,-2.8,.3),
                        UB=c(.98,2.8,2.8,6),
                        Y=Y,bounds=c(-5,5),
                        control=list(rho=1,tol=1e-16,trace=0))

thetaA_new$par
thetaA

thetaB=c(-0.9,-.1,0.5,2)
N=5000
Y=GND_sim(thetaB,N,seed=1,ngrid = 10^6)
plot(Y)
bounds=range(Y)
thetaB_new=Rsolnp::solnp(c(0,0,0,1),
                         fun,
                         ineqfun = ineq2,
                         ineqUB = c(0,Inf,0,Inf),
                         ineqLB = c(-Inf,0,-Inf,0),
                         LB=c(-.98,-2.8,-2.8,.3),
                         UB=c(.98,2.8,2.8,6),
                         Y=Y,bounds=bounds,
                         control=list(rho=1,tol=1e-16,trace=0))

thetaB_new$par
thetaB

# riduzione parametri
# stima un parametro alla volta (prima rho, poi rhp+theta4, poi rho+theta5...)

# theta_new=Rsolnp::solnp(theta_inv,fun,
#                         ineqfun = ineq2,ineqUB = c(0,Inf,0,Inf),ineqLB = c(-Inf,0,-Inf,0),
#                         LB=c(-3.5,-.9,-.9,-5),UB=c(3.5,.9,.9,5),Y=Y)

# theta_new$convergence
# 0 means convergence
# theta_newpars=theta_new$par;theta_newpars
# theta_inv
# thetaA
# theta_newpars[1]=(exp(theta_newpars[1])-1)/(exp(theta_newpars[1])+1)
# theta_newpars[4]=exp(theta_newpars[4])
# theta_newpars
# thetaA
