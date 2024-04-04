## univariate trials

source("GND functions.R")

rho=.5
thetaA=c(rho,1,1,1)
prv=GND_sim(thetaA,T=1000)
thetaB=c(rho,0,0,2)
prv2=GND_sim(thetaB,T=1000)
par(mfrow=c(1,2))
plot(prv[,1],prv[,2],xlab="Y1",ylab="Y2",main="Simulated GND")
plot(prv2[,1],prv2[,2],xlab="Y1",ylab="Y2",main="Simulated GND")

THETA=cbind(thetaA,thetaB)
n=1000
Q=matrix(c(.9,.1,.1,.9),2,2)
init=c(.5,.5)
seed=123

sim=GNDHMM_sim(1000,THETA,Q,init,seed)
par(mfrow=c(1,1))
plot(sim$SimData[,1],sim$SimData[,2],col=sim$states,xlab="Y1",ylab="Y2",
     main="Simulated GND-HMM",pch=19)
legend("topright",legend=c("Regime 1","Regime 2"),col=1:2,pch=19)

y1=sim$SimData[,1]
y2=sim$SimData[,2]


fun = function(thetaa) {
  #THETA=matrix(thetaa,ncol=r)
  
  #theta=compute_theta(theta[1],theta[2],theta[3],theta[4])  
  eta=get_eta(thetaa)
  log_likelihood = -sum(dGND(thetaa,eta,y1,y2,logd=T))
  # for (i in 2:r) {
  #   eta=get_eta(THETA[,i])
  #   log_likelihood = log_likelihood + (-sum(dGND(THETA[,i],eta,y1,y2,logd=T)))
  # }
  
  return(log_likelihood)
}

fun(thetaA)

# OLD
# ineq=function(theta){
#   
#   z1=theta[5]^2-2*theta[6]*theta[1]
#   z2=theta[4]^2-2*theta[6]*theta[2]
#   # z3=theta[1]
#   # z4=theta[2]
#   return(c(z1,z2))
# }
# 
# thetaA=c(.5,1.5,1.5,2)
# y=GND_sim(thetaA,1000)
# 
# 
# fun(thetaA)
# #thetaa = theta_new = stats::optim(par = thetaA, fun, method = "L-BFGS-B")$par
# theta_new=Rsolnp::solnp(c(.5,1,1,2),fun,
#                         #ineqfun = ineq,ineqUB = c(0,0),ineqLB = c(-Inf,-Inf),
#                         LB=c(-.98,-2,-2,0),UB=c(.98,2,2,5))$par

#### Arturo #########
ineq2=function(theta){
  
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
  
  z11=theta5-sqrt(2*theta6*theta1)
  z12=theta5+sqrt(2*theta6*theta1)
  z21=theta4-sqrt(2*theta6*theta2)
  z22=theta4+sqrt(2*theta6*theta2)
  # z3=theta[1]
  # z4=theta[2]
  return(c(z11,z12,z21,z22))
}


theta_new=Rsolnp::solnp(c(.5,1,1,2),fun,
                        ineqfun = ineq2,ineqUB = c(0,Inf,0,Inf),ineqLB = c(-Inf,0,-Inf,0),
                        LB=c(-.98,-2,-2,0),UB=c(.98,2,2,5))$par
###################

# theta1, theta2, and theta3 constraints

rho=seq(-1,1,by=.1)
theta1=1/(1-rho^2)
plot(rho,theta1,type="l",xlab="rho",ylab="theta1",main="theta1 vs rho")

theta3=rho/(1-rho^2)
plot(rho,theta3,type="l",xlab="rho",ylab="theta3",main="theta3 vs rho")

###
source("GND functions.R")
# run ---------------------------------------------------------------------

library(MASS)
y=mvrnorm(100,mu=c(0.1,0.05),Sigma=matrix(c(1,.6,.6,1),2,2))
y1=y[,1]
y2=y[,2]

theta1=compute_theta(0.5,theta4 = 0,theta5=0.5,theta6=2)
theta2=compute_theta(0.9,theta4 = 0,theta5=0.5,theta6=2)

get_eta(y1,y2,theta1)

cdY1=cond_mom(.1,.2,theta1)
cdY2=cond_mom(.1,.2,theta2)

# x=seq(-10,10,by=.1)
# y=gfun_GND(x,theta0,cdY$sig2AB)
# plot(x,y,type='l',lwd=2)

# Set integral boundaries to -5 and 5 for convergence 
int=integrate(f=gfun_GND,lower = -5,upper = 5,theta=theta1
              #,sig2AB=cdY$sig2AB
              )
et=log(int$value)
eta

dGND(theta0,eta,.1,-.1,logd = T)
