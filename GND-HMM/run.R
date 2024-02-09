###
source("GND functions.R")
# run ---------------------------------------------------------------------

library(MASS)
y=rmvnorm(100,mean=c(0.1,0.05),sigma=matrix(c(1,.6,.6,1),2,2))
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
int=integrate(f=gfun_GND,lower = -5,upper = 5,theta=theta0,sig2AB=cdY$sig2AB)
et=log(int$value)
eta

dGND(theta0,eta,.1,-.1,logd = T)
