
CoMoments<-function(theta){
  #CoMoments for the model of the ANOR paper
  
  theta4=theta[1]
  theta5=theta[2]
  theta6=theta[3]
  rho=theta[4]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  integrand_den1=function(y2){
    #returns f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    stheta2=1/sqrt(theta2)
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1)
  }
  
  
  integrand_den2=function(y1){
    #returns f(y1)*exp(eta)
    coef2=theta2+2*theta6*y1^2-2*theta4*y1
    sigma2_21=1/coef2
    sigma_21=sqrt(sigma2_21)
    mu_21=sigma2_21*(theta3*y1+theta5*y1^2)
    #ee2=exp(.5*(theta3*y1+theta5*y1^2)^2*sigma2_21 )
    ee2=exp(.5*mu_21^2/sigma2_21 )
    stheta1=1/sqrt(theta1)
    out2=2*pi*stheta1*ee2*sigma_21*dnorm(y1,sd=stheta1)
    return(out2)
  }
  
  
  integrand_num12=function(y2){
    #returns y2^2*mu_12(y2)*f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    stheta2=1/sqrt(theta2)
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    out=y2^2*mu_12*f2
    return(out)
  }
  
  
  integrand_num21=function(y1){
    #returns y1^2*mu_21(y1)*f(y1)*exp(eta)
    coef2=theta2+2*theta6*y1^2-2*theta4*y1
    sigma2_21=1/coef2
    sigma_21=sqrt(sigma2_21)
    mu_21=sigma2_21*(theta3*y1+theta5*y1^2)
    #ee2=exp(.5*(theta3*y1+theta5*y1^2)^2*sigma2_21 )
    ee2=exp(.5*mu_21^2/sigma2_21 )
    stheta1=1/sqrt(theta1)
    f1=2*pi*stheta1*ee2*sigma_21*dnorm(y1,sd=stheta1)
    out=y1^2*mu_21*f1
    return(out)
  }
  
  
  integrand_num11=function(y2){
    #returns y2*mu_12(y2)*f2(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    stheta2=1/sqrt(theta2)
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    out=y2*mu_12*f2
    return(out)
  }
  
  
  integrand_num22=function(y2){
    #returns y2*( mu_12(y2)^2+ sigma^2_12(y2)  ) *f2(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    stheta2=1/sqrt(theta2)
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    out=y2^2*(mu_12^2 +  sigma2_12  )*f2
    return(out)
  }
  
  # 
  # integrand_num12_2=function(z){
  #   coef2=1+2*theta6*z^2-2*theta5*z
  #   ee2=exp(.5*theta4^2*z^4/coef2)
  #   out=theta4*z^4*ee2*coef2^(-1.5) *dnorm(z)
  #   return(out)
  # }
  # 
  # integrand_num12_3=function(z){
  #   sigma2_21=1/(1+2*theta6*z^2-2*theta4*z)
  #   ee3=exp(.5*theta5^2*z^4*sigma2_21)
  #   out= z*ee3*(  theta5^2*z^4 * sigma2_21^(2.5) +  
  #                   sigma2_21^(1.5))*dnorm(z)
  #   return(out)
  
  NUM12=integrate(integrand_num12,-5,5)$value#exp(eta)*E(Y1*Y2^2)
  NUM21=integrate(integrand_num21,-5,5)$value#exp(eta)*E(Y1^2*Y2)
  NUM11=integrate(integrand_num11,-5,5)$value#exp(eta)*E(Y1*Y2)
  NUM22=integrate(integrand_num22,-5,5)$value#exp(eta)*E(Y1^2*Y2^2)
  DEN1=integrate(integrand_den1,-5,5)$value#exp(eta)
  eta1=log(DEN1)
  Coskew12=NUM12/DEN1 
  Coskew21=NUM21/DEN1
  Covariance=NUM11/DEN1
  Covolatility=NUM22/DEN1
  
  # return(list(eta=c(eta1,eta2),
  #             Coskew=c(Coskew12_1,Coskew12_2,
  #                      Coskew12_3,Coskew12_4)))
  # 
  
  return( list(eta=eta1,Covariance=Covariance,
               Coskew12= Coskew12 , Coskew21= Coskew21,
               Covolatility=Covolatility) ) 
  
}

CoMoments_v2<-function(theta,tol=1e-8){
  #CoMoments for the model of the ANOR paper
  #tol: integral lower and upper limits
  # are set to qnorm(tol,sd=stheta2) and qnorm(1-tol,sd=stheta2)
  # with stheta2=1/sqrt(theta1)=1/sqrt(theta2)
  theta4=theta[1]
  theta5=theta[2]
  theta6=theta[3]
  rho=theta[4]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  stheta1<-stheta2<-1/sqrt(theta1)
  
  integrand_den1=function(y2){
    #returns f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1)
  }
  
  
  integrand_den2=function(y1){
    #returns f(y1)*exp(eta)
    coef2=theta2+2*theta6*y1^2-2*theta4*y1
    sigma2_21=1/coef2
    sigma_21=sqrt(sigma2_21)
    mu_21=sigma2_21*(theta3*y1+theta5*y1^2)
    #ee2=exp(.5*(theta3*y1+theta5*y1^2)^2*sigma2_21 )
    ee2=exp(.5*mu_21^2/sigma2_21 )
    out2=2*pi*stheta1*ee2*sigma_21*dnorm(y1,sd=stheta1)
    return(out2)
  }
  
  
  integrand_num12=function(y2){
    #returns y2^2*mu_12(y2)*f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    out=y2^2*mu_12*f2
    return(out)
  }
  
  
  integrand_num21=function(y1){
    #returns y1^2*mu_21(y1)*f(y1)*exp(eta)
    coef2=theta2+2*theta6*y1^2-2*theta4*y1
    sigma2_21=1/coef2
    sigma_21=sqrt(sigma2_21)
    mu_21=sigma2_21*(theta3*y1+theta5*y1^2)
    #ee2=exp(.5*(theta3*y1+theta5*y1^2)^2*sigma2_21 )
    ee2=exp(.5*mu_21^2/sigma2_21 )
    f1=2*pi*stheta1*ee2*sigma_21*dnorm(y1,sd=stheta1)
    out=y1^2*mu_21*f1
    return(out)
  }
  
  
  integrand_num11=function(y2){
    #returns y2*mu_12(y2)*f2(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    out=y2*mu_12*f2
    return(out)
  }
  
  
  integrand_num22=function(y2){
    #returns y2*( mu_12(y2)^2+ sigma^2_12(y2)  ) *f2(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    out=y2^2*(mu_12^2 +  sigma2_12  )*f2
    return(out)
  }
  
  # 
  # integrand_num12_2=function(z){
  #   coef2=1+2*theta6*z^2-2*theta5*z
  #   ee2=exp(.5*theta4^2*z^4/coef2)
  #   out=theta4*z^4*ee2*coef2^(-1.5) *dnorm(z)
  #   return(out)
  # }
  # 
  # integrand_num12_3=function(z){
  #   sigma2_21=1/(1+2*theta6*z^2-2*theta4*z)
  #   ee3=exp(.5*theta5^2*z^4*sigma2_21)
  #   out= z*ee3*(  theta5^2*z^4 * sigma2_21^(2.5) +  
  #                   sigma2_21^(1.5))*dnorm(z)
  #   return(out)
  
  LowerLim=qnorm(tol,sd=stheta2)
  UpperLim=qnorm(1-tol,sd=stheta2)
  NUM12=integrate(integrand_num12,LowerLim,UpperLim)$value#exp(eta)*E(Y1*Y2^2)
  NUM21=integrate(integrand_num21,LowerLim,UpperLim)$value#exp(eta)*E(Y1^2*Y2)
  NUM11=integrate(integrand_num11,LowerLim,UpperLim)$value#exp(eta)*E(Y1*Y2)
  NUM22=integrate(integrand_num22,LowerLim,UpperLim)$value#exp(eta)*E(Y1^2*Y2^2)
  DEN1=integrate(integrand_den1,LowerLim,UpperLim)$value#exp(eta)
  eta1=log(DEN1)
  Coskew12=NUM12/DEN1 
  Coskew21=NUM21/DEN1
  Covariance=NUM11/DEN1
  Covolatility=NUM22/DEN1
  
  # return(list(eta=c(eta1,eta2),
  #             Coskew=c(Coskew12_1,Coskew12_2,
  #                      Coskew12_3,Coskew12_4)))
  # 
  
  return( list(eta=eta1,Covariance=Covariance,
               Coskew12= Coskew12 , Coskew21= Coskew21,
               Covolatility=Covolatility) ) 
  
}

JointCDF<-function(a,b,theta,tol=1e-8){
  #P(Y1<=a,Y2<=b) for the model of the ANOR paper
  #tol: integral lower and upper limits
  # are set to qnorm(tol,sd=stheta2) and qnorm(1-tol,sd=stheta2)
  # with stheta2=1/sqrt(theta1)=1/sqrt(theta2)
  theta4=theta[1]
  theta5=theta[2]
  theta6=theta[3]
  rho=theta[4]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  stheta1<-stheta2<-1/sqrt(theta1)
  
  integrand_den1=function(y2){
    #returns f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1)
  }
  
  
  integrand_den2=function(y1){
    #returns f(y1)*exp(eta)
    coef2=theta2+2*theta6*y1^2-2*theta4*y1
    sigma2_21=1/coef2
    sigma_21=sqrt(sigma2_21)
    mu_21=sigma2_21*(theta3*y1+theta5*y1^2)
    #ee2=exp(.5*(theta3*y1+theta5*y1^2)^2*sigma2_21 )
    ee2=exp(.5*mu_21^2/sigma2_21 )
    stheta1=1/sqrt(theta1)
    out2=2*pi*stheta1*ee2*sigma_21*dnorm(y1,sd=stheta1)
    return(out2)
  }
  
  
  integrand_num=function(y2){
    #returns Phi( (a-mu_12(y2))/sigma_12(y2) )*f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    A=(a-mu_12)/sigma_12
    out= pnorm(A)*f2
    return(out)
  }
  
  LowerLim=qnorm(tol,sd=stheta2)
  UpperLim=qnorm(1-tol,sd=stheta2)
  DEN=integrate(integrand_den1,LowerLim,UpperLim)$value#exp(eta)
  NUM=integrate(integrand_num,LowerLim,b)$value
  #eta1=log(DEN1)
  FF=NUM/DEN
  
  return( FF ) 
  
}

JointCDF_v2<-function(a,b,theta,tol=1e-8){
  #P(Y1<=a,Y2<=b) for the model of the ANOR paper
  #a and b can be n-vectors
  #The output is an n-vector with components P(Y1<=a[i],Y2<=b[i])
  #i=1,...,n
  #tol: integral lower and upper limits
  # are set to qnorm(tol,sd=stheta2) and qnorm(1-tol,sd=stheta2)
  # with stheta2=1/sqrt(theta1)=1/sqrt(theta2)
  theta4=theta[1]
  theta5=theta[2]
  theta6=theta[3]
  rho=theta[4]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  stheta1<-stheta2<-1/sqrt(theta1)
  
  integrand_den1=function(y2){
    #returns f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1)
  }
  
  
  
  
  integrand_num=function(y2,a){
    #returns Phi( (a-mu_12(y2))/sigma_12(y2) )*f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    A=(a-mu_12)/sigma_12
    out= pnorm(A)*f2
    return(out)
  }
  
  LowerLim=qnorm(tol,sd=stheta2)
  UpperLim=qnorm(1-tol,sd=stheta2)
  NUM=function(a,b){
    integrate(integrand_num,LowerLim,b,a)$value
  }
  DEN=integrate(integrand_den1,LowerLim,UpperLim)$value#exp(eta)
  #eta1=log(DEN1)
  FF=mapply(NUM,a,b)/DEN
  
  return( FF ) 
  
}

JointCDF_v3<-function(a,b,theta,tol=1e-8){
  #P(Y1<=a,Y2<=b) for the model of the ANOR paper
  #a and b can be n-vectors and m-vectors
  #The output is an nxm-matrix with components P(Y1<=a[i],Y2<=b[j])
  #i=1,...,n 
  #j=1,...,m
  #tol: integral lower and upper limits
  # are set to qnorm(tol,sd=stheta2) and qnorm(1-tol,sd=stheta2)
  # with stheta2=1/sqrt(theta1)=1/sqrt(theta2)
  theta4=theta[1]
  theta5=theta[2]
  theta6=theta[3]
  rho=theta[4]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  stheta1<-stheta2<-1/sqrt(theta1)
  
  integrand_den1=function(y2){
    #returns f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1)
  }
  
  
  
  
  integrand_num=function(y2,a){
    #returns Phi( (a-mu_12(y2))/sigma_12(y2) )*f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    A=(a-mu_12)/sigma_12
    out= pnorm(A)*f2
    return(out)
  }
  
  LowerLim=qnorm(tol,sd=stheta2)
  UpperLim=qnorm(1-tol,sd=stheta2)
  NUM=function(a,b){
    integrate(integrand_num,LowerLim,b,a)$value
  }
  DEN=integrate(integrand_den1,LowerLim,UpperLim)$value#exp(eta)
  #eta1=log(DEN1)
  #FF=outer(a,b,NUM)/DEN
  FF=sapply(b, function(y) mapply(NUM,a,y))/DEN
  
  return( FF ) 
  
}


CDF_Y1_Y2<-function(a,b,theta,tol=1e-8){
  #P(Y1<=a) and P(Y2<=b) for the model of the ANOR paper
  #tol: integral lower and upper limits
  # are set to qnorm(tol,sd=stheta2) and qnorm(1-tol,sd=stheta2)
  # with stheta2=1/sqrt(theta1)=1/sqrt(theta2)
  theta4=theta[1]
  theta5=theta[2]
  theta6=theta[3]
  rho=theta[4]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  stheta1<-stheta2<-1/sqrt(theta1)
  
  integrand_den1=function(y1){
    #returns f(y1)*exp(eta)
    coef2=theta2+2*theta6*y1^2-2*theta4*y1
    sigma2_21=1/coef2
    sigma_21=sqrt(sigma2_21)
    mu_21=sigma2_21*(theta3*y1+theta5*y1^2)
    #ee2=exp(.5*(theta3*y1+theta5*y1^2)^2*sigma2_21 )
    ee2=exp(.5*mu_21^2/sigma2_21 )
    out2=2*pi*stheta1*ee2*sigma_21*dnorm(y1,sd=stheta1)
    return(out2)
  }
  
  integrand_den2=function(y2){
    #returns f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1)
  }
  
  
  LowerLim=qnorm(tol,sd=stheta2)
  UpperLim=qnorm(1-tol,sd=stheta2)
  DEN=integrate(integrand_den1,LowerLim,UpperLim)$value#exp(eta)
  
  integrand_num=integrand_den1
  NUM=integrate(integrand_num,LowerLim,a)$value
  F1=NUM/DEN
  
  integrand_num=integrand_den2
  NUM=integrate(integrand_num,LowerLim,b)$value
  F2=NUM/DEN
  
  return( c(F1,F2) ) 
  
}

CDF_Y1_Y2_v2<-function(a,b,theta,tol=1e-8){
  #P(Y1<=a) and P(Y2<=b) for the model of the ANOR paper
  #a and b can be n-vectors
  #The output is an nx2 matrix with i-th row P(Y1<=a[i]), P(Y2<=b[i])
  #i=1,...,n
  #tol: integral lower and upper limits
  # are set to qnorm(tol,sd=stheta2) and qnorm(1-tol,sd=stheta2)
  # with stheta2=1/sqrt(theta1)=1/sqrt(theta2)
  theta4=theta[1]
  theta5=theta[2]
  theta6=theta[3]
  rho=theta[4]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  stheta1<-stheta2<-1/sqrt(theta1)
  
  integrand_den1=function(y1){
    #returns f(y1)*exp(eta)
    coef2=theta2+2*theta6*y1^2-2*theta4*y1
    sigma2_21=1/coef2
    sigma_21=sqrt(sigma2_21)
    mu_21=sigma2_21*(theta3*y1+theta5*y1^2)
    #ee2=exp(.5*(theta3*y1+theta5*y1^2)^2*sigma2_21 )
    ee2=exp(.5*mu_21^2/sigma2_21 )
    out2=2*pi*stheta1*ee2*sigma_21*dnorm(y1,sd=stheta1)
    return(out2)
  }
  
  integrand_den2=function(y2){
    #returns f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1)
  }
  
  
  
  LowerLim=qnorm(tol,sd=stheta2)
  UpperLim=qnorm(1-tol,sd=stheta2)
  DEN=integrate(integrand_den1,LowerLim,UpperLim)$value#exp(eta)
  
  
  
  ff1=function(a){
    #returns \int_-Inf^a f(y1)*exp(eta)
    integrand_num=integrand_den1
    out=integrate(integrand_num,LowerLim,a)$value
    return(out)
  }
  
  NUM=sapply(a,ff1)
  F1=NUM/DEN
  
  
  ff2=function(b){
    #returns \int_-Inf^b f(y2)*exp(eta)
    integrand_num=integrand_den2
    out=integrate(integrand_num,LowerLim,b)$value
    return(out)
  }
  
  NUM=sapply(b,ff2)
  F2=NUM/DEN
  
  return( cbind(F1,F2) ) 
  
}

JointMGF<-function(t1,t2,theta,tol=1e-8){
  #E(exp(t1*Y1+t2*Y2) for the model of the ANOR paper
  #tol: integral lower and upper limits
  # are set to qnorm(tol,sd=stheta2) and qnorm(1-tol,sd=stheta2)
  # with stheta2=1/sqrt(theta1)=1/sqrt(theta2)
  theta4=theta[1]
  theta5=theta[2]
  theta6=theta[3]
  rho=theta[4]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  stheta1<-stheta2<-1/sqrt(theta1)
  
  integrand_den1=function(y2){
    #returns f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1)
  }
  
  
  integrand_den2=function(y1){
    #returns f(y1)*exp(eta)
    coef2=theta2+2*theta6*y1^2-2*theta4*y1
    sigma2_21=1/coef2
    sigma_21=sqrt(sigma2_21)
    mu_21=sigma2_21*(theta3*y1+theta5*y1^2)
    #ee2=exp(.5*(theta3*y1+theta5*y1^2)^2*sigma2_21 )
    ee2=exp(.5*mu_21^2/sigma2_21 )
    out2=2*pi*stheta1*ee2*sigma_21*dnorm(y1,sd=stheta1)
    return(out2)
  }
  
  
  integrand_num=function(y2){
    #returns exp(t2*y2+ t1*mu_12(y2))+t1^2/2*sigma2_12(y2) )*f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    out= exp( t2*y2+t1*mu_12+t1^2/2*sigma2_12 )*f2
    return(out)
  }
  
  LowerLim=qnorm(tol,sd=stheta2)
  UpperLim=qnorm(1-tol,sd=stheta2)
  DEN=integrate(integrand_den1,LowerLim,UpperLim)$value#exp(eta)
  NUM=integrate(integrand_num,LowerLim,UpperLim)$value
  #eta1=log(DEN1)
  MGF=NUM/DEN
  
  return( MGF ) 
  
}

JointMGF_v2<-function(t1,t2,theta,tol=1e-8){
  #E(exp(t1*Y1+t2*Y2) for the model of the ANOR paper
  #t1 and t2 can be n-vectors
  #The output is an n-vector with components E(exp(t1[i]*Y1+t2[i]*Y2)
  #i=1,...,n
  #tol: integral lower and upper limits
  # are set to qnorm(tol,sd=stheta2) and qnorm(1-tol,sd=stheta2)
  # with stheta2=1/sqrt(theta1)=1/sqrt(theta2)
  theta4=theta[1]
  theta5=theta[2]
  theta6=theta[3]
  rho=theta[4]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  stheta1<-stheta2<-1/sqrt(theta1)
  
  integrand_den1=function(y2){
    #returns f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1)
  }
  
  
  integrand_num=function(y2,t1,t2){
    #returns exp(t2*y2+ t1*mu_12(y2))+t1^2/2*sigma2_12(y2) )*f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    out= exp( t2*y2+t1*mu_12+t1^2/2*sigma2_12 )*f2
    return(out)
  }
  
  LowerLim=qnorm(tol,sd=stheta2)
  UpperLim=qnorm(1-tol,sd=stheta2)
  NUM=function(t1,t2){
    integrate(integrand_num,LowerLim,UpperLim,t1,t2)$value
  }
  DEN=integrate(integrand_den1,LowerLim,UpperLim)$value#exp(eta)
  #eta1=log(DEN1)
  MGF=mapply(NUM,t1,t2)/DEN
  
  return( MGF ) 
  
}

JointMGF_v3<-function(t1,t2,theta,tol=1e-8){
  #E(exp(t1*Y1+t2*Y2) for the model of the ANOR paper
  #t1 and t2 can be n-vectors and m-vectors
  #The output is an nxm-matrix with components E(exp(t1[i]*Y1+t2[j]*Y2)
  #i=1,...,n 
  #j=1,...,m
  #tol: integral lower and upper limits
  # are set to qnorm(tol,sd=stheta2) and qnorm(1-tol,sd=stheta2)
  # with stheta2=1/sqrt(theta1)=1/sqrt(theta2)
  theta4=theta[1]
  theta5=theta[2]
  theta6=theta[3]
  rho=theta[4]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  stheta1<-stheta2<-1/sqrt(theta1)
  
  integrand_den1=function(y2){
    #returns f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1)
  }
  
  integrand_num=function(y2,t1,t2){
    #returns exp(t2*y2+ t1*mu_12(y2))+t1^2/2*sigma2_12(y2) )*f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    out= exp( t2*y2+t1*mu_12+t1^2/2*sigma2_12 )*f2
    return(out)
  }
  
  LowerLim=qnorm(tol,sd=stheta2)
  UpperLim=qnorm(1-tol,sd=stheta2)
  NUM=function(t1,t2){
    integrate(integrand_num,LowerLim,UpperLim,t1,t2)$value
  }
  DEN=integrate(integrand_den1,LowerLim,UpperLim)$value#exp(eta)
  #eta1=log(DEN1)
  #MGF=outer(t1,t2,NUM)/DEN
  MGF=sapply(t2, function(y) mapply(NUM,t1,y))/DEN
  
  return( MGF ) 
  
}

JointCF<-function(t1,t2,theta,tol=1e-8){
  #E(exp(i*t1*Y1+i*t2*Y2) for the model of the ANOR paper
  #tol: integral lower and upper limits
  # are set to qnorm(tol,sd=stheta2) and qnorm(1-tol,sd=stheta2)
  # with stheta2=1/sqrt(theta1)=1/sqrt(theta2)
  theta4=theta[1]
  theta5=theta[2]
  theta6=theta[3]
  rho=theta[4]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  stheta1<-stheta2<-1/sqrt(theta1)
  
  integrand_den1=function(y2){
    #returns f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1)
  }
  
  
  integrand_den2=function(y1){
    #returns f(y1)*exp(eta)
    coef2=theta2+2*theta6*y1^2-2*theta4*y1
    sigma2_21=1/coef2
    sigma_21=sqrt(sigma2_21)
    mu_21=sigma2_21*(theta3*y1+theta5*y1^2)
    #ee2=exp(.5*(theta3*y1+theta5*y1^2)^2*sigma2_21 )
    ee2=exp(.5*mu_21^2/sigma2_21 )
    out2=2*pi*stheta1*ee2*sigma_21*dnorm(y1,sd=stheta1)
    return(out2)
  }
  
  
  integrand_num_Re=function(y2){
    #returns exp(t2*y2+ t1*mu_12(y2))+t1^2/2*sigma2_12(y2) )*f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    cc=t2*y2+t1*mu_12
    out= cos(cc)*exp(-t1^2/2*sigma2_12 )*f2
    return(out)
  }
  
  integrand_num_Im=function(y2){
    #returns exp(t2*y2+ t1*mu_12(y2))+t1^2/2*sigma2_12(y2) )*f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    cc=t2*y2+t1*mu_12
    out= sin(cc)*exp(-t1^2/2*sigma2_12 )*f2
    return(out)
  }
  
  LowerLim=qnorm(tol,sd=stheta2)
  UpperLim=qnorm(1-tol,sd=stheta2)
  DEN=integrate(integrand_den1,LowerLim,UpperLim)$value#exp(eta)
  NUM_Re=integrate(integrand_num_Re,LowerLim,UpperLim)$value
  NUM_Im=integrate(integrand_num_Im,LowerLim,UpperLim)$value
  #eta1=log(DEN1)
  CF=(NUM_Re+1i*NUM_Im)/DEN
  
  return( CF ) 
  
}

JointCF_v2<-function(t1,t2,theta,tol=1e-8){
  #E(exp(i*t1*Y1+i*t2*Y2) for the model of the ANOR paper
  #t1 and t2 can be n-vectors
  #The output is an n-vector with components E(exp(1i*t1[i]*Y1+1i*t2[i]*Y2)
  #i=1,...,n
  #tol: integral lower and upper limits
  # are set to qnorm(tol,sd=stheta2) and qnorm(1-tol,sd=stheta2)
  # with stheta2=1/sqrt(theta1)=1/sqrt(theta2)
  theta4=theta[1]
  theta5=theta[2]
  theta6=theta[3]
  rho=theta[4]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  stheta1<-stheta2<-1/sqrt(theta1)
  
  integrand_den1=function(y2){
    #returns f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1)
  }
  
  
  integrand_den2=function(y1){
    #returns f(y1)*exp(eta)
    coef2=theta2+2*theta6*y1^2-2*theta4*y1
    sigma2_21=1/coef2
    sigma_21=sqrt(sigma2_21)
    mu_21=sigma2_21*(theta3*y1+theta5*y1^2)
    #ee2=exp(.5*(theta3*y1+theta5*y1^2)^2*sigma2_21 )
    ee2=exp(.5*mu_21^2/sigma2_21 )
    out2=2*pi*stheta1*ee2*sigma_21*dnorm(y1,sd=stheta1)
    return(out2)
  }
  
  
  integrand_num_Re=function(y2,t1,t2){
    #returns exp(t2*y2+ t1*mu_12(y2))+t1^2/2*sigma2_12(y2) )*f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    cc=t2*y2+t1*mu_12
    out= cos(cc)*exp(-t1^2/2*sigma2_12 )*f2
    return(out)
  }
  
  integrand_num_Im=function(y2,t1,t2){
    #returns exp(t2*y2+ t1*mu_12(y2))+t1^2/2*sigma2_12(y2) )*f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    cc=t2*y2+t1*mu_12
    out= sin(cc)*exp(-t1^2/2*sigma2_12 )*f2
    return(out)
  }
  
  LowerLim=qnorm(tol,sd=stheta2)
  UpperLim=qnorm(1-tol,sd=stheta2)
  NUM_Re=function(t1,t2){
    integrate(integrand_num_Re,LowerLim,UpperLim,t1,t2)$value
  }
  NUM_Im=function(t1,t2){
    integrate(integrand_num_Im,LowerLim,UpperLim,t1,t2)$value
  }
  DEN=integrate(integrand_den1,LowerLim,UpperLim)$value#exp(eta)
  #eta1=log(DEN1)
  CF=(mapply(NUM_Re,t1,t2)+1i*mapply(NUM_Im,t1,t2))/DEN
  
  
  return( CF ) 
  
}

JointCF_v3<-function(t1,t2,theta,tol=1e-8){
  #E(exp(i*t1*Y1+i*t2*Y2) for the model of the ANOR paper
  #t1 and t2 can be n-vectors and m-vectors
  #The output is an nxm-matrix with components E(exp(1i*t1[i]*Y1+1i*t2[j]*Y2)
  #i=1,...,n 
  #j=1,...,m
  #tol: integral lower and upper limits
  # are set to qnorm(tol,sd=stheta2) and qnorm(1-tol,sd=stheta2)
  # with stheta2=1/sqrt(theta1)=1/sqrt(theta2)
  theta4=theta[1]
  theta5=theta[2]
  theta6=theta[3]
  rho=theta[4]
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  
  stheta1<-stheta2<-1/sqrt(theta1)
  
  integrand_den1=function(y2){
    #returns f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1)
  }
  
  
  integrand_den2=function(y1){
    #returns f(y1)*exp(eta)
    coef2=theta2+2*theta6*y1^2-2*theta4*y1
    sigma2_21=1/coef2
    sigma_21=sqrt(sigma2_21)
    mu_21=sigma2_21*(theta3*y1+theta5*y1^2)
    #ee2=exp(.5*(theta3*y1+theta5*y1^2)^2*sigma2_21 )
    ee2=exp(.5*mu_21^2/sigma2_21 )
    out2=2*pi*stheta1*ee2*sigma_21*dnorm(y1,sd=stheta1)
    return(out2)
  }
  
  
  integrand_num_Re=function(y2,t1,t2){
    #returns exp(t2*y2+ t1*mu_12(y2))+t1^2/2*sigma2_12(y2) )*f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    cc=t2*y2+t1*mu_12
    out= cos(cc)*exp(-t1^2/2*sigma2_12 )*f2
    return(out)
  }
  
  integrand_num_Im=function(y2,t1,t2){
    #returns exp(t2*y2+ t1*mu_12(y2))+t1^2/2*sigma2_12(y2) )*f(y2)*exp(eta)
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    f2=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    cc=t2*y2+t1*mu_12
    out= sin(cc)*exp(-t1^2/2*sigma2_12 )*f2
    return(out)
  }
  
  LowerLim=qnorm(tol,sd=stheta2)
  UpperLim=qnorm(1-tol,sd=stheta2)
  NUM_Re=function(t1,t2){
    integrate(integrand_num_Re,LowerLim,UpperLim,t1,t2)$value
  }
  NUM_Im=function(t1,t2){
    integrate(integrand_num_Im,LowerLim,UpperLim,t1,t2)$value
  }
  DEN=integrate(integrand_den1,LowerLim,UpperLim)$value#exp(eta)
  #eta1=log(DEN1)
  CF_Re=sapply(t2, function(y) mapply(NUM_Re,t1,y))/DEN
  CF_Im=sapply(t2, function(y) mapply(NUM_Im,t1,y))/DEN
  CF=CF_Re+1i*CF_Im
  
  return( CF ) 
  
}


Model_Sim<-function(theta,T){
  #MC Simulations for the model of the ANOR paper
  # rr<- 1-rho^2
  # theta1<-theta2<-1/rr
  # theta3<-rho*theta1
  
  theta4=theta[1]
  theta5=theta[2]
  theta6=theta[3]
  rho=theta[4]
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



Conditional_ES_Model_eq1<-function(rho,theta4,theta5,theta6,alpha1,alpha2){
  #Returns E(Y1|Y1<z1,Y2<z2)
  
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  pars=c(theta1,theta2,theta3,theta4,theta5,theta6)
  fun0 = function(y1,y2,pars) { 
    theta1=pars[1]
    theta2=pars[2]
    theta3=pars[3]
    theta4=pars[4]
    theta5=pars[5]
    theta6=pars[6]
    y12=y1^2
    y22=y2^2
    out=-.5*theta1*y12-.5*theta1*y22+theta3*y1*y2+theta4*y1*y22+theta5*y12*y2-theta6*y12*y22
    return(exp(out))
  }
  
  eta_exp=integral2(fun0, -5, 5, -5, 5, sector = FALSE, reltol = 1e-6, abstol = 0, maxlist = 5000, singular = FALSE, vectorized = TRUE, pars)$Q
  
  
  dens = function(y1,y2,pars,eta_exp) { 
    theta1=pars[1]
    theta2=pars[2]
    theta3=pars[3]
    theta4=pars[4]
    theta5=pars[5]
    theta6=pars[6]
    y12=y1^2
    y22=y2^2
    out=-.5*theta1*y12-.5*theta1*y22+theta3*y1*y2+theta4*y1*y22+theta5*y12*y2-theta6*y12*y22
    return(exp(out)/eta_exp)
  }
  
  dens_times_y1 = function(y1,y2,pars,eta_exp) { 
    theta1=pars[1]
    theta2=pars[2]
    theta3=pars[3]
    theta4=pars[4]
    theta5=pars[5]
    theta6=pars[6]
    y12=y1^2
    y22=y2^2
    out=-.5*theta1*y12-.5*theta1*y22+theta3*y1*y2+theta4*y1*y22+theta5*y12*y2-theta6*y12*y22
    return(y1*exp(out)/eta_exp)
  }
  
  
  integrand_den1=function(y1,eta_exp){
    coef2=theta2+2*theta6*y1^2-2*theta4*y1
    sigma2_21=1/coef2
    sigma_21=sqrt(sigma2_21)
    mu_21=sigma2_21*(theta3*y1+theta5*y1^2)
    #ee2=exp(.5*(theta3*y1+theta5*y1^2)^2*sigma2_21 )
    ee2=exp(.5*mu_21^2/sigma2_21 )
    stheta1=1/sqrt(theta1)
    out2=2*pi*stheta1*ee2*sigma_21*dnorm(y1,sd=stheta1)
    return(out2/eta_exp)
  }
  
  cdf1=function(z1){
    I=integrate(integrand_den1,-5,z1,eta_exp)$value
    return(I)
  }
  
  obj1=function(z1,alpha1){
    return(cdf1(z1)-alpha1)
  }
  
  z1=uniroot(obj1,c(-5,5),alpha1)$root
  
  
  integrand_den2=function(y2,eta_exp){
    #marg. density of Y2
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    stheta2=1/sqrt(theta2)
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1/eta_exp)
  }
  
  cdf2=function(z2){
    I=integrate(integrand_den2,-5,z2,eta_exp)$value
    return(I)
  }
  
  obj2=function(z2,alpha2){
    return(cdf2(z2)-alpha2)
  }
  
  z2=uniroot(obj2,c(-5,5),alpha2)$root
  
  
  NUM=integral2(dens_times_y1, -5, z1, -5, z2, sector = FALSE, reltol = 1e-6, abstol = 0, maxlist = 5000, singular = FALSE, vectorized = TRUE, pars,eta_exp)$Q  
  DEN=integral2(dens, -5, z1, -5, z2, sector = FALSE, reltol = 1e-6, abstol = 0, maxlist = 5000, singular = FALSE, vectorized = TRUE, pars,eta_exp)$Q  
  return(-NUM/DEN)  
}


Conditional_SecondMoment_Model_eq1<-function(rho,theta4,theta5,theta6,alpha1,alpha2){
  #Returns E(Y1^2|Y1<z1,Y2<z2)
  
  theta1<-theta2<-1/(1-rho^2)
  theta3<-rho*theta1
  pars=c(theta1,theta2,theta3,theta4,theta5,theta6)
  fun0 = function(y1,y2,pars) { 
    theta1=pars[1]
    theta2=pars[2]
    theta3=pars[3]
    theta4=pars[4]
    theta5=pars[5]
    theta6=pars[6]
    y12=y1^2
    y22=y2^2
    out=-.5*theta1*y12-.5*theta1*y22+theta3*y1*y2+theta4*y1*y22+theta5*y12*y2-theta6*y12*y22
    return(exp(out))
  }
  
  eta_exp=integral2(fun0, -5, 5, -5, 5, sector = FALSE, reltol = 1e-6, abstol = 0, maxlist = 5000, singular = FALSE, vectorized = TRUE, pars)$Q
  
  
  dens = function(y1,y2,pars,eta_exp) { 
    theta1=pars[1]
    theta2=pars[2]
    theta3=pars[3]
    theta4=pars[4]
    theta5=pars[5]
    theta6=pars[6]
    y12=y1^2
    y22=y2^2
    out=-.5*theta1*y12-.5*theta1*y22+theta3*y1*y2+theta4*y1*y22+theta5*y12*y2-theta6*y12*y22
    return(exp(out)/eta_exp)
  }
  
  dens_times_y1_squared = function(y1,y2,pars,eta_exp) { 
    theta1=pars[1]
    theta2=pars[2]
    theta3=pars[3]
    theta4=pars[4]
    theta5=pars[5]
    theta6=pars[6]
    y12=y1^2
    y22=y2^2
    out=-.5*theta1*y12-.5*theta1*y22+theta3*y1*y2+theta4*y1*y22+theta5*y12*y2-theta6*y12*y22
    return(y1^2*exp(out)/eta_exp)
  }
  
  
  integrand_den1=function(y1,eta_exp){
    coef2=theta2+2*theta6*y1^2-2*theta4*y1
    sigma2_21=1/coef2
    sigma_21=sqrt(sigma2_21)
    mu_21=sigma2_21*(theta3*y1+theta5*y1^2)
    #ee2=exp(.5*(theta3*y1+theta5*y1^2)^2*sigma2_21 )
    ee2=exp(.5*mu_21^2/sigma2_21 )
    stheta1=1/sqrt(theta1)
    out2=2*pi*stheta1*ee2*sigma_21*dnorm(y1,sd=stheta1)
    return(out2/eta_exp)
  }
  
  cdf1=function(z1){
    I=integrate(integrand_den1,-5,z1,eta_exp)$value
    return(I)
  }
  
  obj1=function(z1,alpha1){
    return(cdf1(z1)-alpha1)
  }
  
  z1=uniroot(obj1,c(-5,5),alpha1)$root
  
  
  integrand_den2=function(y2,eta_exp){
    #marg. density of Y2
    coef1=theta1+2*theta6*y2^2-2*theta5*y2
    sigma2_12=1/coef1
    sigma_12=sqrt(sigma2_12)
    mu_12=sigma2_12*(theta3*y2+theta4*y2^2)
    #ee1=exp(.5*(theta3*y2+theta4*y2^2)^2*sigma2_12 )
    ee1=exp(.5*mu_12^2/sigma2_12 )
    stheta2=1/sqrt(theta2)
    out1=2*pi*stheta2*ee1*sigma_12*dnorm(y2,sd=stheta2)
    return(out1/eta_exp)
  }
  
  cdf2=function(z2){
    I=integrate(integrand_den2,-5,z2,eta_exp)$value
    return(I)
  }
  
  obj2=function(z2,alpha2){
    return(cdf2(z2)-alpha2)
  }
  
  z2=uniroot(obj2,c(-5,5),alpha2)$root
  
  
  NUM=integral2(dens_times_y1_squared, -5, z1, -5, z2, sector = FALSE, reltol = 1e-6, abstol = 0, maxlist = 5000, singular = FALSE, vectorized = TRUE, pars,eta_exp)$Q  
  DEN=integral2(dens, -5, z1, -5, z2, sector = FALSE, reltol = 1e-6, abstol = 0, maxlist = 5000, singular = FALSE, vectorized = TRUE, pars,eta_exp)$Q  
  return(NUM/DEN)  
}

Conditional_ES_Model_eq1_sim<-function(rho,theta4,theta5,theta6,alpha1,alpha2,T){
  Y=Model_eq1_Sim(rho,theta4,theta5,theta6,T)
  
}






VaR_ES_Model_eq1<-function(rho,theta4,theta5,theta6,omega,alpha,T,nsims){
  #omega: ptf weight, R=omega Y1 + (1-omega) Y2
  #alpha: VaR/ES level
  
  f0<-function(rho,theta4,theta5,theta6,omega,alpha,T){
    YY=Model_eq1_Sim(rho,theta4,theta5,theta6,T)
    R=omega*YY[,1]+(1-omega)*YY[,2]
    #Z1=(YY[,1]-mean(YY[,1]))/sd(YY[,1])
    #Z2=(YY[,2]-mean(YY[,2]))/sd(YY[,2])
    #R=omega*Z1+(1-omega)*Z2
    VaR=-quantile(R,alpha)
    ES=-mean(R[R<=-VaR])
    return(c(VaR,ES))
  }
  
  return(replicate(nsims, 
                   expr = f0(rho,theta4,theta5,theta6,omega,alpha,T), 
                   simplify = TRUE ))
}


VaR_ES_Model_eq1_v2<-function(sigma1,sigma2,rho,theta4,theta5,theta6,omega,alpha,T,nsims){
  #omega: ptf weight, R=omega Y1 + (1-omega) Y2
  #alpha: VaR/ES level
  
  f0<-function(sigma1,sigma2,rho,theta4,theta5,theta6,omega,alpha,T){
    YY=Model_eq1_Sim(rho,theta4,theta5,theta6,T)
    R=omega*sigma1*YY[,1]+(1-omega)*sigma2*YY[,2]
    #Z1=(YY[,1]-mean(YY[,1]))/sd(YY[,1])
    #Z2=(YY[,2]-mean(YY[,2]))/sd(YY[,2])
    #R=omega*Z1+(1-omega)*Z2
    VaR=-quantile(R,alpha)
    ES=-mean(R[R<=-VaR])
    return(c(VaR,ES))
  }
  
  return(replicate(nsims, 
                   expr = f0(sigma1,sigma2,rho,theta4,theta5,theta6,omega,alpha,T), 
                   simplify = TRUE ))
}
