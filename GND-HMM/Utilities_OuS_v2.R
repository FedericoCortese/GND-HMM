options(warn=-1)

garch_fit<-function(data,ugarch.spec){
  
  require(nloptr)
  require(Rsolnp)
  require(rugarch)
  T=dim(data)[1]
  k=dim(data)[2]
  U <- data
  D <- diag(k)
  mfit <- vector(mode="list",length=k)
  mu_fcst<-sigma_fcst<-rep(0,k)
  for(i in 1:k){
    fit = ugarchfit(ugarch.spec,data[,i],solver="hybrid")
    counterA <- 0
    while((fit@fit$convergence!=0) && (counterA<10)){
      print("ugarchfit failed")
      fit=ugarchfit(ugarch.spec,data[,i],solver="hybrid")
      counterA <- counterA+1
    }
    U[,i] <- residuals(fit,standardize=TRUE)
    mfit[[i]] <- fit
    fcst<-ugarchforecast(fit,n.ahead = 1)
    mu_fcst[i]<-fcst@forecast[["seriesFor"]]
    sigma_fcst[i]<-fcst@forecast[["sigmaFor"]]
  }
  
  names(mfit) <- names(data)
  
  U=as.matrix(U)
  
  
  return(list(  U=U, mu_fcst=mu_fcst, sigma_fcst=sigma_fcst ))
}


DCC_fit<-function(data,ugarch.spec){
  require(rmgarch)
  mspec = multispec( replicate(2, ugarch.spec) )
  model_dcc=dccspec(mspec)
  fit.dcc=dccfit(model_dcc,data)
  Rho0=rcor(fit.dcc)
  Rho=rep(0,dim(Rho0)[3])
  for(t in 1:dim(Rho0)[3]){Rho[t]=Rho0[1,2,t]}
  
  U<-garch_fit(data,ugarch.spec)$U
  
  biv_norm_dens<-function(y1,y2,rho,log=FALSE){
    rho2=1-rho^2
    out=-log(2*pi)-.5*log(rho2)-.5*(y1^2-2*rho*y1*y2+y2^2)/rho2
    if(log==FALSE) out=exp(out)
    
    return(out)
  }
  
  ll=sum(biv_norm_dens(U[,1],U[,2],Rho,log=TRUE))
  
  AIC=2*2-2*ll
  BIC=log(dim(U)[1])*2-2*ll
  return(list(LogLik=ll,AIC=AIC,BIC=BIC))
}

Biv_DCC_normal_GAS_unrescaled_score<-function(U,pars0=NULL,X=NULL,theta.hat=NULL){
  #U: bivariate time series of garch residuals
  #pars0 optional argument for the initial parameter values
  #X (lagged) explanatory variables
  #if theta.hat=NULL the model is estimated
  #otherwise the theta.hat provided is used for filtering and prediction
 
  T=dim(U)[1]
  biv_norm_dens<-function(y1,y2,rho,log=FALSE){
    rho2=1-rho^2
    out=-log(2*pi)-.5*log(rho2)-.5*(y1^2-2*rho*y1*y2+y2^2)/rho2
    if(log==FALSE) out=exp(out)
    
    return(out)
  }
  

  score_delta<-function(rho,y1,y2){
    rho2=1-rho^2
    
    out=rho-rho/rho2*(y1^2+y2^2)+(1+rho^2)*y1*y2/rho2
    return(out)
  }
  
  #Y1=data[,1]
  #Y2=data[,2]
  
  obj_fun<-function(pars,U){
    omega=pars[1]
    a=pars[2]
    b=pars[3]
    # nu=2+exp(pars[3])
    # lxi=pars[seq(from=4,length.out=k)]
    # xi=exp(lxi)
    
    
    npars=3 
    XB=rep(0,T)
    if(!is.null(X)){
      kk=dim(X)[2]
      Beta=pars[seq(from=(npars+1),length.out=kk)]
      XB=X%*% matrix(Beta,ncol=1)
    }
    
    
    
    
    t=1
    delta=omega
    rho=tanh(delta)
    logL <-  biv_norm_dens(U[t,1],U[t,2],rho,log=TRUE)
    
    
    for(t in 2:T){
      
      delta=omega+a*score_delta(rho,U[t-1,1],U[t-1,2])+b*(delta-omega)+XB[t-1]
      rho=tanh(delta)
      logL <-  logL+biv_norm_dens(U[t,1],U[t,2],rho,log=TRUE)
    }
    
    
    
    return(-as.numeric(logL))
    
    
  }
  
  
  eval_g_ineq<-function(pars,U,U_bar){
    omega=pars[1]
    a=pars[2]
    b=pars[3]
    return( abs(b)-1  )
  }
  
  opts = list("algorithm"="NLOPT_LN_COBYLA",
              "xtol_rel"=1.0e-15)
  
  lb <- c(-Inf,-Inf,-1)
  ub <- c(Inf,Inf,1)
  
  if(!is.null(X)){
    kk=dim(X)[2]
    lb <- c(lb,rep(-Inf,kk))
    ub <- c(ub,rep(Inf,kk))
  }
  
  #pars0=c(.01,.8,6,rep(1.2,k))
  
if(is.null(theta.hat)){
  
  if(is.null(pars0)){
    #omega=0
    omega=atanh(cor(U[,1],U[,2]))
    a0=.05
    b0=.9
    pars0=c(omega,a0,b0)
    
    if(!is.null(X)){
      kk=dim(X)[2]
      pars0=c(pars0,rep(0,kk))
    }
  }
  
  obj=function(pars){return(obj_fun(pars,U))}
  
  fit= solnp(pars0, obj, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL, 
             ineqUB = NULL, LB = lb, UB = ub,control=list(trace=0) )
  
  
  
  theta.hat=fit$pars

 
  ll=-fit$value[length(fit$value)]
  AIC=2*length(theta.hat)-2*ll
  BIC=log(dim(U)[1])*length(theta.hat)-2*ll
}
  
  
  rho_filtered<-function(pars,U,X){
    omega=pars[1]
    a=pars[2]
    b=pars[3]
    # nu=2+exp(pars[3])
    # lxi=pars[seq(from=4,length.out=k)]
    # xi=exp(lxi)
    
    T=dim(U)[1]
    npars=3 
    XB=rep(0,T)
    if(!is.null(X)){
      kk=dim(X)[2]
      Beta=pars[seq(from=(npars+1),length.out=kk)]
      XB=X%*% matrix(Beta,ncol=1)
    }
    
    
    Rho=rep(0,T)
    
    t=1
    delta=omega
    rho=tanh(delta)
    Rho[t]<-rho
    
    
    for(t in 2:T){
      delta=omega+a*score_delta(rho,U[t-1,1],U[t-1,2])+b*(delta-omega)+XB[t-1]
      rho=tanh(delta)
      Rho[t]<-rho
    }
    
    t=T+1
    delta=omega+a*score_delta(rho,U[t-1,1],U[t-1,2])+b*(delta-omega)+XB[t-1]
    rho=tanh(delta)
    Rho_fcst<-rho
    
    
    
    return(list(Rho=Rho,Rho_fcst=Rho_fcst))
    
    
  }
  
  FILT=rho_filtered(theta.hat,U,X)
  
  if(is.null(theta.hat)){return(list(theta.hat=theta.hat,LogLik=ll,AIC=AIC,BIC=BIC,
              Rho_fcst=FILT$Rho_fcst))}
  else{ return(list(theta.hat=theta.hat,
                Rho_fcst=FILT$Rho_fcst))}
  
  
  
}



Biv_DCC_GAS_t_unrescaled_score<-function(U,pars0=NULL,X=NULL,theta.hat=NULL){
  #U: bivariate time series of garch residuals
  #pars0 optional argument for the initial parameter values
  #X (lagged) explanatory variables
  #if theta.hat=NULL the model is estimated
  #otherwise the theta.hat provided is used for filtering and prediction
  
  
  T=dim(U)[1]
  
  
  biv_t_dens<-function(y1,y2,nu,rho,log=FALSE){
    nu_half=nu/2
    nu2=nu_half+1#(nu+2)/2
    rho2=1-rho^2
    out=log(nu)-log(nu-2)-log(2*pi)-.5*log(rho2)-nu2*log(1+1/(nu-2)*(y1^2-2*rho*y1*y2+y2^2)/rho2  )
    if(log==FALSE) out=exp(out)
    
    return(out)
  }
  
  
  
  
  score_delta<-function(y1,y2,nu,rho){
    nu_half=nu/2
    nu2=nu_half+1#(nu+2)/2
    rho2=1-rho^2
    out=rho -nu2*(  2*rho*(y1^2+y2^2)  -2*(1+rho^2)*y1*y2  )/((nu-2)*rho2+ (y1^2-2*rho*y1*y2+y2^2) )
    return(out)
  }
  
  
  fit_static_t=function(data){
    #data: bivariate time series 
    #mean 0 and var. 1 fo both components
    
    obj=function(pars){
      
      nu=pars[1]
      rho=tanh(pars[2])
      y1=data[,1]
      y2=data[,2]
      
      nu_half=nu/2
      nu2=nu_half+1#(nu+2)/2
      rho2=1-rho^2
      out=lgamma(nu2)-lgamma(nu_half)-log(nu-2)-log(pi)-.5*log(rho2)-nu2*log(1+1/(nu-2)*(y1^2-2*rho*y1*y2+y2^2)/rho2  )
      return(-sum(out))
    }
    
    lb <- c(2,-Inf)
    ub <- c(Inf,Inf)
    
    pars0=c(8,0)
    
    fit= solnp(pars0, obj, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL, 
               ineqUB = NULL, LB = lb, UB = ub  ,control=list(trace=0))
    
    return(fit)
    
  }
  
  
  obj_fun<-function(pars,U){
    nu=pars[1]
    omega=pars[2]
    a=pars[3]
    b=pars[4]
    # nu=2+exp(pars[3])
    # lxi=pars[seq(from=4,length.out=k)]
    # xi=exp(lxi)
    
    
    npars=4
    XB=rep(0,T)
    if(!is.null(X)){
      kk=dim(X)[2]
      Beta=pars[seq(from=(npars+1),length.out=kk)]
      XB=X%*% matrix(Beta,ncol=1)
    }
    
    
    
    
    t=1
    delta=omega
    rho=tanh(delta)
    logL <-  biv_t_dens(U[t,1],U[t,2],nu,rho,log=TRUE)
    
    for(t in 2:T){
      
      delta=omega+a*score_delta(U[t-1,1],U[t-1,2],nu,rho)+b*(delta-omega)+XB[t-1]
      rho=tanh(delta)
      logL <-  logL+biv_t_dens(U[t,1],U[t,2],nu,rho,log=TRUE)
    }
    
    
    
    return(-as.numeric(logL))
    
    
  }
  
  
  
  lb <- c(2,-Inf,-Inf,-1)
  ub <- c(Inf,Inf,Inf,1)
  
  if(!is.null(X)){
    kk=dim(X)[2]
    lb <- c(lb,rep(-Inf,kk))
    ub <- c(ub,rep(Inf,kk))
  }
  
  #pars0=c(.01,.8,6,rep(1.2,k))

if(is.null(theta.hat)){  
  if(is.null(pars0)){
    #omega=0
    
    #omega=atanh(cor(U[,1],U[,2]))
    
    fit0=fit_static_t(U)
    nu0=fit0$pars[1]
    #omega0=atanh(fit0$pars[2])
    omega0=fit0$pars[2]
    a0=.05
    b0=.9
    pars0=c(nu0,omega0,a0,b0)
    
    if(!is.null(X)){
      kk=dim(X)[2]
      pars0=c(pars0,rep(0,kk))
    }
  }
  
  obj=function(pars){return(obj_fun(pars,U))}
  
  fit= solnp(pars0, obj, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL, 
             ineqUB = NULL, LB = lb, UB = ub ,control=list(trace=0))
  
  
  
  theta.hat=fit$pars
  
  ll=-fit$value[length(fit$value)]
  AIC=2*length(theta.hat)-2*ll
  BIC=log(dim(U)[1])*length(theta.hat)-2*ll
}  
  
  
  rho_filtered<-function(pars,U,X){
    nu=pars[1]
    omega=pars[2]
    a=pars[3]
    b=pars[4]
    # nu=2+exp(pars[3])
    # lxi=pars[seq(from=4,length.out=k)]
    # xi=exp(lxi)
    
    T=dim(U)[1]
    npars=4
    XB=rep(0,T)
    if(!is.null(X)){
      kk=dim(X)[2]
      Beta=pars[seq(from=(npars+1),length.out=kk)]
      XB=X%*% matrix(Beta,ncol=1)
    }
    
    
    Rho=rep(0,T)
    
    t=1
    delta=omega
    rho=tanh(delta)
    Rho[t]<-rho
    
    
    for(t in 2:T){
      delta=omega+a*score_delta(U[t-1,1],U[t-1,2],nu,rho)+b*(delta-omega)+XB[t-1]
      rho=tanh(delta)
      Rho[t]<-rho
    }
    
    t=T+1
    delta=omega+a*score_delta(U[t-1,1],U[t-1,2],nu,rho)+b*(delta-omega)+XB[t-1]
    rho=tanh(delta)
    Rho_fcst<-rho
    
    
    
    return(list(Rho=Rho,Rho_fcst=Rho_fcst))
    
    
  }
  
  FILT=rho_filtered(theta.hat,U,X)

if(is.null(theta.hat)){return(list(theta.hat=theta.hat,LogLik=ll,AIC=AIC,BIC=BIC,
              nu_hat=theta.hat[1],Rho_fcst=FILT$Rho_fcst))}
  else{return(list(theta.hat=theta.hat,
                   nu_hat=theta.hat[1],Rho_fcst=FILT$Rho_fcst))}
  
  
}



Biv_DCC_GAS_skew_normal_unrescaled_score<-function(U,pars0=NULL,X=NULL,theta.hat=NULL){
  #U: bivariate time series of garch residuals
  #pars0 optional argument for the initial parameter values
  #X (lagged) explanatory variables
  #if theta.hat=NULL the model is estimated
  #otherwise the theta.hat provided is used for filtering and prediction
  
  
  T=dim(U)[1]
  U1=U[,1]
  U2=U[,2]
 
  biv_std_esnorm_dens<-function(y1,y2,rho,alpha1,alpha2,log=FALSE){
    
    biv_norm_dens<-function(y1,y2,rho,log=FALSE){
      #biv normal
      rho2=1-rho^2
      out=-log(2*pi)-.5*log(rho2)-.5*(y1^2-2*rho*y1*y2+y2^2)/rho2
      if(log==FALSE) out=exp(out)
      
      return(out)
    }
    
    alpha=c(alpha1,alpha2)
    
    alpha2_star=alpha1^2+2*rho*alpha1*alpha2+alpha2^2
    delta=c(alpha1+rho*alpha2,rho*alpha1+alpha2 )/sqrt(1+alpha2_star)
    #tau=0
    #zz1=sn::zeta(1,tau)
    zz1= sqrt(2/pi)
    #zz2=sn::zeta(2,tau)
    zz2=-zz1^2
    
    a1=zz1*delta[1]
    a2=zz1*delta[2]
    
    b1=sqrt(1+zz2*delta[1]^2)
    b2=sqrt(1+zz2*delta[2]^2)
    
    x1=a1+y1*b1
    x2=a2+y2*b2
    
    
    yy=alpha1*x1+alpha2*x2
    #out=biv_norm_dens(x1,x2,rho,log=TRUE)+pnorm(yy,log.p = TRUE)-log(.5)+log(b1)+log(b2)
    
    rho2=1-rho^2
    out=-log(pi)-.5*log(rho2)-.5*(x1^2-2*rho*x1*x2+x2^2)/rho2+pnorm(yy,log.p = TRUE)+log(b1)+log(b2)
    
    
    if(log==FALSE) out=exp(out)
    
    return(out)
  }
  
  
  fit_biv_esnorm<-function(data,hessian=TRUE){
    #Static model
    #data: bivariate time series 
    #mean 0 and var. 1 fo both components
    obj=function(pars){
      
      
      rho=tanh(pars[1])
      #rho=pars[1]
      alpha1=pars[2]
      alpha2=pars[3]
      
      
      y1=data[,1]
      y2=data[,2]
      
      
      out=biv_std_esnorm_dens(y1,y2,rho,alpha1,alpha2,log=TRUE)
      return(-sum(out))
    }
    
    
    
    rho0=cor(data)[1,2]
    alpha1_0=mean(scale(data[,1])^3)
    alpha2_0=mean(scale(data[,2])^3)
    pars0=c(atanh(rho0),alpha1_0,alpha2_0)
    #pars0=c(atanh(rho0),rep(0,2))
    # LB=c(-1,rep(-Inf,2))
    # UB=c(1,rep(Inf,2))
    # fit= solnp(pars0, obj, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL, 
    #            ineqUB = NULL, LB = LB, UB = UB)
    # theta.hat=fit$pars
    
    fit=optim(pars0,obj)
    theta.hat=fit$par
    
    if(hessian){H=optimHess(theta.hat, obj, gr = NULL)}
    
    
    rnam=c("atanh_rho","alpha1","alpha2")
    
    
    
    
    
    #H=fit$hessian
    
    if(hessian){
      H2=solve(H)
      SD=sqrt(diag(H2))
      
      matcoef=matrix(rep(0,length(theta.hat)*4),ncol=4)
      matcoef[,1]=theta.hat
      matcoef[,2]=SD
      matcoef[,3]=theta.hat/SD
      matcoef[,4]=2-2*pnorm(abs(matcoef[,3]))
      dimnames(matcoef)[[2]]=c(" Estimate"," Std. Error"," t value", "Pr(>|t|)" )
      row.names(matcoef) =rnam
      
      ll=-fit$value[length(fit$value)]
      AIC=2*length(theta.hat)-2*ll
      BIC=log(dim(data)[1])*length(theta.hat)-2*ll
      
      
      return(list(fit=fit,
                  conv=fit$convergence,
                  LogLik=ll,AIC=AIC,BIC=BIC,
                  par.hat=theta.hat,
                  SD=SD,
                  matcoef=matcoef) )
    }else{  
      ll=-fit$value[length(fit$value)]
      AIC=2*length(theta.hat)-2*ll
      BIC=log(dim(data)[1])*length(theta.hat)-2*ll
      
      
      return(list(fit=fit,
                  conv=fit$convergence,
                  LogLik=ll,AIC=AIC,BIC=BIC,
                  par.hat=theta.hat) )
    }
    
  }
  
  
  score_rho<-function(y1,y2,rho,alpha1,alpha2){
    
    
    alpha=c(alpha1,alpha2)
    
    alpha2_star=alpha1^2+2*rho*alpha1*alpha2+alpha2^2
    delta=c(alpha1+rho*alpha2,rho*alpha1+alpha2 )/sqrt(1+alpha2_star)
    
    #tau=0
    #zz1=sn::zeta(1,tau)
    zz1= sqrt(2/pi)
    #zz2=sn::zeta(2,tau)
    zz2=-zz1^2
    
    a1=zz1*delta[1]
    a2=zz1*delta[2]
    
    b1=sqrt(1+zz2*delta[1]^2)
    b2=sqrt(1+zz2*delta[2]^2)
    
    x1=a1+y1*b1
    x2=a2+y2*b2
    
    alpha2_star_1=1+alpha2_star
    alpha2_star_1_sqrt=sqrt(alpha2_star_1)
    delta1_prime=alpha2/alpha2_star_1_sqrt -alpha1*alpha2*delta[1]/alpha2_star_1
    delta2_prime=alpha1/alpha2_star_1_sqrt -alpha1*alpha2*delta[2]/alpha2_star_1
    
    a1_prime=zz1*delta1_prime
    a2_prime=zz1*delta2_prime
    
    b1_prime=-a1_prime*a1/b1
    b2_prime=-a2_prime*a2/b2
    
    x1_prime=a1_prime+y1*b1_prime
    x2_prime=a2_prime+y2*b2_prime
    x12_prime=x1_prime*x2+x1*x2_prime
    
    rho2=1-rho^2
    
    
    #r.1pasq <- sqrt(1 + alpha1^2 + 2 * alpha1 * alpha2 * rho + alpha2^2)
    t <- alpha1 * x1 + alpha2 * x2
    t_prime<- alpha1 * x1_prime + alpha2 * x2_prime
    s.rho= rho/rho2+b1_prime/b1+b2_prime/b2-1/(2*rho2)*( 2*x1*x1_prime+2*x2*x2_prime-2*x1*x2-2*rho*x12_prime ) - rho/rho2^2*(x1^2-2*rho*x1*x2+x2^2) + t_prime*dnorm(t)/pnorm(t)
    return(s.rho)
  }
  
  score_delta<-function(y1,y2,rho,alpha1,alpha2){
    s.rho=score_rho(y1,y2,rho,alpha1,alpha2)
    s.delta=s.rho*(1-rho^2)
    return(s.delta)
  }
  
  obj<-function(pars){
    alpha1=pars[1]
    alpha2=pars[2]
    omega=pars[3]
    a=pars[4]
    b=pars[5]
    
    
    npars=5
    XB=rep(0,T)
    if(!is.null(X)){
      kk=dim(X)[2]
      Beta=pars[seq(from=(npars+1),length.out=kk)]
      XB=X%*% matrix(Beta,ncol=1)
    }
    
    
    
    
    t=1
    delta=omega
    rho=tanh(delta)
    logL <-  biv_std_esnorm_dens(U1[t],U2[t],rho,alpha1,alpha2,log=TRUE)
    
    
    for(t in 2:T){
      
      delta=omega+a*score_delta(U1[t-1],U2[t-1],rho,alpha1,alpha2)+b*(delta-omega)+XB[t-1]
      rho=tanh(delta)
      logL <-  logL+biv_std_esnorm_dens(U1[t],U2[t],rho,alpha1,alpha2,log=TRUE)
    }
    
    
    
    return(-as.numeric(logL))
    
    
  }
  
  

if(is.null(theta.hat)){    
  if(is.null(pars0)){
    fit_static<-fit_biv_esnorm(U,hessian=FALSE)
    alpha1_0=fit_static$par.hat[2]
    alpha2_0=fit_static$par.hat[3]
    omega=fit_static$par.hat[1]
    a0=.05
    b0=.9
    pars0=c(alpha1_0,alpha2_0,omega,a0,b0)
    
    if(!is.null(X)){
      kk=dim(X)[2]
      pars0=c(pars0,rep(0,kk))
    }
  }
  
  
  lb=c(rep(-Inf,4),-1)
  ub=c(rep(Inf,4),1)
  if(!is.null(X)){
    kk=dim(X)[2]
    lb <- c(lb,rep(-Inf,kk))
    ub <- c(ub,rep(Inf,kk))
  }
  
  fit= solnp(pars0, obj, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL, 
             ineqUB = NULL, LB = lb, UB = ub,control=list(trace=0))
  
  
  theta.hat=fit$pars
  
  ll=-fit$value[length(fit$value)]
  AIC=2*length(theta.hat)-2*ll
  BIC=log(dim(U)[1])*length(theta.hat)-2*ll
  
  #pars= theta.hat
}  

  
  rho_filtered<-function(pars,U,X){
    alpha1=pars[1]
    alpha2=pars[2]
    omega=pars[3]
    a=pars[4]
    b=pars[5]
    
    ###(rho-zz1^2*delta[1]*delta[2])/(b1*b2)
    ###with zz1= sqrt(2/pi)
    
    T=dim(U)[1]
    npars=5
    XB=rep(0,T)
    if(!is.null(X)){
      kk=dim(X)[2]
      Beta=pars[seq(from=(npars+1),length.out=kk)]
      XB=X%*% matrix(Beta,ncol=1)
    }
    
    
    Rho=rep(0,T)
    
    t=1
    delta=omega
    rho=tanh(delta)
    Rho[t]<-rho
    
    
    
    
    for(t in 2:T){
      delta=omega+a*score_delta(U1[t-1],U2[t-1],rho,alpha1,alpha2)+b*(delta-omega)+XB[t-1]
      rho=tanh(delta)
      Rho[t]<-rho
    }
    
    t=T+1
    delta=omega+a*score_delta(U1[t-1],U2[t-1],rho,alpha1,alpha2)+b*(delta-omega)+XB[t-1]
    rho=tanh(delta)
    Rho_fcst<-rho
    
    
    alpha=c(alpha1,alpha2)
    
    alpha2_star=alpha1^2+2*Rho*alpha1*alpha2+alpha2^2
    
    
    den_delta=sqrt(1+alpha2_star)
    delta1=(alpha1+Rho*alpha2)/den_delta
    delta2=(Rho*alpha1+alpha2 )/den_delta
    #tau=0
    #zz1=sn::zeta(1,tau)
    zz1= sqrt(2/pi)
    #zz2=sn::zeta(2,tau)
    zz2=-zz1^2
    
    
    b1=sqrt(1+zz2*delta1^2)
    b2=sqrt(1+zz2*delta2^2)
    Rho2=(Rho-zz1^2*delta1*delta2)/(b1*b2)#Actual correlation
    
    
    alpha2_star=alpha1^2+2*Rho_fcst*alpha1*alpha2+alpha2^2
    den_delta=sqrt(1+alpha2_star)
    delta1=(alpha1+Rho_fcst*alpha2)/den_delta
    delta2=(Rho_fcst*alpha1+alpha2 )/den_delta
    b1=sqrt(1+zz2*delta1^2)
    b2=sqrt(1+zz2*delta2^2)
    Rho2_fcst=(Rho_fcst-zz1^2*delta1*delta2)/(b1*b2)
    
    return(list(Rho=Rho,Rho_fcst=Rho_fcst,Rho2=Rho2,Rho2_fcst=Rho2_fcst))
    
    
  }
  
  FILT=rho_filtered(theta.hat,U,X)
  
if(is.null(theta.hat)){    
  return(list(theta.hat=theta.hat,LogLik=ll,AIC=AIC,BIC=BIC,
              alpha1_hat=theta.hat[1],alpha2_hat=theta.hat[2],
              Rho_fcst=as.numeric(FILT$Rho_fcst),
              Rho2_fcst=as.numeric(FILT$Rho2_fcst)
  ) )}
  else{  return(list(theta.hat=theta.hat,
                     alpha1_hat=theta.hat[1],alpha2_hat=theta.hat[2],
                     Rho_fcst=as.numeric(FILT$Rho_fcst),
                     Rho2_fcst=as.numeric(FILT$Rho2_fcst)  ) )}
  
  
}




Biv_DCC_GAS_skew_t_unrescaled_score<-function(U,pars0=NULL,X=NULL,theta.hat=NULL){
  #U: bivariate time series of garch residuals
  #pars0 optional argument for the initial parameter values
  #X (lagged) explanatory variables
  #if theta.hat=NULL the model is estimated
  #otherwise the theta.hat provided is used for filtering and prediction
  
  T=dim(U)[1]
  U1=U[,1]
  U2=U[,2]
  
  
  
  biv_std_est_dens<-function(y1,y2,rho,nu,alpha1,alpha2,log=FALSE){
    
    
    
    alpha=c(alpha1,alpha2)
    
    alpha2_star=alpha1^2+2*rho*alpha1*alpha2+alpha2^2
    delta=c(alpha1+rho*alpha2,rho*alpha1+alpha2 )/sqrt(1+alpha2_star)
    #tau=0
    
    zz1= sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
    
    aa=zz1*delta
    a1=aa[1]
    a2=aa[2]
    bb=sqrt(  nu/(nu-2) - aa^2)
    b1=bb[1]
    b2=bb[2]
    
    x1=a1+y1*b1
    x2=a2+y2*b2
    
    rho2=1-rho^2
    Q=(x1^2-2*rho*x1*x2+x2^2)/rho2
    nu2=nu+2
    yy=alpha1*x1+alpha2*x2
    yy2=yy*sqrt( nu2/(nu+Q)   )
    
    
    
    
    out=-log(pi)-.5*log(rho2)-nu2/2* log( 1+ Q/nu )    +pt(yy2,df=nu2,log.p = TRUE)+log(b1)+log(b2)
    
    
    if(log==FALSE) out=exp(out)
    
    return(out)
  }
  
  
  fit_biv_est<-function(data,hessian=TRUE){
    #Static model
    #data: bivariate time series 
    #mean 0 and var. 1 fo both components
    obj=function(pars){
      
      
      rho=tanh(pars[1])
      #rho=pars[1]
      nu=exp(pars[2])+2
      alpha1=pars[3]
      alpha2=pars[4]
      
      
      y1=data[,1]
      y2=data[,2]
      
      
      out=biv_std_est_dens(y1,y2,rho,nu,alpha1,alpha2,log=TRUE)
      return(-sum(out))
    }
    
    
    
    rho0=cor(data)[1,2]
    nu0=8
    alpha1_0=mean(scale(data[,1])^3)
    alpha2_0=mean(scale(data[,2])^3)
    pars0=c(atanh(rho0),log(nu0-2),alpha1_0,alpha2_0)
    #pars0=c(atanh(rho0),rep(0,2))
    # LB=c(-1,rep(-Inf,2))
    # UB=c(1,rep(Inf,2))
    # fit= solnp(pars0, obj, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL, 
    #            ineqUB = NULL, LB = LB, UB = UB)
    # theta.hat=fit$pars
    
    fit=optim(pars0,obj)
    theta.hat=fit$par
    
    if(hessian){H=optimHess(theta.hat, obj, gr = NULL)}
    
    
    rnam=c("atanh_rho","log(nu-2)","alpha1","alpha2")
    
    
    
    
    
    #H=fit$hessian
    
    if(hessian){
      H2=solve(H)
      SD=sqrt(diag(H2))
      
      matcoef=matrix(rep(0,length(theta.hat)*4),ncol=4)
      matcoef[,1]=theta.hat
      matcoef[,2]=SD
      matcoef[,3]=theta.hat/SD
      matcoef[,4]=2-2*pnorm(abs(matcoef[,3]))
      dimnames(matcoef)[[2]]=c(" Estimate"," Std. Error"," t value", "Pr(>|t|)" )
      row.names(matcoef) =rnam
      
      ll=-fit$value[length(fit$value)]
      AIC=2*length(theta.hat)-2*ll
      BIC=log(dim(data)[1])*length(theta.hat)-2*ll
      
      
      return(list(fit=fit,
                  conv=fit$convergence,
                  LogLik=ll,AIC=AIC,BIC=BIC,
                  par.hat=theta.hat,
                  SD=SD,
                  matcoef=matcoef) )
    }else{  
      ll=-fit$value[length(fit$value)]
      AIC=2*length(theta.hat)-2*ll
      BIC=log(dim(data)[1])*length(theta.hat)-2*ll
      
      
      return(list(fit=fit,
                  conv=fit$convergence,
                  LogLik=ll,AIC=AIC,BIC=BIC,
                  par.hat=theta.hat) )
    }
    
  }
  
  
  
  score_rho<-function(y1,y2,rho,nu,alpha1,alpha2){
    
    alpha=c(alpha1,alpha2)
    
    alpha2_star=alpha1^2+2*rho*alpha1*alpha2+alpha2^2
    delta=c(alpha1+rho*alpha2,rho*alpha1+alpha2 )/sqrt(1+alpha2_star)
    #tau=0
    
    zz1= sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
    
    aa=zz1*delta
    a1=aa[1]
    a2=aa[2]
    bb=sqrt(  nu/(nu-2) - aa^2)
    b1=bb[1]
    b2=bb[2]
    
    x1=a1+y1*b1
    x2=a2+y2*b2
    
    rho2=1-rho^2
    Q=(x1^2-2*rho*x1*x2+x2^2)/rho2
    nu2=nu+2
    
    
    
    
    alpha2_star_1=1+alpha2_star
    alpha2_star_1_sqrt=sqrt(alpha2_star_1)
    delta1_prime=alpha2/alpha2_star_1_sqrt -alpha1*alpha2*delta[1]/alpha2_star_1
    delta2_prime=alpha1/alpha2_star_1_sqrt -alpha1*alpha2*delta[2]/alpha2_star_1
    
    a1_prime=zz1*delta1_prime
    a2_prime=zz1*delta2_prime
    
    b1_prime=-a1_prime*a1/b1
    b2_prime=-a2_prime*a2/b2
    
    x1_prime=a1_prime+y1*b1_prime
    x2_prime=a2_prime+y2*b2_prime
    x12_prime=x1_prime*x2+x1*x2_prime
    
    rho2=1-rho^2
    
    
    Q_prime=2/rho2*(  x1*x1_prime+x2*x2_prime-x1*x2-rho*x12_prime +rho*Q     )
    yy=alpha1*x1+alpha2*x2
    yy2=yy*sqrt( nu2/(nu+Q)   )
    yy_prime<- alpha1 * x1_prime + alpha2 * x2_prime
    s.rho= rho/rho2+b1_prime/b1+b2_prime/b2 - (nu2)/(2*nu)*Q_prime/(1+Q/nu)  + sqrt( nu2/(nu+Q)   )*  (yy_prime-yy/2*Q_prime/(nu+Q))*dt(yy2,df=nu2)/pt(yy2,df=nu2)
    return(s.rho)
  }
  
  score_delta<-function(y1,y2,rho,nu,alpha1,alpha2){
    s.rho=score_rho(y1,y2,rho,nu,alpha1,alpha2)
    s.delta=s.rho*(1-rho^2)
    return(s.delta)
  }
  
  
  
  obj<-function(pars){
    nu=pars[1]
    alpha1=pars[2]
    alpha2=pars[3]
    omega=pars[4]
    a=pars[5]
    b=pars[6]
    
    
    npars=6
    XB=rep(0,T)
    if(!is.null(X)){
      kk=dim(X)[2]
      Beta=pars[seq(from=(npars+1),length.out=kk)]
      XB=X%*% matrix(Beta,ncol=1)
    }
    
    
    
    
    t=1
    delta=omega
    rho=tanh(delta)
    logL <-  biv_std_est_dens(U1[t],U2[t],rho,nu,alpha1,alpha2,log=TRUE)
    
    
    for(t in 2:T){
      
      delta=omega+a*score_delta(U1[t-1],U2[t-1],rho,nu,alpha1,alpha2)+b*(delta-omega)+XB[t-1]
      rho=tanh(delta)
      logL <-  logL+biv_std_est_dens(U1[t],U2[t],rho,nu,alpha1,alpha2,log=TRUE)
    }
    
    
    
    return(-as.numeric(logL))
    
    
  }
  
  
  
  
  
if(is.null(theta.hat)){    
  if(is.null(pars0)){
    fit_static<-fit_biv_est(U,hessian=FALSE)
    
    nu0=exp(fit_static$par.hat[2])+2
    alpha1_0=fit_static$par.hat[3]
    alpha2_0=fit_static$par.hat[4]
    omega0=fit_static$par.hat[1]
    a0=.05
    b0=.9
    pars0=c(nu0,alpha1_0,alpha2_0,omega0,a0,b0)
    
    if(!is.null(X)){
      kk=dim(X)[2]
      pars0=c(pars0,rep(0,kk))
    }
  }
  
  
  lb=c(2,rep(-Inf,4),-1)
  ub=c(rep(Inf,5),1)
  if(!is.null(X)){
    kk=dim(X)[2]
    lb <- c(lb,rep(-Inf,kk))
    ub <- c(ub,rep(Inf,kk))
  }
  
  fit= solnp(pars0, obj, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL, 
             ineqUB = NULL, LB = lb, UB = ub ,control=list(trace=0))
  
  theta.hat=fit$pars
  
  ll=-fit$value[length(fit$value)]
  AIC=2*length(theta.hat)-2*ll
  BIC=log(dim(U)[1])*length(theta.hat)-2*ll
  
 
  
  #pars=theta.hat
}

    
  rho_filtered<-function(pars,U,X){
    nu=pars[1]
    alpha1=pars[2]
    alpha2=pars[3]
    omega=pars[4]
    a=pars[5]
    b=pars[6]
    
    
    T=dim(U)[1]
    npars=5
    XB=rep(0,T)
    if(!is.null(X)){
      kk=dim(X)[2]
      Beta=pars[seq(from=(npars+1),length.out=kk)]
      XB=X%*% matrix(Beta,ncol=1)
    }
    
    
    Rho=rep(0,T)
    
    t=1
    delta=omega
    rho=tanh(delta)
    Rho[t]<-rho
    
    
    
    
    for(t in 2:T){
      delta=omega+a*score_delta(U1[t-1],U2[t-1],rho,nu,alpha1,alpha2)+b*(delta-omega)+XB[t-1]
      rho=tanh(delta)
      Rho[t]<-rho
    }
    
    t=T+1
    delta=omega+a*score_delta(U1[t-1],U2[t-1],rho,nu,alpha1,alpha2)+b*(delta-omega)+XB[t-1]
    rho=tanh(delta)
    Rho_fcst<-rho
    
    ###
    ###with zz1= sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
    
    alpha=c(alpha1,alpha2)
    
    alpha2_star=alpha1^2+2*Rho*alpha1*alpha2+alpha2^2
    
    
    delta1=(alpha1+Rho*alpha2)/sqrt(1+alpha2_star)
    delta2=(Rho*alpha1+alpha2)/sqrt(1+alpha2_star)
    
    zz1= sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
    
    
    b1=sqrt(  nu/(nu-2) - (zz1*delta1)^2)
    b2=sqrt(  nu/(nu-2) - (zz1*delta2)^2)
    
    Rho2=(nu/(nu-2)*Rho-zz1^2*delta1*delta2)/(b1*b2)#Actual correlation
    
    
    alpha2_star=alpha1^2+2*Rho_fcst*alpha1*alpha2+alpha2^2
    
    
    delta1=(alpha1+Rho_fcst*alpha2)/sqrt(1+alpha2_star)
    delta2=(Rho_fcst*alpha1+alpha2)/sqrt(1+alpha2_star)
    
    zz1= sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
    
    
    b1=sqrt(  nu/(nu-2) - (zz1*delta1)^2)
    b2=sqrt(  nu/(nu-2) - (zz1*delta2)^2)
    
    Rho2_fcst=(nu/(nu-2)*Rho_fcst-zz1^2*delta1*delta2)/(b1*b2)
    
    return(list(Rho=Rho,Rho_fcst=Rho_fcst,Rho2=Rho2,Rho2_fcst=Rho2_fcst))
    
    
  }
  
  FILT=rho_filtered(theta.hat,U,X)
 
if(is.null(theta.hat)){     
  return(list(theta.hat=theta.hat, LogLik=ll,AIC=AIC,BIC=BIC,
               nu_hat=theta.hat[1],
               alpha1_hat=theta.hat[2],
               alpha2_hat=theta.hat[3],
              Rho_fcst=as.numeric(FILT$Rho_fcst),
              Rho2_fcst=as.numeric(FILT$Rho2_fcst)
  ) )}
  else{  return(list(theta.hat=theta.hat,
                     nu_hat=theta.hat[1],
                     alpha1_hat=theta.hat[2],
                     alpha2_hat=theta.hat[3],
                     Rho_fcst=as.numeric(FILT$Rho_fcst),
                     Rho2_fcst=as.numeric(FILT$Rho2_fcst)
  ) )}
  
  
}


PIT_normal<-function(y1,y2,alpha,pars){
  require(mvtnorm)
  
  mu1<-pars[1]
  mu2<-pars[2]
  sigma1<-pars[3]
  sigma2<-pars[4]
  rho<-pars[5]
  
  u2=pnorm(y2,mean=mu2,sd=sigma2)
  VaR2= mu2+sigma2*qnorm(alpha)
  lower=rep(-Inf,2)
  upper=c(y1,VaR2)
  covariance=rho*sigma1*sigma2
  sigma=matrix(c(sigma1^2,covariance,covariance,sigma2^2),ncol=2)
  u12=pmvnorm(lower, upper, mean=c(mu1,mu2), sigma=sigma)/alpha
  
  

  # N=100000
  # Z1=rnorm(N)
  # Z2=rnorm(N)
  # W2=rho*Z1+sqrt(1-rho^2)*Z2
  # Y1=mu1+sigma1*Z1
  # Y2=mu2+sigma2*W2
  # 
  # VaR2=quantile(Y2,alpha)
  # u2=mean(Y2<= y2)
  # I= Y2<=VaR2
  # Y1_I=Y1[I]
  # u12= mean(Y1_I<=y1)
  
  return(c(u12,u2))
}


PIT_t<-function(y1,y2,alpha,pars,T){
  #T: number of MC simulations 
  mu1<-pars[1]
  mu2<-pars[2]
  sigma1<-pars[3]
  sigma2<-pars[4]
  rho<-pars[5]
  nu<-pars[6]
  
  biv_t_sim<-function(T,rho,nu){
    #T simulations from the density biv_t_dens
    
    
    coef=sqrt((nu-2)/nu)
    q<- rchisq(T,df=nu)/nu
    s_q=sqrt(q)
    Z1=rnorm(T)
    Z2=rnorm(T)
    W2=rho*Z1+sqrt(1-rho^2)*Z2
    Y1=coef*Z1/s_q
    Y2=coef*W2/s_q
    return(cbind(Y1,Y2))
  }
  
  #T=100000
  ZZ=biv_t_sim(T,rho,nu)
  Y1=mu1+sigma1*ZZ[,1]
  Y2=mu2+sigma2*ZZ[,2]
   
  VaR2=quantile(Y2,alpha)
  u2=mean(Y2<= y2)
  I= Y2<=VaR2
  Y1_I=Y1[I]
  u12= mean(Y1_I<=y1)
  
  return(c(u12,u2))
  
  
}


PIT_skew_normal<-function(y1,y2,alpha,pars,T){
  #T: number of MC simulations 
  
  mu1<-pars[1]
  mu2<-pars[2]
  sigma1<-pars[3]
  sigma2<-pars[4]
  rho<-pars[5]
  alpha1=pars[6]
  alpha2=pars[7]
  

  biv_esnorm_sim<-function(T,rho,alpha1,alpha2){
    #T simulations from the density biv_std_esnorm_dens
    #Based on eq. (5.16) of Azzalini's book
    alpha=c(alpha1,alpha2)
    
    alpha2_star=alpha1^2+2*rho*alpha1*alpha2+alpha2^2
    delta=c(alpha1+rho*alpha2,rho*alpha1+alpha2 )/sqrt(1+alpha2_star)
    
    #correlated normal rv
    W1=rnorm(T)
    W2=rnorm(T)
    X0_1=W1
    X0_2=rho*W1+sqrt(1-rho^2)*W2
    
    U1=rnorm(T)
    
    #a_X0=alpha1*X0_1+alpha2*X0_2
    
    X1=(alpha1*X0_1 + alpha2*X0_2-U1)/sqrt(1+alpha1^2+alpha2^2+2*rho*alpha1*alpha2)
    
    #cond=U1> a_X0
    cond=X1>0
    
    Z1=cond*X0_1 - (1-cond)*X0_1
    Z2=cond*X0_2 - (1-cond)*X0_2
    
    zz1= sqrt(2/pi)
    #zz2=sn::zeta(2,tau)
    zz2=-zz1^2
    
    a1=zz1*delta[1]
    a2=zz1*delta[2]
    
    b1=sqrt(1+zz2*delta[1]^2)
    b2=sqrt(1+zz2*delta[2]^2)
    
    Y1=(Z1-a1)/b1
    Y2=(Z2-a2)/b2
    
    return(cbind(Y1,Y2))
  }
  
  
  #T=100000
  ZZ=biv_esnorm_sim(T,rho,alpha1,alpha2)
  Y1=mu1+sigma1*ZZ[,1]
  Y2=mu2+sigma2*ZZ[,2]
  
  VaR2=quantile(Y2,alpha)
  u2=mean(Y2<= y2)
  I= Y2<=VaR2
  Y1_I=Y1[I]
  u12= mean(Y1_I<=y1)
  
  return(c(u12,u2))
  
}


PIT_skew_t<-function(y1,y2,alpha,pars,T){
  #T: number of MC simulations 
  
  mu1<-pars[1]
  mu2<-pars[2]
  sigma1<-pars[3]
  sigma2<-pars[4]
  rho<-pars[5]
  alpha1=pars[6]
  alpha2=pars[7]
  nu<-pars[8]
  
  
  biv_est_sim<-function(T,rho,nu,alpha1,alpha2){
    #T simulations from the density biv_std_esnorm_dens
    #Based on eq. (5.16) of Azzalini's book
    
    alpha=c(alpha1,alpha2)
    
    alpha2_star=alpha1^2+2*rho*alpha1*alpha2+alpha2^2
    delta=c(alpha1+rho*alpha2,rho*alpha1+alpha2 )/sqrt(1+alpha2_star)
    #tau=0
    
    zz1= sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
    
    aa=zz1*delta
    a1=aa[1]
    a2=aa[2]
    bb=sqrt(  nu/(nu-2) - aa^2)
    b1=bb[1]
    b2=bb[2]
    
    #correlated normal rv
    W1=rnorm(T)
    W2=rnorm(T)
    X0_1=W1
    X0_2=rho*W1+sqrt(1-rho^2)*W2
    
    U1=rnorm(T)
    
    #a_X0=alpha1*X0_1+alpha2*X0_2
    
    X1=(alpha1*X0_1 + alpha2*X0_2-U1)/sqrt(1+alpha1^2+alpha2^2+2*rho*alpha1*alpha2)
    
    #cond=U1> a_X0
    cond=X1>0
    
    Z1=cond*X0_1 - (1-cond)*X0_1
    Z2=cond*X0_2 - (1-cond)*X0_2
    
    V=rchisq(T,df=nu)/nu
    VV=V^(-.5)
    
    
    Y1=(VV*Z1-a1)/b1
    Y2=(VV*Z2-a2)/b2
    
    return(cbind(Y1,Y2))
  }
  
  
  #T=100000
  ZZ=biv_est_sim(T,rho,nu,alpha1,alpha2)
  Y1=mu1+sigma1*ZZ[,1]
  Y2=mu2+sigma2*ZZ[,2]
  
  VaR2=quantile(Y2,alpha)
  u2=mean(Y2<= y2)
  I= Y2<=VaR2
  Y1_I=Y1[I]
  u12= mean(Y1_I<=y1)
  
  return(c(u12,u2))
  
}






DuEscanciano<-function(alpha,u12,u2,m){
  
  myBox.test<-function(x,const,m){
    n=length(x)
    y=x-const
    gamma0=var(y)
    gamma=rep(0,m)
    for(j in 1:m){
      gamma[j]=sum( y[(j+1):n]*y[1:(n-j)]  )/(n-j)
    }
    rho=gamma/gamma0
    statistic=n*sum(rho^2)
    p.value=1 - pchisq(statistic, m)
    return(list(statistic=statistic,p.value=p.value))
  }
  
  T=length(u12)
  I=u2<=alpha
  H= (1 - u12)*I
  H_bar=mean(H)
  alpha2=alpha/2
  U=sqrt(T)*(H_bar-alpha2)/sqrt(alpha*(1/3-alpha/4))
  pval.U=2-2*pnorm(abs(U))
  #rho=acf(H2,lag.max=m,plot = FALSE)$ acf 
  #C=T*(T-2)*sum(rho[-1]^2/(T-(1:m)))
  
  nm=length(m)
  C=rep(0,nm)
  pval.C=rep(0,nm)
  for(i in 1:nm){
    #Box=Box.test(H2, lag = m[i], type = "Ljung-Box") 
    Box=myBox.test(H,alpha2,m[i])
    C[i]=Box$ statistic
    pval.C[i]=Box$p.value
  }
  return(list(U=U,pval.U=pval.U,C=C,pval.C=pval.C))
}