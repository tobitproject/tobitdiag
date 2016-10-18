#'Local Influence
#'
#'@description Local influence plots for tobit model.
#'
#'@param model an object of class "tobit" as fitted by tobit.
#'@param tau is the censoring point. The default is zero.
#'@param npoints points to be plotted.
#'
#'
#'@export

diag.tobitn <-function(model,tau=0,npoints=0,perturbation=c("cases","scale","response","explanatory"),l=NULL,ylim=c(0,0.40),plot.ci=FALSE,plot.dmax=FALSE,vecplot = c("theta","beta","sigma"))
{
  X <- model.matrix(model)
  p=ncol(X) #number parameters 
  n=nrow(X) #number observations 
  nu <- model$parms
  y  <- as.numeric(model$y)[1:n]
  c  <- (1*(y>tau))
  muhat <- model$linear.predictors
  sigmahat <- model$scale
  delta <- (y-muhat)/sigmahat
  phi         <- dnorm(delta)
  Phi         <- pnorm(delta)
  Wdelta      <- phi/Phi
  deriv_phi   <- (-delta/sqrt(2*pi))*exp(-(delta^2)/2)
  deriv_Wdelta<- (deriv_phi/(1-Phi))+((phi/(1-Phi))^2) 
  
  
  #############################################################################################
  ##########################################The Hessian Matrix#################################
  #############################################################################################
  
  #Lsigmasigma
  
  Lsigmasigma  <- (1/(sigmahat^2))*(sum((1-c)*(2*Wdelta*delta)))+
    (1/(sigmahat^2))*(sum((1-c)*delta*deriv_Wdelta*delta))+
    (1/(sigmahat^2))*(sum(c*(1-3*(delta^2))))
  
  #Quantidades para obter Lbetabeta
  
  #C?lculo de V
  
  
  parte1       <- (1-c)*(1/(sigmahat^2))*deriv_Wdelta	#Parte com censura
  parte2       <- (1/(sigmahat^2))*c  #Parte sem censura
  vi           <- parte1-parte2
  v            <- as.vector(vi)
  V            <- diag(v)
  
  #Lbetabeta
  
  
  Lbeta        <- t(X)%*%V%*%X
  
  
  #Quantidades para obter Lbetasigma
  
  #C?lculo de h
  
  parte11      <- (1-c)*((1/(sigmahat^2))*(Wdelta+delta*deriv_Wdelta)) #Parte com censura
  parte12      <- -c*(2/(sigmahat^2))*delta  #Parte sem censura
  hi           <- parte11+parte12
  
  #Lbetasigma
  
  Lbetasigma   <- t(X)%*%hi
  
  # Matriz de informa??o observada
  
  L1           <- cbind(Lbeta,Lbetasigma)
  L2           <- cbind(t(Lbetasigma), Lsigmasigma)
  
  
  Hessian <- rbind(L1,L2)
  
  #Obs.Infor <- Hessian
  #Inv.Obs <- solve(-Obs.Infor)
  #sqrt(diag(Inv.Obs))
  # Inversa da matriz de informa??o observada
  
  invobservmatrix <- -solve(Hessian)
  #diag(invobservmatrix)
  
  #Esquema de perturba??o pondera??o de casos
  
  
  #C?lculo da matriz delta
  Lsigma= Hessian[p+1,p+1]
  Lbeta=Hessian[1:p,1:p]
  b11=cbind(matrix(0, p, p), matrix(0, p, 1))
  b12=cbind(matrix(0, 1, p), -Lsigma^(-1))
  B1= rbind(b11, b12)  #parameter beta
  b211 =cbind(-solve(Lbeta), matrix(0, p, 1))
  b212= cbind(matrix(0, 1, p), matrix(0, 1, 1))
  B2=rbind(b211,b212)  # parameter delta
  
  
  ################# New structure of function ################
  if(perturbation == "cases")
  {
    #C?lculo da matriz delta
    bi           <- -(1-c)*(1/sigmahat)*Wdelta+c*(1/sigmahat)*delta
    b            <- as.vector(bi)								   
    fi           <- -(1-c)*(1/sigmahat)*delta*Wdelta +c*(1/sigmahat)*(-1+delta^2)
    f <- as.vector(fi)
    deltabeta    <- t(X)%*%diag(b)
    matrixdelta  <- rbind(deltabeta,t(fi))
    F0            <- t(matrixdelta)%*%(invobservmatrix)%*%matrixdelta #matriz F
    F1 <- t(matrixdelta)%*%(invobservmatrix - B1)%*%matrixdelta #matriz F beta
    F2 <- t(matrixdelta)%*%(invobservmatrix - B2)%*%matrixdelta #matriz F sigma
    
    
    # C?lculo dos autovalores e autovetores da matriz F
    
    Lmax<-eigen(F0,symmetric=TRUE)$val[1]
    dmax<-eigen(F0,symmetric=TRUE)$vec[,1]
    Lmaxb<-eigen(F1,symmetric=TRUE)$val[1]
    dmaxb<-eigen(F1,symmetric=TRUE)$vec[,1]
    Lmaxs<-eigen(F2,symmetric=TRUE)$val[1]
    dmaxs<-eigen(F2,symmetric=TRUE)$vec[,1]
    # Influ?ncia Local Total
    
    
    Ci<- 2*abs(diag(F0))
    Ci<- Ci/sum(Ci)
    
    Cib<- 2*abs(diag(F1))
    Cib<- Cib/sum(Cib)
    
    Cis<- 2*abs(diag(F2))
    Cis<- Cis/sum(Cis)
    #ponto de corte
    corte1 <-2*mean(Ci)
    corte2<-2*mean(Cib)
    corte3<-2*mean(Cis)
    
  } else if(perturbation == "scale"){
    #C?lculo da matriz delta
    
    bi=-(1-c)*(1/(2*sigmahat))*(Wdelta+delta*deriv_Wdelta)+c*(delta/sigmahat)
    b=as.vector(bi)
    fi=-(1-c)*(delta/(2*sigmahat))*(Wdelta+delta*deriv_Wdelta)+c*((delta^2)/sigmahat)
    f=as.vector(fi)
    deltabeta=t(X)%*%diag(b)
    matrixdelta= rbind(deltabeta,t(f))
    
    F0            <- t(matrixdelta)%*%(invobservmatrix)%*%matrixdelta #matriz F
    F1 <- t(matrixdelta)%*%(invobservmatrix - B1)%*%matrixdelta #matriz F beta
    F2 <- t(matrixdelta)%*%(invobservmatrix - B2)%*%matrixdelta #matriz F sigma
    
    # C?lculo dos autovalores e autovetores da matriz F
    
    Lmax<-eigen(F0,symmetric=TRUE)$val[1]
    dmax<-eigen(F0,symmetric=TRUE)$vec[,1]
    Lmaxb<-eigen(F1,symmetric=TRUE)$val[1]
    dmaxb<-eigen(F1,symmetric=TRUE)$vec[,1]
    Lmaxs<-eigen(F2,symmetric=TRUE)$val[1]
    dmaxs<-eigen(F2,symmetric=TRUE)$vec[,1]
    # Influ?ncia Local Total
    
    Ci<- 2*abs(diag(F0))
    Ci<- Ci/sum(Ci)
    
    Cib<- 2*abs(diag(F1))
    Cib<- Cib/sum(Cib)
    
    Cis<- 2*abs(diag(F2))
    Cis<- Cis/sum(Cis)
    #ponto de corte
    corte1<-2*mean(Ci)
    corte2<-2*mean(Cib)
    corte3<-2*mean(Cis)
    
  } else if(perturbation == "response"){
    #C?lculo da matriz delta
    sy <- sd(y)
    bi<- -(1-c)*0+c*(sy/(sigmahat^2))
    b <- as.vector(bi)
    fi <- -(1-c)*0+(c*(2*delta)*sy)/(sigmahat^2)
    f <- as.vector(fi)
    deltabeta <- t(X)%*%diag(b)
    matrixdelta <- rbind(deltabeta,t(f))
    
    
    F0            <- t(matrixdelta)%*%(invobservmatrix)%*%matrixdelta #matriz F
    F1 <- t(matrixdelta)%*%(invobservmatrix - B1)%*%matrixdelta #matriz F beta
    F2 <- t(matrixdelta)%*%(invobservmatrix - B2)%*%matrixdelta #matriz F sigma
    
    # C?lculo dos autovalores e autovetores da matriz F
    
    Lmax<-eigen(F0,symmetric=TRUE)$val[1]
    dmax<-eigen(F0,symmetric=TRUE)$vec[,1]
    Lmaxb<-eigen(F1,symmetric=TRUE)$val[1]
    dmaxb<-eigen(F1,symmetric=TRUE)$vec[,1]
    Lmaxs<-eigen(F2,symmetric=TRUE)$val[1]
    dmaxs<-eigen(F2,symmetric=TRUE)$vec[,1]
    # Influ?ncia Local Total
    
    Ci<- 2*abs(diag(F0))
    Ci<- Ci/sum(Ci)
    
    Cib<- 2*abs(diag(F1))
    Cib<- Cib/sum(Cib)
    
    Cis<- 2*abs(diag(F2))
    Cis<- Cis/sum(Cis)
    #ponto de corte
    corte1<-2*mean(Ci)
    corte2<-2*mean(Cib)
    corte3<-2*mean(Cis)
  } else{
    #C?lculo da matriz delta
    xl <- X[,l]
    sxl <- sd(xl)
    betal <- model$coefficients[l]
    fi  <- ((1-c)/(sigmahat^2))*(betal*sxl)*(delta*deriv_Wdelta+Wdelta)- (2*c*betal*sxl*delta)/(sigmahat^2)
    f <- as.vector(fi)
    xaux <- matrix(0,nrow=nrow(X),ncol=ncol(X))
    xaux[,l] <- ones(nrow(X),1)
    x0 <- xaux
    bi <- (1-c)*betal*sxl*deriv_Wdelta/(sigmahat*sigmahat) - c*sxl*betal/(sigmahat*sigmahat)
    b <- as.vector(bi)	
    bijt <-  -(1-c)*sxl*Wdelta/sigmahat + c*delta*sxl/sigmahat
    deltabeta <- crossprod(X,diag(b)) + crossprod(x0,diag(bijt))
    matrixdelta  <- rbind(deltabeta,t(f))
    F0            <- t(matrixdelta)%*%(invobservmatrix)%*%matrixdelta #matriz F
    F1 <- t(matrixdelta)%*%(invobservmatrix - B1)%*%matrixdelta #matriz F beta
    F2 <- t(matrixdelta)%*%(invobservmatrix - B2)%*%matrixdelta #matriz F sigma
    
    # C?lculo dos autovalores e autovetores da matriz F
    
    Lmax<-eigen(F0,symmetric=TRUE)$val[1]
    dmax<-eigen(F0,symmetric=TRUE)$vec[,1]
    Lmaxb<-eigen(F1,symmetric=TRUE)$val[1]
    dmaxb<-eigen(F1,symmetric=TRUE)$vec[,1]
    Lmaxs<-eigen(F2,symmetric=TRUE)$val[1]
    dmaxs<-eigen(F2,symmetric=TRUE)$vec[,1]
    
    Ci<- 2*abs(diag(F0))
    Ci<- Ci/sum(Ci)
    
    Cib<- 2*abs(diag(F1))
    Cib<- Cib/sum(Cib)
    
    Cis<- 2*abs(diag(F2))
    Cis<- Cis/sum(Cis)
    #ponto de corte
    corte1<-2*mean(Ci)
    corte2<-2*mean(Cib)
    corte3<-2*mean(Cis)
    
  }
  
  if(plot.ci == TRUE && perturbation != "explanatory" )
  {
    
    if(vecplot=="theta")
    {
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(Ci, ylab="C",type="h", cex=0.5, ylim=ylim, pch=16,xlab="I")
      abline(h=corte1,lty=2)
      if(npoints != 0) identify(Ci,n=npoints)
    } else if(vecplot=="beta"){
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(Cib, ylab="Cb",type="h", cex=0.5, ylim=ylim, pch=16,xlab="I")
      abline(h=corte1,lty=2)
      if(npoints != 0) identify(Cib,n=npoints)
    }
    else{
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(Cis, ylab="Cs",type="h", cex=0.5, ylim=ylim, pch=16,xlab="I")
      abline(h=corte1,lty=2)
      if(npoints != 0) identify(Cis,n=npoints)
    }
  }
  else if(plot.dmax == TRUE && perturbation != "explanatory")
  {
    if(vecplot=="theta")
    {
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(abs(dmax), ylab="L",type="h", ylim=ylim, cex=0.5, pch=16,xlab="I")
      if(npoints != 0) identify(abs(dmax),n=npoints)
    } else  if(vecplot=="beta")
    {
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(abs(dmaxb), ylab="Lb",type="h", ylim=ylim, cex=0.5, pch=16,xlab="I")
      if(npoints != 0) identify(abs(dmaxb),n=npoints)
    }else{
      
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(abs(dmaxs), ylab="Ls",type="h", ylim=ylim, cex=0.5, pch=16,xlab="I")
      if(npoints != 0) identify(abs(dmaxs),n=npoints)
    }
  }
  else if(plot.ci == TRUE && perturbation == "explanatory"){
    if(vecplot=="theta")
    {
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(xl, Ci, ylab="C",type="h", cex=0.5, ylim=ylim, pch=16,xlab="I")
      abline(h=corte1,lty=2)
      if(npoints != 0) identify(xl,Ci,n=npoints)
    } else  if(vecplot=="beta")
    {
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(xl,Cib, ylab="Cb",type="h", cex=0.5, ylim=ylim, pch=16,xlab="I")
      abline(h=corte1,lty=2)
      if(npoints != 0) identify(xl,Cib,n=npoints)
    }
    else{
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(xl,Cis, ylab="Cs",type="h", cex=0.5, ylim=ylim, pch=16,xlab="I")
      abline(h=corte1,lty=2)
      if(npoints != 0) identify(xl,Cis,n=npoints)
    }
  }
  else if(plot.dmax == TRUE && perturbation == "explanatory")
  {
    if(vecplot=="theta")
    {
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(xl,abs(dmax), ylab="L",type="h", ylim=ylim, cex=0.5, pch=16,xlab="I")
      if(npoints != 0) identify(xl,abs(dmax),n=npoints)
    } else  if(vecplot=="beta")
    {
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(xl,abs(dmaxb), ylab="Lb",type="h", ylim=ylim, cex=0.5, pch=16,xlab="I")
      if(npoints != 0) identify(xl,abs(dmaxb),n=npoints)
    } else{
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(xl,abs(dmaxs), ylab="Ls",type="h", ylim=ylim, cex=0.5, pch=16,xlab="I")
      if(npoints != 0) identify(xl,abs(dmaxs),n=npoints) 
    }
  } else{
    out <- list(ci = list(ctheta=Ci,cbeta=Cib,csigma=Cis),dmax = list(dmaxtheta=abs(dmax),dmaxbeta=abs(dmaxb),dmaxsigma=abs(dmaxs)),corte=list(ctheta=corte1,cbeta=corte2,csigma=corte3)) 
    return(out)
  }
  
}


#'Local Influence
#'
#'@description Local influence plots for tobit-t model.
#'
#'@param model an object of class "tobit" as fitted by tobit.
#'@param tau is the censoring point. The default is zero.
#'@param npoints points to be plotted.
#'
#'
#'@export

diag.tobitt <-function(model,tau=0,npoints=0,perturbation=c("cases","scale","response","explanatory"),l=NULL,ylim=c(0,0.40),plot.ci=FALSE,plot.dmax=FALSE,vecplot = c("theta","beta","sigma"))
{
  X <-model.matrix(model)
  p=ncol(X) #number parameters 
  n=nrow(X) #number observations 
  nu <- model$parms
  y  <- as.numeric(model$y)[1:n]
  c  <- (1*(y>tau))
  muhat <- model$linear.predictors
  sigmahat <- model$scale
  deltahat <-(y-muhat)/sigmahat
  phi         <- dt(deltahat,df=nu)
  Phi         <- pt(deltahat,df=nu)
  Wdelta      <- phi/Phi
  const <- gamma((nu+1)/2)/((gamma(nu/2)*sqrt(pi*nu)))
  deriv_phi <- -const*((nu+1)/nu)*deltahat*((1+((deltahat^2)/nu))^(-(nu+3)/2))
  deriv_Wdelta <- (deriv_phi/Phi)+((phi/Phi)^2)
  
  ##############################################Lsigmasigma####################################
  parte1 <- (1/(sigmahat^2))*sum((1-c)*(2*Wdelta*deltahat + (deltahat^2)*deriv_Wdelta))
  parte2 <- (1/(sigmahat^2))*sum(c*(1 - (nu+1)*((3*(deltahat^2)*nu + (deltahat^4))/(((deltahat^2)+nu)^2)))) 
  Lsigmasigma  <- (1/(sigmahat^2))*(sum((1-c)*(2*Wdelta*deltahat))) + (1/(sigmahat^2))*(sum((1-c)*deltahat*deriv_Wdelta*deltahat))+
    (1/(sigmahat^2))*(sum(c*(1-3*(((nu+1)/nu)*(deltahat^2)/(1+((deltahat^2)/nu)))+
                               2*((nu+1)/nu^2)*(deltahat^4)/(1+(deltahat^2)/nu)^2)))
  
  #################################################Lbetabeta ##################################
  parte1       <- (1-c)*(1/(sigmahat^2))*deriv_Wdelta	#Parte com censura
  parte2       <- c*(1/(sigmahat^2))*((((nu+1)/nu)/(1+((deltahat^2)/nu))) +
                                        2*((nu+1)/nu^2)*(deltahat^2)/(1+(deltahat^2)/nu)^2) #Parte sem censura
  vi           <- parte1 - parte2
  v            <- as.vector(vi)
  V            <- diag(v)
  Lbeta        <- t(X)%*%V%*%X
  
  ######################################################Lbetasigma############################### 
  parte11      <- (1-c)*((1/(sigmahat^2))*(Wdelta+deltahat*deriv_Wdelta)) #Parte com censura
  parte12      <- c*(2/(sigmahat^2))*((-(nu+1)/nu)*(deltahat/(1+((deltahat^2)/nu))) 
                                      +((nu+1)/nu^2)*(deltahat^3)/(1+(deltahat^2)/nu)^2)#Parte sem censura
  hi           <- parte11 + parte12
  Lbetasigma   <- t(X)%*%hi
  
  L1           <- cbind(Lbeta,Lbetasigma)
  L2           <- cbind(t(Lbetasigma), Lsigmasigma)
  
  Hessian <- rbind(L1,L2)
  
  #Obs.Infor <- Hessian
  #Inv.Obs <- solve(-Obs.Infor)
  #sqrt(diag(Inv.Obs))
  # Inversa da matriz de informa??o observada
  
  invobservmatrix <- -solve(Hessian)
  #diag(invobservmatrix)
  
  #Esquema de perturba??o pondera??o de casos
  Lsigma= Hessian[p+1,p+1]
  Lbeta=Hessian[1:p,1:p]
  b11=cbind(matrix(0, p, p), matrix(0, p, 1))
  b12=cbind(matrix(0, 1, p), -Lsigma^(-1))
  B1= rbind(b11, b12)  #parameter beta
  b211 =cbind(-solve(Lbeta), matrix(0, p, 1))
  b212= cbind(matrix(0, 1, p), matrix(0, 1, 1))
  B2=rbind(b211,b212)  # parameter delta
  
  b311 =cbind(matrix(0, p, p), matrix(0, p, 1))
  b312= cbind(matrix(0, 1, p), matrix(0, 1, 1))
  B3=rbind(b311,b312)  # parameter theta
  
  
  if(perturbation == "case")
  {
    #C?lculo da matriz delta
    bi           <- -(1-c)*(1/sigmahat)*Wdelta+ c*(1/sigmahat)*((nu+1)/nu)*deltahat/(1+((deltahat^2)/nu))
    b            <- as.vector(bi)								   
    fi           <- -(1-c)*(1/sigmahat)*deltahat*Wdelta +c*(1/sigmahat)*(-1+((nu+1)/nu)*(deltahat^2)/(1+((deltahat^2)/nu)))
    deltabeta    <- t(X)%*%diag(b)
    matrixdelta  <- rbind(deltabeta,t(fi))
    F0            <- t(matrixdelta)%*%(invobservmatrix)%*%matrixdelta #matriz F
    F1 <- t(matrixdelta)%*%(invobservmatrix - B1)%*%matrixdelta #matriz F beta
    F2 <- t(matrixdelta)%*%(invobservmatrix - B2)%*%matrixdelta #matriz F sigma
    
    # C?lculo dos autovalores e autovetores da matriz F
    
    Lmax<-eigen(F0,symmetric=TRUE)$val[1]
    dmax<-eigen(F0,symmetric=TRUE)$vec[,1]
    Lmaxb<-eigen(F1,symmetric=TRUE)$val[1]
    dmaxb<-eigen(F1,symmetric=TRUE)$vec[,1]
    Lmaxs<-eigen(F2,symmetric=TRUE)$val[1]
    dmaxs<-eigen(F2,symmetric=TRUE)$vec[,1]
    # Influ?ncia Local Total
    
    Ci<- 2*abs(diag(F0))
    Ci<- Ci/sum(Ci)
    
    Cib<- 2*abs(diag(F1))
    Cib<- Cib/sum(Cib)
    
    Cis<- 2*abs(diag(F2))
    Cis<- Cis/sum(Cis)
    #ponto de corte
    corte1<-2*mean(Ci)
    corte2<-2*mean(Cib)
    corte3<-2*mean(Cis)
    
    
  } else if(perturbation == "scale"){
    #C?lculo da matriz delta
    bi           <- -(1-c)*(1/(2*sigmahat))*(Wdelta+deltahat*deriv_Wdelta) + c*((nu+1)/(nu*sigmahat))*( deltahat/( 1+( (deltahat^2)/nu) ) - ((deltahat^3)/(nu*(1+(deltahat^2)/nu)^2) ))
    b            <- as.vector(bi)								   
    fi           <- -(1-c)*(deltahat/(2*sigmahat))*(Wdelta+deltahat*deriv_Wdelta) + c*((nu+1)/(nu*sigmahat))*deltahat*( deltahat/( 1+( (deltahat^2)/nu) ) - ((deltahat^3)/(nu*(1+(deltahat^2)/nu)^2)))
    deltabeta    <- t(X)%*%diag(b)
    matrixdelta  <- rbind(deltabeta,t(fi))
    F0            <- t(matrixdelta)%*%(invobservmatrix)%*%matrixdelta #matriz F
    F1 <- t(matrixdelta)%*%(invobservmatrix - B1)%*%matrixdelta #matriz F beta
    F2 <- t(matrixdelta)%*%(invobservmatrix - B2)%*%matrixdelta #matriz F sigma
    
    # C?lculo dos autovalores e autovetores da matriz F
    
    Lmax<-eigen(F0,symmetric=TRUE)$val[1]
    dmax<-eigen(F0,symmetric=TRUE)$vec[,1]
    Lmaxb<-eigen(F1,symmetric=TRUE)$val[1]
    dmaxb<-eigen(F1,symmetric=TRUE)$vec[,1]
    Lmaxs<-eigen(F2,symmetric=TRUE)$val[1]
    dmaxs<-eigen(F2,symmetric=TRUE)$vec[,1]
    # Influ?ncia Local Total
    
    
    Ci<- 2*abs(diag(F0))
    Ci<- Ci/sum(Ci)
    
    Cib<- 2*abs(diag(F1))
    Cib<- Cib/sum(Cib)
    
    Cis<- 2*abs(diag(F2))
    Cis<- Cis/sum(Cis)
    #ponto de corte
    corte1<-2*mean(Ci)
    corte2<-2*mean(Cib)
    corte3<-2*mean(Cis)
    
    
  } else if(perturbation == "response"){
    #C?lculo da matriz delta
    sy<- sd(y)
    bi<- -(1-c)*0+c*((sy*(nu+1))/(nu*(sigmahat^2)))*(1/(1+( (deltahat^2)/nu))-(2*(deltahat^2))/(nu*(1+((deltahat^2)/nu))^2))
    b <- as.vector(bi)
    firesp<- -(1-c)*0+c*2*sy*((nu+1)/(nu*sigmahat*sigmahat))*(deltahat/(1+(deltahat^2)/nu)- (deltahat^3)/(nu*((1+(deltahat^2)/nu)^2) )  )
    fi<- as.vector(firesp)
    deltabeta    <- t(X)%*%diag(b)
    matrixdelta  <- rbind(deltabeta,t(fi))
    F0            <- t(matrixdelta)%*%(invobservmatrix)%*%matrixdelta #matriz F
    F1 <- t(matrixdelta)%*%(invobservmatrix - B1)%*%matrixdelta #matriz F beta
    F2 <- t(matrixdelta)%*%(invobservmatrix - B2)%*%matrixdelta #matriz F sigma
    
    # C?lculo dos autovalores e autovetores da matriz F
    
    Lmax<-eigen(F0,symmetric=TRUE)$val[1]
    dmax<-eigen(F0,symmetric=TRUE)$vec[,1]
    Lmaxb<-eigen(F1,symmetric=TRUE)$val[1]
    dmaxb<-eigen(F1,symmetric=TRUE)$vec[,1]
    Lmaxs<-eigen(F2,symmetric=TRUE)$val[1]
    dmaxs<-eigen(F2,symmetric=TRUE)$vec[,1]
    # Influ?ncia Local Total
    
    Ci<- 2*abs(diag(F0))
    Ci<- Ci/sum(Ci)
    
    Cib<- 2*abs(diag(F1))
    Cib<- Cib/sum(Cib)
    
    Cis<- 2*abs(diag(F2))
    Cis<- Cis/sum(Cis)
    #ponto de corte
    corte1<-2*mean(Ci)
    corte2<-2*mean(Cib)
    corte3<-2*mean(Cis)
    
    
  } else{
    #C?lculo da matriz delta
    xl <- X[,l]
    sxl <- sd(xl)
    betal <- model$coefficients[l]
    fi <- (1-c)*(betal*sxl/(sigmahat^2))*(deltahat*deriv_Wdelta+Wdelta) - c*2*betal*sxl*((nu+1)/(nu*sigmahat*sigmahat))*( deltahat/(1+(deltahat^2)/nu) - (deltahat*deltahat*deltahat)/((nu*(1+ (deltahat^2)/nu)^2) ))
    xaux <- matrix(0,nrow=nrow(X),ncol=ncol(X))
    xaux[,l] <- ones(nrow(X),1)
    x0 <- xaux
    bi <- (1-c)*betal*sxl*deriv_Wdelta/(sigmahat*sigmahat) - c*sxl*betal*((nu+1)/(nu*sigmahat))*( 1/(1+(deltahat^2)/nu) - (2*deltahat*deltahat)/((nu*(1+ (deltahat^2)/nu)^2) ))
    b <- as.vector(bi)	
    bijt <-  -(1-c)*Wdelta/sqrt(sigmahat) + c*((nu+1)/(nu*sqrt(sigmahat)))*(deltahat/(1+(deltahat^2)/nu))
    deltabeta <- crossprod(X,diag(b)) + sxl*crossprod(x0,diag(bijt))
    matrixdelta  <- rbind(deltabeta,t(fi))
    F0            <- t(matrixdelta)%*%(invobservmatrix)%*%matrixdelta #matriz F
    F1 <- t(matrixdelta)%*%(invobservmatrix - B1)%*%matrixdelta #matriz F beta
    F2 <- t(matrixdelta)%*%(invobservmatrix - B2)%*%matrixdelta #matriz F sigma
    
    # C?lculo dos autovalores e autovetores da matriz F
    
    Lmax<-eigen(F0,symmetric=TRUE)$val[1]
    dmax<-eigen(F0,symmetric=TRUE)$vec[,1]
    Lmaxb<-eigen(F1,symmetric=TRUE)$val[1]
    dmaxb<-eigen(F1,symmetric=TRUE)$vec[,1]
    Lmaxs<-eigen(F2,symmetric=TRUE)$val[1]
    dmaxs<-eigen(F2,symmetric=TRUE)$vec[,1]
    
    Ci<- 2*abs(diag(F0))
    Ci<- Ci/sum(Ci)
    
    Cib<- 2*abs(diag(F1))
    Cib<- Cib/sum(Cib)
    
    Cis<- 2*abs(diag(F2))
    Cis<- Cis/sum(Cis)
    #ponto de corte
    corte1<-2*mean(Ci)
    corte2<-2*mean(Cib)
    corte3<-2*mean(Cis)
    
  }
  
  if(plot.ci == TRUE && perturbation != "explanatory" )
  {
    
    if(vecplot=="theta")
    {
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(Ci, ylab="C",type="h", cex=0.5, ylim=ylim, pch=16,xlab="I")
      abline(h=corte1,lty=2)
      if(npoints != 0) identify(Ci,n=npoints)
    } else if(vecplot=="beta"){
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(Cib, ylab="Cb",type="h", cex=0.5, ylim=ylim, pch=16,xlab="I")
      abline(h=corte1,lty=2)
      if(npoints != 0) identify(Cib,n=npoints)
    }
    else{
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(Cis, ylab="Cs",type="h", cex=0.5, ylim=ylim, pch=16,xlab="I")
      abline(h=corte1,lty=2)
      if(npoints != 0) identify(Cis,n=npoints)
    }
  }
  else if(plot.dmax == TRUE && perturbation != "explanatory")
  {
    if(vecplot=="theta")
    {
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(abs(dmax), ylab="L",type="h", ylim=ylim, cex=0.5, pch=16,xlab="I")
      if(npoints != 0) identify(abs(dmax),n=npoints)
    } else  if(vecplot=="beta")
    {
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(abs(dmaxb), ylab="Lb",type="h", ylim=ylim, cex=0.5, pch=16,xlab="I")
      if(npoints != 0) identify(abs(dmaxb),n=npoints)
    }else{
      
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(abs(dmaxs), ylab="Ls",type="h", ylim=ylim, cex=0.5, pch=16,xlab="I")
      if(npoints != 0) identify(abs(dmaxs),n=npoints)
    }
  }
  else if(plot.ci == TRUE && perturbation == "explanatory"){
    if(vecplot=="theta")
    {
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(xl, Ci, ylab="C",type="h", cex=0.5, ylim=ylim, pch=16,xlab="I")
      abline(h=corte1,lty=2)
      if(npoints != 0) identify(xl,Ci,n=npoints)
    } else  if(vecplot=="beta")
    {
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(xl,Cib, ylab="Cb",type="h", cex=0.5, ylim=ylim, pch=16,xlab="I")
      abline(h=corte1,lty=2)
      if(npoints != 0) identify(xl,Cib,n=npoints)
    }
    else{
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(xl,Cis, ylab="Cs",type="h", cex=0.5, ylim=ylim, pch=16,xlab="I")
      abline(h=corte1,lty=2)
      if(npoints != 0) identify(xl,Cis,n=npoints)
    }
  }
  else if(plot.dmax == TRUE && perturbation == "explanatory")
  {
    if(vecplot=="theta")
    {
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(xl,abs(dmax), ylab="L",type="h", ylim=ylim, cex=0.5, pch=16,xlab="I")
      if(npoints != 0) identify(xl,abs(dmax),n=npoints)
    } else  if(vecplot=="beta")
    {
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(xl,abs(dmaxb), ylab="Lb",type="h", ylim=ylim, cex=0.5, pch=16,xlab="I")
      if(npoints != 0) identify(xl,abs(dmaxb),n=npoints)
    } else{
      par(mar=c(4.0,4.0,0.1,0.1))
      plot(xl,abs(dmaxs), ylab="Ls",type="h", ylim=ylim, cex=0.5, pch=16,xlab="I")
      if(npoints != 0) identify(xl,abs(dmaxs),n=npoints) 
    }
  } else{
    out <- list(ci = list(ctheta=Ci,cbeta=Cib,csigma=Cis),dmax = list(dmaxtheta=abs(dmax),dmaxbeta=abs(dmaxb),dmaxsigma=abs(dmaxs)),corte=list(ctheta=corte1,cbeta=corte2,csigma=corte3)) 
    return(out)
  }
  
}


#' Cook's distance
#' 
#' @description A measure that combines the information of leverage and residual of the observation
#' 
#' @param model x
#' @param tau x
#' @param npoints x
#' @param plot x
#' @param dist x
#' 
#' @export

<<<<<<< HEAD
cooks.distn <- function(formula, tau=0, npoints=0, dist="t", plot=FALSE,xlab=NULL,ylab=NULL,ylim=NULL,xlim=NULL,pch=19,cex=0.5,type="p",col=NULL)
{
  
  if(dist=="t")
  {
    model <-  tobit(formula, dist="t")
    theta <- as.vector(c(model$coef,model$scale))
    X    <- model.matrix(model) 
    n    <-dim(X)[1]
    p    <-dim(X)[2]
    y  <- as.numeric(model$y)[1:n]
    c  <- (1*(y>tau))
    #Lsigmasigma
    nu <- model$parms
    muhat <- model$linear.predictors
    sigmahat <- model$scale
    deltahat <-(y-muhat)/sigmahat
    phi         <- dt(deltahat,model$parms)
    Phi         <- pt(deltahat,model$parms)
    Wdelta      <- phi/Phi
    
    
    #############################################################################################
    ##########################################The Hessian Matrix#################################
    #############################################################################################
    
    const <- gamma((nu+1)/2)/((gamma(nu/2)*sqrt(pi*nu)))
    deriv_phi <- -const*((nu+1)/nu)*deltahat*((1+((deltahat^2)/nu))^(-(nu+3)/2))
    deriv_Wdelta <- (deriv_phi/Phi)+((phi/Phi)^2)
    
    ##############################################Lsigmasigma####################################
    parte1 <- (1/(sigmahat^2))*sum((1-c)*(2*Wdelta*deltahat + (deltahat^2)*deriv_Wdelta))
    parte2 <- (1/(sigmahat^2))*sum(c*(1 - (nu+1)*((3*(deltahat^2)*nu + (deltahat^4))/(((deltahat^2)+nu)^2)))) 
    Lsigmasigma  <- (1/(sigmahat^2))*(sum((1-c)*(2*Wdelta*deltahat))) + (1/(sigmahat^2))*(sum((1-c)*deltahat*deriv_Wdelta*deltahat))+
      (1/(sigmahat^2))*(sum(c*(1-3*(((nu+1)/nu)*(deltahat^2)/(1+((deltahat^2)/nu)))+
                                 2*((nu+1)/nu^2)*(deltahat^4)/(1+(deltahat^2)/nu)^2)))
    
    #################################################Lbetabeta ##################################
    parte1       <- (1-c)*(1/(sigmahat^2))*deriv_Wdelta	#Parte com censura
    parte2       <- c*(1/(sigmahat^2))*((((nu+1)/nu)/(1+((deltahat^2)/nu))) +
                                          2*((nu+1)/nu^2)*(deltahat^2)/(1+(deltahat^2)/nu)^2) #Parte sem censura
    vi           <- parte1 - parte2
    v            <- as.vector(vi)
    V            <- diag(v)
    Lbeta        <- t(X)%*%V%*%X
    
    ######################################################Lbetasigma############################### 
    parte11      <- (1-c)*((1/(sigmahat^2))*(Wdelta+deltahat*deriv_Wdelta)) #Parte com censura
    parte12      <- c*(2/(sigmahat^2))*((-(nu+1)/nu)*(deltahat/(1+((deltahat^2)/nu))) 
                                        +((nu+1)/nu^2)*(deltahat^3)/(1+(deltahat^2)/nu)^2)#Parte sem censura
    hi           <- parte11 + parte12
    Lbetasigma   <- t(X)%*%hi
    
    L1           <- cbind(Lbeta,Lbetasigma)
    L2           <- cbind(t(Lbetasigma), Lsigmasigma)
    observmatrix <- rbind(L1,L2)
    
    values <- vector()
    for(i in 1:n)
    {
      fit <- tobit(formula,dist="t",subset = -c(i) )
      thetai<-as.vector(c(fit$coef,fit$scale))
      values[i]<-(1/(p+1))*t(theta-thetai)%*%(-observmatrix)%*%(theta-thetai)
    }
  } 
  else{
=======
cooks.distn <- function(formula, tau=0, npoints=NULL, dist="t", plot=FALSE,xlab=NULL,ylab=NULL,ylim=NULL,xlim=NULL,pch=19,cex=0.5,type="p",col=NULL)
  {
  
  if(dist=="t")
  {
  modt <-  tobit(formula, dist="t")
  thetas <- as.vector(c(modt$coef,modt$scale))
  X    <- model.matrix(modt) 
  n    <-dim(X)[1]
  p    <-dim(X)[2]
  y  <- as.numeric(model$y)[1:n]
  c  <- (1*(y>tau))
  #Lsigmasigma
  nu <- modt$parms
  muhat <- modt$linear.predictors
  sigmahat <- modt$scale
  deltahat <-(y-muhat)/sigmahat
  phi         <- dt(deltahat,modt$parms)
  Phi         <- pt(deltahat,modt$parms)
  Wdelta      <- phi/Phi
  
  
  #############################################################################################
  ##########################################The Hessian Matrix#################################
  #############################################################################################
  
  const <- gamma((nu+1)/2)/((gamma(nu/2)*sqrt(pi*nu)))
  deriv_phi <- -const*((nu+1)/nu)*deltahat*((1+((deltahat^2)/nu))^(-(nu+3)/2))
  deriv_Wdelta <- (deriv_phi/Phi)+((phi/Phi)^2)
  
  ##############################################Lsigmasigma####################################
  parte1 <- (1/(sigmahat^2))*sum((1-c)*(2*Wdelta*deltahat + (deltahat^2)*deriv_Wdelta))
  parte2 <- (1/(sigmahat^2))*sum(c*(1 - (nu+1)*((3*(deltahat^2)*nu + (deltahat^4))/(((deltahat^2)+nu)^2)))) 
  Lsigmasigma  <- (1/(sigmahat^2))*(sum((1-c)*(2*Wdelta*deltahat))) + (1/(sigmahat^2))*(sum((1-c)*deltahat*deriv_Wdelta*deltahat))+
    (1/(sigmahat^2))*(sum(c*(1-3*(((nu+1)/nu)*(deltahat^2)/(1+((deltahat^2)/nu)))+
                               2*((nu+1)/nu^2)*(deltahat^4)/(1+(deltahat^2)/nu)^2)))
  
  #################################################Lbetabeta ##################################
  parte1       <- (1-c)*(1/(sigmahat^2))*deriv_Wdelta	#Parte com censura
  parte2       <- c*(1/(sigmahat^2))*((((nu+1)/nu)/(1+((deltahat^2)/nu))) +
                                        2*((nu+1)/nu^2)*(deltahat^2)/(1+(deltahat^2)/nu)^2) #Parte sem censura
  vi           <- parte1 - parte2
  v            <- as.vector(vi)
  V            <- diag(v)
  Lbeta        <- t(X)%*%V%*%X
  
  ######################################################Lbetasigma############################### 
  parte11      <- (1-c)*((1/(sigmahat^2))*(Wdelta+deltahat*deriv_Wdelta)) #Parte com censura
  parte12      <- c*(2/(sigmahat^2))*((-(nu+1)/nu)*(deltahat/(1+((deltahat^2)/nu))) 
                                      +((nu+1)/nu^2)*(deltahat^3)/(1+(deltahat^2)/nu)^2)#Parte sem censura
  hi           <- parte11 + parte12
  Lbetasigma   <- t(X)%*%hi
  
  L1           <- cbind(Lbeta,Lbetasigma)
  L2           <- cbind(t(Lbetasigma), Lsigmasigma)
  observmatrix <- rbind(L1,L2)
  
  values <- vector()
   for(i in 1:n)
    {
    observ <- i
    fit <- tobit(formula,dist="t",subset = -observ)
    thetai<-as.vector(c(fit$coef,fit$scale))
    values[i]<-(1/(p+1))*t(theta-thetai)%*%(-observmatrix)%*%(theta-thetai)
    }
  } 
    else{
>>>>>>> 3927da0b4283228af32c64cc74388102d73c101e
    model <- tobit(formula)
    theta <-c(model$coef,model$scale)
    theta <-as.vector(theta)
    X    <-model.matrix(model) 
    n    <-dim(X)[1]
    p    <-dim(X)[2]
    y  <- as.numeric(model$y)[1:n]
    c  <- (1*(y>tau))
    muhat <- model$linear.predictors
    sigmahat <- model$scale
    delta <- (y-muhat)/sigmahat
    phi         <- dnorm(delta)
    Phi         <- pnorm(delta)
    Wdelta      <- phi/Phi
    deriv_phi   <- (-delta/sqrt(2*pi))*exp(-(delta^2)/2)
    deriv_Wdelta<- (deriv_phi/(1-Phi))+((phi/(1-Phi))^2) 
    
    #############################################################################################
    ##########################################The Hessian Matrix#################################
    #############################################################################################
    #Lsigmasigma
    Lsigmasigma  <- (1/(sigmahat^2))*(sum((1-c)*(2*Wdelta*delta)))+
      (1/(sigmahat^2))*(sum((1-c)*delta*deriv_Wdelta*delta))+
      (1/(sigmahat^2))*(sum(c*(1-3*(delta^2))))
    #Quantidades para obter Lbetabeta
    #C?lculo de V
    parte1       <- (1-c)*(1/(sigmahat^2))*deriv_Wdelta	#Parte com censura
    parte2       <- (1/(sigmahat^2))*c  #Parte sem censura
    vi           <- parte1-parte2
    v            <- as.vector(vi)
    V            <- diag(v)
    #Lbetabeta
    Lbeta        <- t(X)%*%V%*%X
    #Quantidades para obter Lbetasigma
    #C?lculo de h
    parte11      <- (1-c)*((1/(sigmahat^2))*(Wdelta+delta*deriv_Wdelta)) #Parte com censura
    parte12      <- -c*(2/(sigmahat^2))*delta  #Parte sem censura
    hi           <- parte11+parte12
    #Lbetasigma
    Lbetasigma   <- t(X)%*%hi
    # Matriz de informa??o observada
    L1           <- cbind(Lbeta,Lbetasigma)
    L2           <- cbind(t(Lbetasigma), Lsigmasigma)
    observmatrix <- rbind(L1, L2)
    
<<<<<<< HEAD
    values <- vector()
    for(i in 1:n)
    {
      fit <- tobit(formula,subset = -c(i) )
      thetai<-as.vector(c(fit$coef,fit$scale))
      values[i]<-(1/(p+1))*t(theta-thetai)%*%(-observmatrix)%*%(theta-thetai)     }
  }
  
  if(plot==TRUE)
  {  
    par(mar=c(4.0,4.0,0.1,0.1)) #graphical parameter of the margin.
    plot(values, ylab=ylab,type=type,xlim = xlim, ylim=ylim, cex=cex, pch=pch,xlab=xlab,col=col,main = NULL,sub = NULL)
    if(npoints!=0)identify(values,n=npoints)
  } else return(values)
=======
     values <- vector()
     for(i in 1:n)
     {
     observ <- i 
     print(observ)
     fit <- tobit(formula,subset = -observ)
     thetai<-as.vector(c(fit$coef,fit$scale))
     print(thetai)
     values[i]<-(1/(p+1))*t(theta-thetai)%*%(-observmatrix)%*%(theta-thetai)
     print(values)
     }
    }
  
if(plot==TRUE)
{
if(is.null(npoints))
{  
par(mar=c(4.0,4.0,0.1,0.1)) #graphical parameter of the margin.
plot(values, ylab=ylab,type=type, ylim=ylim, cex=cex, pch=pch,xlab=xlab,col=col,main = NULL,sub = NULL)
} 
  else{
  par(mar=c(4.0,4.0,0.1,0.1)) #graphical parameter of the margin.
  plot(values, ylab=ylab,type=type, ylim=ylim, cex=cex, pch=pch,xlab=xlab,col=col,main = NULL,sub = NULL)
  identify(values,n=npoints)
  }
} else return(values)
>>>>>>> 3927da0b4283228af32c64cc74388102d73c101e
  
}



###################################################
#information criterions
####################################################
#' Akaike's An Information Criterion
#' 
#' AICc is AIC with a correction for finite sample sizes.
#' 
#' @param object a fitted model object for which there exists a logLik method to extract the corresponding log-likelihood, or an object inheriting from class logLik
#' @param object1 optionally more fitted model objects
#' 
#' @return 
#' If just one object is provided, a numeric value with the corresponding AICc.
#' 
#' If two objects are provided, a data.frame with rows corresponding to the objects and columns 
#' representing the number of 
#' parameters in the model (df) and the AICc.
#' 
#' @export

AICc <- function (object,object1=NULL) 
{
  if(is.null(object1))
  {
    aic <- AIC(object)
    if (!is.numeric(aic)) 
      stop("Cannot calculate AIC!")
    k <- object$df
    n <- length(object$y)
    aicc <- aic + ((2 * k * (k + 1))/(n - k - 1))
    return(aicc) 
  } 
  else{
    aic <- AIC(object);aic1 <- AIC(object1)
    if (!is.numeric(aic)&&!is.numeric(aic1)) 
      stop("Cannot calculate AIC!")
    k <- object$df ; k1 <- object$df
    n <- length(object$y)
    aicc <- aic + ((2 * k * (k + 1))/(n - k - 1))
    aicc1 <- aic1 + ((2 * k1 * (k1 + 1))/(n - k1 - 1))
    
    out <- cbind(c(k,k1),c(aicc,aicc1))    
    colnames(out) <- c("df","AICc")    
    return(out[order(-out[2,],out[1,]),])
  }
}


