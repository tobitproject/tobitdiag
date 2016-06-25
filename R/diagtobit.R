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

diag.tobitn <-function(model,tau=0,npoints=0)
{
  X <- model.matrix(model)
  p=ncol(X) #number parameters 
  n=nrow(X) #number observations 
  nu <- model$parms
  y  <- as.numeric(model$y)[1:n]
  c  <- (1*(y>tau))
  muhat <- model$linear.predictors
  sigmahat <- model$scale
  delta <-(y-muhat)/sigmahat
  phi         <- dnorm(delta)
  Phi         <- pnorm(delta)
  Wdelta      <- phi/Phi
  deriv_phi   <- (-delta/sqrt(2*pi))*exp(-(delta^2)/2)
  deriv_Wdelta<- (deriv_phi/(1-Phi))+((phi/(1-Phi))^2) #estou aqui!
  
  
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
  
  
  #C?lculo da matriz delta
  
  
  bi           <- -(1-c)*(1/sigmahat)*Wdelta+c*(1/sigmahat)*delta
  b            <- as.vector(bi)								   
  fi           <- -(1-c)*(1/sigmahat)*delta*Wdelta +c*(1/sigmahat)*(-1+delta^2)
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
  
  dev.new()
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(Ci, ylab="C",type="h", cex=0.5, ylim=c(0,0.40), pch=16,xlab="I")
  abline(h=corte1,lty=2)
  identify(Ci,n=npoints)
  
  dev.new()
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(Cib, ylab="Cb",type="h", cex=0.5, ylim=c(0,0.40), pch=16,xlab="I")
  abline(h=corte1,lty=2)
  identify(Cib,n=npoints)
  
  dev.new()
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(Cis, ylab="Cs",type="h", cex=0.5, ylim=c(0,0.40), pch=16,xlab="I")
  abline(h=corte1,lty=2)
  identify(Cis,n=npoints)
  
  dev.new()
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(abs(dmax), ylab="L",type="h", ylim=c(0,1), cex=0.5, pch=16,xlab="I")
  identify(abs(dmax),n=npoints)
  
  dev.new()
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(abs(dmaxb), ylab="Lb",type="h", ylim=c(0,1), cex=0.5, pch=16,xlab="I")
  identify(abs(dmaxb),n=npoints)
  
  
  dev.new()
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(abs(dmaxs), ylab="Ls",type="h", ylim=c(0,1), cex=0.5, pch=16,xlab="I")
  identify(abs(dmaxs),n=npoints)
  
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

diag.tobitt <-function(model,tau=0,npoints=0)
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
  
  dev.new()
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
  
  dev.new()
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(Ci, ylab="C",type="h", cex=0.5, ylim=c(0,0.40), pch=16,xlab="I")
  abline(h=corte1,lty=2)
  identify(Ci,n=npoints)
  
  dev.new()
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(Cib, ylab="Cb",type="h", cex=0.5, ylim=c(0,0.40), pch=16,xlab="I")
  abline(h=corte1,lty=2)
  identify(Cib,n=npoints)
  
  dev.new()
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(Cis, ylab="Cs",type="h", cex=0.5, ylim=c(0,0.40), pch=16,xlab="I")
  abline(h=corte1,lty=2)
  identify(Cis,n=npoints)
  
  dev.new()
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(abs(dmax), ylab="L",type="h", ylim=c(0,1), cex=0.5, pch=16,xlab="I")
  identify(abs(dmax),n=npoints)
  
  dev.new()
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(abs(dmaxb), ylab="Lb",type="h", ylim=c(0,1), cex=0.5, pch=16,xlab="I")
  identify(abs(dmaxb),n=npoints)
  
  
  dev.new()
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(abs(dmaxs), ylab="Ls",type="h", ylim=c(0,1), cex=0.5, pch=16,xlab="I")
  identify(abs(dmaxs),n=npoints)
  
}
