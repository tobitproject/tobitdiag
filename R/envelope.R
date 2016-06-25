#'Residuals
#'
#'@description Calculates martingale, deviance, martingale type residuals for tobit-t model.
#'
#'@param model an object of class "tobit" as fitted by tobit.
#'@param type what type of residuals should be used. 
#'Possible type are: "deviance", "martingale" and "martingale type". 
#'@param tau is the censoring point. The default is zero.
#'
#'@return The returned object is a vector with one element for each subject. 
#'
#'@export

residuals <- function(model,type,tau=0)
{ # begin function
  n <- summary(model)$n[1] 
  y  <- as.numeric(model$y)[1:n]
  c  <- (1*(y>tau))
  nu <- model$parms
  muhat = model$linear.predictors
  sigmahat <- model$scale
  deltahat <-(y-muhat)/sigmahat
  S <- 1-pt(deltahat,nu)
  rM <- c+log(S) # Martingale residual
  rtM <- (sign(y-muhat)*(-2*(rM+c*log(c-rM)))^(1/2))#Martingale-type residual
  rDC <- sign(y-muhat)*sqrt(-2*log(S))
  
  switch(type,
         deviance = rDC,
         martingale = rM,
         martingalet = rtM
  )
  
} #end function


#'Envelope
#'
#'@description Computes simulation envelopes of a tobit-t model.
#'
#'@param model an object of class "tobit" as fitted by tobit.
#'@param res character string indicating the type of residual 
#'desired. Possible values are "deviance", "martingale" and "martingale type".
#'@param nboot Number of simulated point patterns to be generated when computing the envelopes.
#'@param alpha the confidence level required. The default is to find 95 confidence envelopes.
#'@param intercept logical. Should an intercept be included in the null model?
#'@param td assumed distribution for the dependent variable y. 
#'
#'
#'@export

envelope <- function(model,res="martingalet",nboot = 19,alpha=0.05,tau=0,intercept = "TRUE",td="t")
{
  n  <- summary(model)$n[1]
  y  <- as.numeric(model$y)[1:n]
  nu <- model$parms
  muhat <- model$linear.predictors
  sigmahat <- model$scale
  deltahat <-(y-muhat)/sigmahat
  X <- model.matrix(model)
  var.explic <- X[,-1]
  
  rD <- residuals(model,res,tau=tau)#Martingale-type residual
  
  alpha1 <- ceiling(nboot*alpha)
  alpha2 <- ceiling(nboot*(1-alpha))
  e <- matrix(0,n,nboot)
  
  for(i in 1:nboot){
    
    ygerado  <- sigmahat*rt(n,nu) + muhat
    n1  <- summary(model)$n[2] #n?mero de obs. cens.
    pc       <- n1/n  #propor??o de obs. cens.
    tau1      <- sort(ygerado)[pc*n]
    yestrela <- ifelse(ygerado>tau1,ygerado,0)
    
    if(intercept == "FALSE") form = yestrela ~ var.explic - 1 else  form = yestrela ~ var.explic
    
    model1 = tobit(form,dist=td)
    
    rD2 <- residuals(model1,res,tau=tau)
    e[,i] <- sort(rD2)
  }
  
  e1<- numeric(n)
  e2<- numeric(n)	
  
  for(j in 1:n){
    
    eo    <- sort(e[j,])
    e1[j] <- eo[alpha1]
    e2[j] <- eo[alpha2]
  }
  a<-  qqnorm(e1,plot.it=FALSE)$x
  a1<-  qqnorm(e1,plot.it=FALSE)$y
  b<-  qqnorm(e2,plot.it=FALSE)$x
  b1<-  qqnorm(e2,plot.it=FALSE)$y
  r<-  qqnorm(rD,plot.it=FALSE)$x
  r1<-  qqnorm(rD,plot.it=FALSE)$y
  
  xx <- c(a,rev(b))
  yy <- c(a1,rev(b1))
  med   <- apply(e,1,mean)
  faixa <- range(rD,e1,e2,med)
  par(mar=c(4.0,4.0,0.1,0.1))   
  plot(r,r1,type = "n", ylim=faixa,axes=FALSE,xlab="",ylab="")
  par(new=T)
  polygon(xx,yy,col="gray80",border=NA)
  par(new=T)
  qqnorm(rD,main="", ylim=faixa, ylab="R",xlab="Q", cex=0.7, pch=18)
  par(new=T)
  qqnorm(e1,axes=F,type="l",main="",ylim=faixa,lty=1, xlab="", ylab="",col="gray")
  par(new=T)
  qqnorm(e2,axes=F,type="l",main="",ylim=faixa,lty=1, xlab="", ylab="",col="gray")
  par(new=T)
  qqnorm(med,axes=F,type="l",main="",ylim=faixa,lty=2, xlab="", ylab="",col="black")
  
}

#' Envelope 
#'
#'@description Computes simulation envelopes of a tobit model.
#'
#'@param model an object of class "tobit" as fitted by tobit.
#'@param res character string indicating the type of residual 
#'desired. Possible values are "deviance", "martingale" and "martingale type".
#'@param nboot Number of simulated point patterns to be generated when computing the envelopes.
#'@param alpha the confidence level required. The default is to find 95 confidence envelopes.
#'@param intercept logical. Should an intercept be included in the null model?
#'
#' @export
#' 

envelope.normal <- function(model,nboot = 19,alpha=0.05,tau=0,intercept = "TRUE")
{
  
  n  <- summary(model)$n[1]
  y  <- as.numeric(model$y)[1:n]
  c <-  (1*(y>tau))
  muhat <- model$linear.predictors
  sigmahat <- model$scale
  deltahat <-(y-muhat)/sigmahat
  X <- model.matrix(model)
  var.explic <- X[,-1]
  
  S <- 1-pnorm(deltahat)
  rM <- c+log(S)
  
  # Res?duo componente do desvio Martingal 
  
  rD<- (sign((y-muhat))*(-2*(rM+c*log(c-rM)))^(1/2))
  
  
  alpha1 <- ceiling(nboot*alpha)
  alpha2 <- ceiling(nboot*(1-alpha))
  e <- matrix(0,n,nboot)
  
  for(i in 1:nboot){
    
    ygerado  <- sigmahat*rnorm(n,0,1)+ muhat
    n1  <- summary(model)$n[2] #n?mero de obs. cens.
    pc       <- n1/n  #propor??o de obs. cens.
    tau1      <- sort(ygerado)[pc*n]
    yestrela <- ifelse(ygerado>tau1,ygerado,0)
    status   <- (1*(ygerado>tau1))
    
    if(intercept == "FALSE") form = yestrela ~ var.explic - 1 else  form = yestrela ~ var.explic;
    
    model1 <- tobit(form)
    muhat1 <- model1$linear.predictors
    sigmahat1 <- model1$scale
    deltahat1 <- (yestrela - muhat1)/sigmahat1
    S1 <-1- pnorm(deltahat1)
    rM1 <- status + log(S1)
    rD2<- (sign((yestrela-muhat1))*(-2*(rM1+status*log(status-rM1)))^(1/2)) # Res?duo componente do desvio Martingal 
    
    e[,i]    <- sort(rD2)
  }
  
  e1<- numeric(n)
  e2<- numeric(n)	
  
  for(j in 1:n){
    
    eo    <- sort(e[j,])
    e1[j] <- eo[alpha1]
    e2[j] <- eo[alpha2]
  }
  
  med   <- apply(e,1,mean)
  faixa <- range(rD,e1,e2,med)
  a<-  qqnorm(e1,plot.it=FALSE)$x
  a1<-  qqnorm(e1,plot.it=FALSE)$y
  b<-  qqnorm(e2,plot.it=FALSE)$x
  b1<-  qqnorm(e2,plot.it=FALSE)$y
  r<-  qqnorm(rD,plot.it=FALSE)$x
  r1<-  qqnorm(rD,plot.it=FALSE)$y
  
  xx <- c(a,rev(b))
  yy <- c(a1,rev(b1))
  med   <- apply(e,1,mean)
  faixa <- range(rD,e1,e2,med)
  par(mar=c(4.0,4.0,0.1,0.1))   
  plot(r,r1,type = "n", ylim=faixa,axes=FALSE,xlab="",ylab="")
  par(new=T)
  polygon(xx,yy,col="gray80",border=NA)
  par(new=T)
  qqnorm(rD,main="", ylim=faixa, ylab="R",
         xlab="Q", cex=0.7, pch=18)
  par(new=T)
  qqnorm(e1,axes=F,type="l",main="",ylim=faixa,lty=1, xlab="", ylab="",col="gray")
  par(new=T)
  qqnorm(e2,axes=F,type="l",main="",ylim=faixa,lty=1, xlab="", ylab="",col="gray")
  par(new=T)
  qqnorm(med,axes=F,type="l",main="",ylim=faixa,lty=2, xlab="", ylab="",col="black")
}
