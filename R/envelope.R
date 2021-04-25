#'@name tobitdiag
#'
#'@aliases residuals
#'
#'@description Calculates martingale, deviance, martingale type residuals for normal tobit and tobit-t models.
#'
#'@usage
#'residuals(model,type,tau=0,dist ="t")
#'
#'@param model an object of class "tobit" as fitted by tobit. The formula should be a symbolic description of a regression model of type \code{y ~ x1 + x2 + ... + xp}.
#'@param type what type of residuals should be used. 
#'Possible type are: "deviance", "martingale" and "martingale type". 
#'@param tau is the censoring point. The default is zero.
#'@param dist assumed distribution for the dependent variable y. The default is student-t distribution.
#'
#'@return 
#'The returned object is a vector with one element for each subject. 
#'
#'@references
#'Barros, M., Galea, M., Leiva, V. and Santos-Neto, M. Generalized tobit models: Diagnostics and application in econometrics. (submitted manuscript)
#'
#'@examples
#'data("PSID1976")
#'y<-PSID1976$wage
#'x1<- PSID1976$age
#'x2<- PSID1976$education
#'x3<- PSID1976$youngkids
#'x4<- PSID1976$oldkids
#'x5<- PSID1976$experience
#'
#'form <- (y ~ x1+x2+x3+x4+x5)
#'mt<- tobit(form, dist="t") #tobit-t model
#'mnormal <- tobit(form) # normal tobit model
#'
#'#Analysis under the normal tobit model
#'rest <- residuals(mnormal,"martingalet",dist="normal")
#'plot(rest,pch=16,cex=0.7, ylim=c(-3,6), ylab="R",xlab="I")
#'abline(h=c(-3,3),lty=2)
#'identify(rest,n=7,cex=1.0)
#'
#'#Analysis under the tobit t model
#'rest <- residuals(mt,"martingalet",dist="t")
#'plot(rest,pch=16,cex=0.7, ylim=c(-3,6), ylab="R",xlab="I")
#'abline(h=c(-3,3),lty=2)
#'
#'
#'@export

residuals <- function(model, type, tau = 0, dist ="t")
{ # begin function
  
  if (summary(model)$n[2] != 0) cens <- 'left' else cens <- 'right'
  
  if (dist == "t") {
    
    n <- summary(model)$n[1] 
    y  <- as.numeric(model$y)[1:n]
    if (cens == 'left') c  <- (1*(y > tau)) else c  <- (1*(y <= tau)) 
    nu <- model$parms
    muhat = model$linear.predictors
    sigmahat <- model$scale
    deltahat <- (y - muhat)/sigmahat
    S <- 1 - pt(deltahat,nu)
    rM <- c + log(S) # Martingale residual
    rtM <- (sign(y - muhat)*(-2*(rM + c*log(c - rM)))^(1/2))#Martingale-type residual
    rDC <- sign(y - muhat)*sqrt(-2*log(S))
  } else{
    n  <- summary(model)$n[1]
    y  <- as.numeric(model$y)[1:n]
    if (cens == 'left') c  <- (1*(y > tau)) else c  <- (1*(y <= tau)) 
    muhat <- model$linear.predictors
    sigmahat <- model$scale
    deltahat <- (y - muhat)/sigmahat
    S <- 1 - pnorm(deltahat)
    rM <- c + log(S)
    rtM <- (sign(y - muhat)*(-2*(rM + c*log(c - rM)))^(1/2))#Martingale-type residual
    rDC <- sign(y - muhat)*sqrt(-2*log(S))
  }
  switch(type,
         deviance = rDC,
         martingale = rM,
         martingalet = rtM
  )
  
} #end function


#'@name tobitdiag
#'
#'@aliases envelope
#'
#'@description Computes simulation envelopes of a normal tobit and tobit-t model.
#'
#'@usage
#'envelope(model,res="martingalet",nboot = 19,alpha=0.05,tau=0,
#'intercept = "TRUE",dist="t",col="gray80",xlab="Q",ylab="R")
#'
#'@param model an object of class "tobit" as fitted by tobit. The formula should be a symbolic description of a regression model of type \code{y ~ x1 + x2 + ... + xp}.
#'@param res character string indicating the type of residual 
#'desired. Possible values are "deviance", "martingale" and "martingale type".
#'@param nboot Number of simulated point patterns to be generated when computing the envelopes.
#'@param alpha the confidence level required. The default is to find 95 confidence envelopes.
#'@param intercept logical. Should an intercept be included in the null model?
#'@param dist assumed distribution for the dependent variable y. 
#'@param xlab a label for the x axis, defaults to a description of x.
#'@param ylab a label for the y axis, defaults to a description of y.
#'@param col The colors for points. The default is gray80.
#'
#'@references
#'Barros, M., Galea, M., Leiva, V. and Santos-Neto, M. Generalized tobit models: Diagnostics and application in econometrics. (submitted manuscript)
#'
#'@examples
#'data("PSID1976")
#'y<-PSID1976$wage
#'x1<- PSID1976$age
#'x2<- PSID1976$education
#'x3<- PSID1976$youngkids
#'x4<- PSID1976$oldkids
#'x5<- PSID1976$experience
#'
#'form <- (y ~ x1+x2+x3+x4+x5)
#'mt<- tobit(form, dist="t") #tobit-t model
#'mnormal <- tobit(form) # tobit normal model
#'
#'#Analysis under the normal tobit model
#'envelope(mnormal,dist = "normal")
#'
#'#Analysis under the tobit-t model
#'envelope(mt,dist = "t")
#'@importFrom AER tobit
#'@importFrom grDevices dev.new
#'@importFrom graphics abline identify par polygon
#'@importFrom stats AIC dnorm dt formula model.matrix pnorm pt qqnorm rnorm rt sd
#'@export

envelope <- function(model, res = "martingalet", nboot = 19, alpha = 0.05, tau = 0, intercept = "TRUE", dist = "t", col = "gray80", xlab = "Q", ylab = "R")
{
  
  if (summary(model)$n[2] != 0) cens <- 'left' else cens <- 'right'
  
  
  if (dist == "t") {
    n  <- summary(model)$n[1]
    y  <- as.numeric(model$y)[1:n]
    nu <- model$parms
    muhat <- model$linear.predictors
    sigmahat <- model$scale
    deltahat <- (y - muhat)/sigmahat
    X <- model.matrix(model)
    var.explic <- X[,-1]
    
    rD <- residuals(model, res, tau = tau, dist = "t")
    
    alpha1 <- ceiling(nboot*alpha)
    alpha2 <- ceiling(nboot*(1 - alpha))
    e <- matrix(0,n,nboot)
    
    for (i in 1:nboot) {
      
      ygerado  <- sigmahat*rt(n,nu) + muhat
      if (cens == 'left') n1  <- summary(model)$n[2]  else n1  <- summary(model)$n[4]
      pc       <- n1/n  
      tau1      <- sort(ygerado)[pc*n]
      if (cens == 'left') yestrela <- ifelse(ygerado > tau1, ygerado,tau) else yestrela <- ifelse(ygerado <= tau1, ygerado,tau) 
      if (intercept == "FALSE") form <- yestrela ~ var.explic - 1 else  form <- yestrela ~ var.explic
      if (cens == 'left') model1 <- tobit(form, left = tau, dist = 't') else model1 <- tobit(form, right = tau, dist = 't')
      rD2   <- residuals(model1, res, tau = tau)
      e[,i] <- sort(rD2)
    }
    
    e1 <- numeric(n)
    e2 <- numeric(n)	
    
    for (j in 1:n) {
      eo    <- sort(e[j,])
      e1[j] <- eo[alpha1]
      e2[j] <- eo[alpha2]
    }
    
    a  <-  qqnorm(e1, plot.it = FALSE)$x
    a1 <-  qqnorm(e1, plot.it = FALSE)$y
    b  <-  qqnorm(e2, plot.it = FALSE)$x
    b1 <-  qqnorm(e2, plot.it = FALSE)$y
    r  <-  qqnorm(rD, plot.it = FALSE)$x
    r1 <-  qqnorm(rD, plot.it = FALSE)$y
    
    xx    <- c(a,rev(b))
    yy    <- c(a1,rev(b1))
    med   <- apply(e,1,mean)
    faixa <- range(rD,e1,e2,med)
    par(mar = c(4.0,4.0,0.1,0.1))   
    plot(r,r1,type = "n", ylim = faixa,axes = FALSE, xlab = "", ylab = "")
    par(new = T)
    polygon(xx, yy, col = col, border = NA)
    par(new = T)
    qqnorm(rD, main = "", ylim = faixa, ylab = ylab,xlab = xlab, cex = 0.7, pch = 18)
    par(new = T)
    qqnorm(e1, axes = F,type = "l",main = "",ylim = faixa, lty = 1, xlab = "", ylab = "", col = "gray")
    par(new = T)
    qqnorm(e2, axes = F, type = "l", main = "",ylim = faixa, lty = 1, xlab = "", ylab = "", col = "gray")
    par(new = T)
    qqnorm(med, axes = F, type = "l", main = "", ylim = faixa, lty = 2, xlab = "", ylab = "", col = "black")
  } else{
    n  <- summary(model)$n[1]
    y  <- as.numeric(model$y)[1:n]
    muhat <- model$linear.predictors
    sigmahat <- model$scale
    deltahat <- (y - muhat)/sigmahat
    X <- model.matrix(model)
    var.explic <- X[,-1]
    rD <- residuals(model,res,tau = tau, dist = "normal")#Martingale-type residual
    alpha1 <- ceiling(nboot*alpha)
    alpha2 <- ceiling(nboot*(1 - alpha))
    e <- matrix(0,n,nboot)
    
    for (i in 1:nboot) {
      ygerado  <- sigmahat*rnorm(n,0,1) + muhat
      if (cens == 'left') n1  <- summary(model)$n[2]  else n1  <- summary(model)$n[4]
      pc       <- n1/n
      tau1      <- sort(ygerado)[pc*n]
      
      if (cens == 'left') yestrela <- ifelse(ygerado > tau1, ygerado,tau) else yestrela <- ifelse(ygerado <= tau1, ygerado,tau) 
      if (intercept == "FALSE") form <- yestrela ~ var.explic - 1 else  form <- yestrela ~ var.explic
      if (cens == 'left') model1 <- tobit(form, left = tau) else model1 <- tobit(form, right = tau)
      rD2   <- residuals(model1, res, tau = tau, dist = "normal")
      e[,i]    <- sort(rD2)
    }
    
    e1 <- numeric(n)
    e2 <- numeric(n)	
    
    for(j in 1:n){
      
      eo    <- sort(e[j,])
      e1[j] <- eo[alpha1]
      e2[j] <- eo[alpha2]
    }
    
    med   <- apply(e,1,mean)
    faixa <- range(rD,e1,e2,med)
    a     <-  qqnorm(e1,plot.it = FALSE)$x
    a1    <-  qqnorm(e1,plot.it = FALSE)$y
    b     <-  qqnorm(e2,plot.it = FALSE)$x
    b1    <-  qqnorm(e2,plot.it = FALSE)$y
    r     <-  qqnorm(rD,plot.it = FALSE)$x
    r1    <-  qqnorm(rD,plot.it = FALSE)$y
    
    xx <- c(a,rev(b))
    yy <- c(a1,rev(b1))
    med   <- apply(e,1,mean)
    faixa <- range(rD,e1,e2,med)
    par(mar = c(4.0,4.0,0.1,0.1))   
    plot(r,r1,type = "n", ylim = faixa, axes = FALSE, xlab = "", ylab = "")
    par(new = T)
    polygon(xx,yy,col = "gray80",border = NA)
    par(new = T)
    qqnorm(rD,main = "", ylim = faixa, ylab = ylab,
           xlab = xlab, cex = 0.7, pch = 18)
    par(new = T)
    qqnorm(e1,axes = F,type = "l",main = "",ylim = faixa,lty = 1, xlab = "", ylab = "", col = "gray")
    par(new = T)
    qqnorm(e2, axes = F,type = "l", main = "", ylim = faixa, lty = 1, xlab = "", ylab = "", col = "gray")
    par(new = T)
    qqnorm(med, axes = F, type = "l", main = "", ylim = faixa, lty = 2, xlab = "", ylab = "", col = "black") 
    
    
  }
  
}


