% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagtobit.R
\name{AICc}
\alias{AICc}
\title{Akaike's An Information Criterion}
\usage{
AICc(object, object1 = NULL)
}
\arguments{
\item{object}{a fitted model object for which there exists a logLik method to extract the corresponding log-likelihood, or an object inheriting from class logLik.}

\item{object1}{optionally more fitted model objects.}
}
\value{
If just one object is provided, a numeric value with the corresponding AICc.

If two objects are provided, a data.frame with rows corresponding to the objects and columns 
representing the number of 
parameters in the model (df) and the AICc.
}
\description{
AICc is AIC with a correction for finite sample sizes.
}
\examples{
 
data("PSID1976")
y<-PSID1976$wage
x1<- PSID1976$age
x2<- PSID1976$education
x3<- PSID1976$youngkids
x4<- PSID1976$oldkids
x5<- PSID1976$experience

form <- (y ~ x1+x2+x3+x4+x5)
mt<- tobit(form, dist="t") #tobit-t model
mnormal <- tobit(form) # normal tobit model

AICc(mnormal,mt)

}

