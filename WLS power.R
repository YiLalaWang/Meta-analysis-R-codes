###########################################
## Fixed- & Mixed-model meta-WLS procedures
## This is for a priori power estimation 
## when raw data are available
###########################################

## Raw data (in .csv) must be organized in column in the following order: 
# column 1: sample size (n) or variance (v),  
# column 2: intercept (b0, note that all b0 are filled with 1), 
# column 3 and on: values of all predictors (starting from b1).
rawdata<-read.table(file.choose(),header=T,sep=",") 
alpha=0.05 # Please specify the desired alpha level for the test
p=3 # Please specify the number of predictors (betas)
## Please specify the expected or actual values of all predictor betas (starting from b1).
# If a priori power, use the expected values
# if post hoc power, use the actual values
beta.pred=c(0.25, 0.25, 0) 
b <- 1 # Please specify the beta to be tested for power
desire.power<-0.80 # Please specify the desired power to estimate the desired beta value

#### This is to calculate the lumda value for the power of test for the fit of the fixed-effect regression model
# Specify the amount of heterogeneity (tau2) across effect sizes in fixed effect
hetero.fit = "small" #Choose one out of: small, moderate, large

##### This is for the power of test for residual variance tau2=0 
# If tau2=0, then it should be a fixed effect model rather than a mixed-effect model
## Specify the amount of heterogeneity in mixed effect (residual variance)
## If tau2 is known, remove the first "#" sign in the line below and specify the tau2 value 
#tau2=.036 
## If tau2 is unknown, it can be calculated using Hedges & Pigott's (2004) formulas 
# The calculation is based on an informed assumption of the magnitude of heterogeneity
hetero.res = "moderate" # Choose one out of: small, moderate, large, tau2 (if tau2 is known)

######### No specification/changes are needed after this point #############
v<-if(is.null(rawdata$n)) rawdata$v else (1/(rawdata$n-3))
k<-nrow(rawdata) 
tau2<-if(hetero.res=="small"){
  rep(mean(v)/3,k)
} else if(hetero.res=="moderate"){
  rep(2*mean(v)/3,k)
} else if(hetero.res=="large"){
  rep(mean(v),k)
} else if(hetero.res=="tau2"){
  rep(tau2,k)
}
V<-diag(v+tau2) 
Var<-diag(v)
w<-1/v
X<-data.matrix(rawdata[1:k, 2:(2+p)],rownames.force=NA)
beta.pred<-as.matrix(beta.pred)

################ The power of omnibus test for all betas ###################
c.omni<-qchisq(1-alpha,p,ncp=0,lower.tail=T, log.p=F)

#### For fixed-effects (FE) model
Sum.fe<-solve(t(X) %*% solve(Var) %*% X) # Matrix for all betas
Sum.fe.omni<-Sum.fe[2:(1+p),2:(1+p)] # Matrix for all predictor betas (except beta0)
lumda.fe<-t(beta.pred) %*% solve(Sum.fe.omni) %*% beta.pred
phi.fe.omni<-pchisq(c.omni,p,ncp=lumda.fe,lower.tail=T,log.p=F) 
# h is the probability to falsely determine 
# that all predictor betas equal 0 when not all of them are 0
power.fe.omni<-1-phi.fe.omni  
# This is the power to reject that all predictor betas equal 0

#### For mixed-effects (ME) model
Sum.me<-solve(t(X) %*% solve(V) %*% X) # Matrix for all betas
Sum.me.omni<-Sum.me[2:(1+p),2:(1+p)] # Matrix for all predictor betas (except beta0)
lumda.me<-t(beta.pred) %*% solve(Sum.me.omni) %*% beta.pred
phi.me.omni<-pchisq(c.omni, df=p, ncp=lumda.me,lower.tail=T,log.p=F) 
# h is the probability to falsely determine 
# that all predictor betas equal 0 when not all of them are 0
power.me.omni<-1-phi.me.omni  
############################################################################


############ Power of test for individual regression coefficients ##########
c.ind.1tail<-qnorm(1-alpha,mean=0,sd=1,lower.tail=T, log.p=F)
c.ind.2tail<-qnorm(1-alpha/2,mean=0,sd=1,lower.tail=T, log.p=F)

#### For fixed-effects (FE) model
var.fe.b<-Sum.fe[1+b,1+b]
# One-sided power
pred.b<-beta.pred[b,1]
phi.b.1tail.fe<-pnorm(c.ind.1tail-abs(pred.b)/sqrt(var.fe.b),mean=0,sd=1,lower.tail=T, log.p=F)
power.b.1side.fe<-1-phi.b.1tail.fe
# Two-sided power
phi.b.2tail.fe<-pnorm(c.ind.2tail-abs(pred.b)/sqrt(var.fe.b),mean=0,sd=1,lower.tail=T, log.p=F)
phi.b.2tail.neg.fe<-pnorm(-c.ind.2tail-abs(pred.b)/sqrt(var.fe.b),mean=0,sd=1,lower.tail=T, log.p=F)
power.b.2side.fe<-1-phi.b.2tail.fe+phi.b.2tail.neg.fe

#### For mixed-effects (ME) model
var.me.b<-Sum.me[1+b,1+b]
# One-sided power
pred.b<-beta.pred[b,1]
phi.b.1tail.me<-pnorm(c.ind.1tail-abs(pred.b)/sqrt(var.me.b),mean=0,sd=1,lower.tail=T, log.p=F)
power.b.1side.me<-1-phi.b.1tail.me
# Two-sided power
phi.b.2tail.me<-pnorm(c.ind.2tail-abs(pred.b)/sqrt(var.me.b),mean=0,sd=1,lower.tail=T, log.p=F)
phi.b.2tail.neg.me<-pnorm(-c.ind.2tail-abs(pred.b)/sqrt(var.me.b),mean=0,sd=1,lower.tail=T, log.p=F)
power.b.2side.me<-1-phi.b.2tail.me+phi.b.2tail.neg.me

## Desired beta value to reach desired power 
phi.bpower<-qnorm(1-desire.power,mean=0,sd=1,lower.tail=T,log.p=F)
desire.b.1tail.fe<-sqrt(var.fe.b)*(c.ind.1tail-phi.bpower)
desire.b.2tail.fe<-sqrt(var.fe.b)*(c.ind.2tail-phi.bpower)
desire.b.1tail.me<-sqrt(var.me.b)*(c.ind.1tail-phi.bpower)
desire.b.2tail.me<-sqrt(var.me.b)*(c.ind.2tail-phi.bpower)
#################################################################################


############## Power of test for Goodness of fit for the FE model ##############
c.fit<-qchisq(1-alpha,k-p-1,ncp=0,lower.tail=T, log.p=F)
lumda.fit<-if(hetero.fit=="small") {
  (k-p-1)/3 
  } else if(hetero.fit=="moderate") {
    (2*(k-p-1)/3) 
    }else {
      (k-p-1)}
phi.fit<-pchisq(c.fit,k-p-1,ncp=lumda.fit,lower.tail=T,log.p=F)
power.fit<-1-phi.fit
#################################################################################


############# Power of test for residual variance (tau2) ########################
############# This is a test of assumption for the ME model #####################
############# i.e., whether the magnitude of tau2 matches the assumtption) ######
M<-solve(Var) %*% X %*% solve(t(X) %*% solve(Var) %*% X) %*% t(X) %*% solve(Var)
Var.2<-solve(Var) %*% solve(Var)
V2<-V %*% V
V.2<-solve(V) %*% solve(V)
Var.Q<-2*(sum(diag(Var.2 %*% V2))-2*sum(diag(M %*% solve(Var) %*% V2))+sum(diag(M %*% V %*% M %*% V)))
a<-sum(w)-sum(diag(solve(t(X) %*% solve(Var) %*% X) %*% t(X) %*% Var.2 %*% X))
Miu.Q<-a*mean(tau2)+k-p-1
r<-Var.Q/(2*Miu.Q)
s<-(2*(Miu.Q^2))/Var.Q
# Power of the test to reject that tau2=0
c.tau2<-qchisq(1-alpha, k-p-1, ncp=0,lower.tail=T, log.p=F)
phi.tau2<-pchisq(c.tau2/r,s,ncp=0,lower.tail=T,log.p=F)
power.tau2<-1-phi.tau2
#################################################################################

######################## Round all results ######################################
power.fe.omni<-round(power.fe.omni,digits=2) 
power.me.omni<-round(power.me.omni,digits=2)
beta.1tail.fe<-paste("power:",round(power.b.1side.fe,digits=2),"/desired beta:",round(desire.b.1tail.fe,digits=2))
beta.2tail.fe<-paste("power:",round(power.b.2side.fe,digits=2),"/desired beta:",round(desire.b.2tail.fe,digits=2))
beta.1tail.me<-paste("power:",round(power.b.1side.me,digits=2),"/desired beta:",round(desire.b.1tail.me,digits=2))
beta.2tail.me<-paste("power:",round(power.b.2side.me,digits=2),"/desired beta:",round(desire.b.2tail.me,digits=2))
power.fe.fit<-round(power.fit,digits=2)
power.me.tau2<-round(power.tau2,digits=2)
#################################################################################

######################## Organizing all results #################################
power.fe.omni
power.me.omni
power.fe.fit
power.me.tau2
beta.1tail.fe
beta.2tail.fe
beta.1tail.me
beta.2tail.me
#################################################################################