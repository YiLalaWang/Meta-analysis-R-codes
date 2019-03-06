##########################################################################
## Fixed- & Mixed-model meta-WLS procedures (Multiple Regression analogues) 
## based on Hedges & Pigott (2004) 
## This is for a priori and post hoc power estimation 
## when raw data (i.e., sample sizes) are available
## Please contact Yi Wang (wangyilalala@gmail.com) for questions
###########################################################################
setwd("G:/SkyDrive/methods&stats/Meta-analysis/Calculation/My own files/R code power for meta moderation")

WLSpower = function(x, b, tau2, beta.pred){
  rawdata = as.data.frame(x)
  alpha = 0.05 # Desired alpha is set as 0.5 by default
  desire.power = 0.80 # Desired power is set as 0.80 by default
  v=rawdata$vi 
  k=nrow(rawdata) 
  w=1/v
  beta.pred=as.matrix(beta.pred)
  p = nrow(beta.pred)
  q = ncol(rawdata)
  X=data.matrix(rawdata[1:k, (q-p):q],rownames.force=NA)
  theta=data.matrix(rawdata$es)
  Var=diag(v)
  Var.2=solve(Var) %*% solve(Var)
  a=sum(w)-sum(diag(solve(t(X) %*% solve(Var) %*% X) %*% t(X) %*% Var.2 %*% X))
  
  tau2=rep(tau2,k)
  V=diag(v+tau2) 
  V2=V %*% V
  V.2=solve(V) %*% solve(V)
  
  ################ The power of omnibus test for all betas ###################
  c.omni=qchisq(1-alpha,p,ncp=0,lower.tail=T, log.p=F)
  
  #### For fixed-effects (FE) model
  Sum.fe=solve(t(X) %*% solve(Var) %*% X) # Matrix for all betas
  Sum.fe.omni=Sum.fe[2:(1+p),2:(1+p)] # Matrix for all predictor betas (except beta0)
  lumda.fe=t(beta.pred) %*% solve(Sum.fe.omni) %*% beta.pred
  phi.fe.omni=pchisq(c.omni,p,ncp=lumda.fe,lower.tail=T,log.p=F) 
  # h is the probability to falsely determine 
  # that all predictor betas equal 0 when not all of them are 0
  power.fe.omni=1-phi.fe.omni  
  # This is the power to reject that all predictor betas equal 0
  
  #### For mixed-effects (ME) model
  Sum.me=solve(t(X) %*% solve(V) %*% X) # Matrix for all betas
  Sum.me.omni=Sum.me[2:(1+p),2:(1+p)] # Matrix for all predictor betas (except beta0)
  lumda.me=t(beta.pred) %*% solve(Sum.me.omni) %*% beta.pred
  phi.me.omni=pchisq(c.omni, df=p, ncp=lumda.me,lower.tail=T,log.p=F) 
  # h is the probability to falsely determine 
  # that all predictor betas equal 0 when not all of them are 0
  power.me.omni=1-phi.me.omni  
  ############################################################################
  
  
  ############ Power of test for individual regression coefficients ##########
  c.ind.1tail=qnorm(1-alpha,mean=0,sd=1,lower.tail=T, log.p=F)
  c.ind.2tail=qnorm(1-alpha/2,mean=0,sd=1,lower.tail=T, log.p=F)
  
  #### For fixed-effects (FE) model
  var.fe.b=Sum.fe[1+b,1+b]
  # One-sided power
  pred.b=beta.pred[b,1]
  phi.b.1tail.fe=pnorm(c.ind.1tail-abs(pred.b)/sqrt(var.fe.b),mean=0,sd=1,lower.tail=T, log.p=F)
  power.b.1side.fe=1-phi.b.1tail.fe
  # Two-sided power
  phi.b.2tail.fe=pnorm(c.ind.2tail-abs(pred.b)/sqrt(var.fe.b),mean=0,sd=1,lower.tail=T, log.p=F)
  phi.b.2tail.neg.fe=pnorm(-c.ind.2tail-abs(pred.b)/sqrt(var.fe.b),mean=0,sd=1,lower.tail=T, log.p=F)
  power.b.2side.fe=1-phi.b.2tail.fe+phi.b.2tail.neg.fe
  
  #### For mixed-effects (ME) model
  var.me.b=Sum.me[1+b,1+b]
  # One-sided power
  pred.b=beta.pred[b,1]
  phi.b.1tail.me=pnorm(c.ind.1tail-abs(pred.b)/sqrt(var.me.b),mean=0,sd=1,lower.tail=T, log.p=F)
  power.b.1side.me=1-phi.b.1tail.me
  # Two-sided power
  phi.b.2tail.me=pnorm(c.ind.2tail-abs(pred.b)/sqrt(var.me.b),mean=0,sd=1,lower.tail=T, log.p=F)
  phi.b.2tail.neg.me=pnorm(-c.ind.2tail-abs(pred.b)/sqrt(var.me.b),mean=0,sd=1,lower.tail=T, log.p=F)
  power.b.2side.me=1-phi.b.2tail.me+phi.b.2tail.neg.me
  
  ## Desired beta value to reach desired power 
  phi.bpower=qnorm(1-desire.power,mean=0,sd=1,lower.tail=T,log.p=F)
  desire.b.1tail.fe=sqrt(var.fe.b)*(c.ind.1tail-phi.bpower)
  desire.b.2tail.fe=sqrt(var.fe.b)*(c.ind.2tail-phi.bpower)
  desire.b.1tail.me=sqrt(var.me.b)*(c.ind.1tail-phi.bpower)
  desire.b.2tail.me=sqrt(var.me.b)*(c.ind.2tail-phi.bpower)
  #################################################################################
  
  
  ############## Power of test for Goodness of fit for the FE model ##############
  c.fit=qchisq(1-alpha,k-p-1,ncp=0,lower.tail=T, log.p=F)
  lumda.fit=t(theta)%*%(solve(Var)-solve(Var) %*% X %*% solve(t(X) %*% solve(Var) %*% X) %*% t(X) %*% solve(Var))%*%theta
  phi.fit=pchisq(c.fit,k-p-1,ncp=lumda.fit,lower.tail=T,log.p=F)
  power.fit=1-phi.fit
  #################################################################################
  
  
  ############# Power of test for residual variance (tau2) ########################
  ############# This is a test of assumption for the ME model #####################
  ############# i.e., whether the magnitude of tau2 matches the assumtption) ######
  M=solve(Var) %*% X %*% solve(t(X) %*% solve(Var) %*% X) %*% t(X) %*% solve(Var)
  Var.Q=2*(sum(diag(Var.2 %*% V2))-2*sum(diag(M %*% solve(Var) %*% V2))+sum(diag(M %*% V %*% M %*% V)))
  Miu.Q=a*mean(tau2)+k-p-1
  r=Var.Q/(2*Miu.Q)
  s=(2*(Miu.Q^2))/Var.Q
  # Power of the test to reject that tau2=0
  c.tau2=qchisq(1-alpha, k-p-1, ncp=0,lower.tail=T, log.p=F)
  phi.tau2=pchisq(c.tau2/r,s,ncp=0,lower.tail=T,log.p=F)
  power.tau2=1-phi.tau2
  #################################################################################
  
  ######################## Round all results ######################################
  power.fe.omni=round(power.fe.omni,digits=2) 
  power.me.omni=round(power.me.omni,digits=2)
  power.fe.fit=round(power.fit,digits=2)
  power.me.tau2=round(power.tau2,digits=2)
  
  power.b.1side.fe=round(power.b.1side.fe,digits=2)
  desire.b.1tail.fe=round(desire.b.1tail.fe,digits=2)
  power.b.2side.fe=round(power.b.2side.fe,digits=2)
  desire.b.2tail.fe=round(desire.b.2tail.fe,digits=2)
  
  power.b.1side.me=round(power.b.1side.me,digits=2)
  desire.b.1tail.me=round(desire.b.1tail.me,digits=2)
  power.b.2side.me=round(power.b.2side.me,digits=2)
  desire.b.2tail.me=round(desire.b.2tail.me,digits=2)
  #################################################################################
  
  ################## Organizing all results #############################
  WLSpower.table=data.frame(matrix(nrow = 6,ncol = 2))
  colnames(WLSpower.table)=c("fixed effects","mixed effects")
  rownames(WLSpower.table)=c("omnibus test (QM)", "Model fit (QE)",
                             "actual power for beta (one-tail test)",
                             "beta needed for desired power (one-tail test)",
                             "actual power for beta (two-tail test)",
                             "beta needed for desired power (two-tail test)")
  WLSpower.table["fixed effects"]=c(power.fe.omni, power.fe.fit,
                                    power.b.1side.fe, desire.b.1tail.fe,
                                    power.b.2side.fe, desire.b.2tail.fe)
  WLSpower.table["mixed effects"]=c(power.me.omni, power.me.tau2,
                                    power.b.1side.me, desire.b.1tail.me,
                                    power.b.2side.me, desire.b.2tail.me)
  return(WLSpower.table)
}

#################### Instructions for specifying values in the WLSpower function ##########################

# use x to specify the name for the raw data file
### Raw data must be in .csv format, and contains the following columns: 
### column 1, entitled as vi: variance of effect sizes,  
### column 2, entitled as b0: intercept of meta-regression model, filled with 1 for all rows 
### column 3, entitled as b1: moderator values corresponding with n (or v)
### column 4 and on: similar with column 3, entitled as b2, b3 and so on

# use b to specify which beta you want to analyze for power (1 if testing for beta1, etc.)

# use tau2 to specify the amount of unaccounted heterogeneity in mixed effect
### tau2 is R metafor default output
### Some software may term tau2 as "random effect variance component v"

# beta.pred needs to be specified before using the WLSpower function for power calculation
### If a priori power, fill in the parentheses using the expected beta estimates from meta-regression
### if post hoc power, fill in the parentheses using the actual beta estimates from meta-regression
beta.pred=c(b1, b2, b3) 

# For an example using Table 2 data in Hedges & Pigott (2004), see below: 

rawdata=read.table(file="Table 2.csv",header=TRUE,sep=",")
res=rma(yi=es,vi,data=rawdata,level=95,mods=~b1+b2+b3)
# b1,b2,b3 are columns 2-4 in the "X matrix" in Table 2 (assuming that column 1 represents b0)
beta.pred=res$b[c("b1","b2","b3"),1] 
WLSpower(x=rawdata, b=1, tau2=res$tau2, beta.pred)
