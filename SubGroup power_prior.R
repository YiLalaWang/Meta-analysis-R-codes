##########################################################################
## Fixed- & Mixed-model sub-group analyses (ANOVA analogues) based on Hedges & Pigott (2004) 
## This is for a priori and post hoc power estimation 
## when raw data (i.e., sample sizes) are available
## Please contact Yi Wang (wangyilalala@gmail.com) for questions
###########################################################################

SubGrouppower.prior=function(x,es.g,hetero){
  rawdata = as.data.frame(x)
  require(metafor)
  rawdata$w=1/rawdata$vi
  k=nrow(rawdata) 
  r.g=as.matrix(es.g)
  p=nrow(r.g)
  df_p=p-1
  alpha=0.05
  m=round(k/p,0) # assuming that all groups have the same K
  
  res.lumda=data.frame(matrix(nrow = p,ncol=4))
  colnames(res.lumda)=c("vi","rg","rg_rm.diff","vi_me")
  res.lumda[,"vi"]=1/(sum(rawdata$w)/3) 
  # assuming the same effect sizes thus same between-study within-group variance
  res.lumda[,"rg"]=r.g
  est.rm=sum(res.lumda$rg/res.lumda$vi)/sum(1/res.lumda$vi)
  res.lumda[,"rg_rm.diff"]=res.lumda$rg-est.rm
  res.lumda[,"a"]=sum(rawdata$w)/3-sum(rawdata$w^2)/sum(rawdata$w)
  a_all=sum(res.lumda$a)
  
  tau2=if(hetero=="small"){
    mean(rawdata$vi)/3
  } else if(hetero=="medium"){
    2*mean(rawdata$vi)/3
  } else if(hetero=="large"){
    mean(rawdata$vi)
  } 
  rawdata$vi_me=rawdata$vi+tau2
  rawdata$w_me=1/rawdata$vi_me
  res.lumda[,"vi_me"]=1/(sum(rawdata$w_me)/3)

  ############# Power of omnibus test #############
  c.alpha.omni=qchisq(1-alpha,df_p,ncp=0,lower.tail=T, log.p=F)
  
  ## fixed-effects
  ratio.R=var(r.g)/mean(rawdata$vi)
  lumda.FE.omni=m*(p-1)*ratio.R
  power.FE.omni=1-pchisq(c.alpha.omni, df=df_p, ncp=lumda.FE.omni,lower.tail=T,log.p=F)
  
  ## mixed effects
  lumda.ME.omni=if(hetero=="small"){
    3*ratio.R*m*(p-1)/4
  } else if (hetero=="medium"){
    3*ratio.R*m*(p-1)/5
  } else if (hetero=="large"){
    ratio.R*m*(p-1)/2
  } 
  power.ME.omni=1-pchisq(c.alpha.omni, df=df_p, ncp=lumda.ME.omni, lower.tail=T,log.p=F)
  
  ############# power of group contrasts tests ############
  require(Surrogate)
  c=RandVec(a=-1,b=1,s=0,n=p,m=1,Seed=1)$RandVecOutput # generate a random vector with row=p & sum=0
  res.lumda=cbind(res.lumda,c)
  names(res.lumda)[ncol(res.lumda)]="c"
  gamma=abs(sum(res.lumda$c*res.lumda$rg))
  
  #### Power for planned (a priori) comparison
  c.contrast.1tail=qnorm(1-alpha,mean=0,sd=1,lower.tail=T, log.p=F)
  c.contrast.2tail=qnorm(1-alpha/2,mean=0,sd=1,lower.tail=T, log.p=F)
  
  ## Fixed effects model
  # One-sided power
  v.contrast.FE=sum((res.lumda$c^2)*res.lumda$vi)
  phi.contrast.1tail.FE=pnorm(c.contrast.1tail-gamma/sqrt(v.contrast.FE),mean=0,sd=1,lower.tail=T, log.p=F)
  power.contrast.1side.FE=1-phi.contrast.1tail.FE
  
  # Two-sided power
  phi.contrast.2tail.FE=pnorm(c.contrast.2tail-gamma/sqrt(v.contrast.FE),mean=0,sd=1,lower.tail=T, log.p=F)
  phi.contrast.2tail.FE.neg=pnorm(-c.contrast.2tail-gamma/sqrt(v.contrast.FE),mean=0,sd=1,lower.tail=T, log.p=F)
  power.contrast.2side.FE=1-phi.contrast.2tail.FE+phi.contrast.2tail.FE.neg
  
  ## mixed effects model
  # One-sided power
  v.contrast.ME=sum((res.lumda$c^2)*res.lumda$vi_me)
  phi.contrast.1tail.ME=pnorm(c.contrast.1tail-gamma/sqrt(v.contrast.ME),mean=0,sd=1,lower.tail=T, log.p=F)
  power.contrast.1side.ME=1-phi.contrast.1tail.ME
  # Two-sided power
  phi.contrast.2tail.ME=pnorm(c.contrast.2tail-gamma/sqrt(v.contrast.ME),mean=0,sd=1,lower.tail=T, log.p=F)
  phi.contrast.2tail.ME.neg=pnorm(-c.contrast.2tail-gamma/sqrt(v.contrast.ME),mean=0,sd=1,lower.tail=T, log.p=F)
  power.contrast.2side.ME=1-phi.contrast.2tail.ME+phi.contrast.2tail.ME.neg
  
  ############# Power of within-group heterogeneity (QE) #############
  c.fit=qchisq(1-alpha,k-p,ncp=0,lower.tail=T, log.p=F)
  lumda.fit=if(hetero=="small") {
    (k-p)/3 
  } else if(hetero=="medium") {
    (2*(k-p)/3) 
  } else if(hetero=="large"){
    (k-p)
  } 
  phi.fit=pchisq(c.fit,k-p,ncp=lumda.fit,lower.tail=T,log.p=F)
  power.fit=1-phi.fit
  
  ############# Power of the Between-Studies variance component (tau2) ########################
  Miu.Q=a_all*tau2+k-p
  res.lumda[,"S"]=sum(rawdata$w)/3
  res.lumda[,"S2"]=sum(rawdata$w^2)/3
  res.lumda[,"S3"]=sum(rawdata$w^3)/3
  res.lumda[,"b"]=sum(m-1+2*(res.lumda[,"S"]-res.lumda[,"S2"]/res.lumda[,"S"])*tau2+
                              (res.lumda[,"S2"]-2*res.lumda[,"S3"]/res.lumda[,"S"]+
                                 res.lumda[,"S2"]^2/(res.lumda[,"S"]^2))*(tau2^2))
  
  b_all=sum(res.lumda$b)
  sigma2.Q=2*b_all
  r=sigma2.Q/(2*Miu.Q)
  s=2*(Miu.Q^2)/sigma2.Q
  # Power of the test to reject that tau2=0
  c.tau2=qchisq(1-alpha, k-p, ncp=0,lower.tail=T, log.p=F)
  phi.tau2=pchisq(c.tau2/r,s,ncp=0,lower.tail=T,log.p=F)
  power.tau2=1-phi.tau2
  
  #################### round all results ##########################
  
  power.fe.omni=round(power.FE.omni,digits=2) 
  power.fe.fit=round(power.fit,digits=2)
  power.contrast.1side.FE=if(is.null(power.contrast.1side.FE)){NULL}else{round(power.contrast.1side.FE,digits=2)}
  power.contrast.2side.FE=if(is.null(power.contrast.2side.FE)){NULL}else{round(power.contrast.2side.FE,digits=2)}

  power.me.omni=round(power.ME.omni,digits=2)
  power.me.tau2=round(power.tau2,digits=2)
  power.contrast.1side.ME=if(is.null(power.contrast.1side.ME)){NULL}else{round(power.contrast.1side.ME,digits=2)}
  power.contrast.2side.ME=if(is.null(power.contrast.2side.ME)){NULL}else{round(power.contrast.2side.ME,digits=2)}

  ################## Organizing all results #############################
  SubGrouppower.table.prior=data.frame(matrix(nrow = 4,ncol = 2))
  colnames(SubGrouppower.table.prior)=c("fixed effects","mixed effects")
  rownames(SubGrouppower.table.prior)=c("omnibus test (QM)", "Model fit (QE)",
                                        "group contrast (one-tail test)",
                                        "group contrast (two-tail test)")
  SubGrouppower.table.prior["fixed effects"]=c(power.fe.omni, power.fe.fit,
                                               power.contrast.1side.FE, power.contrast.2side.FE)
  SubGrouppower.table.prior["mixed effects"]=c(power.me.omni, power.me.tau2,
                                               power.contrast.1side.ME, power.contrast.2side.ME)
  return(SubGrouppower.table.prior)
}

############### Instructions for specifying values in the SubGrouppower function ##########################

# use x to specify the name for the raw data file
### Raw data must be organized with the following columns: 
### vi: variance of effect sizes
### es: effect sizes 

# use "es.g" to specify the expected mean effect size for each group 
### e.g., es.g=c(0.5,0.3,0.3) means that there are 3 groups, each with expected mean effect size of 0.5, 0.3, and 0.3

# use "hetero" to specify the amount of heterogeneity (tau2) across effect sizes 
### To specify hetero, choose one out the options (including the quotes in your code): 
### "small", "moderate", "large"

# For an example using Table 1 data in Hedges & Pigott (2004), see below: 

rawdata=read.table(file="Table 1.csv",header=TRUE,sep=",") 
rawdata$yi=(1-3/(4*(rawdata$n-2)-1))*rawdata$es 
rawdata$vi=2/(rawdata$n/2)+rawdata$yi^2/(2*rawdata$n)
# yi & vi converted using standard conversion for standardized mean difference
SubGrouppower.prior(rawdata,es.g=c(.48,.6,.6),hetero="large")
