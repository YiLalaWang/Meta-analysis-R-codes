
# This is the code that contains the function for: 
## (1) doing meta-analysis based on the metafor package, 
## (2) calculating the 95%CI for tau2, and 
## (3) organizing all meta-analytic results into table format. 

# rawdata should at least include r(effect sizes) and n(sample sizes) for all studies
## sample/study info (citation) is strongly recommended

resMeta=function(rawdata){
  require(metafor)
  rawdata=data.frame(rawdata)
  ## meta-analysis using RMEL estimator for heterogeneity
  rawdata.rescale=escalc(measure="ZCOR",ri=r,ni=n,data=rawdata)
  x=rma(measure="ZCOR",yi,vi,data=rawdata.rescale,level=95,slab=rawdata.rescale$citation)
  ci=predict(x,level=95,transf=transf.ztor) # 95% CI
  cv=predict(x,level=90,transf=transf.ztor) # 90% CV
  rm=round(ci$pred,2)
  ci.l=round(ci$ci.lb,2)
  ci.u=round(ci$ci.ub,2)
  ci=paste(ci.l,",",ci.u)
  ci.diff=ci.u-ci.l
  cv.l=round(cv$cr.lb,2)
  cv.u=round(cv$cr.ub,2)
  cv=paste(cv.l,",",cv.u)
  cv.diff=cv.u-cv.l
  k=as.numeric(nrow(rawdata))
  N=sum(rawdata$n)
  
  ## organizing results into table format
  res.table=data.frame(matrix(nrow=10, ncol=1))
  colnames(res.table)="res.meta"
  rownames(res.table)=c("N","K","rm","95%CI","CI width",
                        "90%CV","CV width","I2","Q(df,p)","tau2(95%CI)")
  res.table[c("N","K","rm","I2"),1]=c(N,k,rm,round(x$I2,2))
  res.table[c("95%CI","CI width"),1]=c(ci,ci.diff)
  res.table[c("90%CV","CV width"),1]=c(cv,cv.diff)
  res.table["Q(df,p)",1]=paste(round(x$QE,2),"(",x$k-x$m,",",round(x$QEp,2),")")
  ## 95%CI of tau2 formulas see Borenstein et al. (2009) chapter16
  if(x$QE>(x$k-x$m+1)){
    B=0.5*((log(x$QE)-log(x$k-x$m))/(sqrt(2*x$QE)-sqrt(2*(x$k-x$m)-1))) #formula (16.14)
  } else if(x$QE<=(x$k-x$m+1)){
    B=sqrt(1/(2*(x$k-x$m-1)*(1-(1/(3*(x$k-x$m-1)^2))))) #formula(16.15)
  }
  C=(x$QE-x$k-x$m-1)/x$tau2 #formula(16.5)
  L=exp(0.5*log(x$QE/(x$k-x$m))-1.96*B) #formula(16.16)
  tau2.CIL=(x$k-x$m)*(L^2-1)/C #formula(16.18)
  tau2.CIL=round(tau2.CIL,2)
  U=exp(0.5*log(x$QE/(x$k-x$m))+1.96*B) #formula(16.17)
  tau2.CIU=(x$k-x$m)*(U^2-1)/C #formula(16.19)
  tau2.CIU=round(tau2.CIU,2)
  res.table["tau2(95%CI)",1]=paste(round(x$tau2,2),"(",tau2.CIL,",",tau2.CIU,")")
  return(res.table)
}


