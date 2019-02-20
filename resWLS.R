
# This is the code for organizing WLS regression results
## The code also includes the calculation for 95%CI of tau2
## input x should be the WLS results from the rma function based on metafor package

resWLS.table=function(x){
  res.table=data.frame(matrix(nrow = x$p,ncol=2))
  colnames(res.table)=c("beta (SE)","95%CI (LL,UL)")
  rownames(res.table)=row.names(x$beta)
  res.table["beta (SE)"]=paste(round(x$beta,2)," (",round(x$se,2),")")
  res.table["95%CI (LL,UL)"]=paste("(",round(x$ci.lb,2),",", round(x$ci.ub,2),")")
  res.table["# of cultures",1]=x$p
  res.table[c("K","I2"),1]=c(x$k, round(x$I2,2))
  res.table["QM (df,p)",1]=paste(round(x$QM,2),"(",x$m,",",round(x$QMp,2),")")
  res.table["R2",1]=round(x$R2/100,2)
  res.table["QE (df,p)",1]=paste(round(x$QE,2),"(",x$k-x$m-1,",",round(x$QEp,2),")")
  ## 95%CI of tau2 formulas see Borenstein et al. (2009) chapter16
  if(x$QE>(x$k-x$m)){
    B=0.5*((log(x$QE)-log(x$k-x$m-1))/(sqrt(2*x$QE)-sqrt(2*(x$k-x$m-1)-1))) #formula (16.14)
  } else if(x$QE<=(x$k-x$m)){
    B=sqrt(1/(2*(x$k-x$m-1-1)*(1-(1/(3*(x$k-x$m-1-1)^2))))) #formula(16.15)
  }
  C=(x$QE-x$k-x$m-1)/x$tau2 #formula(16.5)
  L=exp(0.5*log(x$QE/(x$k-x$m-1))-1.96*B) #formula(16.16)
  tau2.CIL=(x$k-x$m-1)*(L^2-1)/C #formula(16.18)
  tau2.CIL=round(tau2.CIL,2)
  U=exp(0.5*log(x$QE/(x$k-x$m-1))+1.96*B) #formula(16.17)
  tau2.CIU=(x$k-x$m-1)*(U^2-1)/C #formula(16.19)
  tau2.CIU=round(tau2.CIU,2)
  res.table["tau2 (LL,UL)",1]=paste(round(x$tau2,2),"(",tau2.CIL,",",tau2.CIU,")")
  return(res.table)
}

