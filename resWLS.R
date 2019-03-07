
# This is the code for organizing WLS regression results
## The code also includes the calculation for 95%CI of tau2
## input x should be the dataset, "moderators" should specify the moderator(s) to be tested
## As an example with rawdata as dataset and a & b as moderators: 
res.WLS(rawdata,c("a","b"))

res.WLS=function(x,moderators){
  rawdata=as.data.frame(x)
  require(metafor)
  res.mod=rma(measure="COR",ri=r,ni=N,data=rawdata,level=95,mods=as.formula(paste("~",paste(moderators,collapse="+"))))

  res.table=data.frame(matrix(nrow = res.mod$p,ncol=1))
  colnames(res.table)=c("95%CI [UL, LL]")
  rownames(res.table)=c("intercept",moderators)
  res.table["95%CI [UL, LL]"]=paste(round(res.mod$beta,2)," (",round(res.mod$ci.lb,2),",", round(res.mod$ci.ub,2),")")
  #res.table["# of cultures",1]=res.mod$p
  res.table[c("K","I2"),1]=c(res.mod$k, round(res.mod$I2,2))
  res.table["QM (df,p)",1]=paste(round(res.mod$QM,2),"(",res.mod$m,",",round(res.mod$QMp,2),")")
  res.table["R2",1]=round(res.mod$R2/100,2)
  res.table["QE (df,p)",1]=paste(round(res.mod$QE,2),"(",res.mod$k-res.mod$m-1,",",round(res.mod$QEp,2),")")
  ## 95%CI of tau2 formulas see Borenstein et al. (2009) chapter16
  if(res.mod$QE>(res.mod$k-res.mod$m)){
    B=0.5*((log(res.mod$QE)-log(res.mod$k-res.mod$m-1))/(sqrt(2*res.mod$QE)-sqrt(2*(res.mod$k-res.mod$m-1)-1))) #formula (16.14)
  } else if(res.mod$QE<=(res.mod$k-res.mod$m)){
    B=sqrt(1/(2*(res.mod$k-res.mod$m-1-1)*(1-(1/(3*(res.mod$k-res.mod$m-1-1)^2))))) #formula(16.15)
  }
  C=(res.mod$QE-res.mod$k-res.mod$m-1)/res.mod$tau2 #formula(16.5)
  L=exp(0.5*log(res.mod$QE/(res.mod$k-res.mod$m-1))-1.96*B) #formula(16.16)
  tau2.CIL=(res.mod$k-res.mod$m-1)*(L^2-1)/C #formula(16.18)
  tau2.CIL=round(tau2.CIL,2)
  U=exp(0.5*log(res.mod$QE/(res.mod$k-res.mod$m-1))+1.96*B) #formula(16.17)
  tau2.CIU=(res.mod$k-res.mod$m-1)*(U^2-1)/C #formula(16.19)
  tau2.CIU=round(tau2.CIU,2)
  res.table["T2 (95%CI[LL,UL])",1]=paste(round(res.mod$tau2,2),"(",tau2.CIL,",",tau2.CIU,")")
  return(res.table)
}

