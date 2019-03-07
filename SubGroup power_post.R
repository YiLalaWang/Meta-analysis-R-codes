##########################################################################
## Fixed- & Mixed-model sub-group analyses (ANOVA analogues) based on Hedges & Pigott (2004) 
## This is for post hoc power estimation 
## when raw data are available
## Please contact Yi Wang (wangyilalala@gmail.com) for questions
###########################################################################

SubGrouppower.post=function(x,group){
  rawdata = as.data.frame(x)
  require(dplyr)
  rawdata=rawdata%>%na.omit(rawdata[group])
  require(metafor)
  rawdata$w=1/rawdata$vi
  k=nrow(rawdata) 
  colnames(rawdata)[colnames(rawdata)==group]="group"
  res.mod=rma(es,vi,data=rawdata,level=95,mods=~group-1)
  r.g=res.mod$b
  p=nrow(r.g)
  df_p=p-1
  alpha=0.05
  
  names=gsub("group","",row.names(r.g))
  rownames(r.g)=names
  res.lumda=data.frame(matrix(nrow = p,ncol=6))
  colnames(res.lumda)=c("vi","rg","rg_rm.diff","QE","m","vi_me")
  row.names(res.lumda)=names
  for (g in names){
    name=g
    res.lumda[name,"vi"]=1/sum(subset(rawdata,group==g)["w"])
    res.lumda[name,"rg"]=r.g[g,1]
  }
  est.rm=sum(res.lumda$rg/res.lumda$vi)/sum(1/res.lumda$vi)
  for (g in names){
    name=g
    res.lumda[name,"rg_rm.diff"]=res.lumda[g,"rg"]-est.rm
    res.lumda[name,"QE"]=rma(yi,vi,data=subset(rawdata,group==g))$QE
  }
  QE_all=sum(res.lumda$QE)
  for (g in names){
    name=g
    res.lumda[name,"a"]=sum(subset(rawdata,group==g)$w)-sum((subset(rawdata,group==g)$w)^2)/sum(subset(rawdata,group==g)$w)
    res.lumda[name,"m"]=nrow(subset(rawdata,group==g))
  }
  for (g in names){
    name=g
    rawdata[ which(rawdata$group==g),"m"]=res.lumda[g,"m"]
    rawdata[ which(rawdata$group==g),"rg"]=res.lumda[g,"rg"]
  }
  
  a_all=sum(res.lumda$a)
  QE_all=sum(res.lumda$QE)
  
  tau2=(QE_all-(k-p))/a_all
  
  rawdata$vi_me=rawdata$vi+tau2
  rawdata$w_me=1/rawdata$vi_me
  for (g in names) {
    name=g
    res.lumda[name,"vi_me"]=1/sum(subset(rawdata,group==g)["w_me"])
  }
  
  ############# Power of omnibus test #############
  c.alpha.omni=qchisq(1-alpha,df_p,ncp=0,lower.tail=T, log.p=F)
  
  ## fixed-effects
  ratio.R=var(res.lumda$rg)/mean(rawdata$vi)
  lumda.FE.omni=sum(res.lumda$rg_rm.diff^2/res.lumda$vi)
  power.FE.omni=1-pchisq(c.alpha.omni, df=df_p, ncp=lumda.FE.omni,lower.tail=T,log.p=F)
  
  ## mixed effects
  lumda.ME.omni=sum(res.lumda$rg_rm.diff^2/res.lumda$vi_me)
  power.ME.omni=1-pchisq(c.alpha.omni, df=df_p, ncp=lumda.ME.omni, lower.tail=T,log.p=F)
  
  ############# power of group contrasts tests ############
  require(Surrogate)
  c=RandVec(a=-1,b=1,s=0,n=p,m=1,Seed=1)$RandVecOutput # generate a random vector with row=p & sum=0
  res.lumda=cbind(res.lumda,c)
  names(res.lumda)[ncol(res.lumda)]="c"
  gamma=abs(sum(res.lumda$c*res.lumda$rg))
  v.contrast.FE=sum((res.lumda$c^2)*res.lumda$vi)
  v.contrast.ME=sum((res.lumda$c^2)*res.lumda$vi_me)

  #### power of two-tail unplanned (post hoc) comparison
  ## using Bonferroni procedures
  l.posthoc=2
  alpha.contrast.Bonfer=alpha/(l.posthoc*2)
  c.contrast.Bonfer=qnorm(1-alpha.contrast.Bonfer,mean=0,sd=1,lower.tail=T, log.p=F)
  # Fixed-effects model
  phi.contrast.FE_Bonfer=pnorm(c.contrast.Bonfer-gamma/sqrt(v.contrast.FE),mean=0,sd=1,lower.tail=T, log.p=F)
  phi.contrast.FE.neg_Bonfer=pnorm(-c.contrast.Bonfer-gamma/sqrt(v.contrast.FE),mean=0,sd=1,lower.tail=T, log.p=F)
  power.contrast.FE_Bonfer=1-phi.contrast.FE_Bonfer+phi.contrast.FE.neg_Bonfer
  # mixed-effects model
  phi.contrast.ME_Bonfer=pnorm(c.contrast.Bonfer-gamma/sqrt(v.contrast.ME),mean=0,sd=1,lower.tail=T, log.p=F)
  phi.contrast.ME.neg_Bonfer=pnorm(-c.contrast.Bonfer-gamma/sqrt(v.contrast.ME),mean=0,sd=1,lower.tail=T, log.p=F)
  power.contrast.ME_Bonfer=1-phi.contrast.ME_Bonfer+phi.contrast.ME.neg_Bonfer
  
  ## Using Scheffe procedures
  c.contrast.Scheffe=sqrt(qchisq(1-alpha,df_p,ncp=0,lower.tail=T, log.p=F))
  # fixed effects model
  phi.contrast.FE_Scheffe=pnorm(c.contrast.Scheffe-gamma/sqrt(v.contrast.FE),mean=0,sd=1,lower.tail=T, log.p=F)
  phi.contrast.FE.neg_Scheffe=pnorm(-c.contrast.Scheffe-gamma/sqrt(v.contrast.FE),mean=0,sd=1,lower.tail=T, log.p=F)
  power.contrast.FE_Scheffe=1-phi.contrast.FE_Scheffe+phi.contrast.FE.neg_Scheffe
  # mixed effects model
  phi.contrast.ME_Scheffe=pnorm(c.contrast.Scheffe-gamma/sqrt(v.contrast.ME),mean=0,sd=1,lower.tail=T, log.p=F)
  phi.contrast.ME.neg_Scheffe=pnorm(-c.contrast.Scheffe-gamma/sqrt(v.contrast.ME),mean=0,sd=1,lower.tail=T, log.p=F)
  power.contrast.ME_Scheffe=1-phi.contrast.ME_Scheffe+phi.contrast.ME.neg_Scheffe
  
  ############# Power of within-group heterogeneity (QE) #############
  c.fit=qchisq(1-alpha,k-p,ncp=0,lower.tail=T, log.p=F)
  lumda.fit=sum((rawdata$es-rawdata$rg)^2*rawdata$w,na.rm=TRUE)
  phi.fit=pchisq(c.fit,k-p,ncp=lumda.fit,lower.tail=T,log.p=F)
  power.fit=1-phi.fit
  
  ############# Power of the Between-Studies variance component (tau2) ########################
  Miu.Q=a_all*tau2+k-p
  for (g in names) {
    name=g
    res.lumda[name,"S"]=sum(subset(rawdata,group==g)$w)
    res.lumda[name,"S2"]=sum((subset(rawdata,group==g)$w)^2)
    res.lumda[name,"S3"]=sum((subset(rawdata,group==g)$w)^3)
    res.lumda[name,"b"]=sum(res.lumda[g,"m"]-1+2*(res.lumda[g,"S"]-res.lumda[g,"S2"]/res.lumda[g,"S"])*tau2+
                              (res.lumda[g,"S2"]-2*res.lumda[g,"S3"]/res.lumda[g,"S"]+
                                 res.lumda[g,"S2"]^2/(res.lumda[g,"S"]^2))*(tau2^2))
  }
  
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
  power.contrast.FE_Bonfer=round(power.contrast.FE_Bonfer,digits=2)
  power.contrast.FE_Scheffe=round(power.contrast.FE_Scheffe,digits=2)
  
  power.me.omni=round(power.ME.omni,digits=2)
  power.me.tau2=round(power.tau2,digits=2)
  power.contrast.ME_Bonfer=round(power.contrast.ME_Bonfer,digits=2)
  power.contrast.ME_Scheffe=round(power.contrast.ME_Scheffe,digits=2)
  
  ################## Organizing all results #############################
  SubGrouppower.table.post=data.frame(matrix(nrow = 4,ncol = 2))
  colnames(SubGrouppower.table.post)=c("fixed effects","mixed effects")
  rownames(SubGrouppower.table.post)=c("omnibus test (QM)", "Model fit (QE)",
                                       "group contrast (Bonferroni two tail)",
                                       "group contrast (Scheffe two tail)")
  SubGrouppower.table.post["fixed effects"]=c(power.fe.omni, power.fe.fit,
                                              power.contrast.FE_Bonfer, power.contrast.FE_Scheffe)
  SubGrouppower.table.post["mixed effects"]=c(power.me.omni, power.me.tau2,
                                              power.contrast.ME_Bonfer, power.contrast.ME_Scheffe)
  return(SubGrouppower.table.post)
}

############### Instructions for specifying values in the SubGrouppower.post function ##########################

# use x to specify the name for the raw data file
### Raw data must be organized with the following columns: 
### vi: variance of effect sizes
### es: effect sizes 
### groups for comparison

# use "group" to specify the grouping variable (including the quotes in your code, e.g., group="group")

# For an example using Table 1 data in Hedges & Pigott (2004), see below: 
rawdata=read.table(file="Table 1.csv",header=TRUE,sep=",") 
rawdata$yi=(1-3/(4*(rawdata$n-2)-1))*rawdata$es 
rawdata$vi=2/(rawdata$n/2)+rawdata$yi^2/(2*rawdata$n)
# yi & vi converted using standard conversion for standardized mean difference
SubGrouppower.post(rawdata,group="SES")
