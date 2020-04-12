# This is the function that conducts publication bias tests based on:
## (1) Egger's test:
### column "Egger's test" produces the estimate and 95% CI of the test;

## (2) PET-PEESE method: 
### column "PET" produces the PET estimate and its 95% CI, 
### column "PEESE" produces the PEESE estimate to compare with rm; 

## (3) cumulative meta-analysis based on precision: 
### column "cum.dir" produces the direction of the potential drift of the forest plot, 
#### "Neg" = negative drift as precision decrease, "Pos" = positive drift,
### column "cum.drift(percent) produces the magnitude and % of the drift of rm between the 
#### top 20% most precise (top 20% largest N) studies and the top 20% least precise 
#### (top 20% smallest N) studies;

## (4) fixed-random trim and fill acompanied with contour-enhanced funnel plot: 
### column "FE k fill (side)" produces the number of filled studies and at which side they are filled in the funnel plot
### column "FE k fill non-sig" produces the number of filled studies that are in the non-significant area of the contour-enhanced funnel plot
### column "FE rm fill (95%CI)" produces the rm and the 95% CI after studies are filled in
### column "FE diff.fill (diff.perc)" produces the difference between the original and filled rm and the % of difference from rm.

# rawdata should at least include r(effect sizes) and N(sample sizes) for all studies
# note that N is in capital letter
## sample/study info (e.g., citation) is strongly recommended

PubBias=function(rawdata){
  require(metafor)
  rawdata=data.frame(rawdata)
  rawdata.rescale=escalc(measure="ZCOR",ri=r,ni=N,data=rawdata)
  rawdata.rescale$sei=sqrt(rawdata.rescale$vi)
  
  ###### original meta-analysis using RMEL estimator for heterogeneity #######
  x=rma(measure="ZCOR",yi,vi,data=rawdata.rescale,level=95,slab=rawdata.rescale$SampleID)
  ci=predict(x,level=95,transf=transf.ztor) # 95% CI
  rm=round(ci$pred,2)
  ci.l=round(ci$ci.lb,2)
  ci.u=round(ci$ci.ub,2)
  ci=paste(ci.l,",",ci.u)
  k=as.numeric(nrow(rawdata))
  N=sum(rawdata$N)
  se=(ci.u-rm)/1.96
  
  ###### Egger's regression test ######
  egger.x=regtest(x,model="lm",predictor="sei")
  egger.a=round(coef(egger.x$fit)["Xsei","Estimate"],2) ## slope of regression (alpha)
  egger.a.se=coef(egger.x$fit)["Xsei","Std. Error"]
  egger.a.LL=round(egger.a-egger.a.se*1.96,2) ## 95%CI to test if alpha=0
  egger.a.UL=round(egger.a+egger.a.se*1.96,2)
  
  ##########PET-PEESE ANALYSIS ##########
  # based on equations in Stanley & Doucouliagos (2014) Res.Syn.Met.
  PET.gamma = round(coef(egger.x$fit)["Xintrcpt","Estimate"],2) # gamma1 (formula 6)
  PET.gamma.se=coef(egger.x$fit)["Xintrcpt","Std. Error"]
  PET.gamma.LL=round(PET.gamma-PET.gamma.se*1.96,2) ## 95%CI to test if gamma1=0
  PET.gamma.UL=round(PET.gamma+PET.gamma.se*1.96,2)
  PET.gamma.p=coef(egger.x$fit)["Xintrcpt","Pr(>|t|)"]
  
  reg.PEESE = regtest(x,model="lm",predictor="vi") # formula 8
  PEESE.beta = round(coef(reg.PEESE$fit)["Xintrcpt","Estimate"],2) # beta1 (formula 8)
  
  ###### (fixed-random) trim & fill + contour-enhanced funnel plot######
  x.fe=rma(measure="ZCOR",yi,vi,data=rawdata.rescale,level=95,method="FE")
  tnf.fe.x=trimfill(x.fe,estimator="L0")
  fill.fe.x=tnf.fe.x$fill
  res.tnf.fe.x=rma(tnf.fe.x$yi,tnf.fe.x$vi,level=95)
  res.tnf.fe.x$fill=fill.fe.x
  class(res.tnf.fe.x)=c("rma.uni.trimfill","rma")
  k.fe.fill = tnf.fe.x$k0
  side.fe.fill = tnf.fe.x$side
  rm.fe.fill = round(predict(res.tnf.fe.x,level=95,transf=transf.ztor)$pred, 2)
  rm.fe.LL.fill = round(predict(res.tnf.fe.x,level=95,transf=transf.ztor)$ci.lb, 2)
  rm.fe.UL.fill = round(predict(res.tnf.fe.x,level=95,transf=transf.ztor)$ci.ub, 2)
  ## comparing diff in rm before & after fill using percentage
  diff.fe.fill = round(rm.fe.fill-rm,2) 
  perc.diff.fe.fill = round((abs(diff.fe.fill)/abs(rm))*100,2)
  
  dat.fil=as.data.frame(cbind(res.tnf.fe.x$yi.f,res.tnf.fe.x$vi.f,res.tnf.fe.x$fill))
  colnames(dat.fil)=c("yi","vi","fill")
  dat.fil.t=subset(dat.fil,fill==1)
  dat.fil.t$vi.95=1.96*sqrt(dat.fil.t$vi)
  dat.fil.t$yi.abs=abs(dat.fil.t$yi)
  dat.fil.t$non.sig=ifelse(dat.fil.t$yi.abs>dat.fil.t$vi.95,0,1)
  ### if abs(yi) < 1.96*sqrt(vi) then it's within the non-significant area
  k.fil.fe.ng=sum(dat.fil.t$non.sig)
  
  ## for visual inspection only 
  # funnel(res.tnf.fe.x,level=c(95, 99), shade=c("white", "gray55"), refline=0,
  #         atransf=transf.ztor, yaxis="vinv",  # vinv is used so y-axis looks more like actual N
  #         legend=TRUE,back="grey75",hlines="grey75",
  #         digits=c(2,0),font=6, font.lab=6) #(font=6 turns fonts into Times New Roman)
  
  ########## cumulative meta-analysis by precision ################
  #Sort by sei
  rawdata.rescale = rawdata.rescale[order(rawdata.rescale$sei),] 
  #Calculate cumulative N values
  rawdata.rescale$cumsumN = round(cumsum(rawdata.rescale$N), digits = 2)
  rawdata.rescale$cumsumN = paste0(rawdata.rescale$N," (", rawdata.rescale$cumsumN,")")
  
  res.cma = cumul(x, order=order(rawdata.rescale$sei))
  #forest(res.cma, transf=transf.ztor) # for visual inspection only
  dat.cum=as.data.frame(cbind(res.cma$estimate,sort(rawdata.rescale$N,decreasing=TRUE)))
  colnames(dat.cum)=c("res.cul","N")
  trend=lm(formula=res.cul~N,data=dat.cum) # a positive drift = the larger N, the larger r
  dir.drift=if((as.list(trend$coefficients)$N<0)){"Neg"} else {"Pos"}
  
  ### compute & compare the mean r with the 20% most/least precise studies 
  select = round(nrow(rawdata.rescale)*.20,0)
  dat.first=rawdata.rescale[seq(from = 1, to = select),]
  dat.bottom=rawdata.rescale[seq(from=nrow(rawdata.rescale)-select+1,to=nrow(rawdata.rescale)),]
  x.first=rma(measure="ZCOR",yi,vi,data=dat.first,level=95)
  ci.first=predict(x.first,level=95,transf=transf.ztor) 
  rm.first=ci.first$pred
  
  x.bottom=rma(measure="ZCOR",yi,vi,data=dat.bottom,level=95)
  ci.bottom=predict(x.bottom,level=95,transf=transf.ztor) 
  rm.bottom=ci.bottom$pred
  
  diff.drift = round(rm.first - rm.bottom,2) 
  perc.diff.drift = round((abs(diff.drift)/abs(rm.bottom))*100,2)
  
  ####### organizing results into table format #######
  res.table=data.frame(matrix(nrow=11, ncol=1))
  colnames(res.table)="results"
  rownames(res.table)=c("K","rm (95%CI)","Egger's test", "PET", "PEESE",
                        "cum.dir", "cum.drift(percent)",  
                        "FE k fill (side)","FE k fill non-sig",
                        "FE rm fill (95%CI)", "FE diff.fill (diff.perc)")
  res.table["K",1]=k
  res.table["rm (95%CI)",1] = paste(rm,"(",ci,")")
  res.table["Egger's test",1]=paste(egger.a,"(",egger.a.LL,",",egger.a.UL,")")
  res.table["PET",1]=paste(PET.gamma,"(",PET.gamma.LL,",",PET.gamma.UL,")")
  if(PET.gamma.p < .05){
    res.table["PEESE",1]=PEESE.beta
  } else {is.null(res.table["PEESE",1])}
  
  res.table["cum.dir",1]=dir.drift
  res.table["cum.drift(percent)",1]=paste(diff.drift,"(",perc.diff.drift,")")
  
  res.table["FE k fill (side)",1]=if(k.fe.fill==0){0} else {paste(k.fe.fill, "(", side.fe.fill,")")}
  if(k.fe.fill == 0) {
    is.null(res.table["FE rm fill (95%CI)",1])
    is.null(res.table["FE k fill non-sig",1])
    is.null(res.table["FE diff.fill (diff.perc)",1])
  } else {
    res.table["FE rm fill (95%CI)",1] = paste(rm.fe.fill,"(", rm.fe.LL.fill, ",", rm.fe.UL.fill, ")")
    res.table["FE k fill non-sig",1]=k.fil.fe.ng
    res.table["FE diff.fill (diff.perc)",1] = paste(diff.fe.fill,"(",perc.diff.fe.fill,")")
  }
  return(res.table)
}
