######## R code for composite calculation ###########
## This function produces a new dataset where multiple effect sizes and multiple reliabilities 
      # in the same sample are combined into the same composite effect size and composite reliability.

## The calculation requires preinstalling the psychmeta package. 
## "rawdata_comp" is the original data files for composite that should include at least the following: 
# (1) columns for sample identifiers, entitled "author", "year", and "sampleID";
# (1) a column for all effect sizes (r between Xs and Ys) that go into the composite, entitled "r";
# (2) a column for sample size n, entitled "n";
    # sample size has to be the same for all effect sizes in the same composite. 
# (3) two columns for reliabilities of all Xs and Ys, entitled as "rxx" and "ryy".
# (4) columns "x1" to "x7", which contain the correlations among all Xs;
# (5) columns "y1" to "y7", which contain the correlations among all Ys. 

## The new dataset contains three columns which helps trouble-shooting: 
# (1) "k_vars_x": Number of Xs that go into the composite;
# (2) "k_vars_y": Number of Yx that go into the composite; 
# (3) "k_xy": Number of X-Y correlations that go into the composite (k_xy = k_vars_x * k_vars_y).
# All these columns should produce integers. 
# If not, then there may be miscalculation or mispecification in rawdata_comp. 

composite_w_rel=function(rawdata_comp){
  require(psychmeta)
  rawdata_comp$ID=apply(as.data.frame(rawdata_comp[c("author","year","sampleID")]),
                        1,paste,collapse=".")
  ## subset without duplicated samples
  rawdata=rawdata_comp[!(duplicated(rawdata_comp$ID)|duplicated(rawdata_comp$ID,fromLast = TRUE)),]
  ## subset with duplicated samples
  rawdata_comp=rawdata_comp[duplicated(rawdata_comp$ID)|duplicated(rawdata_comp$ID,fromLast = TRUE),]
  comp=unique(rawdata_comp$ID)
  ##### Calculating composite r for each duplicated samples #####
  for (id in comp) {
    rxy=rawdata_comp$r[rawdata_comp$ID==id]
    intercor_x=subset(rawdata_comp,ID==id)[c("x1","x2","x3","x4","x5","x6","x7")]
    intercor_x=as.matrix(intercor_x)
    intercor_y=subset(rawdata_comp,ID==id)[c("y1","y2","y3","y4","y5","y6","y7")]
    intercor_y=as.matrix(intercor_y)
    rel_x=subset(rawdata_comp,ID==id)[c("rxx")]
    rel_x=as.matrix(rel_x)
    rel_y=subset(rawdata_comp,ID==id)[c("ryy")]
    rel_y=as.matrix(rel_y)
    ##### composite function based on Ghiselli et al. (1981) ####
    composite_r=function(rxy,intercor_x,intercor_y,rel_x,rel_y){
      k_xy=sum(!is.na(rxy))
      k_x=sum(!is.na(intercor_x))
      k_y=sum(!is.na(intercor_y))
      # based on an nonlinear simultaneous equation system: 
      # (1)k_vars_x=k_xy/k_vars_y, (2)k_y=(k_vars_y*(k_vars_y-1))/2, (3)k_x=(k_vars_x*(k_vars_x-1))/2 
      k_vars_y = (k_xy^2-2*2*(k_y*k_x))/(2*k_x+k_xy)
      k_vars_x = k_xy/k_vars_y
      
      mean_rxy=mean(rxy)
      mean_intercor_x=if(k_vars_x==1){NULL} else {mean(intercor_x,na.rm=TRUE)}
      mean_intercor_y=if(k_vars_y==1){NULL} else {mean(intercor_y,na.rm=TRUE)}
      k_vars_x = if(k_vars_x==1){NULL} else {k_vars_x}
      k_vars_y = if(k_vars_y==1){NULL} else {k_vars_y}
      r_comp=composite_r_scalar(mean_rxy,k_vars_x,mean_intercor_x,k_vars_y,mean_intercor_y)
      r_comp=round(r_comp,2)
      
      if(length(rel_x)==0){rel_x_comp=NA} else if(is.null(k_vars_x)){rel_x_comp=mean(rel_x,na.rm=TRUE)
      } else {
        mean_rel_x=mean(rel_x,na.rm=TRUE)
        rel_x_comp=composite_rel_scalar(mean_rel_x,mean_intercor_x,k_vars_x)}
      rel_x_comp=if(is.na(rel_x_comp)){NA} else {round(rel_x_comp,2)}
      
      if(length(rel_y)==0){rel_y_comp=NA} else if(is.null(k_vars_y)){rel_y_comp=mean(rel_y,na.rm=TRUE)
      } else {
        mean_rel_y=mean(rel_y,na.rm=TRUE)
        rel_y_comp=composite_rel_scalar(mean_rel_y,mean_intercor_y,k_vars_y)}
      rel_y_comp=if(is.na(rel_x_comp)){NA} else {round(rel_y_comp,2)}
      
      k_vars_x = if(is.null(k_vars_x)){1} else {k_vars_x}
      k_vars_y = if(is.null(k_vars_y)){1} else {k_vars_y}
      list(r_comp=r_comp,k_xy=k_xy,rel_x_comp=rel_x_comp,rel_y_comp=rel_y_comp,
           k_vars_x=k_vars_x,k_vars_y=k_vars_y)
    }
    r_composite=composite_r(rxy,intercor_x,intercor_y,rel_x,rel_y)
    rawdata_comp$r_comp[rawdata_comp$ID==id]=r_composite$r_comp
    rawdata_comp$rel_x_comp[rawdata_comp$ID==id]=r_composite$rel_x_comp
    rawdata_comp$rel_y_comp[rawdata_comp$ID==id]=r_composite$rel_y_comp
    rawdata_comp$k_xy[rawdata_comp$ID==id]=r_composite$k_xy
    rawdata_comp$k_vars_x[rawdata_comp$ID==id]=r_composite$k_vars_x
    rawdata_comp$k_vars_y[rawdata_comp$ID==id]=r_composite$k_vars_y
  }
  #### Combining both datasets
  rawdata_comp$r=rawdata_comp$r_comp
  rawdata_comp$rxx=rawdata_comp$rel_x_comp
  rawdata_comp$ryy=rawdata_comp$rel_y_comp
  rawdata_comp$r_comp<-rawdata_comp[c("x1","x2","x3","x4","x5","x6","x7")]<-NULL
  rawdata_comp$rel_x_comp<-rawdata_comp$rel_y_comp<-NULL
  rawdata_comp[c("y1","y2","y3","y4","y5","y6","y7")]<-NULL
  rawdata_comp=rawdata_comp[!duplicated(rawdata_comp$ID),]
  rawdata[c("x1","x2","x3","x4","x5","x6","x7")]<-rawdata[c("y1","y2","y3","y4","y5","y6","y7")]<-NULL
  rawdata$k_vars_x=rawdata$k_vars_y=rawdata$k_xy=1
  rawdata=rbind(rawdata,rawdata_comp)
  rawdata
}