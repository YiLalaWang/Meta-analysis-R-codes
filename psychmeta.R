######################################################################
########### R code for psychometric meta-analysis main effect ########

## R code that contains the function for: 
# (1) doing meta-analysis corrected for measurement errors in X & Y and weighted by sample size; 
        # (a) Individual corrections are used; 
        # (b) Missing reliabilities are imputed with the sample-size weighted mean of the full distribution of artifacts 
              # Although this imputation method is not recommended by default, 
              # it saves lots of computation compared to 
              # other methods such as those based on bootstrapping or simulation. 
        # (c) 95% CI and 90% CV are produced for the meta-analysis. 
# (2) organizing all meta-analytic results into table format.

## The calculation requires preinstalling the psychmeta and metafor packages. 
## "rawdata" is the original data files for meta-analysis that should include at least the following: 
# (1) a column for effect size r, entitled "r";
# (2) a column for sample size n, entitled "n";
# (3) two columns for reliabilities of the correlated constructs, entitled as "rxx" and "ryy".

## At least 3 effect sizes (i.e., K>=3) is required to produce results. 
## Function would return null if K<3. 

psychmeta=function(rawdata){
  rawdata=data.frame(rawdata)
  require(psychmeta)
  if (nrow(rawdata)<3) {return(NULL)
  } else {
    res.table=ma_r(data=rawdata,ma_method="ic",rxyi=r,n=n,rxx=rxx, ryy=ryy,
                   wt_type="sample_size", correct_rxx=TRUE, correct_ryy=TRUE,
                   control=control_psychmeta(error_type="sample",conf_level=.95, cred_level=.9,
                                             conf_method="norm",cred_method="norm", var_unbiased=TRUE,
                                             pairwise_ads=FALSE, residual_ads=FALSE,impute_artifacts=TRUE,
                                             impute_method="wt_mean_full"))$meta_tables$`analysis_id: 1`$individual_correction$true_score
    res.table=round(res.table,2)
    res.table[,c(1:2)]=round(res.table[,c(1:2)],0)
    
    require(metafor)
    rawdata.rescale=escalc(measure="ZCOR",ri=r,ni=n,data=rawdata)
    x=rma(measure="ZCOR",yi,vi,data=rawdata.rescale,level=95)
    res.table[1,"I2"]=round(x$I2,2)
    res.table[1,"Q(df,p)"]=paste(round(x$QE,2),"(",x$k-x$m,",",round(x$QEp,2),")")
    
    return(res.table)}
}