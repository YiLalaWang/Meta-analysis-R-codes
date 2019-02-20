# Below is the function for producing correlation matrix with p-values
## This function formats the correlation results into APA format

#### Note that "x" here refers to the data input to be calculated for correlation

corrtest = function(x){
  require(psych)
  x=as.data.frame(x)
  corr=corr.test(x,use="pairwise",method="pearson")
  r=round(corr$r,2)
  p=round(corr$p,2)
  n=corr$n
  corr.p=matrix(paste(r," (",n,", ",p,")"),ncol = ncol(x))
  diag(corr.p)=paste(diag(r)," ",sep="")
  rownames(corr.p)=colnames(x)
  colnames(corr.p)=paste(colnames(x),"",sep = "")
  corr.p=as.matrix(corr.p)
  corr.p[upper.tri(corr.p,diag = TRUE)]=""
  corr.p=as.data.frame(corr.p)
  return(corr.p)
}

