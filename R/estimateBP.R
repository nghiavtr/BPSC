#' Estimate parameters of a beta-Poisson model for a data vector
#'
#' @param x A vector of input data points
#' @param para.num Mode of beta-Poisson model: 3, 4 (default) or 5 parameters
#' @param tbreak.num Number of breaks for binning
#' @param break.thres A parameter setting of \code{\link{getTbreak}} function
#' @param useExt A parameter setting of \code{\link{getTbreak}} function that allows to extend the last bin to infinity or not
#' @param param0 Initial parameters for the model. If it is not set, default initial parameters will be obtained from \code{\link{getInitParam}} function
#' @return An optimal model to the input data including optimal parameters (par), X2 test results (X2 and PVAL), etc..
#' @export
#' @examples
#' set.seed(2015)
#' #Create a simulated gene expression by randomly generating 100 data points 
#' #from a beta-poisson model
#' alp=0.6;bet=1.5;lam1=20;lam2=0.05
#' par0=c(alp,bet,lam1,lam2)
#' bp.vec=rBP(100,par0)
#' #Estimate parameters of the four-parameter beta-Poisson model from the data points
#' res=estimateBP(bp.vec,para.num=4)
#' #Print the goodness-of-fit of the model and the optimal parameters
#' res$X2
#' res$par
estimateBP<-function(x,para.num=4,tbreak.num=10,break.thres=10,useExt=FALSE,param0=NULL){
  tbreak=getTbreak(x,tbreak.num=tbreak.num,break.thres=break.thres,useExt=useExt)
  y = table(cut(x, tbreak))
  ### if bins are too sparse, get bins from expressed values
  expBin.ratio=0.4
  if (sum(y>0)<tbreak.num*expBin.ratio){
    x2=x[x>0];
    tbreak=getTbreak(x2,tbreak.num=tbreak.num,break.thres=break.thres,useExt=useExt)
    y = table(cut(x, tbreak))
    ## if it is still there
    if (sum(y>0)<tbreak.num*expBin.ratio){
      tbreak=getTbreak(x,tbreak.num=tbreak.num,break.thres=break.thres,useExt=useExt,useQuantile=TRUE)
      y = table(cut(x, tbreak))
    }
  }

  ### estimate initial parameters
  ff3 = function(param,y,n,tbreak) {
    if (sum(param<0)>0) return (1e+10)
    alpha=param[1]
    beta=param[2]
    lambda1=param[3]
#    lambda2=param[4]
    if (alpha > 1e+3)  return (1e+10)
    if (beta > 1e+3)  return (1e+10)  
    bp.y=NULL;
    for (i in 1:(length(tbreak)-1)) {
      freq=n*pBPi(tbreak[i],tbreak[i+1],alpha,beta,lambda1)
      bp.y=c(bp.y,freq)
    }
    bp.y=abs(bp.y)
    o.y=y
    res=1e+10
    naNum=sum(is.na(bp.y)) # to deal with tobp high beta or alpha
    if (naNum==0 & sum(bp.y==0)!=length(bp.y)){
      esp.val.E=1e-06
      keepIds=which(bp.y>=0)
      #  res=sum((o.y[keepIds]-bp.y[keepIds])^2/(bp.y[keepIds]+esp.val.E),na.rm = TRUE)  
      res=2*sum(o.y[keepIds]*log(o.y[keepIds]/(bp.y[keepIds]+esp.val.E)),na.rm = TRUE)
    }
    res
  }
  
  ff4 = function(param,y,n,tbreak) {
  if (sum(param<0)>0) return (1e+10)
  alpha=param[1]
  beta=param[2]
  lambda1=param[3]
  lambda2=param[4]
  if (alpha > 1e+3)  return (1e+10)
  if (beta > 1e+3)  return (1e+10)  
  bp.y=NULL;
  for (i in 1:(length(tbreak)-1)) {
    freq=n*pBPi(tbreak[i],tbreak[i+1],alpha,beta,lambda1,lambda2)
    bp.y=c(bp.y,freq)
  }
  bp.y=abs(bp.y)
  o.y=y
  res=1e+10
  naNum=sum(is.na(bp.y)) # to deal with tobp high beta or alpha
  if (naNum==0 & sum(bp.y==0)!=length(bp.y)){
    esp.val.E=1e-06
    keepIds=which(bp.y>=0)
  #  res=sum((o.y[keepIds]-bp.y[keepIds])^2/(bp.y[keepIds]+esp.val.E),na.rm = TRUE)  
    res=2*sum(o.y[keepIds]*log(o.y[keepIds]/(bp.y[keepIds]+esp.val.E)),na.rm = TRUE)
  }
  res
  }

  ff5 = function(param,y,n,tbreak) {
  if (sum(param<0)>0) return (1e+10)
  alpha=param[1]
  beta=param[2]
  lambda1=param[3]
  lambda2=param[4]
  prob0=param[5]  
  if (alpha > 1e+3)  return (1e+10)
  if (beta > 1e+3)  return (1e+10)
  if (prob0 > 1) return (1e+10)  
  bp.y=NULL;
  for (i in 1:(length(tbreak)-1)) {
    dens=pBPi(tbreak[i],tbreak[i+1],alpha,beta,lambda1,lambda2)
    bp.y=c(bp.y,dens)
  }
  bp.y=(n-round(prob0*n))*bp.y
  bp.y[1]=bp.y[1]+round(prob0*n)
  bp.y=abs(bp.y)
  o.y=y  
  res=1e+10
  naNum=sum(is.na(bp.y)) # to deal with tobp high beta or alpha
  if (naNum==0 & sum(bp.y==0)!=length(bp.y)){
    esp.val.E=1e-06
    keepIds=which(bp.y>=0)
    #  res=sum((o.y[keepIds]-bp.y[keepIds])^2/(bp.y[keepIds]+esp.val.E),na.rm = TRUE)  
    res=2*sum(o.y[keepIds]*log(o.y[keepIds]/(bp.y[keepIds]+esp.val.E)),na.rm = TRUE)
  }
  res
  }

  obp=list()
  try({
    if (is.null(param0)) param0=getInitParam(x,para.num)  
    if (para.num==3) obp = optim(param0,ff3,y=y,n=sum(y),tbreak=tbreak)
    if (para.num==4) obp = optim(param0,ff4,y=y,n=sum(y),tbreak=tbreak)
    if (para.num==5) obp = optim(param0,ff5,y=y,n=sum(y),tbreak=tbreak)
  }, silent=TRUE) # keep silent if errors occur  

  if (length(obp)!=0){
    obp$n=sum(y)
    obp$tbreak=tbreak
    obp$y=y
    obp$X2=obp$value
    obp$df=length(y)-1-length(param0)
    obp$PVAL=pchisq(obp$X2,obp$df,lower.tail = FALSE)  
    ### keep other information and compute chisquared test
    obp$AIC=2*length(param0) + obp$X2
    obp$param0=param0
    obp$prop0=sum(x==0)/length(x)
  } else { # if optimiztion can not obtain    
    obp$PVAL=NA
    obp$X2=NA
    obp$df=NA
    obp$AIC=NA
    obp$par=NA    
    obp$param0=NA

    obp$n=sum(y)
    obp$tbreak=tbreak
    obp$y=y
    obp$prop0=sum(x==0)/length(x)
  }

  return(obp)
}