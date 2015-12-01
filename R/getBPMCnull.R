#' Generate a Monte-Carlo null distribution for a beta-Poisson model
#'
#' @param param A vector of parameters of the beta-Poisson model
#' @param n The number of data points for each simulation
#' @param tbreak A vector of break points
#' @param tbreak.num Number of breaks
#' @param sim.num A number of simulations
#' @param ran.num The number of data points generated from the beta-Poisson model to approximate the theoretical model
#' @param E.esp An small value added to expected value when computing X2
#' @param useDebug A parameter just used for debug and checking, so useDebug=FALSE by default
#' @return A list containning the Monte-Carlo null distribution and other values
#' @export
#' @examples
#' set.seed(2015)
#' par=c(0.6,1.5,20,0.05)
#' MCnull.res=getBPMCnull(par,n=100,sim.num=100)
#' hist(MCnull.res$X2,xlab="goodness-of-fit statistic", 
#' main="Monte-Carlo null distribution",breaks=20, prob=TRUE)
getBPMCnull = function(param,n=200,tbreak=NULL,tbreak.num=20,sim.num=1000,ran.num=100000,E.esp=0.0,useDebug=FALSE){
  alpha=param[1]; beta=param[2]; lambda1=param[3]; lambda2=1; if (length(param)> 3) lambda2=param[4]; prob0=0;if (length(param)> 4) prob0=param[5]
  # theoretical beta poisson distribution
  bp.count=rBP(n=ran.num,alp=alpha,bet=beta,lam1=lambda1,lam2=lambda2,prob0=prob0)
  X2=PVAL=df=NULL;
  for (i in 1:sim.num){
    ### generate n data points
    bp.sub.count=rBP(n=n,alp=alpha,bet=beta,lam1=lambda1,lam2=lambda2,prob0=prob0)
    rr = range(c(bp.sub.count,bp.count), na.rm = TRUE)
    if (is.null(tbreak)){
      tbreak = seq(rr[1], rr[2], length = tbreak.num)
      tbreak[1] = tbreak[1]-1.0e-06
      tbreak[length(tbreak)] = tbreak[length(tbreak)]+1.0e-06
      tbreak=unique(tbreak) 
    }
    y = table(cut(bp.sub.count, tbreak))    
    bp.y=table(cut(bp.count, tbreak))
    bp.y=bp.y*n/ran.num
    obp=NULL
    obp$df=length(y)-1-length(param)
    if (obp$df > 0){      
      X2Val=2*sum(y*log(y/(bp.y+E.esp)),na.rm = TRUE)
      if (X2Val < 0) X2Val=0 # too strict, it should change to X2Val=Inf or X2Val=max(X2)
      obp$X2=X2Val
      if (useDebug){ 
        obp$PVAL=pchisq(obp$X2,obp$df,lower.tail = FALSE)# else obp$PVAL=NA
        PVAL=c(PVAL,obp$PVAL)
        df=c(df,obp$df)        
      }
      X2=c(X2,obp$X2)
    }
  }
  if (!useDebug) PVAL=df=NA
  return(list(PVAL=PVAL,X2=X2,df=df,param=param))
}
