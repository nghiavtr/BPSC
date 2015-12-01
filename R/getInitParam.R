#' Estimate initial parameters of the beta-Poisson model for a data vector
#'
#' @param x A vector of input data points
#' @param para.num Mode of beta-Poisson model: 3, 4 (default) or 5 parameters
#' @return Initial parameters for BP models of the data points
#' @export
#' @examples
#' set.seed(2015)
#' bp.vec=rBP(100,alp=0.6,bet=1.5,lam1=20,lam2=0.05)
#' init.par=getInitParam(bp.vec,para.num=4)
getInitParam<-function(x,para.num=4){  
  ### estimate initial parameters
  if (sum(x)>0){
    lambda2=1; 
    if (para.num!=3){
      #scale.thres=length(x)*0.75/10
      #if (quantile(x[x>0],prob=0.50)<=scale.thres) lambda2=0.1
      # if isoform expression is very small (<1.0)
      medianNonzero=median(x[x>0]); if (medianNonzero < 1.0) lambda2=medianNonzero*0.1
    }
    
    Ey=mean(x/lambda2); Vy=var(x/lambda2);
    lambda1=max(x)/lambda2
    if (lambda1-Ey < 1e-6) lambda1=2*Ey; #avoid zeros of alpha and beta
    alpha=(Ey*(lambda1-Ey)/(Vy-Ey)-1)*Ey/lambda1
    beta= alpha*(lambda1-Ey)/Ey
    if (alpha==0) alpha=1;if (beta==0) beta=1;
    prob0=0
    } else {
      # extreme case: all values = 0
      alpha=beta=lambda1=lambda2=0
      prob0=1
    }
  
  if (para.num==3){ 
    param0 = c(alpha,beta,lambda1)
  }
  if (para.num==4){ 
     param0 = c(alpha,beta,lambda1,lambda2)
  }
  if (para.num==5){ 
    param0 = c(alpha,beta,lambda1,lambda2,prob0)
  }
  param0=abs(param0);
  return(param0)
}