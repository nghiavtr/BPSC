#' Get breaks for binning of data points
#'
#' @param x A vector of input data points
#' @param tbreak.num Number of breaks
#' @param break.thres A threshold to decide whether use breaks from logs cale or not. We decide to use breaks from log scale if 75 percent of input data less than this threshold
#' @param useExt An option to extend the last bin to infinity or not
#' @param useQuantile An option to bin the data points by equal-size quantile values
#' @return A vector of break values
#' @export
#' @examples
#' set.seed(2015)
#' bp.vec=rBP(100,alp=0.6,bet=1.5,lam1=20,lam2=0.05)
#' tbreak=getTbreak(bp.vec,tbreak.num=10)
getTbreak=function(x,tbreak.num=10,break.thres=10,useExt=FALSE,useQuantile=FALSE){
  if (useExt) x=x[x <= quantile(x[x>0],prob=0.95)]

  if (useQuantile){
    x2=x[x>0];
    qtvec=seq(0,1,length = tbreak.num-1)
    tbreak=quantile(x2,prob=qtvec)
  }else{    
    if (quantile(x,prob=0.75)>break.thres){
      rr = range(x, na.rm = TRUE)
      tbreak = seq(rr[1], rr[2], length = tbreak.num)
      tbreak=sort(unique(tbreak))
    }else{ 
      rr = range(log2(x+1), na.rm = TRUE)
      tbreak = seq(rr[1], rr[2], length = tbreak.num)
      tbreak=2^tbreak-1
      tbreak=sort(unique(tbreak))
    }
  }
  tbreak[1] = tbreak[1]-1.0e-06
  tbreak[length(tbreak)] = tbreak[length(tbreak)]+1.0e-06
  if (useExt) tbreak[length(tbreak)]=Inf
  tbreak=unique(tbreak);
  if (anyNA(tbreak) || length(tbreak)==1) tbreak=tbreak.num
  return (tbreak)
}
