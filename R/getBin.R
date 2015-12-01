#' Get bins of data points
#'
#' @param x A vector of input data points
#' @param tbreak A vector of break points
#' @param tbreak.num Number of breaks
#' @param minVal Under bound of the binning range
#' @param maxVal Upper bound of the binning range
#' @return A vector of bin values
#' @export
#' @examples
#' set.seed(2015)
#' bp.vec=rBP(100,alp=0.6,bet=1.5,lam1=20,lam2=0.05)
#' # get 10 bins for the data
#' getBin(bp.vec,tbreak.num=10)
#' # get bins from a predefined tbreak
#' tbreak=getTbreak(bp.vec,tbreak.num=10,break.thres=10)
#' getBin(bp.vec,tbreak=tbreak)
getBin=function(x,tbreak=NULL,tbreak.num=NULL,minVal=0,maxVal=NULL){
  if (is.null(tbreak)){
    if (is.null(maxVal)) rr = range(x, na.rm = TRUE) else rr=c(minVal,maxVal)
    tbreak = seq(rr[1], rr[2], length = tbreak.num)
    tbreak[1]=minVal-1.0e-06
  }
  y = table(cut(x, tbreak))
  return (y)
}
