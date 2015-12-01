#' Generate random values from a beta-Poisson model
#'
#' @param n The number of random data points
#' @param alp A non-negative value, the parameter of Beta function (alpha)
#' @param bet A non-negative value, the parameter of Beta function (beta)
#' @param lam1 A non-negative value, the parameter for scaling
#' @param lam2 A non-negative value, the parameter for smoothing, used for BP4 or BP5
#' @param prob0 A non-negative value in [0;1), the parameter for the portion of fixed density of zero in the distribution, used for BP5
#' @return A vector of the generated data points
#' @export
#' @importFrom stats ppois pnorm rpois rbeta dpois dnorm
#' @importFrom statmod gauss.quad
#' @examples
#' set.seed(2015)
#' #BP3
#' bp.vec=rBP(100,0.2,1.5,200)
#' hist(bp.vec, prob=TRUE)
#' #BP4
#' bp.vec=rBP(100,0.6,1.5,20,0.05)
#' hist(bp.vec, prob=TRUE)
#' #BP5
#' bp.vec=rBP(100,0.2,1.5,200,0.1,0.05)
#' hist(bp.vec, prob=TRUE)
rBP = function(n,alp,bet,lam1=1,lam2=1,prob0=0) {
	if (missing(bet)){
  	par=alp;
  	alp=par[1]; bet=par[2]; lam1=par[3]; lam2=1; if (length(par)> 3) lam2=par[4]; prob0=0;if (length(par)> 4) prob0=par[5]
    }

    bp.count=lam2*rpois(n-round(prob0*n),lambda=lam1*rbeta(n,shape1=alp,shape2=bet))
    bp.count=c(bp.count,rep(0,round(prob0*n)))# add zeros
    return(bp.count)
}