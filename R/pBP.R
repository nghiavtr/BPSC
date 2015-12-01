#' Compute the probability distribution function (pdf) of a single quantile from a four-parameter beta-Poisson model
#'
#' @param x A non-negative quantile
#' @param alp A non-negative value, the parameter of Beta function (alpha)
#' @param bet A non-negative value, the parameter of Beta function (beta)
#' @param lam1 A non-negative value, the parameter for scaling
#' @param lam2 A non-negative value, the parameter for smoothing, used for BP4 or BP5
#' @return pdf of the input value
#' @export
#' @importFrom stats ppois pnorm rpois rbeta dpois dnorm
#' @importFrom statmod gauss.quad
#' @examples
#' x=0.3
#' pBP(x,alp=0.6,bet=1.5,lam1=20,lam2=0.05)
pBP = function(x,alp,bet,lam1=1,lam2=1) {
  if (missing(bet)){
    par=alp;
    alp=par[1]; bet=par[2]; lam1=par[3]; lam2=1; if (length(par)> 3) lam2=par[4];
  }

  x=x/lam2
  ff = function(k,m){
    if (max(m) < 100000) res= dpois(k,m) else res = dnorm(k,m,sqrt(m))
    return(res)
  }
  w = gauss.quad(10,'jacobi', alpha=bet-1, beta=alp-1)
  gs = sum(w$weight*ff(x,m=lam1*(1+w$node)/2))
  prob = 1/beta(alp,bet)*2^(-alp-bet+1)*gs   # get P(X=x)
  prob=prob
  return(prob)
}