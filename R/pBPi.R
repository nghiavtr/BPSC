#' Compute the probability distribution function (pdf) of an interval from a four-parameter beta-Poisson model
#'
#' @param x1 Under bound of the interval
#' @param x2 Upper bound of the interval
#' @param alp A non-negative value, the parameter of Beta function (alpha)
#' @param bet A non-negative value, the parameter of Beta function (beta)
#' @param lam1 A non-negative value, the parameter for scaling
#' @param lam2 A non-negative value, the parameter for smoothing, used for BP4 or BP5
#' @return pdf of the interval
#' @export
#' @importFrom stats ppois pnorm rpois rbeta dpois dnorm
#' @importFrom statmod gauss.quad
#' @examples
#' x1=0.2
#' x2=0.4
#' pBPi(x1,x2,alp=0.6,bet=1.5,lam1=20,lam2=0.05)
pBPi = function(x1,x2,alp,bet,lam1=1,lam2=1) {
  if (missing(bet)){
    par=alp;
    alp=par[1]; bet=par[2]; lam1=par[3]; lam2=1; if (length(par)> 3) lam2=par[4];
  }

  ## get interval prob: P(x1 <X <= x2)
  x1=x1/lam2
  x2=x2/lam2
  fn3 = function(x1,x2,m){
    if (max(m) < 100000) res=ppois(x2,m) -ppois(x1,m) else res=pnorm(x2,m,sqrt(m)) -pnorm(x1,m,sqrt(m))
    return(res)
  }
  w = gauss.quad(10,'jacobi', alpha=bet-1, beta=alp-1)
  gs = sum(w$weight*fn3(x1,x2, m=lam1*(1+w$node)/2))
  prob = 1/beta(alp,bet)*2^(-alp-bet+1)*gs
  prob=prob
  return(prob)
}
