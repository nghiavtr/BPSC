#' Get mean of a beta-Poisson model
#'
#' @param alp A non-negative value, the pareter of Beta function (alpha)
#' @param bet A non-negative value, the pareter of Beta function (beta)
#' @param lam1 A non-negative value, the pareter for scaling
#' @param lam2 A non-negative value, the pareter for smoothing, used for BP4 or BP5
#' @param prob0 A non-negative value in [0;1), the pareter for the portion of fixed density of zero in the distribution, used for BP5
#' @return The mean of the beta-Poisson model
#' @export
#' @examples
#' meanBP(alp=0.6,bet=1.5,lam1=20,lam2=0.05)
meanBP = function(alp,bet,lam1=1,lam2=1,prob0=0){
  if (missing(bet)){
  	par=alp;
  	alp=par[1]; bet=par[2]; lam1=par[3]; lam2=1; if (length(par)> 3) lam2=par[4]; prob0=0;if (length(par)> 4) prob0=par[5]
  }
  m4=lam2*lam1*alp/((alp+bet))
  m5=m4*(1-prob0)
  res=m5
  return(res)
}