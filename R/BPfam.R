#' Family function for beta-Poisson generalized linear model
#'
#' @param alp A non-negative value, the parameter of Beta function (alpha)
#' @param bet A non-negative value, the parameter of Beta function (beta)
#' @param lam1 A non-negative value, the parameter for scaling
#' @param lam2 A non-negative value, the parameter for smoothing, used for BP4 or BP5
#' @param link The link function, as a character string, similar to quasi-poisson family, typically log link
#' @return An object of class "family" for function glm() to fit a beta-Poisson generalized linear model
#' @export
#' @examples
#' set.seed(2015)
#' control.dat=rBP(100,alp=0.6,bet=1.5,lam1=20,lam2=0.05)
#' treated.dat=rBP(100,alp=0.2,bet=5,lam1=200,lam2=0.1)
#' x=c(control.dat,treated.dat)
#' y=c(rep(1,length(control.dat)),rep(2,length(treated.dat)))
#' fdat=data.frame(x=x,group=y)
#' fam0=do.call("BPfam", list(alp=0.6,bet=1.5,lam1=20,lam2=0.05, link = "log"))
#' fit=glm(x~.,data=fdat,family=fam0)
#' fit
#' summary(fit)$coefficients
BPfam <-function (alp, bet, lam1, lam2, link = "log"){
    #### begin: codes modified from quasipoisson function of stats package
    linktemp <- substitute(link)
    if (!is.character(linktemp)) 
        linktemp <- deparse(linktemp)
    okLinks <- c("log", "identity", "sqrt")
    if (linktemp %in% okLinks) 
        stats <- make.link(linktemp)
    else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    }
    else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name)) 
                linktemp <- stats$name
        }
        else {
            stop(gettextf("link \"%s\" not available for betaPoisson family; available links are %s", 
                linktemp, paste(sQuote(okLinks), collapse = ", ")), 
                domain = NA)
        }
    }
    #### end: codes modified from quasipoisson function of stats package
######################################################
    #### codes for beta_Poisson model
#    phi1=alp/(alp+bet)
    phi2=bet/(alp*(alp+bet+1))
#    mymu=lam2*lam1*phi1
#    myVar=mymu*lam2+mymu^2*phi2
    variance <- function(mu){ mu*lam2 + mu^2 * phi2}

    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the 'betaPoisson' family")
        n <- rep.int(1, nobs)
        mustart <- y + (y == 0)/1e9

    })
###################################################### 
    #### begin: codes modified from quasipoisson function of stats package
    validmu <- function(mu) all(is.finite(mu)) && all(mu > 0)
    aic <- function(y, n, mu, wt, dev) NA
    dev.resids <- function(y, mu, wt) {
        r <- mu * wt
        p <- which(y > 0)
        r[p] <- (wt * (y * log(y/mu) - (y - mu)))[p]
        2 * r
    }
    stats <- make.link(link)
    #### end: codes modified from quasipoisson function of stats package

    famname <- paste("Beta Poisson(", format(round(alp,4)), format(round(bet,4)),format(round(lam1,4)),format(round(lam2,4)), ")", sep = "")
    structure(list(family = famname, link = link, linkfun = stats$linkfun, aic = aic,dev.resids = dev.resids,
        linkinv = stats$linkinv, variance = variance,  mu.eta = stats$mu.eta, validmu = validmu,initialize = initialize,
         valideta = stats$valideta),
        class = "family")
}
