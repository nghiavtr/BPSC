% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/varBP.R
\name{varBP}
\alias{varBP}
\title{Get variance of a beta-Poisson model}
\usage{
varBP(alp, bet, lam1 = 1, lam2 = 1, prob0 = 0)
}
\arguments{
\item{alp}{A non-negative value, the parameter of Beta function (alpha)}

\item{bet}{A non-negative value, the parameter of Beta function (beta)}

\item{lam1}{A non-negative value, the parameter for scaling}

\item{lam2}{A non-negative value, the parameter for smoothing, used for BP4 or BP5}

\item{prob0}{A non-negative value in [0;1), the parameter for the portion of fixed density of zero in the distribution, used for BP5}
}
\value{
The variance of a beta-Poisson model
}
\description{
Get variance of a beta-Poisson model
}
\examples{
varBP(alp=0.6,bet=1.5,lam1=20,lam2=0.05)
}

