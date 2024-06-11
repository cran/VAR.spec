\name{simulate.VAR}
\alias{simulate.VAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates a bivariate series from a bivariate VAR (Vector autoregression) model.
}
\description{
Simulates a Gaussian bivariate series from a VAR model, described by an object of class \code{var}.
}
\usage{
simulate.VAR(a.var, sample.size, burn.in = 1000)
}
\arguments{
  \item{a.var}{
An object of class \code{var}, obtained from \code{calculate.VAR}. Actually only the attributes \code{order} and \code{ar.list} of \code{a.var} are needed.
}
  \item{sample.size}{
Integer. The size of the bivariate sample to be generated.
}
  \item{burn.in}{
Integer. The number of initial observations to be discarded. The bivariate series is simulated by iterating the VAR recursion, starting with white noise (with the covariance structure of the innovations, given in \code{ar.list$var.pred} ). This iteration converges to the correct distribution after a \code{burn.in} period, provided the roots of \eqn{det} are outside the unit circle, i.e. their inverses have a modulus less than 1. The closer roots of \eqn{det} are to the unit circle, the longer the \code{burn.in} period should be chosen.
}
}

\value{
A numeric array of dimensions \code{sample.size}x2, the bivariate sample series simulated from the VAR model.

}



\seealso{
\code{\link{VAR.spec-package}}, \code{\link{Init.var}}, \code{\link{calculate.VAR}}, 
\code{\link{plot_VAR.spectra}}, \code{\link{plot_VAR.Phase.details}},
\code{\link{calc.VAR.spec.from.coefs}}
}
\examples{

my.var <- Init.var(grid=501, order.max.init=10, inv.roots.def=NULL)
my.var$inv.roots[2,]<- c(0.98,0.017261,2,3,1,1,2, rep(0,8))
my.var$inv.roots[3,]<- c(0.92,0.897598,2,1,1,1,2, rep(0,8))
my.var$inv.roots[4,]<- c(0.98,1.795196,1,1,0,1,1, rep(0,8))
my.var <- calculate.VAR(a.var=my.var,calc.method="from.det.cross",
                        suppr.spec.check.warn=TRUE)
print(my.var$validity.msg)
my.sample <-simulate.VAR(a.var=my.var, sample.size=250, burn.in = 500)

}

