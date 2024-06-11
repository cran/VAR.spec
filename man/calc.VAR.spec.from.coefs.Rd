\name{calc.VAR.spec.from.coefs}
\alias{calc.VAR.spec.from.coefs}
\title{
Calculates the spectral matrix of a multivariate VAR (Vector autoregression) model.
}
\description{
Calculates the spectral matrix for all grid points for any (multivariate) VAR model, represented by an object of class \code{var}, with specified coefficients (in attribute \code{ar.list}). 
}
\usage{
calc.VAR.spec.from.coefs(a.var)}

\arguments{
  \item{a.var}{
an object of class \code{var} for which the spectral matrix will be calculated. Actually only the attributes \code{grid, order} and \code{ar.list} of \code{a.var} are needed.}
}
\details{
Is being called in \code{calculate.VAR}, so there is no need to call it again after \code{calculate.VAR} has been called. Can also be used to calculate the spectral matrix of a VAR model which has been fitted to data, e.g. using function \code{ar}. See Brockwell and Davis, 1990, Example 11.8.1, Hannan, 1970, Chapter II, Section 5.iv or Brillinger, 2001, Sections 2.8 and 2.9.
}
\value{

Returns the object of class \code{var} after having calculated and set the following attributes:
\item{freq}{A one dimensional array of length \code{grid}, containing the grid-points at which \code{spec} has been calculated. }
\item{spec}{A three dimensional array (\code{grid}xkxk), where k is the number of series of the VAR model. Contains for each grid-point the spectral matrix, as calculated at the basis of \code{ar.list}. That is, \code{spec[j,1,1]} is the spectrum of series 1 at grid-point j, \code{spec[j,1,2]} is the cross-spectrum between series 1 and 2 at grid-point j, etc.}
}
\references{
Brillinger, D.R. (2001). Time Series: Data Analysis and Theory, second ed. Classics in Applied Mathematics, vol. 36. SIAM, Philadelphia.

Brockwell, P. J., Davis, R. A. (1990). Time Series: Theory and Methods, second ed. Springer, New York.

Hannan, E. J. (1970). Multiple Time Series. Wiley, New York.

}



\seealso{
\code{\link{VAR.spec-package}},\code{\link{Init.var}}, \code{\link{calculate.VAR}}, \code{\link{simulate.VAR}},

\code{\link{plot_VAR.spectra}}, \code{\link{plot_VAR.Phase.details}}
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
est.var<-list()
est.var$ar.outp <- ar(x=my.sample, aic=FALSE, order.max=6, demean=FALSE, method="ols", 
		series=c("ser1","ser2"))
est.var$label <-"MY VAR(6)"
est.var$order<-dim(est.var$ar.outp$ar)[1]
est.var$ar.list$ar <- est.var$ar.outp$ar
est.var$ar.list$var.pred <- est.var$ar.outp$var.pred
est.var$grid <- 501 
est.var <- calc.VAR.spec.from.coefs(est.var)
plot_VAR.spectra(a.var=est.var,both=FALSE)


}

