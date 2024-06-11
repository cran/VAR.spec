\name{plot_VAR.Phase.details}
\alias{plot_VAR.Phase.details}
\title{
Plots details related to the Phase spectrum of a bivariate VAR (Vector autoregression) model.
}
\description{
Plots in a 2x2 layout the Phase spectrum of a VAR model for a bivariate series, represented by an object of class \code{var}, the Phase spectrum divided by frequency -which gives the lag-/lead-structure between the two series in units of time --, the derivative of the Phase spectrum with respect to frequency -also known as the group-delay--, and the Coherency spectrum of the VAR model.}
\usage{
plot_VAR.Phase.details(a.var)
}
\arguments{
  \item{a.var}{
An object of class \code{var} obtained from \code{calculate.VAR} for which details of the Phase spectrum are to be plotted.}
}
\value{
No return value, called for side effects.
}


\seealso{
\code{\link{VAR.spec-package}},\code{\link{Init.var}}, \code{\link{calculate.VAR}},

\code{\link{simulate.VAR}}, \code{\link{plot_VAR.spectra}}, 

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
plot_VAR.Phase.details (a.var=my.var)

}

