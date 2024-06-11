\name{plot_VAR.spectra}
\alias{plot_VAR.spectra}
\title{
Plots spectra of a bivariate VAR (Vector autoregression) model.
}
\description{
Plots in a 2x2 layout the marginal spectra, the squared Coherency spectrum and the Phase spectrum of a VAR model, represented by an object of class \code{var}. If \code{both}=TRUE, each one is calculated (and plotted) in two ways: once directly from \eqn{det, cross, chi.1, chi.2}, and once from the coefficients in \code{ar.list}. For each quantity the values from the two ways of calculation should coincide, if the calculation of the VAR model undertaken in \code{calculate.VAR}  was correct. If \code{both}=FALSE,  only the calculation based on the coefficients in attribute \code{ar.list} is plotted.
}
\usage{
plot_VAR.spectra(a.var, both = TRUE)
}
\arguments{
  \item{a.var}{
An object of class \code{var} for which spectra are to be plotted.
}
  \item{both}{
Logical. If TRUE, each quantity will be plotted calculated in two ways: once directly from \eqn{det, cross, chi.1, chi.2} and once from the coefficients in \code{ar.list}. If FALSE, each quantity will be plotted calculated from the coefficients in \code{ar.list} only.
}
}
\details{
If \code{both}=TRUE the var should have been obtained from \code{calculate.VAR}. If \code{both}=FALSE only the attributes \code{grid, freq} and \code{spec} are required. Thus, it can also be used to plot spectra of a VAR model fitted to a bivariate series, e.g. using function \code{ar}, the spectra of which have previously been calculated by a call to \code{calc.VAR.spec.from.coefs}.
}
\value{
No return value, called for side effects.
}


\seealso{
\code{\link{VAR.spec-package}},\code{\link{Init.var}}, \code{\link{calculate.VAR}},

\code{\link{simulate.VAR}}, \code{\link{plot_VAR.Phase.details}},

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
my.var <- calc.VAR.spec.from.coefs(a.var=my.var)
plot_VAR.spectra(a.var=my.var)

}

