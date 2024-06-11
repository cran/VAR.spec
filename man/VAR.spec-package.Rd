\name{VAR.spec-package}
\alias{VAR.spec-package}
\alias{VAR.spec}
\docType{package}
\title{
\packageTitle{VAR.spec}
}
\description{
\packageDescription{VAR.spec}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{VAR.spec}
\packageIndices{VAR.spec}

The specification of the VAR (Vector autoregression) model is based on the following fact (see Ioannidis (2007))

<doi:10.1016/j.jspi.2005.12.013>: 

For any four complex polynomials \eqn{det(z), cross(z)} of degree \eqn{2p} and \eqn{chi.1(z), chi.2(z)} of degree \eqn{p}, satisfying
  
  EQ(*) \deqn{|det(z)|^2+|cross(z)|^2 = |chi.1(z)|^2 * |chi.2(z)|^2} on \eqn{|z|=1}, there exists a bivariate  VAR(p) with marginal spectra \deqn{f.i (w) = (1/2\pi) |chi.i(z)|^2 /|det(z)|^2,} and cross-spectrum \deqn{f.1.2 (w) = (1/2\pi) z^{(-p)}cross(z) /|det(z)|^2,} where \eqn{z=exp(-iw)}. The squared Coherency is then given by \deqn{|cross(z)|^2 /(|det(z)|^2+|cross(z)|^2).}
  
  The idea in this package is to define the necessary polynomials by specifying their roots, instead of their coefficients. This is an idea  borrowed from the Time Series Library of R. Dahlhaus,  where it is used for defining ARMA models for univariate time series. 
  
  Moreover, the package allows the user to specify only some of the roots of \eqn{det(z), cross(z)} and \eqn{chi.1(z), chi.2(z)}, while it attempts to find further non-specified roots in a way such that EQ(*) is satisfied. 
  
  By specifying certain roots of \eqn{det(z), cross(z)} and \eqn{chi.1(z), chi.2(z)} one can induce desired features in the spectra and in the series. For example,
  \itemize{
  \item if one wishes a strong presence of frequency \eqn{w_0} components in series 1,
  there should be a root with modulus close to 1 and angle equal to \eqn{w_0}, the multiplicity of which for \eqn{det(z)} is higher than its multiplicity for \eqn{chi.1(z)}, forcing the spectrum of series 1 to have a peak at \eqn{w_0}. 
  \item If the frequency-\eqn{w_0}-components should be strongly correlated between series 1 and series 2, the multiplicity of this root for \eqn{det(z)} should be higher than its multiplicity for \eqn{cross(z)}, forcing the Coherency at \eqn{w_0} close to 1.
  

  }
    Required polynomials (see argument \code{calc.method} of function \code{calculate.VAR} and related \code{Details}) are first passed to function \code{Init.var}, via a \code{data.frame} or a text file, which contains the multiplicities of desired roots (rows) for each required polynomial (columns).
  
  After calling \code{Init.var}, function \code{calculate.VAR} must be called, which attempts to find a VAR model which is compatible with the polynomials specified in its attribute \code{inv.roots}, by making necessary adjustments, calculates its spectra, coefficients and order, makes the necessary checks and plots it's spectra. Then, one can simulate from the specified model calling \code{simulate.VAR}.
}
\author{
\packageAuthor{VAR.spec}

Maintainer: \packageMaintainer{VAR.spec}
}
\references{
Ioannidis, E. E. (2007). Spectra of bivariate VAR(p) models. Journal of Statistical Planning and Inference 137(2), 554-566. 

Ioannidis, E. E. and Chronis, G. A. (2005). Extreme spectra of VAR models and orders of near-cointegration. J. Time Ser. Anal. 26, 399-421.}

\keyword{ package }
\seealso{
\code{\link{Init.var}}, \code{\link{calculate.VAR}}, \code{\link{simulate.VAR}},

\code{\link{plot_VAR.spectra}}, \code{\link{plot_VAR.Phase.details}},

\code{\link{calc.VAR.spec.from.coefs}}
}
\examples{
my.var <- Init.var(grid=501, order.max.init=10, inv.roots.def=NULL)
my.var$inv.roots[2,]<- c(0.98,0.017261,2,3,1,1,2, rep(0,8))
my.var$inv.roots[3,]<- c(0.92,0.897598,2,1,1,1,2, rep(0,8))
my.var$inv.roots[4,]<- c(0.98,1.795196,1,1,0,1,1, rep(0,8))
my.var <- calculate.VAR (a.var = my.var, calc.method="from.det.cross", 
            plot.spectra=TRUE,suppr.spec.check.warn=TRUE)
print(my.var$validity.msg)

}
