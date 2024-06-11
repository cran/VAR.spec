\name{calculate.VAR}
\alias{calculate.VAR}
\title{
Attempts to define a bivariate VAR (Vector autoregression) model.
}
\description{
Attempts to find a VAR model which is compatible with the polynomials specified in its attribute \code{inv.roots}, makes necessary adjustments, calculates its spectra, coefficients and order, makes the necessary checks and plots the spectra. Finding a VAR model involves finding complex polynomials \eqn{det(z), cross(z), chi.1(z)} and \eqn{chi.2(z)}, satisfying 
EQ(*): \deqn{|det(z)|^2+|cross(z)|^2 = |chi.1(z)|^2 |chi.2(z)|^2} on \eqn{z|=1}. If \code{calculate.VAR} is successful in finding such polynomials, there exists a bivariate VAR(p) model with marginal spectra \deqn{f.i (w) = (1/2\pi)* |chi.i(z)|^2 /|det(z)|^2,} and cross-spectrum \deqn{f.1.2 (w) = (1/2\pi)* z^{(-p)}*cross(z) /|det(z)|^2,} where \eqn{z=exp(-iw)}. The squared Coherency is then given by \deqn{|cross(z)|^2 /(|det(z)|^2+|cross(z)|^2).}
}
\usage{
calculate.VAR(a.var, calc.method = "from.det.cross", M.fact = 1.1, 
              plot.spectra = TRUE,suppr.spec.check.warn=FALSE)
}
\arguments{
  \item{a.var}{
An object of class \code{var}, previously initialized by \code{Init.var}.
}
  \item{calc.method}{
One of 'from.det.cross' or 'from.eta.ksi.zeta'. See details for further explanations.
}
  \item{M.fact}{
Numeric. A factor >1, by which the \eqn{min (|ksi.c(z) /[eta.1(z)eta.2(z)zeta(z)]|^2 )} is multiplied, in order to ensure that \eqn{|chi.1(z)*chi.2(z)|^2 - |det(z)|^2} is positive on \eqn{|z|=1}. Only relevant when \code{calc.method} = 'from.eta.ksi.zeta'. Affects the range of values of the Coherency spectrum.
}
  \item{plot.spectra}{
Logical. If TRUE, the marginal spectra, Coherency and Phase spectra will be plotted.
}
\item{suppr.spec.check.warn}{
Logical. If TRUE, warnings on whether spectra obtained directly from \eqn{cross, det, chi.1} and \eqn{chi.2}  (in \code{spec.1, spec.2, Coher, Phase}) differ by more than \code{eps.for.spectra} from the spectra (in \code{spec}) obtained from the coefficients of the VAR model (see \code{Note}, below)  are suppressed.
}
}
\details{
In order to define a VAR model, one seeks to define \eqn{det(z), cross(z)} and \eqn{chi.1(z), chi.2(z)}, such that EQ(*) holds. There are two ways implemented for doing so, which are specified by argument \code{calc.method}:
\itemize{
\item \code{calc.method} = 'from.det.cross'. (See Ioannidis (2007), Proposition 2). One must first specify the roots for \eqn{det(z)} and \eqn{cross(z)}. They fix these two polynomials, which will not be modified during the calculation. One may also specify some roots for \eqn{chi.1(z)} and \eqn{chi.2(z)}. If these two latter are compatible with EQ(*) they will be kept. If not, they may be modified during the calculation. Further roots may be added by \code{compute.VAR} during the calculation to \eqn{chi.2} to ensure that EQ(*) holds. All changes are recorded in attribute \code{inv.roots}.  If the user wishes some of the roots of \eqn{chi.2} to be rather added to \eqn{chi.1}, he may do so by editing attribute \code{inv.roots} after running the function once, by increasing the multiplicity of the root for \eqn{chi.1} and decreasing it by the same amount for \eqn{chi.2}, and then, re-running the function. With the option \code{calc.method} = 'from.det.cross' one has more direct control over peaks for the marginal spectra (via \eqn{chi.i} and \eqn{det}) and the Coherency being close to 0 or 1 (via \eqn{det} and \eqn{cross}), as well as the Phase spectrum, i.e., the lag-/lead- structure between the two series (via \eqn{cross}). 
\item \code{calc.method} = 'from.eta.ksi.zeta'. (See Ioannidis (2007), Proposition 1). One must specify the roots of \eqn{eta.1, eta.2, ksi.1, ksi.2} and \eqn{zeta}. These fix \deqn{chi.i(z) = M* eta.i(z) zeta(z) ksi.j(z)/ksi.c(z))} and \deqn{det(z)=M *zeta(z) ksi.1(z) ksi.2(z)/ksi.c(z),} where \eqn{M} equals \code{M.fact}*\eqn{min (|ksi.c(z) /[eta.1(z) eta.2(z) zeta(z)]|^2 )} and \eqn{ksi.c} is the maximal common divisor of \eqn{ksi.1(z)} and \eqn{ksi.2(z)}. Then, \eqn{cross(z)} is automatically determined so as that EQ(*) is satisfied. The result is that marginal spectra are given by \deqn{f.i (w) = |eta.i(z)|^2 /|ksi.i(z)|^2} and the squared Coherency is given by \deqn{1-|ksi.c(z)|^2 /|M* eta.1(z) eta.2(z) zeta(z)|^2.} With the option \code{calc.method} = 'from.eta.ksi.zeta' one has more direct control over peaks (via \eqn{ksi}) and zeroes (via \eqn{eta}) of the marginal spectra and, partially, over the Coherency (via \eqn{ksi.c, eta} and \eqn{zeta}), but not of the Phase spectrum. However, for certain aspects of the Coherency, a more direct control over the relation of \eqn{cross} and \eqn{det} would be required, which is not available with this option. For example, with this option, it is not possible to set the multiplicity of a root to be higher for \eqn{cross} than for \eqn{det}.
}
}
\value{
An object of class var, basically a list including elements:
\item{grid }{The number of equidistant gridpoints in \eqn{[0,\pi]}, used for plotting and numerical approximations.}
\item{p.max.init}{An initial limitation on the order of the VAR model.}
\item{inv.roots}{A \code{data.frame} containing the roots of certain polynomials, used to specify the VAR model.}
\item{niter}{An integer specifying the maximal number of iterations for the Innovations algorithm.}
\item{eps.max.for.UIA}{Numeric close to 0. Specifies the threshold for the relative increase in precision achieved at a step of the Innovations algorithm, for accepting that the algorithm has converged.}
\item{eps.for.roots}{Numeric close to 0. A threshold for considering two roots as identical, if a) the inverses of their radii and b) their angles (in their polar representation) differ by less than this threshold.}
 \item{eps.for.spectra}{
Numeric close to 0. A threshold for checking the validity of the calculation of the VAR model. It's spectra are calculated on the basis of the polynomials used to specify the VAR model in \code{inv.roots.def} should coincide with the spectra calculated on the basis of the VAR model's coefficients. If they differ by more than \code{a.eps.for.spectra} warning messages will be written out to the console.}
\item{validity.msg}{A message containing more precise information on a) the convergence of the Innovations Algorithm and b) on the maximal difference between the spectra of the VAR model as calculated on the basis of the polynomials used to specify it and the spectra obtained on the basis of the VAR model's coefficients. The message should be printed after each call to \code{calculate.var}.}
\item{order}{(Integer) The final order of the VAR model.}
\item{ar.list}{A list with entries (calculated solving the Yule-Walker equations): 
\itemize{
  \item \code{order} (Integer) The final order of the VAR model. 
  \item	\code{ar} A three dimensional array (orderx2x2), containing the coefficients of the VAR model. 
  \item	\code{var.pred} A 2x2 matrix containing the covariance matrix of the innovations.}
  }
\item{cov.1, cov.2, cov.cross, cov.cross.neg}{Three one dimensional arrays of length \code{order}+1 containing the auto- and the cross- covariances (for lags and leads respectively) of the two series for lags=0,...,\code{order}. They are being calculated as Dirichlet-approximations to the Fourier transforms of the spectra, which have been calculated from \eqn{det, cross, chi.1} and \eqn{chi.2}.}
\item{freq}{A one dimensional array of length \code{grid}. It contains the grid-points at which spectra are being calculated.} 
\item{spec}{A three dimensional array (\code{grid}x2x2). Contains for each grid-point the spectral matrix, as calculated \bold{at the basis of the VAR model's coefficients}, given in attribute \code{ar.list}. That is, \code{spec}[j,1,1] is the spectrum of series 1 at grid-point j, \code{spec}[j,1,2] is the cross-spectrum between the series 1 and 2 at grid-point j, etc.}
\item{spec.1, spec.2}{Two one dimensional arrays of length \code{grid} containing the \bold{log} of the values of the spectra of series 1 and 2, respectively, \bold{as calculated from \eqn{chi, det}}, i.e.  \eqn{\code{spec.1}[j]=log[(1/2\pi)  |chi.i(z)|^2 /|det(z)|^2]}, where \eqn{z=exp(-ij\pi/(\code{grid}-1))}.}
\item{Coher}{A one dimensional array of length \code{grid} containing the values of the log of the squared Coherency spectrum between series 1 and 2, \bold{as calculated from \eqn{cross, det}}, i.e.  \eqn{\code{Coher}[j] = log(|cross(z)|^2 /(|det(z)|^2+|cross(z)|^2])}, \eqn{z=exp(-ij\pi/(\code{grid}-1))}.}
\item{Phase}{A one dimensional array of length \code{grid} containing the values of the Phase spectrum between series 1 and 2, \bold{as calculated from \eqn{cross}}, i.e., \eqn{\code{Phase} (w) = arg[ z^{(-p)}*cross(z)]} , where \eqn{z=exp(-ij\pi/(\code{grid}-1))}.}
\item{Phase.div.freq, group.delay }{Two one dimensional arrays of length \code{grid} containing the values of \code{Phase}(w)/w, which gives the lead/lag in units of time, and the derivative of \code{Phase}(w) with respect to w, respectively.}
\item{det, cross, chi.1, chi.2, ksi.1, ksi.2, eta.1, eta.2, zeta}{Contain detailed information on the respective polynomials (the last five only if \code{calc.method} = 'from.eta.ksi.zeta' was used), as modified during the calculation, each one with attributes:
\itemize{
\item \code{const} the constant of the polynomial;
\item \code{inv.roots} the modulus and angle for the inverses of the roots of the polynomial (for conjugate pairs only the one with angle in \eqn{[0,\pi]});
\item	\code{inv.roots.number}, the number of the roots of the polynomial (pairs of complex roots are counted only once);
\item	\code{order}, the degree of the polynomial;
\item	\code{coefs}, an array of length \eqn{2( \code{inv.roots.number} +1)} containing the coefficients of the polynomial;
\item	\code{fourier.coefs}, an array of length \eqn{2\code{pmax.init}+1} containing the Fourier transform of the \eqn{|polynomial(z)|^2} on \eqn{|z|=1};
\item	\code{values$spec}, an array of length \code{grid}, containing the values of the squared modulus of the polynomial, i.e. \eqn{|polynomial(z)|^2} on \eqn{|z|=1};
\item \code{Inv.values$spec}, an array of length \code{grid}, containing the \bold{log of the inverses} of the values of the squared modulus of the polynomial, i.e. \eqn{-log [ |polynomial(z)|^2 ]} on \eqn{|z|=1}.
}
}

}
\references{
Ioannidis, E. E. (2007). Spectra of bivariate VAR(p) models. Journal of Statistical Planning and Inference 137(2), 554-566. 

Ioannidis, E. E. and Chronis, G. A. (2005). Extreme spectra of VAR models and orders of near-cointegration. J. Time Ser. Anal. 26, 399-421.

}

\note{
\code{calculate.VAR} performs a check on the validity of the output and stores relevant information in attribute \code{validity.msg}. There are certain reasons for which the output may  not be valid: First, the calculation involves numerical approximations; the relevant algorithms may have not achieved convergence to the desired precision. Second, there may be features in the spectra, which cannot be well represented by a VAR model, resulting in polynomials \eqn{cross, det, chi.1} and \eqn{chi.2} not satisfying EQ(*). Thus, it is imperative to check a) the validity of EQ(*). It is also important to check b) that the spectra obtained directly from \eqn{cross, det, chi.1} and \eqn{chi.2}  (in \code{spec.1, spec.2, Coher, Phase}) coincide with the spectra (in \code{spec}) obtained from the coefficients of the VAR model (in \code{ar.list}). The checks performed by \code{calculate.VAR} concern exactly the checks under a) and b): the maximum absolute difference between the relevant quantities (in log scale, where appropriate) is passed in \code{validity.msg}. Moreover, warnings are written out to the console if the Innovations algorithm did not converge or if any of the quantities in 'b' exceeds \code{eps.for.spectra}, provided that \code{suppr.spec.check.warn}=FALSE. The same check is conducted visually if \code{plot.spectra}=TRUE: each spectrum is calculated and plotted according to both methods. The two should coincide.

The option \code{method}='from.det.cross' seems generally numerically more stable than the option

\code{method}='from.eta.ksi.zeta'.

\bold{Suggestions on how to obtain desired spectral characteristics:}

Let \eqn{z_0:= r_0^{-1}exp(iw_0)}, with \eqn{r_0} close to 1. 
For a polynomial \eqn{p(z)} denote by \eqn{[p(z_0)]} the multiplicity of \eqn{z_0} as root of \eqn{p(z)}.

To obtain
\itemize{

\item a peak in \eqn{f.1(w_0)}
\itemize{
\item with \code{method}='from.det.cross': choose \eqn{[det(z_0)] > [chi.1(z_0)]}.
\item with \code{method}='from.eta.ksi.zeta': choose \eqn{[ksi.1(z_0)] > [eta.1(z_0)]}.
}


\item a trough in \eqn{f.1(w_0)}
\itemize{
\item with \code{method}='from.eta.ksi.zeta': choose \eqn{[ksi.1(z_0)] < [eta.1(z_0)]}.
}


\item a \eqn{Coh(w_0)} close to 1
\itemize{
\item with \code{method}='from.det.cross': choose \eqn{[det(z_0)] > [cross(z_0)]}.
\item with \code{method}='from.eta.ksi.zeta': choose \eqn{min([ksi.1(z_0)],[ksi.2(z_0)]) > [eta.1(z_0)]+[eta.2(z_0)]+[zeta(z_0)]}.
}

\item a \eqn{Coh(w_0)} close to 0
\itemize{
\item with \code{method}='from.det.cross': choose \eqn{[det(z_0)] < [cross(z_0)]}.
\item with \code{method}='from.eta.ksi.zeta': choose \eqn{min([ksi.1(z_0)],[ksi.2(z_0)]) < [eta.1(z_0)]+[eta.2(z_0)]+[zeta(z_0)]}.
}


}
}


\seealso{
\code{\link{VAR.spec-package}},\code{\link{Init.var}}, \code{\link{simulate.VAR}},

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

