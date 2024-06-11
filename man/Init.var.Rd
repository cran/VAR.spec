\name{Init.var}
\alias{Init.var}
\title{Initializes an object of class \code{var}
}
\description{Initializes an object of class \code{var}, representing a VAR (Vector autoregression) model, sets some of its parameters and reads  the specification of certain polynomials which determine the VAR model from a \code{data.frame} or from a tab-delimited text file (, which is converted to a \code{data.frame}). In this \code{data.frame} each polynomial is specified via it's (complex) roots and constant. If no such text file or \code{data.frame} is specified, an appropriate default \code{data.frame} is created, which can then be modified by the user. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Init.var(grid = 1001, order.max.init = 10, inv.roots.def = NULL,
              a.niter = 5000, a.eps.max.for.UIA = 1e-10, 
              a.eps.for.roots = 1e-05,a.eps.for.spectra=1E-4)
              }
\arguments{
  \item{grid}{
An integer specifying the number of equidistant gridpoints in \eqn{[0,\pi]} used for plotting and numerical approximations. They are given by \eqn{j\pi /((grid-1)), j=0,...,grid}.
}
  \item{order.max.init}{
An integer providing an initial limitation on the order of the VAR model. The degrees of the polynomials specified in the \code{inv.roots.def} should not exceed 2*\code{order.max.init} and the lines of the \code{inv.roots.def} should be at least equal to 6*\code{order.max.init}+2 (if this not the case, appropriate lines will be automatically added). 
}
  \item{inv.roots.def}{
Either NULL or a \code{data.frame} or a tab-delimited text file (, which is converted to a \code{data.frame}). This \code{data.frame} contains the roots (rows) of the polynomials (columns) used to specify the VAR model. See \code{Details} for the specific structure of the text file. If set to NULL, a \code{data.frame} of the appropriate structure will be created, which contains no roots, and corresponds, thus, to a bivariate white noise series. The \code{data.frame} in any of the above cases may be modified after the call of \code{init.var}, via e.g.  \code{var.name$inv.roots <- edit(var.name$inv.roots)}.
}
  \item{a.niter}{
An integer specifying the maximal number of iterations for the Innovations algorithm.
}
  \item{a.eps.max.for.UIA}{
Numeric close to 0. Specifies the threshold for the relative increase in precision achieved at a step of the Innovations algorithm, for accepting that the algorithm has converged. 
}
  \item{a.eps.for.roots}{
Numeric close to 0. A threshold for considering two roots as identical, if a) the inverses of their radii and b) their angles (in their polar representation) differ by less than this threshold.
}
  \item{a.eps.for.spectra}{
Numeric close to 0. A threshold for checking the validity of the calculation of the VAR model. It's spectra are calculated on the basis of the polynomials used to specify the VAR model in \code{inv.roots.def} should coincide with the spectra calculated on the basis of the VAR model's coefficients. If they differ by more than \code{a.eps.for.spectra} warning messages will be written out to the console.
}
}
\details{
The text file [/ \code{data.frame}] passed in argument \code{inv.roots.def} should have a specific structure:
\itemize{
\item The first row of the file (headers) [/ names of columns of the \code{data.frame}] should be 'radius angle det cross chi.1 chi.2 chi.1.prod.2 ma.1 ma.2 eta.1 eta.2 ksi.1 ksi.2 ksi.c zeta'.
\item The second row of the file [/ first row-'observation' in the data.frame] contains the constants of the polynomials, usually set initially to 1, while the first two entries, corresponding to columns radius and angle, are set to '#NA' [/ 'NA'].
\item Subsequent rows, containing each one the specification of a root and the multiplicity of this root for each relevant polynomial, should have the following entries: The first entry is the inverse of its modulus (radius). For all polynomials, except for \code{cross}, this should be in \eqn{[0,1)}. The second entry is it's argument (angle). This should be a figure in \eqn{[0,\pi]}.  Complex conjugates for each non-real root will be automatically added, in order to ensure real coefficients of the polynomials. Subsequent entries in the row are non-negative integers determining the multiplicity of the given root for the polynomial corresponding to the column of the entry. Not all polynomials need to be specified (see parameter \code{calc.method} in function \code{calculate.VAR}). Non relevant polynomials should have '0' as entry for their multiplicity of the root. 
\item Subsequent rows containing no specified roots should have the following structure: the first two entries, corresponding to radius and angle, should equal '#NA' [/ 'NA']. The rest of the row should be filled with '0'. The total number of rows should be at least 6*\code{order.max.init}+2 for the file (, including headers) [/ or 6*\code{order.max.init}+1 for the \code{data.frame}].
}
Two example text files of the appropriate structure, \code{VAR.inv.roots.from.det.cross.txt} and \code{VAR.inv.roots.from.eta.ksi.zeta.txt} may be found in the \code{extdata} folder of the package. A path to them may be obtained via 

\code{fpath<-system.file("extdata","text.file",package="VAR.spec").} 

The corresponding example \code{data.frames} are also directly accessible to the user as 

\code{VAR.inv.roots.from.det.cross} and \code{VAR.inv.roots.from.eta.ksi.zeta}.
}

\value{
An object of class var, basically a list including elements:
\item{grid }{The integer contained in argument \code{grid}. The number of equidistant gridpoints in \eqn{[0,\pi]}, used for plotting and numerical approximations.}
\item{p.max.init}{The integer contained in argument \code{order.max.init}. An initial limitation on the order of the VAR model.}
\item{inv.roots}{A \code{data.frame} containing the roots of certain polynomials, required to specify the VAR model.}
\item{niter}{The integer contained in argument \code{a.niter}. An integer specifying the maximal number of iterations for the Innovations algorithm.}
\item{eps.max.for.UIA}{Numeric close to 0. Specifies the threshold for the relative increase in precision achieved at a step of the Innovations algorithm, for accepting that the algorithm has converged.}
\item{eps.for.roots}{Numeric close to 0. A threshold for considering two roots as identical, if a) the inverses of their radii and b) their angles (in their polar representation) differ by less than this threshold.}
 \item{eps.for.spectra}{
Numeric close to 0. A threshold for checking the validity of the calculation of the VAR model. It's spectra are calculated on the basis of the polynomials used to specify the VAR model in \code{inv.roots.def} should coincide with the spectra calculated on the basis of the VAR model's coefficients. If they differ by more than \code{a.eps.for.spectra} warning messages will be written out to the console.}

}


\seealso{
\code{\link{VAR.spec-package}}, \code{\link{calculate.VAR}}, \code{\link{simulate.VAR}},

\code{\link{plot_VAR.spectra}}, \code{\link{plot_VAR.Phase.details}},

\code{\link{calc.VAR.spec.from.coefs}}
}
\examples{
my.var <- Init.var(grid=201, order.max.init=10, inv.roots.def=NULL)
my.var$inv.roots[2,]<- c(0.98,0.017261,2,3,1,1,2, rep(0,8))
my.var$inv.roots[3,]<- c(0.92,0.897598,2,1,1,1,2, rep(0,8))
my.var$inv.roots[4,]<- c(0.98,1.795196,1,1,0,1,1, rep(0,8))
my.var$inv.roots

}

