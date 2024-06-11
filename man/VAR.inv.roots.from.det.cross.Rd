\name{VAR.inv.roots.from.det.cross}
\alias{VAR.inv.roots.from.det.cross}
\title{
An example \code{data.frame} defining a VAR (Vector autoregression) model.
}
\description{
An example \code{data.frame} to be used for argument \code{inv.roots.def} of function \code{Init.var}, when subsequently argument \code{method} of function \code{calculate.VAR} is set to "from.det.cross".
}
\usage{
VAR.inv.roots.from.det.cross
}
\format{A \code{data.frame} containing 3 roots (rows) and their multiplicities for the polynomials \eqn{det(z), cross(z)}, \eqn{chi.1(z)} and \eqn{chi.2(z)}, as well as their constants in the first row.}


\seealso{
\code{Init.var}, \code{calculate.VAR}
}
\examples{

my.var <- Init.var(grid=501, order.max.init=10, inv.roots.def=VAR.inv.roots.from.det.cross)
my.var <- calculate.VAR(a.var=my.var,calc.method="from.det.cross",
                        suppr.spec.check.warn=TRUE)
print(my.var$validity.msg)
}

