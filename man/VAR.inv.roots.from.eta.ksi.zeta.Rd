\name{VAR.inv.roots.from.eta.ksi.zeta}
\alias{VAR.inv.roots.from.eta.ksi.zeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An example \code{data.frame} defining a VAR model.
}
\description{
An example \code{data.frame} to be used for argument \code{inv.roots.def} of function \code{Init.var}, when subsequently argument \code{method} of function \code{calculate.VAR} is set to "from.eta.ksi.zeta".
}
\usage{
VAR.inv.roots.from.eta.ksi.zeta
}
%- maybe also 'usage' for other objects documented here.
\format{A \code{data.frame} containing 3 roots (rows) and their multiplicities for the polynomials \eqn{eta.1(z), eta.2(z)}, \eqn{ksi.1(z), ksi.2(z)} and \eqn{zeta(z)}, as well as their constants in the first row.}


\seealso{
\code{Init.var}, \code{calculate.VAR}
}
\examples{

my.var <- Init.var(grid=501, order.max.init=10, 
                        inv.roots.def=VAR.inv.roots.from.eta.ksi.zeta)
my.var <- calculate.VAR(a.var=my.var,calc.method="from.eta.ksi.zeta")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
