% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pep8.R, R/run.R
\name{pepfit}
\alias{pepfit}
\alias{print.pepfit}
\alias{plot.pepfit}
\alias{fitted.pepfit}
\alias{makeFit}
\alias{getFit}
\title{Fit a model}
\usage{
pepfit(TimePoint, RelAb, Channel, data, Elements, Abundance = NULL,
  se = FALSE, time = exists("proc.time"), maxit = 1000,
  setup = parsetup.default, time.unit = NA, name = NULL,
  method = c("binomial", "regression"), useRelAb = FALSE, ...)

makeFit(seq, relAb, isotope.method, ...)

getFit(fit)
}
\arguments{
\item{TimePoint}{The timepoint the observation was taken at}

\item{RelAb}{The relative abundance of that observation}

\item{Channel}{The number of heavy elements for that the observation}

\item{data}{An optional data frame to take the TimePoint, RelAb, and Channel variables
from.  If names for any of the three previous variables are missing, they are
looking for in the data frame using these specific nams. Also, can be \code{RelAbTimes}
object from \code{relAbForTimes}, in which the three previous
variables are ignored.}

\item{Elements}{A list of elements known to be in this peptide and their quantity.  The first
element is the one that will be allowed to change relative abundances
over time.}

\item{Abundance}{A list of elements and their natural abundance.  C, H, O, N, and S
need not be included, as those natural abundances are included by default.}

\item{se}{Boolean, set TRUE to compute standard errors using the Hessian from optim}

\item{time}{Boolean, set FALSE to not measure elapsed time}

\item{maxit}{Maximum iterations to allow optim}

\item{setup}{The desired setup function}

\item{time.unit}{For use in plots and other output}

\item{name}{For use in plots and other output}

\item{method}{Experimental}

\item{useRelAb}{Use relative abundance instead of proportions}

\item{\dots}{Additional parameters passed to the setup function}

\item{seq}{the sequence to get the relative abundance for}

\item{relAb}{the relative abundance to use}

\item{isotope.method}{desired method, passed as type to pepfit}

\item{fit}{a pepfit}
}
\value{
A list of S3 class pepfit, with elements
  \item{par}{Fitted value of desired parameters}
  \item{value}{Fitted likelihood}
  \item{counts}{count data from optim function}
  \item{convergence}{convergence information from optim function}
  \item{message}{message from optim function}
  \item{Day}{The TimePoint input values}
  \item{RelAb}{RelAb, from input}
  \item{Count}{Count, from input}
  \item{Element}{The name, count, and natural abundance of the changing element}
  \item{Elements}{Elements, from input}
  \item{Abundance}{the natural abundance of each element used}
  \item{NatAb}{the combined natural abundance over non-changing elements}
  \item{p}{the output from the setup function}
  \item{setup}{the setup function used}
  \item{parmatrix}{the values of all parameters, both fitted and specified, as in matrix form}
  \item{time}{Elapsed time}
}
\description{
Fit a model of how mass changes over time of a specific element
}
\examples{
data(isodata)
a <- pepfit(TimePoint, RelAb, Channel, data=isodata, Elements=list(N=12,C=45,H=73,O=15))
summary(a)
plot(a)
 
data(isoincorp)
b <- pepfit(TimePoint, RelAb, Channel, data=isoincorp,
type="incorporation", Elements=list(N=12,C=45,H=73,O=15))
summary(b)
plot(b)
}
\seealso{
\code{\link{relAbForTimes}}
}
