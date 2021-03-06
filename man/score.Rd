% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lackfit.R
\name{score}
\alias{score}
\alias{score.visual}
\alias{score.dev}
\alias{score.N}
\title{Get a score for a given fit}
\usage{
score.visual(obj)

score.dev(obj, level = 0.05)

score.N(obj, level = 0.05)
}
\arguments{
\item{obj}{a pepfit object}

\item{level}{significance level to use to get associated sample size}
}
\value{
\item{overall}{Overall Score}
 \item{byday}{Score for each day separately}
}
\description{
Get a score for a given fit
}
\details{
The score is a value that in some way measures how
 well the model fits the actual data.

score.N returns simply the sample size associated with the given
 significance level

score.dev converts this sample size to a value between 0 and 100.
 1000 corresponds to roughly 80, 5000 to roughly 90, 20000 to roughly
 95.

score.visual is the ratio of the total difference between the fitted
 and observed values to the total observed values.
}
\examples{
data(isodata)
a0 <- pepfit(data=isodata, Elements = list(N = 12,
     C = 45, H = 73, O = 15))
score.N(a0)
score.dev(a0)
score.visual(a0)
}
\seealso{
\code{\link{pepfit}}
}
