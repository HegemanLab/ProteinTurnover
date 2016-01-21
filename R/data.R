#' Sample Count Data
#'
#' A sample of the kind of count data used by this package.
#'
#' Count data like this can be converted to relative abundances using relAbFromCounts.
#'
#' @name isocounts
#' @format   A data frame with 1224 observations on the following 4 variables.
#'  \describe{
#'    \item{\code{TimePoint}}{Day the observation was made on}
#'    \item{\code{Channel}}{Channel the observation was on (that is, how
#'      many ``heavy'' elements there were)}
#'    \item{\code{RT}}{Retention time the observation was made at}
#'    \item{\code{Count}}{Number of items observed}
#'  }
#' @examples
#' data(isocounts)
#' head(isocounts)
#' counts <- within(isocounts,{
#'   Channel <- factor(Channel)
#'   TimePoint <- factor(TimePoint)
#'   logCount <- log10(Count)
#' })
#' coplot(logCount~RT|Channel*TimePoint, data=counts, show.given=FALSE)
#' @keywords datasets
NULL

#' Sample Peptide/Isotope Data
#'
#'   A sample of the kind of data used by this package.
#'
#' @name isodata
#' @format   A data frame with 68 observations on the following 3 variables.
#'  \describe{
#'    \item{\code{TimePoint}}{The day the observation was made}
#'    \item{\code{Channel}}{Number of ``heavy'' elements}
#'    \item{\code{RelAb}}{The relative abundance of that count on that
#'      day, compared to the other counts on that day}
#'  }
#' @examples
#' data(isodata)
#' head(isodata)
#' coplot(RelAb~Channel|factor(TimePoint), data=isodata)
#' pepfit(TimePoint, RelAb, Channel, data=isodata, Elements=list(N=12,C=45,H=73,O=15))
#' @keywords datasets
NULL

#' Sample Peptide/Isotope Data
#'
#'   A sample of the incorporation data used by this package.
#'
#' @name isoincorp
#' @format   A data frame with 73 observations on the following 3 variables.
#'  \describe{
#'    \item{\code{TimePoint}}{The day the observation was made}
#'    \item{\code{Channel}}{Number of ``heavy'' elements}
#'    \item{\code{RelAb}}{The relative abundance of that count on that
#'      day, compared to the other counts on that day}
#'  }
#' @examples
#' data(isoincorp)
#' head(isoincorp)
#' coplot(RelAb~Channel|factor(TimePoint), data=isoincorp)
#' pepfit(TimePoint, RelAb, Channel, data=isoincorp,
#' Elements=list(N=12,C=45,H=82,O=12), type="incorporation")
#' @keywords datasets
NULL
