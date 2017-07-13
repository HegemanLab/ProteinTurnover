###################################################
### fitting exponential curve functions
###################################################

expcurve<-function(x, intercept, limit, halflife, vary=TRUE) {
  p <- exp((log(1/2)*x)/halflife)
  (intercept-limit)*p+limit
}

#' @export
fitted.pepcurve <- function(object, TimePoint=sort(unique(object$data$TimePoint)), ...) {
  out <- data.frame(TimePoint=TimePoint)
  for(i in 1:nrow(object$par)) {
    nn <- rownames(object$par)[i]
    out[[nn]] <- expcurve(out$TimePoint,
                          object$par[nn,"intercept"],
                          object$par[nn,"limit"],
                          object$par[nn,"halflife"])
  }
  out
}

#' @export
plot.pepcurve <- function(x, ylab="proportion", ylim=c(0,1), ...) {
  maxtime <- ceiling(max(x$data$TimePoint)*1.2)
  times <- seq(0, maxtime, length=101)
  fff <- stats::fitted(x, TimePoint=times)
  d <- x$data
  if(! "group" %in% names(d)) { d[["group"]] <- 1 }
  group <- 1 ## also have a global group to prevent note in check
  nn <- rownames(x$par)
  names(nn) <- nn
  ds <- lapply(nn, function(n) {
    di <- data.frame(var=n, d[,c("TimePoint",n,"group")])
    names(di)[3] <- "val"
    di
  })
  ddd <- do.call(rbind, ds)
  lattice::xyplot(val~TimePoint|var, ylim=ylim, ylab=ylab, ..., 
       group=group, data=ddd, type=c("p","a"),
       col="black",
       col.line="gray",
       fun=function(x) {mean(x, na.rm=TRUE)},
       panel=function(x,y,subscripts=NULL,groups=NULL,...) {
         lattice::panel.superpose(x,y,subscripts=subscripts,groups=groups,...)
         var <- as.character(ddd$var[subscripts[1]])
         lattice::panel.lines(fff$TimePoint, fff[[var]], col="black")
       }
      )
}

#' @export
print.pepcurve<-function(x,...) {
  print(x$par)
}

nacurve <- function(x, fillNA) {
  if(!fillNA) {
    NULL
  } else {
    x <- getparmatrix(x)
    nn <- c("pi","r","alpha")
    nn <- nn[sapply(nn, function(v) { v %in% colnames(x) && any(!is.na(x[[v]])) })]
    if(length(nn)==0) { nn <- "NA" }
    nNA <- rep(NA, length(nn))
    structure(list(par=data.frame(intercept=nNA, limit=nNA, halflife=nNA, row.names=nn),
                   data=x),
              class = c("pepcurve", "list"))
  }
}

#' Fit an exponential curve to fitted pi, r, and alpha values
#'
#' Fit an exponential curve to fitted pi, r, and alpha values
#'
#' @param x     A pepfit object, or a data frame with columns TimePoint, and at
#'    least one of pi, r, and alpha
#' @param \dots used to pass nab to the appropriate function
#' @return   A list of class "pepcurve".
#' 
#'  The first element has four values. The first three are from the pi fit:
#' the fitted value at day 0, the lower limit, and the halflife.  The last
#' is the halflife from the r fit, for which the fitted value at day 0 is
#' set to 1 and the lower limit is set to 0.
#' 
#' Remaining elements contain the pi and r values from the
#' pepfit model and the fitted curve objects for pi and r, I think.  This
#' should be checked to make sure it's still true.
#' @examples
#' data(isodata)
#' a0 <- pepfit(TimePoint, RelAb, Channel, data=isodata, Elements = list(N = 12,
#'      C = 45, H = 73, O = 15))
#'  
#' c0 <- pepcurve(a0)
#' c0
#' plot(c0)
#' @export
pepcurve <- function(x, ...) {
  UseMethod("pepcurve")
}

#' @rdname pepcurve
#' @param nab   for pepcurve.default, sets the natural abundance for use in either
#'  intercept or the limit of pi; NA finds the best fit instead
#'  for pepcurve.pepfit, can be numeric, or TRUE to use the known natural
#'  abundance from the pepfit, or FALSE to find the best fit
#' @export
pepcurve.pepfit <- function(x, nab=TRUE, fillNA=FALSE, ...) {
  if(length(unique(x$data$TimePoint))==1) {return(nacurve(x, fillNA))}
  if(isTRUE(nab)) {
    nab <- x$Element$nab
  }
  d <- parmatrix(x, default=FALSE, as.df=TRUE)
  pepcurve.default(d, nab=nab, fillNA=fillNA, ...)
}

#' @rdname pepcurve
#' @param fillNA should output have NA?
#' @param intercept set intercept values for pi, alpha, or r
#' @export
pepcurve.default <- function(x, nab=NA, fillNA=FALSE, intercept=c(pi=NA, alpha=NA, r=NA), ...) {
  if(is.null(x)) {
    stop("No data available.")
  }
  getint <- function(x, y) {
    ok <- !is.na(x) & !is.na(y)
    ifelse(stats::coef(stats::lsfit(x[ok], y[ok]))[2]>0, 0, 1)
  }
  out <- list()
  pi0 <- intercept[["pi"]]
  alpha0 <- intercept[["alpha"]]
  r0 <- intercept[["r"]]
  if("pi" %in% names(x) && length(unique(pis <- x$pi[!is.na(x$pi)]))>1) {
    if(is.na(pi0)) {pi0 <- getint(x$TimePoint, x$pi)}
    if(pi0 == 0) {
      # incorporation
      out$pi <- fitexpcurve(x$TimePoint, x$pi, intercept=nab, limit=1)
    } else {
      # dilution
      out$pi <- fitexpcurve(x$TimePoint, x$pi, limit=nab)
    }
  }
  if("alpha" %in% names(x) && length(unique(alphas <- x$alpha[!is.na(x$alpha)]))>1) {
    if(is.na(alpha0)) {alpha0 <- getint(x$TimePoint, x$alpha) }
    out$alpha <- fitexpcurve(x$TimePoint, x$alpha, intercept=alpha0, limit=1-alpha0)
  }
  if("r" %in% names(x) && length(unique(rs <- x$r[!is.na(x$r)]))>1) {
    if(is.na(r0)) {r0 <- getint(x$TimePoint, x$r) }
    out$r <- fitexpcurve(x$TimePoint, x$r, intercept=r0, limit=1-r0)
  }
  out <- do.call(rbind, out)
  out <- list(par=out, data=x)
  out <- structure(out, class=c("pepcurve","list") )
  if(is.null(out$par)) {
    nacurve(x, fillNA=fillNA)
  } else {
    out
  }
}

getparmatrix <- function(x, ...) {
  UseMethod("getparmatrix")
}
getparmatrix.pepfit <- function(x, ...) {
   parmatrix(x, default=FALSE, as.df=TRUE)
}
getparmatrix.default <- function(x, ...) {x}
