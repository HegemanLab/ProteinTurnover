#'   Relative Abundance from Counts, for multiple time points
#'
#'   Gets relative abundance of each channel from counts taken at several
#'  retention times, and over multiple time points
#'
#'   The data from each time point is fit separately, using
#'  relAbFromCounts.  Plots are created using the lattice library.
#' 
#' @param Count     Observed counts
#' @param Channel     Channels corresponding to the observed counts
#' @param RT     Retention times corresponding to the observed counts
#' @param TimePoint     Time points corresponding to the observed counts
#' @param data     An optional data frame to take the preceding variables from
#' @param nboot Experimental, for bootstrapping
#' @return   A list is returned with the following elements.  It is of class
#'  "regRelAbTimes" if a regression method is used and "logRelAbtimes" if a log
#'  model is used.
#'  \item{data.long}{The long data frames from each fit, bound together
#'    into one data frame.}
#'  \item{norm_channel}{The chosen baseline channel for each timepoint}
#'  \item{coefs}{For the regression model, a list with the fitted coefficients at each
#'  time points.}
#'  \item{relab}{The fitted relative abundance at each channel and timepoint}
#'  \item{proportion}{The fitted relative abundance, converted to a
#'    proportion}
#'  \item{method}{The method used}
#'  \item{threshold}{The threshold used}
#' @seealso   \code{\link{relAbFromCounts}}
#' @examples
#' data(isocounts)
#' @export
relAbForTimes <- function(Count, Channel, RT, TimePoint, data, nboot=0, ...) {
  if (missing(data) || !is.data.frame(data) || nrow(data) <= 0) {
        data <- NULL
  }
  if(missing(Count)) Count <- as.name('Count')
  if(missing(Channel)) Channel <- as.name('Channel')
  if(missing(RT)) RT <- as.name('RT')
  if(missing(TimePoint)) TimePoint <- as.name('TimePoint')
  df.call <- call("data.frame",
                  Count = substitute(Count),
                  Channel = substitute(Channel),
                  RT = substitute(RT),
                  TimePoint = substitute(TimePoint))
  d.all <- eval(df.call, data, parent.frame())
  
  d.byt <- split(d.all[,c("Count", "Channel", "RT")], d.all$TimePoint)
  times <- as.numeric(names(d.byt))
  out <- lapply(seq_along(d.byt), function(i)
                relAbFromCounts(data=d.byt[[i]], TimePoint=times[i], nboot=nboot, ...))
  names(out) <- times
  
  d.long <- do.call(rbind, lapply(out, function(x) x$data.long))
  rownames(d.long) <- 1:nrow(d.long)
  data <- do.call(rbind, lapply(out, function(x) x$data))
  rownames(data) <- 1:nrow(data)
  
  structure(list(data.long=d.long,
                 data=data,
                 norm_channel=sapply(out, function(x) x$norm_channel),
                 coefs=lapply(out, function(x) x$coefs),          
                 Rsq=lapply(out, function(x) x$Rsq),
                 boot=do.call(rbind, lapply(out, function(x) x$boot)),
                 method=out[[1]]$method, threshold=out[[1]]$threshold, cutoff=out[[1]]$cutoff),
            class=c(paste(class(out[[1]])[1:2],"Times",sep=""),"list") )
}

#' @rdname relAbForTimes
#' @param x A regRelAbTimes or logRelAbTimes object to be plotted
#' @param \dots Additional parameters to be passed
#' @export
regressionPlot <- function (x, ...) {
  UseMethod("regressionPlot")
}

#' @export
regressionPlot.regRelAbTimes <- function(x, ...) {
  obj <- x
  d.long <- obj$data.long
  d.long$Channel <- factor(d.long$Channel)
  d.long$TimePoint <- factor(d.long$TimePoint)
  bcn <- obj$norm_channel
  ubcn <- unique(bcn)
  if(length(ubcn)==1) {
    xlab <- paste("Base Count, Channel", ubcn)
  } else {
    xlab <- paste("Base Count, Channels ", paste(bcn, collapse=", "), ", respectively", sep="")
  }
  p <- lattice::xyplot(Count~BaseCount|Channel*TimePoint, data=d.long, as.table=TRUE,
         xlab=xlab,
         panel=function(x,y,subscripts=NULL) {
           if(!all(is.na(y))) {
             lattice::panel.points(x,y, pch=c(4,1)[d.long$use[subscripts]+1])
             timename <- d.long$TimePoint[subscripts[1]]
             chname <- d.long$Channel[subscripts[1]]
             chrow <- which(rownames(obj$coefs[[timename]])==chname)
             lattice::panel.abline(obj$coefs[[timename]][chrow,])
           } 
         }, ...)
  latticeExtra::useOuterStrips(p)
}

#' @export
regressionPlot.logRelAbTimes <- function(x, ...) {
  d.long <- x$data.long
  d.long$Channel <- factor(d.long$Channel)
  d.long$TimePoint <- factor(d.long$TimePoint)
  p <- lattice::xyplot(logCount ~ RT | Channel*TimePoint, data=d.long, as.table=TRUE,
         panel=function(x,y,subscripts=NULL) {
           lattice::panel.points(x,y,pch=c(4,1)[d.long$use[subscripts]+1])
           lattice::panel.lines(x,d.long$predicted[subscripts])
         })
  latticeExtra::useOuterStrips(p)
}

#' @rdname relAbForTimes
#' @param type desired type of plot
#' @export
plot.RelAbTimes <- function(x, ..., type=c("relAb", "regression")) {
    type <- match.arg(type)
    if(type=="regression") {
        regressionPlot(x, ...)
    } else {
        d.long <- x$data.long
        d.long$Channel <- factor(d.long$Channel)
        d.long$TimePoint <- factor(d.long$TimePoint)
        p <- lattice::xyplot(Count ~ RT | Channel*TimePoint, data=d.long, as.table=TRUE,
                    panel=function(x,y,subscripts=NULL) {
                        lattice::panel.points(x,y,pch=c(4,1)[d.long$use[subscripts]+1])
                    })
        latticeExtra::useOuterStrips(p)
    }
}
