#' Relative Abundance from Counts
#'
#' Gets relative abundance of each channel from counts taken at several retention times
#'
#' The "lm", "rlm", "lqs", and "rq" methods plots the counts for each channel against the
#' baseline channel and computes a regression; the slope is the relative
#' abundance.
#' 
#' The "sum" method simply adds the counts and gets the proportion in
#' each channel.
#' 
#' The "log" method takes the logs of the counts and computes a two-way
#' ANOVA on RT and Channel; the coefficients for the channels are then
#' converted back to relative abundances.
#'  
#' Plots are created using the lattice library.
#'
#' @name relAbFromCounts
#' @param Count     Observed counts
#' @param Channel     Channels corresponding to the observed counts
#' @param RT     Retention times corresponding to the observed counts
#' @param data     An optional data frame to take the preceding variables from.
#' @param TimePoint         The time point that this data is from, to be optionally added to the
#'    output data frames.
#' @param method Options are "lm", "rlm", "lqs", "rq", "sum", and "log".
#' @param norm_channel     The channel to use as baseline for the regression methods, in the
#'    output of the relative abundance, and when computing the threshold.
#'    By default, the channel with the highest total counts is used.
#' @param threshold     A proportion used to remove values corresponding to retention times
#'    for which the value of the norm_channel is less than this threshold
#'    times the max of the norm_channel.
#' @param cutoff     All counts less than or equal to this are set to missing.
#' @param nboot Experimental, for bootstrapping
#' @return   A list is returned with the following elements.  It is of class
#'  "regRelAb" if a regression method is used and "logRelAb" if a log
#'  model is used.
#'  \item{data.wide}{A table of counts, by Channel, and RT}
#'  \item{data.long}{The counts in long data frame format}
#'  \item{data}{The fitted relative abundance and proportions at each channel.}
#'  \item{norm_channel}{The chosen baseline channel}
#'  \item{coefs}{For the regression model, the fitted coefficients}
#'  \item{method}{The method used}
#'  \item{threshold}{The threshold used}
#' @examples
#' data(isocounts)
NULL

#' @rdname relAbFromCounts
#' @export
relAbFromCounts <- function(Count, Channel, RT, data, TimePoint,
                            method=c("lm","rlm", "lqs", "rq", "sum","log"),
                            norm_channel,
                            threshold,
                            cutoff,
                            nboot=0) {
  method <- match.arg(method)

  if (method=="rq" && !requireNamespace("quantreg", quietly = TRUE)) {
    stop("quantreg package needed to specify rq method. Please install it.",
      call. = FALSE)
  }

  if (method %in% c("rlm", "lqs") && !requireNamespace("MASS", quietly = TRUE)) {
    stop("quantreg package needed to specify ", method, " method. Please install it.",
      call. = FALSE)
  }

  if (missing(data) || !is.data.frame(data) || nrow(data) <= 0) {
        data <- NULL
  }
  if(missing(Count)) Count <- as.name('Count')
  if(missing(Channel)) Channel <- as.name('Channel')
  if(missing(RT)) RT <- as.name('RT')
  df.call <- call("data.frame",
                  Count = substitute(Count),
                  Channel = substitute(Channel),
                  RT = substitute(RT))
  d.long <- eval(df.call, data, parent.frame())
  d.long$use <- 1

  RTs <- sort(unique(d.long$RT))
  d.wide <- split(d.long[,c("Count","RT")], d.long$Channel)
  d.wide <- do.call(cbind, lapply(d.wide, function(x) x$Count[match(RTs, x$RT)]))
  rownames(d.wide) <- RTs

  if(nrow(d.wide)==1 & method!="sum") {
    warning("only one RT value found; method switched to `sum` method")
    method <- "sum"
  }
  
  if(missing(norm_channel)) {
    norm_channel <- which.max(colSums(d.wide, na.rm=TRUE))
  } else if(! norm_channel %in% colnames(d.wide)) {
    norm_channel <- which.max(colSums(d.wide, na.rm=TRUE))
    warning("Specified norm_channel not in data set; setting to max.")
  } else {
    norm_channel <- which(colnames(d.wide)==norm_channel)
  }

  if (!missing(threshold) && !is.na(threshold) && threshold > 0) {
    intensity_limit <- threshold * max(d.wide[,norm_channel], na.rm=TRUE)
    skip <- d.wide[,norm_channel] < intensity_limit
    d.wide <- d.wide[!skip,]
    d.long$use <- (d.long$RT %in% rownames(d.wide))*1
  } else {
    threshold <- NA
  }
  if(!missing(cutoff) && !is.na(cutoff)) {
    d.wide[d.wide<=cutoff] <- NA
    d.long$use[d.long$Count<=cutoff] <- NA
  } else {
    cutoff <- NA
  }
  d.base <- subset(d.long, Channel==colnames(d.wide)[norm_channel])[,c("Count", "RT")]
  names(d.base)[1] <- "BaseCount"
  d.long <- merge(d.long, d.base, all=TRUE)
  d.long <- d.long[order(d.long$Channel, d.long$RT),c("Channel","RT","Count","BaseCount","use")]
  rownames(d.long) <- 1:nrow(d.long)
  Rsq <- rep(NA, ncol(d.wide))
  relab.est <- relab.se <- NA
  if(method %in% c("lm","rlm", "lqs", "rq")) {
    regression.function <-
      if(method=="lm") { function(x,y) stats::lm(y~x) }
      else if(method=="rlm") { function(x,y) MASS::rlm(y~x, maxit=50) }
      else if(method=="lqs") { function(x,y) MASS::lqs(y~x,method='lms') }
      else if(method=="rq") { function(x,y) quantreg::rq(y~x, tau=0.5) }
    
    channels <- 1:ncol(d.wide)
    names(channels) <- colnames(d.wide)
    coefs <- se <- matrix(nrow=length(channels), ncol=2, dimnames=list(names(channels), c("intercept", "slope")))
    if(nboot>0) {
        boot <- matrix(nrow=length(channels), ncol=nboot, dimnames=list(names(channels), paste0("boot", 1:nboot)))
    } else {
        boot <- NULL
    }
    min.relab.se <- min(d.wide[d.wide>0])/max(d.wide[,norm_channel])
    for(idx in channels) {
      if(norm_channel==idx) {
        coefs[idx,] <- c(0, 1)
        se[idx,] <- c(0, 0)
        Rsq[idx] <- 1
      } else {
        x <- d.wide[,norm_channel]
        y <- d.wide[,idx]
        ok <- !is.na(x) & !is.na(y)
        if(sum(ok) >=2) {
          x <- x[ok]
          y <- y[ok]
          m <- regression.function(x,y)
          coefs[idx,] <- stats::coef(m)
          se[idx,] <- stats::coef(summary(m))[,2]
          Rsq[idx] <- summary(m)$r.squared
        } else {
          coefs[idx,] <- se[idx,] <- c(NA, NA)
        }
      }
      if(nboot>0) {
          boot[idx,] <- stats::rnorm(nboot, mean=coefs[idx,2], sd=se[idx,2])
      }
    }
    if(nboot>0) {
        boot[is.na(boot)] <- 0
        boot[boot<0] <- 0
    }
    relab.obs <- as.numeric(coefs[,2])
    relab.se <- as.numeric(se[,2])
    relab.se <- pmax(relab.se, min.relab.se)
    relab <- relab.obs
    relab[is.na(relab)] <- 0
    relab[relab<0] <- 0
    out.class <- "regRelAb"
  } else if(method=="sum") {
    relab <- as.numeric(colSums(d.wide, na.rm=TRUE))
    relab <- relab/sum(relab)
    relab <- relab/relab[norm_channel]
    relab[is.na(relab)] <- 0
    coefs <- cbind(intercept=0, slope=relab)
    rownames(coefs) <- colnames(d.wide)
    out.class <- "regRelAb"
  } else if(method=="log") {
    d.long$logCount <- log(d.long$Count)
    dx <- within(d.long[d.long$use>0,], {
      RT <- factor(RT)
      Channel <- factor(Channel)
    })
    m <- stats::lm(logCount ~ RT + Channel, data=dx)
    coefs <- stats::coef(m)
    d.long$predicted[d.long$use>0] <- stats::predict(m, newdata=dx)
    relab <- as.numeric(c(1,exp(coefs[grepl("Channel", names(coefs))])))
    relab <- relab/sum(relab)
    relab <- relab/relab[norm_channel]
    out.class <- "logRelAb"
  }
  out <- data.frame(Channel=as.numeric(colnames(d.wide)),
                    RelAb=relab,
                    RelAb.obs=relab.obs,
                    RelAb.se=relab.se,
                    proportion=relab/sum(relab),
                    row.names=colnames(d.wide))
  if(!missing(TimePoint)) {
    out <- data.frame(TimePoint=TimePoint, out)
    d.long <- data.frame(TimePoint=TimePoint, d.long)
  }
  structure(list(data.wide=d.wide, data.long=d.long,
                 data=out,
                 norm_channel=out$Channel[norm_channel],
                 coefs=coefs,
                 se=se,
                 boot=boot,
                 Rsq=Rsq,
                 method=method, threshold=threshold),
            class=c(out.class,"RelAb","list"))
} 

#' @export
plot.regRelAb <- function(x, ...) {
  obj <- x
  d.long <- obj$data.long
  d.long$Channel <- factor(d.long$Channel)
  lattice::xyplot(Count~BaseCount|Channel, data=d.long, as.table=TRUE,
         xlab=paste("Base Count, Channel",colnames(obj$data.wide)[obj$base_channel]),
         panel=function(x,y,subscripts=NULL) {
           if(!all(is.na(y))) {
             lattice::panel.points(x,y, pch=c(4,1)[d.long$use[subscripts]+1])
             chname <- d.long$Channel[subscripts[1]]
             chrow <- which(rownames(obj$coefs)==chname)
             lattice::panel.abline(obj$coefs[chrow,])
           } 
         })
}

#' @export
plot.logRelAb <- function(x, ...) {
  d.long <- x$data.long
  d.long$Channel <- factor(d.long$Channel)
  lattice::xyplot(logCount ~ RT | Channel, data=d.long, as.table=TRUE,
         panel=function(x,y,subscripts=NULL) {
           lattice::panel.points(x,y,pch=c(4,1)[d.long$use[subscripts]+1])
           lattice::panel.lines(x,d.long$predicted[subscripts])
         })
}
