# 2/23, pep7a (renamed from pep7)
# Attached is a new version of the peptide code (pep7.R) that has two
# new features:  it can change which element is changing, and it can
# take into the natural abundance of multiple other elements.  It now
# takes an input called Elements which is a list of the elements that it
# should take account of and how many there are of each; the first in
# the list is the one that will change over time.  It knows about the
# natural abundance of N, C, O, H, and S; others can be added within the
# call, if desired.  This version of the code is also considerably
# faster than the last one I sent you.

# 2/25, pep7b
# Wen Ping: I found a minor bug. The fitting
# curves (green line) for the old peptide distribution on days 2, 6 and 12 are
# kind of off a little bit when beta-binomial model is used.
# Aaron: It was using the M from the betabinomial to get
# the lines for all three distributions instead of just for the one that
# was changing.

# 3/5, pep7c
# WenPing: Did you get chance to check why the non-linear regression code won't
# work if there were only three data points for pi? It works correctly if there
# are at least four time points.
# Aaron: It needs to have the plataeu set in order to work; otherwise the
# number of parameters is equal to the number of data points.  When
# fitting the curve to a full model, it was supposed to set this to the
# natural abundance by default, but there was a typo.  It's fixed in the
# attached version.

# 4/1, pep7d
# adds an error if Relative Abundance detected when
# calculated probability is zero
# also calls isotopes7d

# 4/24, pep7e
# add option to use FFT convolution method
# calls isotopes7e

# 7/24, pep7f
# use options to switch convolution methods rather than passing a variable
# also include faster exact function as default
#   (though still not quite as fast as FFT)
# calls isotopes7f

# 8/13, pep7g
# adds code to model random noise

# 9/3, pep7h
# adds code to perform a lack of fit test

# 10/13, pep8
# add C code for convolution

# 12/10, pep9
# revisit parameter code
###################################################
### probability code
###################################################

dbetabinom <- function(k, N, mu, M) {
  if(is.na(M)) {
    dbinom(k,N,mu)
  } else {
    alpha <- M*mu
    beta  <- M*(1-mu)
    not0 <- !(k<0 | k>N)
    usek <- k[not0]
    out <- numeric(length(k))
    out[not0]<-choose(N,usek)*exp(lbeta(alpha+usek, beta+N-usek) - lbeta(alpha, beta))
    out
  }
}


###################################################
### helper code for all days
###################################################


estCorr <- function(vcov) {
  sdi <- diag(1/sqrt(diag(vcov)))
  corr <- sdi %*% vcov %*% sdi
  rownames(corr) <- colnames(corr) <- colnames(vcov)
  corr
}

###################################################
### Day j code
###################################################

probjk <- function(k, N, nab, NatAb,
                   parmatrix, row, verbose=FALSE) {
  alpha <- parmatrix[row,1]
  r     <- parmatrix[row,2]
  pi0   <- parmatrix[1,3]
  piJ   <- parmatrix[row,3]
  M     <- parmatrix[row,4]
  M0     <- parmatrix[1,4]
  k0 <- 0:max(k)
  delta <- dbetabinom(k0, N, nab, NA)
  if(alpha==1) {
    out<-delta
    gamma0<-gammaJ<-0
  } else {
    gamma0 <- dbetabinom(k0, N, pi0, M0)
    if(r==0) {
      gammaJ <- 0
    } else {
      gammaJ <- dbetabinom(k0, N, piJ, M)
    }
    out <- alpha*delta + (1-alpha) * ((1-r)*gamma0 + r*gammaJ)    
  }
  if(verbose) {
    out <- data.frame(delta =alpha*convC(delta, NatAb, scale=FALSE)[k+1],
                      gamma0=(1-alpha)*(1-r)*convC(gamma0, NatAb, scale=FALSE)[k+1],
                      gammaJ=(1-alpha)*r*convC(gammaJ, NatAb, scale=FALSE)[k+1],
                      total = convC(out, NatAb, scale=FALSE)[k+1])
    out
  } else {
      convC(out, NatAb, scale=FALSE)[k+1]  
  }
}

getprob <- function(parmatrix, Day, Count, N, nab, NatAb) {
    days <- sort(unique(Day))
    p <- numeric(length(Count))
    for(j in seq_along(days)) {
        use <- Day == days[j]
        p[use] <- probjk(k=Count[use], N=N, nab=nab, NatAb=NatAb, parmatrix=parmatrix, row=j)
    }
    p
}

getRelAb <- function(parmatrix, Day, Count, N, nab, NatAb, norm_channel) {
    days <- sort(unique(Day))
    relab <- numeric(length(Count))
    for(j in seq_along(days)) {
        use <- Day == days[j]
        tmp <- probjk(k=Count[use], N=N, nab=nab, NatAb=NatAb, parmatrix=parmatrix, row=j)
        nc <- which(Count[use] == norm_channel[j])
        tmp <- tmp/tmp[nc]
        tmp[nc] <- NA
        relab[use] <- tmp
    }
    relab
}

fitLogLik2 <- function(par, to.parmatrix, Day, Count, N, nab, NatAb, relab.obs, relab.se, norm_channel) {
    relab.est <- getRelAb(parmatrix=to.parmatrix(par), Day=Day, Count=Count,
                          N=N, nab=nab, NatAb=NatAb, norm_channel=norm_channel)    
    sum(dnorm(relab.est, mean=relab.obs, sd=relab.se, log=TRUE), na.rm=TRUE)
    
}

LogLik1 <- function(p.est, p.obs) {
  if(any(p.obs > 0 & p.est == 0)) {
    stop("Observed proportion non-zero where estimated probability is zero.\n")
  }
  ok <- p.obs > 0
  sum(p.obs[ok] * log(p.est[ok]))
}

fitLoglik<-function(par, to.parmatrix, Day, Count, N, nab, NatAb, LogLik, ...){
  p.est <- getprob(parmatrix=to.parmatrix(par), Day=Day, Count=Count, N=N, nab=nab, NatAb=NatAb)
  LogLik(p.est, ...)
}

fitN <- function(p.obs, p.est) {
    d2 <- function(logn, fi, pi) {
        n <- exp(logn)
        ai <- fi*n + 1
        digamma(sum(ai))*sum(fi)*n - sum(fi*digamma(ai))*n + sum(n*fi*log(pi))
    }
    d3 <- function(logn, fi, pi) {
        n <- exp(logn)
        ai <- fi*n + 1
        digamma(sum(ai))*sum(fi)*n + trigamma(sum(ai))*sum(fi)^2*n^2 -
            sum(fi*digamma(ai))*n - sum(fi^2*trigamma(ai))*n^2 +
                sum(n*fi*log(pi))
    }
    logN <- uniroot(d2, log(c(10, 1e8)), fi=p.obs, pi=p.est)$root
    se <- sqrt(-1/d3(logN, p.obs, p.est))
    c(logN=logN, se=se)
}

#' Fit a model
#'
#' Fit a model of how mass changes over time of a specific element
#'
#' @aliases print.pepfit plot.pepfit fitted.pepfit
#' @param TimePoint The timepoint the observation was taken at
#' @param RelAb The relative abundance of that observation
#' @param Channel The number of heavy elements for that the observation
#' @param data An optional data frame to take the TimePoint, RelAb, and Channel variables
#'   from.  If names for any of the three previous variables are missing, they are
#'   looking for in the data frame using these specific nams. Also, can be \code{RelAbTimes}
#'   object from \code{relAbForTimes}, in which the three previous
#'   variables are ignored. 
#' @param Elements A list of elements known to be in this peptide and their quantity.  The first
#'   element is the one that will be allowed to change relative abundances
#'   over time.
#' @param Abundance A list of elements and their natural abundance.  C, H, O, N, and S
#'   need not be included, as those natural abundances are included by default.
#' @param se Boolean, set TRUE to compute standard errors using the Hessian from optim
#' @param time Boolean, set FALSE to not measure elapsed time
#' @param maxit Maximum iterations to allow optim
#' @param setup The desired setup function
#' @param time.unit For use in plots and other output
#' @param name For use in plots and other output
#' @param method Experimental
#' @param useRelAb Use relative abundance instead of proportions
#' @param \dots Additional parameters passed to the setup function
#' @return A list of S3 class pepfit, with elements
#'   \item{par}{Fitted value of desired parameters}
#'   \item{value}{Fitted likelihood}
#'   \item{counts}{count data from optim function}
#'   \item{convergence}{convergence information from optim function}
#'   \item{message}{message from optim function}
#'   \item{Day}{The TimePoint input values}
#'   \item{RelAb}{RelAb, from input}
#'   \item{Count}{Count, from input}
#'   \item{Element}{The name, count, and natural abundance of the changing element}
#'   \item{Elements}{Elements, from input}
#'   \item{Abundance}{the natural abundance of each element used}
#'   \item{NatAb}{the combined natural abundance over non-changing elements}
#'   \item{p}{the output from the setup function}
#'   \item{setup}{the setup function used}
#'   \item{parmatrix}{the values of all parameters, both fitted and specified, as in matrix form}
#'   \item{time}{Elapsed time}
#' @seealso \code{\link{relAbForTimes}}
#' @examples
#' data(isodata)
#' a <- pepfit(TimePoint, RelAb, Channel, data=isodata, Elements=list(N=12,C=45,H=73,O=15))
#' summary(a)
#' plot(a)
#'  
#' data(isoincorp)
#' b <- pepfit(TimePoint, RelAb, Channel, data=isoincorp,
#' type="incorporation", Elements=list(N=12,C=45,H=73,O=15))
#' summary(b)
#' plot(b)
#' @export

pepfit <- function(TimePoint, RelAb, Channel, data, 
                   Elements, Abundance=NULL,
                   se=FALSE,
                   time=exists("proc.time"), maxit=1000,
                   setup=parsetup.default,
                   time.unit=NA,
                   name=NULL, method=c("binomial", "regression"),
                   useRelAb=FALSE,
                   ...) {
  method <- match.arg(method)
  if(time) {start.time <- proc.time()}
  if(!missing(data) && inherits(data, "RelAbTimes")) {
    d <- data$data
  } else {
    if(missing(data) || !is.data.frame(data) || nrow(data) <= 0 ) {
      data <- NULL
    }
    if(missing(Channel)) Channel <- as.name('Channel')
    if(missing(RelAb)) RelAb <- as.name('RelAb')
    if(missing(TimePoint)) TimePoint <- as.name('TimePoint')  
    df.call <- call("data.frame",
                    TimePoint=substitute(TimePoint),
                    RelAb=substitute(RelAb),
                    Channel=substitute(Channel))
    d <- eval(df.call, data, parent.frame())
  }
  Abundance<-getabundance(Abundance)
  nameEl<-names(Elements)[1]
  if(! nameEl %in% names(Abundance)) {
    stop("No abundance given for element", nameEl, ".")
  }
  abundEl<-Abundance[[nameEl]]
  if(length(abundEl)>2) {
    stop("Only elements with the possibility of one extra abundance are allowed.")
  }
  p <- setup(d$TimePoint,...)     
  N   <- Elements[[1]]
  nab <- abundEl[2]/sum(abundEl)     
  NatAb <- if (length(Elements)>1) calcEnv(Elements[-1], abundance=Abundance) else NULL
  start.par<-p$par$start
  names(start.par)<-rownames(p$par)
  ts <- as.numeric(factor(d$TimePoint))
  d$proportion <- d$RelAb/tapply(d$RelAb, ts, sum, na.rm=TRUE)[ts]
  usecol <- if(useRelAb) {"RelAb"} else {"proportion"}
  if(method=="binomial") {
      ans1 <- optim(start.par, fitLoglik, to.parmatrix=p$parmatrix,
                    Day=d$TimePoint, Count=d$Channel, N=N, nab=nab, NatAb=NatAb,
                    LogLik=LogLik1, p.obs=d[[usecol]], 
                    method="L-BFGS-B", control=list(fnscale=-1, maxit=maxit),
                    lower=p$par$lower, upper=p$par$upper, hessian=se)
  } else if(method=="regression") {
      ans1 <- optim(start.par, fitLogLik2, to.parmatrix=p$parmatrix,
                    Day=d$TimePoint, Count=d$Channel, N=N, nab=nab, NatAb=NatAb,
                    relab.obs=d$RelAb.obs, relab.se=d$RelAb.se, norm_channel=data$norm_channel,
                    method="L-BFGS-B", control=list(fnscale=-1, maxit=maxit),
                    lower=p$par$lower, upper=p$par$upper, hessian=se)
  }
  if(se) {
      ans1$vcov <- tryCatch(solve(-ans1$hessian), error=function(e) {out <- ans1$hessian; out[] <- NA; out})
  }
  #p.est <- getprob(parmatrix=p$parmatrix(ans1$par), Day=d$TimePoint, Count=d$Channel,
  #                 N=N, nab=nab, NatAb=NatAb)
  #ans1$N <- fitN(p.obs=d$proportion, p.est=p.est)
  ans1$data <- d
  ans1$Day <- d$TimePoint
  ans1$RelAb <- d$RelAb
  ans1$Count <- d$Channel
  ans1$Element <- list(name=names(Elements)[1], n=N, nab=nab)
  ans1$Elements  <- Elements
  ans1$Abundance <- Abundance
  ans1$NatAb <- NatAb
  ans1$p<-p
  ans1$setup<-setup
  ans1$dots <- list(...)
  ans1$name <- name
  ans1$parmatrix <- p$parmatrix(ans1$par)
  ans1$TimeUnit <- time.unit
  class(ans1)<- c("pepfit", class(ans1))
  if(time) {ans1$time <- proc.time() - start.time}
  ans1
}

#' Gets the fitted parameters from a pepfit object
#'
#' Gets the fitted parameters from a pepfit object
#'
#' @param pepfit A pepfit object
#' @param default Whether or not parameters given as defaults should be included or not
#' @param as.df Whether it should return a data frame or a matrix
#' @return The fitted parameters from the pepfit object, in the desired format
#' @export
parmatrix <- function(pepfit, default=TRUE, as.df=FALSE) {
  out <- pepfit$p$parmatrix(pepfit$par, default=default)
  if(as.df) {
    out <- data.frame(TimePoint=as.numeric(rownames(out)),
                      out)
  }
  out
}

#' Summarize and print results from a pepfit
#'
#' Summarize and print results from a pepfit
#'
#' @param object A pepfit object
#' @param digits How many digits to print
#' @param \dots Additional arguments passed to \code{\link{pepfittest}}
#' @return   list of class summary.pepfit with elements
#'  \item{par}{Fitted parameters}
#'  \item{hessian}{Fitted hessian matrix, if in original pepfit}
#'  \item{time}{Elapsed time}
#'  \item{test}{the results from \code{\link{pepfittest}}}
#' @seealso  \code{\link{pepfittest}}
#' @examples
#' data(isodata)
#' a <- pepfit(data=isodata,
#' Elements=list(N=12,C=45,H=73,O=15))
#' summary(a)
#' @export
summary.pepfit <- function(object, digits=3, ...) {
  out <- list(par=object$par,
              vcov=object$vcov,
              time=object$time,
              test=pepfittest(object,...),
              digits=digits)
  structure(out, class = c("summary.pepfit", class(out)))
}

#' @export
print.pepfit <- function(x, digits=3, ...) {
  print(summary(x, digits=digits, ...))
}

#' @export
print.summary.pepfit <- function(x, digits=x$digits, ...) {
  obj <- x
  cat("\nParameter Estimates:\n")
  est<-cbind(Estimate=obj$par)
  if(!is.null(obj$vcov)) {
    est<-cbind(est, `Std. Error`=sqrt(diag(obj$vcov)))
  }
  print(est)
  if(!is.null(obj$vcov)) {
    cat("\nEstimated Correlations:\n")
    print(estCorr(obj$vcov))
  }
  cat("\n")
  print(obj$test)
  if (!is.null(obj$time)) {
    cat("\nTime Elapsed:", obj$time[3],"sec\n\n")
  }
}

###################################################
### plotting functions
###################################################

getFitFor <- function(obj,daynum,range) {
  probjk(k=range, N=obj$Element$n, nab=obj$Element$nab, NatAb=obj$NatAb,
              parmatrix=obj$parmatrix,
              row=daynum, verbose=TRUE)
}

#' @export
fitted.pepfit <- function(object, range=0:max(object$data$Channel), ...) {
  times <- sort(unique(object$data$TimePoint))
  n <- seq_along(times)
  d <- lapply(n, function(i) {
    out <- probjk(k=range, N=object$Element$n, nab=object$Element$nab,
                  NatAb=object$NatAb, parmatrix=object$parmatrix,
                  row=i, verbose=TRUE)
    data.frame(TimePoint=times[i], Channel=range, out)
  })
  out <- merge(object$data, do.call(rbind, d))
  out <- out[order(out$TimePoint, out$Channel),]
  rownames(out) <- 1:nrow(out)
  out
}

#' @export
plot.pepfit <- function(x, as.table=TRUE, main=x$name,
                        sub=paste("Isotope:", x$Element$name),
                        xlab="Isotope Channel", ylim, ...) {
  f <- fitted(x)
  if(!is.na(x$TimeUnit)) {
    tp.levels <- paste(x$TimeUnit, sort(unique(f$TimePoint)))
    tp <- paste(x$TimeUnit, f$TimePoint)
    f$TimePoint <- factor(tp, levels=tp.levels)
  } else {
    f$TimePoint <- factor(f$TimePoint)
  }
  if(missing(ylim)) {
    ymax <- max(f$proportion, f$total, na.rm=TRUE)
    ylim <- ymax*c(-0.01, 1.05)
  }
  panel.linesifnz <- function(x, y, ...) {
    if(!all(y==0)) {lattice::panel.lines(x,y, ...)}
  }
  lattice::xyplot(proportion ~ Channel|TimePoint, data=f,
         as.table=as.table, ylim=ylim,
         panel=function(x,y,subscripts=NULL) {
           lattice::panel.points(x,y,type="h", col="gray")
           panel.linesifnz(x,f$delta[subscripts], col="red")
           panel.linesifnz(x,f$gamma0[subscripts], col="green")
           panel.linesifnz(x,f$gammaJ[subscripts], col="blue")
           panel.linesifnz(x,f$total[subscripts], col="black", lty=2)
         }, main=main, sub=sub, ...)
}
