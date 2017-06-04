
logit<-function(p) {log(p/(1-p))}
ilogit<-function(o) {if(is.null(o)) {1} else {exp(o)/(1+exp(o))}}


topar<-function(p,days) {
  if(inherits(p, "function")) {
    p(days)
  } else {
    ndays<-length(days)
    if(length(p)==1) {
      rep(p, ndays)
    } else if(length(p)==2) {
      c(p[1], rep(p[2],ndays-1))[1:ndays]
    } else if(length(p)==ndays) {
      p
    } else {
      stop("length of p not 1, 2, or the number of days")
    }
  }
}

#' Create a function that will set up a parmatrix
#'
#'   Creates a function that takes a simple vector of parameters and creates a
#'   matrix, repeating parameters and substituting default values as required.
#'
#'   Suppose a function takes a matrix of predictors, and you want to
#'  optimize those predictors, but with some of them always set equal to
#'  each other and others always set to some default value.  Since the optimize function
#'  uses not the full matrix as input and output but instead, a simple
#'  vector of the parameters that are allowed to change.  This function
#'  creates a function that will convert this simple vector to the full
#'  matrix, repeating values as necessary and filling in default values.
#'
#' @param ix     A matrix that defines which value from the simple vector go
#' in that location.  If a default value is needed instead, the value should be
#' set to NA.
#' @param dpar   A matrix that defines the default values for locations in the ix
#'  matrix set to NA.  Values in other locations are ignored.
#' @return   A function that takes a simple vector of parameters and creates a
#'  matrix, repeating parameters and substituting default values as required.
#' @examples
#' (par.index <- matrix(c(1,1,2,NA), nrow=2))
#'  
#' (par.default <- matrix(c(NA,NA,NA,100), nrow=2))
#'  
#' par.function <- parmatrix.byix(par.index, par.default)
#' par.function(c(1.1, 1.2))
#' par.function(c(21, 22))
#' @export
parmatrix.byix <- function(ix, dpar) {
  ix<-ix
  dpar<-dpar
  function(par, default=TRUE) {
    p<-ix
    p[!is.na(ix)]<-par[ix[!is.na(ix)]]
    if(default) {
      p[is.na(ix)]<-dpar[is.na(ix)]
    } 
    p
  }
}


#' Parameter setup functions
#'
#'   Functions to assist in setting up parameters. To be used in creating
#'  new arrangements of parameters.
#'
#' @name parsetup
#' @param Day     A vector of the days that the data were taken.  It's processed using
#' unique, so duplicates can either be included or not.
#' @param type Sets default function to either parsetup.dilution or parsetup.incorporation
#' @param alpha   For parsetup.full, a list of additional parameters for alpha,
#' suitable for passing to make.par. See new.par below.  Defaults are type = "one",
#' start = 0.2, lower = 0.001, and upper = 0.999.
#' For parsetup.default, only the desired type can be set.
#' @param r   For parsetup.full, a list of additional parameters for r, suitable for
#' passing to make.par. See
#'  new.par below.  Defaults are type = "many", start = a function to give
#'  a suitable starting value for each day separately, lower = 0.001,
#'  upper = 0.999, and default0 = 0. 
#'  Cannot be changed using parsetup.default.
#' @param pi   For parsetup.full, a list of additional parameters for pi,
#' suitable for passing to make.par. See
#'  new.par below.  Defaults are type = "many", start = a function to give
#'  a suitable starting value for each day separately, lower = 0.001, and 
#'  upper = 0.999. 
#'  Cannot be changed using parsetup.default.
#' @param M   For parsetup.full, a list of additional parameters for pi,
#' suitable for passing to make.par. See
#'  new.par below.  Defaults are type = "none", start = 4, lower = 0.001, and 
#'  upper = Inf.
#'  For parsetup.default, only the desired type can be set.
#' @param noise     For parsetup.full, a list of additional parameters for noise,
#' suitable for passing to make.par. See
#'  new.par below.  Defaults are type = "many", start = a function to give
#'  a suitable starting value for each day separately, lower = 0.001, and 
#'  upper = 0.999.
#'  Cannot be changed using parsetup.default. 
#'  This parameter is currently ignored in the model fitting.
#' @param new.par   A list of additional parameters used to set up the column in the
#'  parameter matrix for each parameter.
#'  \describe{
#'    \item{type:}{
#'      Takes the values "one", "many", or "none".  The value "one" sets
#'      this parameter equal for all days.  The value "many" allows it to
#'      be different for all days.  The value "none" sets it to NA.
#'    }
#'    \item{start:}{
#'      The starting value for the optimization.  Can be set either as a
#'      numeric value or as a function which takes the day as a single parameter.
#'    }
#'    \item{default:}{
#'      The default value for the parameter, for all days (except possibly
#'      the first day, see default0).  Used when the
#'      parameter should not be optimized.  Leave missing if the parameter
#'      should be optimized over.
#'    }
#'    \item{default0:}{
#'      The default value for the parameter for the first day (usually day
#'      0).  If missing, the default as set in default is used.
#'      Set equal to null to set to NA.
#'    }
#'    \item{lower:}{
#'      The lower bound for the optimization.  Can be set either as a
#'      numeric value or as a function which takes the day as a single parameter.
#'    }		
#'    \item{upper:}{
#'      The upper bound for the optimization.  Can be set either as a
#'      numeric value or as a function which takes the day as a single parameter.
#'    }
#'  }
#' @param default.par Default values for the parameters as described in
#'   new.par
#' @param \dots Additional items to be passed along.
#' @return   A list with the following elements:
#'  \item{par}{A dataframe with the starting value, upper bound, and lower
#'  bound for each parameter to be optimized over.}
#'  \item{days}{A vector containing the unique days that were input, sorted in increasing order.}
#'  \item{parmatrix}{A function to convert a vector of parameters to a
#'    matrix suitable for passing to the fitting function.}
#' @examples
#' # With all default parameters
#' p <- parsetup.default(0:3)
#' p$parmatrix(p$par$start)
#'  
#' # With many alphas and one M 
#' p <- parsetup.default(0:3, alpha="many", M="one")
#' p$parmatrix(p$par$start)
#'  
#' # Forcing r to be 1 for all days
#' p <- parsetup.full(0:3, r=list(default=1, default0=1))
#' p$parmatrix(p$par$start)
#'
NULL

#' @rdname parsetup
#' @export
make.par<-function(Day, new.par, default.par) {
  days <- sort(unique(Day))
  p<-default.par
  if(!missing(new.par)) {
    p[names(new.par)]<-new.par
  }
  type<-match.arg(p$type, c("one", "many","none"))
  n<-length(days)
  dpar<-rep(NA, n)
  if(type=="one") {
    ix<-rep(1,n)
    ns<-rep(p$name,n)
  } else if(type=="many") {
    ix<-1:n
    ns<-paste(p$name,days,sep="")
  } else if(type=="none") {
    ix<-rep(NA,n)
    ns<-rep(NA,n)
  }
  ix0<-ix
  ns0<-ns
  if("default" %in% names(p)) {
    if(!is.null(p$default)) {
      ix<-ns<-rep(NA,n)
      dpar<-rep(p$default,n)
    }
  }
  if("default0" %in% names(p)) {
    if(is.null(p$default0)) {
      ix[1]<-ix0[1]
      ns[1]<-ns0[1]
      dpar[1]<-NA
    } else {
      ix[1]<-ns[1]<-NA;
      dpar[1]<-p$default0
    }
  }
  if(!all(is.na(ix))) {
    ix.min<-min(ix, na.rm=TRUE)
    ix <- ix-ix.min+1
  }
  list(ix=ix, dpar=dpar, start=topar(p$start,days), lower=topar(p$lower,days), upper=topar(p$upper,days), n=ns)
}

#' @rdname parsetup
#' @export
parsetup.full <- function(Day, alpha, r, pi, M, noise) {
  alpha.par<-list(name="alpha",
                  type="one",
                  start=0.2,
                  lower=0.01,
                  upper=0.99)
  r.par    <-list(name="r",
                  type="many",
                  start=function(day) {1-round((1.2+day/2)^(-1/2),2)},
                  lower=0.01,
                  upper=0.99,
                  default0=0)
  pi.par <-  list(name="pi",
                  type="many",
                  start=function(day) {round((1.2+day/2)^(-1/2),2)},
                  lower=0.01,
                  upper=0.99)
  M.par <-   list(name="M",
                  type="none",
                  start=4,
                  lower=0.01,
                  upper=Inf)
  noise.par<-list(name="noise",
                  type="none",
                  start=logit(0.01),
                  lower=-Inf,
                  upper=Inf)
  
  days  <- sort(unique(Day))
  p<-list()
  p$alpha <- make.par(days, alpha, alpha.par)
  p$r <- make.par(days, r, r.par)
  p$pi <- make.par(days, pi, pi.par)
  p$M <- make.par(days, M, M.par)
  p$noise <- make.par(days, noise, noise.par)
  ns<-names(p[[1]])
  p<-lapply(as.list(ns), function(ni) sapply(p, function(pi) pi[[ni]]))
  names(p)<-ns
  toadd<-cumsum(c(0,apply(p$ix,2,function(x) {if(all(is.na(x))) 0 else max(x, na.rm=TRUE)}))[1:ncol(p$ix)])
  p$ix <- sweep(p$ix,2,toadd,`+`)
  rownames(p$ix)<-rownames(p$dpar)<-days
  foo<-lapply(p[-c(1:2)], function(x) x[match(1:max(p$ix, na.rm=TRUE), p$ix)])
  dd <- do.call(data.frame, foo[1:3])
#  dd <-as.data.frame(do.call(rbind, lapply(foo[1:3], function(x) x)))
  rownames(dd) <- foo$n
  out<-list()
  out$par<-dd
  out$days<-days
  out$parmatrix<-parmatrix.byix(p$ix, p$dpar)
  out
}

#' @rdname parsetup
#' @export
parsetup.default <- function(Day, type=c("dilution", "incorporation", "full"), ...) {
  type <- match.arg(type)
  if(type=="dilution") {
    parsetup.dilution(Day, ...)
  } else if(type=="incorporation") {
    parsetup.incorporation(Day, ...)
  } else {
    parsetup.full(Day, ...)
  }
}


#' @rdname parsetup
#' @export
parsetup.dilution<-function(Day, alpha=c("one","many"), M=c("none","one","many"), noise=c("none","one","many")) {
  alpha<-match.arg(alpha)
  M<-match.arg(M)
  noise <- match.arg(noise)
  parsetup.full(Day, 
           alpha=list(type=alpha),
           M=list(type=M, default0=NA),
           noise=list(type=noise))
}

#' @rdname parsetup
#' @export
parsetup.incorporation <- function(Day, alpha=c("many", "log2k", "log2kplateau", "k", "kplateau"), M=c("none","one","many"), noise=c("none","one","many")) {
    M <- match.arg(M)
    alpha <- match.arg(alpha)
    noise <- match.arg(noise)
    if(alpha=="many") {
        parsetup.full(Day, 
                      alpha=list(type="many", default0=1),
                      pi=list(type="many", default0=0),
                      r=list(default=1, default0=1),
                      M=list(type=M, default0=NA),
                      noise=list(type=noise))
    } else {
        p <- parsetup.full(Day, 
                           alpha=list(type="none", default=NA),
                           pi=list(type="many", default0=0),
                           r=list(default=1, default0=1),
                           M=list(type=M, default0=NA),
                           noise=list(type=noise))
        if(alpha %in% c("kplateau", "log2kplateau")) {
            p$par <- rbind(plateau=c(0.02, 0.01, Inf), p$par)
        }
        if(alpha %in% c("k", "kplateau")) {
            p$par <- rbind(k=c(log(2)/mean(p$days), 0.01, Inf), p$par)
        } else {
            p$par <- rbind(log2k=c(log2(log(2)/mean(p$days)), -Inf, Inf), p$par)
        }
        newp <- function(days, origsetup, plateau=FALSE, logk=TRUE) {
            days <- days           
            origsetup <- origsetup
            plateau <- plateau
            logk <- logk
            function(par, default=TRUE) {
                if(plateau) {
                    plat <- par[2]
                    k <- 1:2
                } else {
                    plat <- 0
                    k <- 1
                }                
                foo <- origsetup(par[-k], default=default)
                if(logk) {
                    foo[,"alpha"] <- plat + (1-plat)*exp(-2^(par[1])*days)                    
                } else {
                    foo[,"alpha"] <- plat + (1-plat)*exp(-par[1]*days)
                }
                foo
            }
        }
        p$parmatrix <- newp(p$days, p$parmatrix,
                            plateau="plateau" %in% rownames(p$par),
                            logk="log2k" %in% rownames(p$par))
        p
    } 
}

