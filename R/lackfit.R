##############################
# code for lack of fit test
# this is a Likelihood Ratio Test between the fitted model
# and the full model with proportions equal to the observed
# so asympotically, -2loglikehood ratio is chi-squared with df equal to
# the change in the number of model parameters
# for the full model, this is the sum of the number of possible states
#   on each day minus the number of days

#' Test a pepfit
#'
#' Test a pepfit.  This is usually not run alone, as it is included in the summary of a pepfit.
#'
#' @param object a pepfit object
#' @return   A list of class pepfittest, with elements
#'  \item{deviance.fit}{deviance of the fitted model}
#'  \item{deviance.null}{deviance of the null model}
#'  \item{X2}{Pearson X^2 of the fitted model}
#'  \item{nobs}{number of observations in data set}
#'  \item{df.fit}{degrees of freedom for the fitted model}
#'  \item{df.null}{degrees of freedom for the null model}
#'  \item{sigma2}{the quasi-multinomial scaling paramter}
#'  \item{R2}{the Nagelkerke R^2}
#'  \item{data}{observed, expected, and deviance values for each observation}
#' @seealso   \code{\link{pepfit}}
#' @examples
#' data(isodata)
#' a0 <- pepfit(data=isodata, Elements = list(N = 12,
#'      C = 45, H = 73, O = 15))
#' pepfittest(a0)
#' @export
pepfittest<-function(object) {
  days <- sort(unique(object$Day))
  ndays <- length(days)
  data <- NULL
  for( j in seq_along(days)) {
    day <- days[j]
    sel <- object$Day==day
    x <- object$Count[sel]
    y <- object$RelAb[sel]    
    obs <- y/sum(y)
    exp <- getFitFor(object, j, x)[,"total"]
    null <- 1/length(x)
    data <- rbind(data, data.frame(Day=day, x=x, y=y, obs=obs, exp=exp, null=null))
  }
  data <- within(data, {
    deviance.full <- -2*ifelse(obs==0, 0, obs*log(obs))
    deviance.fit  <- -2*obs*log(exp) - deviance.full
    deviance.null <- -2*obs*log(null) - deviance.full
    Pearson.r <- (obs - exp)/sqrt(exp)
    X2 <- Pearson.r^2
    rm(deviance.full)
  })
  nobs <- nrow(data)
  out <- as.list(colSums(data[c("deviance.fit", "deviance.null", "X2")]))
  out <- c(out, list(nobs=nobs, 
                     df.fit = nobs - length(object$par),
                     df.null = nobs - ndays))
  out$sigma2 <- with(out, X2/df.fit)
  out$R2 <- with(out, (1-exp((deviance.fit-deviance.null)/nobs))/(1-exp(-deviance.null/nobs)))
  out$data <- data
  class(out)<-c("pepfittest", class(out))
  out
}

#' @rdname pepfittest
#' @param x the pepfittest object
#' @param sig.digits desired number of significant digits
#' @param \dots to match the generic print function
#' @export
print.pepfittest<-function(x, sig.digits=5, ...) {
  with(x, {
    cat("Fitted Deviance: ", signif(deviance.fit, sig.digits), ", with ", df.fit, " df\n", sep="")
    cat("  Null Deviance: ", signif(deviance.null, sig.digits), ", with ", df.null, " df\n", sep="")
    cat(sprintf("X2 = %s, sigma2 = %s, R2 = %s\n", 
            signif(X2, sig.digits), signif(sigma2, sig.digits), signif(R2, sig.digits)))
  })
}

############################################
#' Get a score for a given fit
#' 
#' Get a score for a given fit
#'
#'   The score is a value that in some way measures how
#'  well the model fits the actual data.
#' 
#'  score.N returns simply the sample size associated with the given
#'  significance level
#' 
#'  score.dev converts this sample size to a value between 0 and 100.
#'  1000 corresponds to roughly 80, 5000 to roughly 90, 20000 to roughly
#'  95.
#' 
#'  score.visual is the ratio of the total difference between the fitted
#'  and observed values to the total observed values.
#'
#' @name score
#' @param obj a pepfit object
#' @param level   significance level to use to get associated sample size
#' @return  \item{overall}{Overall Score}
#'  \item{byday}{Score for each day separately}
#' @seealso   \code{\link{pepfit}}
#' @examples
#' data(isodata)
#' a0 <- pepfit(data=isodata, Elements = list(N = 12,
#'      C = 45, H = 73, O = 15))
#' score.N(a0)
#' score.dev(a0)
#' score.visual(a0)
NULL

#' @rdname score
#' @export
score.visual <- function(obj) {
  days <- sort(unique(obj$Day))
  d <- NULL
  for( j in seq_along(days)) {
    day <- days[j]
    sel <- obj$Day==day
    x <- obj$Count[sel]
    y <- obj$RelAb[sel]    
    obs <- y/sum(y)
    exp <- getFitFor(obj, j, x)[,"total"]
    d <- rbind(d, data.frame(Day=day, obs=obs, exp=exp))
  }
  di<-sapply(split(abs(d$exp-d$obs), d$Day), sum)
  to<-sapply(split(d$obs,d$Day), sum)
  byday <- pmax(100-di/to*100,0)
  out<-list(overall=mean(byday), time=byday)
  class(out)<-c("pepscore", class(out))
  out
}

#' @rdname score
#' @export
score.dev<-function(obj, level=0.05) {
  obj<-pepfittest(obj)
  ndays <- length(unique(obj$data$Day))
  deviance.byday <- with(obj$data, tapply(deviance.fit, Day, sum))
  h2<-function(NN, b=7, c=0.5) {100*exp(-b/NN^c)}
  N1 <- qchisq(level, obj$df.fit, lower.tail=FALSE)/obj$deviance.fit
  N2 <- qchisq(level, obj$df.fit/ndays, lower.tail=FALSE)/deviance.byday
  out<-list(overall=h2(N1), time=h2(N2))
  class(out)<-c("pepscore",class(out))
  out
}

#' @rdname score
#' @export
score.N <- function (obj, level=0.05) {
    obj<-pepfittest(obj)
    ndays <- length(unique(obj$data$Day))
    deviance.byday <- with(obj$data, tapply(deviance.fit, Day, sum))
    N1 <- qchisq(level, obj$df.fit, lower.tail = FALSE)/obj$deviance.fit
    N2 <- qchisq(level, obj$df.fit/ndays, lower.tail = FALSE)/deviance.byday
    out <- list(overall = N1, time = N2)
    class(out) <- c("pepscore", class(out))
    out
}

#' @export
print.pepscore <- function(x, ...) {
  score <- x
  samelength<-function(x) {
    x<-as.character(x)
    n<-nchar(x)
    a<-sapply(max(n)-n, function(i) paste(rep(" ", i), collapse=""))
    paste(a, x, sep="")
  }
  bd <- score$time
  names(bd) <- paste(" Time ", samelength(names(bd)), ":", sep="")
  scs <- c(score$overall, bd)
  names(scs)[1]  <- "Overall:"
  scs <- format(round(scs,1),nsmall=1)
  scs <- paste(names(scs), scs)
  cat(paste(scs,collapse="\n"), "\n", sep="")
  }
#########################################

#' Create a set of simulated data
#'
#' Create a set of simulated data, to check if the change in deviance
#' between two models is really chi-squared.
#'
#' @param object a pepfit object
#' @param Nobs the sample size to be used for the simulation
#' @param Nrep how many simulations to run
#' @param x a pepsim object
#' @param n how many simulations to plot
#' @param \dots currently ignored
#' @return   A list of class pepsim containing:
#'  \item{sim}{a data frame of simulated relative abundances for each day
#'    and count value}
#'  \item{Nobs}{Nobs, as input}
#'  \item{Nrep}{Nrep, as input}
#'  \item{parent}{the original pepfit object}
#' @examples
#' data(isodata)
#' a0 <- pepfit(data=isodata, Elements = list(N = 12,
#'     C = 45, H = 73, O = 15))
#' a0.sim <- pepsim(a0, Nrep = 50)
#' plot(plot(a0.sim))
#' @export
pepsim<-function(object, Nobs=1/pepfittest(object)$sigma2, Nrep=10) {
  Nobs <- floor(Nobs)
  days <- sort(unique(object$data$TimePoint))
  out <- NULL
  for(j in seq_along(days)) {
    day <- days[j]
    sel <- object$data$TimePoint==day
    x <- object$data$Channel[sel]
    fit <- getFitFor(object, j, x)[,"total"]
    fit.p <- fit/sum(fit)
    foo<-matrix(sample(x, Nobs*Nrep, replace=TRUE, prob=fit.p), nrow=Nobs, ncol=Nrep)
    sim.count<-apply(foo, 2, function(a) xtabs(~factor(a, labels=x, levels=x)))
    out<-rbind(out, data.frame(TimePoint=day, Channel=x, RelAb=sim.count))
  }
  out<-list(sim=out, Nobs=Nobs, Nrep=Nrep, parent=object)
  class(out)<-c("pepsim", class(out))
  out
}

#' @rdname pepsim
#' @export
print.pepsim<-function(x, ...){
  with(x, {
    cat(Nrep, "simulations of", Nobs, "observations each.\n")
  })
}

#' @rdname pepsim
#' @export
plot.pepsim<-function(x, n=30, ...) {
  object <- x
  n <- min(n, object$Nrep)
  newdata<-object$sim[1:(n+2)]
  obs <- object$parent$data[,c("TimePoint", "Channel", "proportion")]
  obs$proportion <- obs$proportion * object$Nobs
  newdata <- merge(newdata, obs)
  newdata <- newdata[order(newdata$TimePoint, newdata$Channel),]
  f<-as.formula(paste(paste(paste("RelAb",1:n, sep="."), collapse="+"),"+ proportion ~ Channel|factor(TimePoint)"))
  usecol <- c(rep("grey", n), "red")
  lattice::xyplot(f, data=newdata, type="l", ylab="Abundance",
         as.table=TRUE,
         par.settings=list(superpose.line=list(col=usecol)))
}

#' Test a set of simulated data
#'
#' Fit a separate model to each simulation and find the deviance, to see if
#' the distribution of the deviance is really chi-squared.
#'
#' @param object a pepsim object
#' @param ask   Since this may take a while, it will ask for permission before
#'   starting.  Disable with ask=FALSE.
#' @param n   The number of simulations to fit models to.
#' @param model A pepfit object with the model that you want to fit the simulated data to.
#' @return   A list of class "pepsimtest" containing:
#'  \item{sim.statistic}{The deviances from the fitted models for each
#'    simulation}
#'  \item{obs.test}{The test of the original fitted model}
#'  \item{sim}{The original pepsim object}
#' @seealso  \code{\link{pepsim}}
#' @examples
#' \dontrun{
#' data(isodata)
#' a0 <- pepfit(data=isodata, Elements = list(N = 12,
#'     C = 45, H = 73, O = 15)))
#' a0.sim <- pepsim(a0, Nrep = 50)
#' a0.test.sim <- pepsimtest(a0.sim, n = 10, ask = FALSE)
#' a0.test.sim
#' plot(a0.test.sim)
#' }
#' @export
pepsimtest<-function(object, ask=TRUE, n=object$Nrep, model=object$parent) {
  n <- min(n, object$Nrep)
  if(ask) {
    if (!is.null(object$parent$time)) {
      cat("Expected Time for", n, "simulations:", model$time[3]*n,"sec\n")
    } else {
      cat("This may take a very long; there are", n, "simulations to run.\n")
    }
    ans<-readline("Do you want to proceed? [yes/no] ")
    if(ans!="yes") {
      cat("Simulation Test Aborted.\n")
      return(invisible(NULL))
    }
  }
  outm <- do.call(rbind, lapply(1:n, function(i) {
    cat("starting", i, "of", n, "\n")
    s1<-pepfit(TimePoint=object$sim$TimePoint,
               RelAb=object$sim[,i+2],
               Channel=object$sim$Channel,
               Elements=object$parent$Elements,
               Abundance=object$parent$Abundance,
               setup=function(...) {model$p} )
    unlist(pepfittest(s1)[1:8])
  }))
  out <- list(sim.test=as.data.frame(outm), 
              obs.test=as.data.frame(pepfittest(object$parent)[1:8]), 
              sim.data=object)
  class(out) <- c("pepsimtest", class(out))
  out
}

#' @rdname pepsimtest
#' @param \dots to match the generic print function
#' @export
print.pepsimtest<-function(x, ...) {
  object <- x
  cat("mean simulated deviance:", mean(object$sim.test$deviance.fit),"\n")
  cat("observed deviance:", object$obs.test$deviance.fit, "\n")
  ts <- mean(object$sim.test$deviance.fit > object$obs.test$deviance.fit)
  if(ts==0) {ts<-paste("<", 1/length(object$sim.statistic))}
  cat("simulated p-value:", ts, "\n")
}

#' @rdname pepsimtest
#' @param x the pepsimtest object
#' @param ylim the vertical limits of the plot
#' @export
plot.pepsimtest<-function(x, ylim, ...) {
  object <- x
  sim <- object$sim.test$deviance.fit/object$sim.test$sigma2
  obs <- object$obs.test$deviance.fit/object$obs.test$sigma2
  df.fit <- object$obs.test$df.fit
  df.sim <- mean(sim)
  xmax <- max(c(sim, obs)*1.1)
  xx <- 0.001
  xrange <- range(sim*1.05,
                  qchisq(c(xx,1-xx), df.fit),
                  qchisq(c(xx,1-xx), df.sim))

  d <- density(sim, from=xrange[1], to=xrange[2])
  dd <- data.frame(x=d$x, density=d$y)
  dd$red <- dchisq(dd$x, df.fit)
  dd$red2 <- dchisq(dd$x, df.sim)

  if(missing(ylim))
    ylim <- c(0, max(dd$density, dd$red, dd$red2)*1.05)
  
  lattice::xyplot(density~x, data=dd, ylim=ylim, type="l", 
         panel=function(x,y, subscripts=NULL) {
           lattice::panel.lines(x,y, col="red")
           lattice::panel.lines(x,dd$red[subscripts], col="black")
           lattice::panel.lines(x,dd$red2[subscripts], col="red", lty=2)
           lattice::panel.abline(v=obs)
           })
}

#' Compare two fitted models
#'
#' Compare two models, fitted using different parameter options.
#'
#'  The models are compared by computing the change in deviance, but the
#' sample size is unknown, so a p-value cannot be computed without
#'   assuming a given sample size.
#'  
#'   Also, this assumes the probability distribution of deviance is
#'   chi-squared; this can be checked with the pepsim function.
#' 
#' @param pepfit1 The first fitted model
#' @param pepfit2 The second fitted model
#' @param sig.digits  Number of significant digits to report.
#' @param x A pepcompare object
#' @param \dots Ignored, here for consistency with generic print function
#' @return   A list with class "pepcompare" and the following elements:
#'  \item{deviance}{A data frame with deviance and degrees of freedom for
#'    the two models, and the difference between them.}
#'  \item{N}{The input value "N"}
#'  \item{p.value}{The p value corresponding to the input value "N"}
#'  \item{level}{The input value "level"}
#'  \item{N.value}{The sample size corresponding to the input value
#'    "level"}
#'  \item{sig.digits}{Number of significant digits to report.}
#' @examples
#'   data(isodata)
#'  a0 <- pepfit(TimePoint, RelAb, Channel, data=isodata, Elements = list(N = 12,
#'     C = 45, H = 73, O = 15))
#'  a2 <- pepfit(TimePoint, RelAb, Channel, data=isodata, Elements = list(N = 12,
#'     C = 45, H = 73, O = 15), M = "one")
#'  pepcompare(a0, a2)
#' @export
pepcompare<-function(pepfit1, pepfit2, sig.digits=5) {
  t1<-m1<-pepfittest(pepfit1)
  t2<-m2<-pepfittest(pepfit2)
  if(m1$df.fit < m2$df.fit) {
    t1<-m2
    t2<-m1
  }
  out <- getF(t1, t2)
  out$sig.digits <- sig.digits
  structure(out, class=c("pepcompare", class(out)))
}

getF <- function(pa, pb) {
  d.dev <- (pa$deviance.fit - pb$deviance.fit)
  d.df <- (pa$df.fit - pb$df.fit)
  FF <- (d.dev/d.df)/pb$sigma2
  pp <- pf(FF, d.df, pb$df.fit, lower.tail=FALSE)
  data.frame(F=FF, ndf=d.df, ddf=pb$df.fit, p.value=pp, 
             deviance1=pa$deviance.fit, df1=pa$df.fit,
             deviance2=pb$deviance.fit, df2=pb$df.fit,
             sigma2=pb$sigma2)
}

#' @rdname pepcompare
#' @export
print.pepcompare <- function(x, sig.digits=x$sig.digits, ...) {
  cat("Model 1 Deviance: ", signif(x$deviance1, sig.digits), ", with ", x$df1, " df\n", sep="")
  cat("Model 2 Deviance: ", signif(x$deviance2, sig.digits), ", with ", x$df2, " df\n", sep="")
  cat("          Sigma2: ", signif(x$sigma2, sig.digits), ", with ", x$df2, " df\n", sep="")
  cat("     F Statistic: ", signif(x$F, sig.digits), ", with ", x$ndf, "/", x$ddf, " df\n", sep="")
  cat("         p-value: ", signif(x$p.value, sig.digits), "\n", sep="")
}
