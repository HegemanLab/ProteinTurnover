fitexpcurve <- function(x, y, intercept=NA, limit=NA, halflife=NA, halflifemax=max(x, na.rm=TRUE)*10) {
  npar <- sum(c(is.na(intercept), is.na(limit), is.na(halflife)))
  if(npar==0) {
    return(c(intercept=intercept, limit=limit, halflife=halflife))
  }
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  nx <- length(x)
  ny <- length(y)
  stopifnot(nx==ny)
  if(nx <= npar) {
    return(c(intercept=NA, limit=NA, halflife=NA))
  }
  fun <- if(is.na(intercept) & is.na(limit)) {
           function(p) findintlim(p, y)
         } else if(is.na(intercept)) {
           function(p) findint(p, y, limit)
         } else if(is.na(limit)) {
           function(p) findlim(p, y, intercept)
         } else {
           function(p) list(int=intercept, lim=limit, ss=ssil(p, y, intercept, limit))
         }
  if(is.na(halflife)) {
    fop <- function(halflife) {
      p <- exp((log(1/2) * x)/halflife)
      fun(p)$ss
    }
    o <- optimize(fop, c(0, halflifemax))
    halflife <- o$minimum
  }
  p <- exp((log(1/2) * x)/halflife)
  res <- fun(p)
  out <- c(res$int, res$lim, halflife)  
  # don't trust the original names as could have extra characters
  # from previous names
  names(out) <- c("intercept", "limit", "halflife")
  out
}

ssil <- function(p, y, intercept, limit) {
  p <- (intercept - limit) * p + limit
  sum((y-p)^2)
}

findintlim <- function(p, y) {
  m <- lsfit(p, y)
  ss <- sum(m$residuals^2)
  cf <- as.numeric(m$coef)
  int <- cf[1] + cf[2]
  lim <- cf[1]
  if(int < 0 | int > 1 | lim < 0 | lim > 1) {
    a <- rbind(i0=findlim(p, y, int=0),
               i1=findlim(p, y, int=1),
               l0=findint(p, y, lim=0),
               l1=findint(p, y, lim=1))
    a[which.min(a[,"ss"]),]                 
  } else {
    list(int=int, lim=lim, ss=ss)
  }
}
  
findlim <- function(p, y, int) {
  m <- lsfit(p-1, y-int, intercept=FALSE)
  lim <- as.numeric(int-m$coef)
  ss <- sum(m$residuals^2)
  if(lim<0 | lim>1) {
    ss.l0 <- ssil(p, y, int, 0)
    ss.l1 <- ssil(p, y, int, 1)
    if (ss.l0 < ss.l1) {
      lim <- 0
      ss <- ss.l0
    } else {
      lim <- 1
      ss <- ss.l1
    }
  }
  list(int=int, lim=lim, ss=ss)
}

findint <- function(p, y, lim) {
  m <- lsfit(p, y-lim, intercept=FALSE)
  int <- as.numeric(lim+m$coef)
  ss <- sum(m$residuals^2)
  if(int<0 | int>1) {
    ss.i0 <- ssil(p, y, 0, lim)
    ss.i1 <- ssil(p, y, 1, lim)
    if (ss.i0 < ss.i1) {
      int <- 0
      ss <- ss.i0
    } else {
      int <- 1
      ss <- ss.i1
    }
  }
  list(int=int, lim=lim, ss=ss)
}
