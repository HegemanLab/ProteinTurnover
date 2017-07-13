sim.fit <- function(fit, N=50) {
    x <- stats::fitted(fit)
    tt <- x$total
    pp <- unname(x$proportion)
    tp <- (tt+pp)/2
    tpq <- sqrt(tp*(1-tp))
    qq3 <- (tt-pp)/tpq
    qs <- split(qq3, x$TimePoint)
    out <- as.data.frame(t(sapply(qs, function(qi) {
        ## r <- as.vector(ar.ols(qi, order.max=1, demean=FALSE, aic=FALSE)$ar)
        r <- NA
        c(sd=sqrt(mean(qi^2)), r=r, n=length(qi))
    })))    
    xx <- replicate(N, {
        xx <- unlist(lapply(1:nrow(out), function(i) {
            ## as.vector(arima.sim(list(ar=out$r[i]), out$n[i])*out$sd[i])
            stats::rnorm(out$n[i], sd=out$sd[i])
        }))
        pmax(pmin(tt+xx*tpq,1),0)
    })
    xx
}

boot.fit <- function(seq, relAb, fit, N, M, alpha, isotope.method) {
    xx <- sim.fit(fit, N)
    relAbtmp <- relAb
    parboot <- do.call(rbind, lapply(seq_len(ncol(xx)), function(i) {
        relAbtmp$data$RelAb <- xx[,i]
        fiti <- makeFit(seq=seq, relAb=relAbtmp, M=M , alpha=alpha,
                        isotope.method=isotope.method, se=FALSE)
        fiti$par
    }))
    rbind(fit=fit$par, boot.fit=colMeans(parboot), boot.se=apply(parboot, 2, stats::sd))
}

## plot(tt, type="l")
## points(pp, type="h")
## points(xx[,i])
## for(i in 1:ncol(xx)) points(xx[,i])
