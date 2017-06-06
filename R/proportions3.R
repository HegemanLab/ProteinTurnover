isotopeProportion <- function(RelAb, Elements, Abundance=NULL) {
  RelAb <- RelAb/sum(RelAb)
  NatAb <- calcEnv(Elements[-1], abundance=Abundance)
  N <- Elements[[1]]
  k <- 0:(length(RelAb)-1)
  nll <- function(pi) {
    pjk <- convC(stats::dbinom(k,N,pi), NatAb)[k+1]
    ra0 <- RelAb==0
    -sum(RelAb[!ra0]*log(pjk[!ra0]))
  }
  stats::optimize(nll, c(0,1))$minimum
}

isotopeProportions <- function(d, ncol.about) {
  m <- which(colnames(d)=="M0")
  sapply(1:nrow(d), function(i) {
    El <- as.list(d[i,(ncol.about+1):(m-1)])
    ok <- El>0
    ok[1] <- TRUE
    El <- El[ok]
    x <- as.numeric(d[i,m:ncol(d)])
    x <- x[!is.na(x)]
    isotopeProportion(x, El)*100
  })
}

getLabeledProportions <- function(d, ncol.about, digits=5) {
  est <- round(isotopeProportions(d, ncol.about),digits)
  cbind(d[,1:ncol.about,drop=FALSE],Estimate=est, d[,-c(1:ncol.about),drop=FALSE])
}

#########################################
# based on the original version
#getp0 <- function(x, Elements) {
#  d <- data.frame(TimePoint=0, RelAb=x, Channel=seq_along(x)-1)
#  #with(d, {graphics::plot(Channel, RelAb); lines(Channel, RelAb)})
#  p <- pepfit(data=d, Elements=Elements, setup=parsetup.X, Abundance=ab)
#  unname(p$par)
#}
# 
#parsetup.X <- function(Day, alpha, r, pi, M, noise) {
#  alpha.par<-list(name="alpha",
#                  type="one",
#                  start=0.2,
#                  lower=0.001,
#                  upper=0.999,
#                  default0=0)
#  r.par    <-list(name="r",
#                  type="many",
#                  start=function(day) {1-round((1.2+day/2)^(-1/2),2)},
#                  lower=0.001,
#                  upper=0.999,
#                  default0=0)
#  pi.par <-  list(name="pi",
#                  type="many",
#                  start=function(day) {round((1.2+day/2)^(-1/2),2)},
#                  lower=0.0001,
#                  upper=0.9999)
#  M.par <-   list(name="M",
#                  type="none",
#                  start=4,
#                  lower=0.001,
#                  upper=Inf)
#  noise.par<-list(name="noise",
#                  type="none",
#                  start=NA,#logit(0.001),
#                  lower=-Inf,
#                  upper=Inf)
#  
#  days  <- sort(unique(Day))
#  p<-list()
#  p$alpha <- make.par(days, alpha, alpha.par)
#  p$r <- make.par(days, r, r.par)
#  p$pi <- make.par(days, pi, pi.par)
#  p$M <- make.par(days, M, M.par)
#  p$noise <- make.par(days, noise, noise.par)
#  ns<-names(p[[1]])
#  p<-lapply(as.list(ns), function(ni) do.call(cbind, lapply(p, function(pi) pi[[ni]]))) 
#  names(p)<-ns
#  toadd<-cumsum(c(0,apply(p$ix,2,function(x) {if(all(is.na(x))) 0 else max(x, na.rm=TRUE)}))[1:ncol(p$ix)])
#  p$ix <- sweep(p$ix,2,toadd,`+`)
#  rownames(p$ix)<-rownames(p$dpar)<-days
#  foo<-lapply(p[-c(1:2)], function(x) x[match(1:max(p$ix, na.rm=TRUE), p$ix)])
#  dd <-as.data.frame(do.call(cbind, lapply(foo[1:3], function(x) x)))
#  rownames(dd)<-foo$n
#  out<-list()
#  out$par<-dd
#  out$days<-days
#  out$parmatrix<-parmatrix.byix(p$ix, p$dpar)
#  out
#}
