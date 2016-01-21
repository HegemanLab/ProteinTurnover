#' Make data frame from a sequence
#' 
#' Make data frame from a sequence
#'
#' @param seq A sequence
#' @export
makeData <- function(seq) {
    k <- which(!is.na(seq$Times$data))
    do.call(rbind, lapply(k, function(i) {
        di <- seq$data[[seq$Times$data[i]]]
        chs <- seq_len(ncol(di$Count)) - 1
        rts <- di$RT                 
        data.frame(TimePoint=seq$Times$TimePoint[i],
                   Channel=rep(chs, each=length(rts)),
                   RT=rep(rts, times=length(chs)),
                   Count=as.vector(di$Count))
    }))
}

#' @rdname relAbForTimes
#' @param seq the sequence to get the relative abundance for
#' @param regression.model passed as method to relAbFromCounts
#' @export
makeRelAb <- function(seq, regression.model, nboot=0) {
    tryCatch({
        if(!is.null(seq$error)) { stop(seq$error) }
        relAbForTimes(data=makeData(seq), method=regression.model, nboot=nboot)
    }, error=identity)
}

#' @rdname pepfit
#' @param seq the sequence to get the relative abundance for
#' @param relAb the relative abundance to use
#' @param isotope.method desired method, passed as type to pepfit
#' @export
makeFit <- function(seq, relAb, isotope.method, ...) {
    tryCatch({
        if(!is.null(seq$error)) { stop(seq$error) }
        if(inherits(relAb, "error")) { stop(relAb$message) }
        pepfit(data=relAb, Elements=as.list(seq$Elements), type=isotope.method,
               name=seq$Sequence, time.unit=seq$pars$time.unit, ...)
    }, error=identity)
}

#' @rdname pepcurve
#' @param fit fit from makeFit
#' @param isotope.method used to set up intercepts correctly for incorporation experiments
#' @param plim sent to pepcurve as nab parameter
#' @export
makeCurve <- function(fit, plim=NA, isotope.method=c("other", "incorporation"), ...) {
    tryCatch({
        if(inherits(fit, "error")) { stop(fit$message) }
        ## set intercepts properly if doing incorporation
        isotope.method <- match.arg(isotope.method)
        if(isotope.method == "incorporation") { 
            int <- c(pi=0, alpha=1, r=NA)
        } else {
            int <- c(pi=NA, alpha=NA, r=NA)
        }
        pepcurve(fit, nab=plim, fillNA=TRUE, intercept=int, ...)
    }, error=identity)
}

#' @rdname pepcurve
#' @param curve a fitted curve
#' @export
getCurve <- function(curve) {
    if(inherits(curve, "error")) return(curve)
    out <- data.frame(parameter=rownames(curve$par), curve$par, stringsAsFactors=FALSE)
    rownames(out) <- NULL
    out
}

#' @rdname pepfit
#' @param fit a pepfit
#' @export
getFit <- function(fit) {
    if(inherits(fit, "error")) return(fit)
    out <- parmatrix(fit, default=FALSE, as.df=TRUE)
    keep <- c("pi","r", "alpha")
    keep <- keep[keep %in% colnames(out)]
    out <- out[,c("TimePoint", keep)]
    rownames(out) <- NULL
    out
}

#' @rdname makeSequences
#' @param seq a sequence
#' @export
getData <- function(seq) {
    x <- seq$Times
    used <- ifelse(is.na(x$data), "No", "Yes")
    x$data <- NULL
    data.frame(Used=used, x, stringsAsFactors=FALSE)
}

#' Get output from a list of completed fits
#'
#' Get output from a list of completed fits. 
#' @param seqs A list of completed fits. Each element should be a list with three elements:\itemize{
#'   \item{seq}{the original sequence information}
#'   \item{relAb}{the relative abundance}
#'   \item{fit}{the fit}
#' }
#' @param score the function to compute the score with
#' @return a data frame with the output
#' @export
getPars <- function(seqs, score=score.visual) {
    ns <- unique(unlist(lapply(seqs, function(s) { names(s$fit$par) } )))
    ns <- ns[order(gsub("[0-9]", "", ns), as.numeric(gsub("[^0-9]", "", ns)))]
    nr <- unique(unlist(lapply(seqs, function(x) rownames(x$curve$par))))
    nc <- unique(unlist(lapply(seqs, function(x) colnames(x$curve$par))))
    ns2 <- as.vector(sapply(nr, paste, nc, sep="."))    
    fits <- as.data.frame(matrix(nrow=length(seqs), ncol=2*length(ns)+length(ns2) + 4))
    names(fits) <- c("ID", "Sequence", "Score", ns, paste0("se.", ns), ns2, "Error")
    for(i in seq_along(seqs)) {
        s <- seqs[[i]]
        fits$ID[i] <- s$seq$ID
        fits$Sequence[i] <- s$seq$Sequence
        if(inherits(s$fit, "error")) {
            fits$Error[i] <- s$fit$message
        } else {
            fits$Score[i] <- score(s$fit)[[1]]
            fits[i,names(s$fit$par)] <- s$fit$par
            if(!is.null(s$fit$vcov)) {
                par.se <- sqrt(diag(s$fit$vcov))
                names(par.se) <- paste0("se.", names(par.se))
                fits[i, names(par.se)] <- par.se                
            }
            if(!is.null(s$curve$par)) {
                nr <- rownames(s$curve$par)
                nc <- colnames(s$curve$par)
                ns2 <- as.vector(sapply(nr, paste, nc, sep="."))
                fits[i, ns2] <- as.vector(t(s$curve$par))
            }
        }
    }
    fits
}

#####################################################################
makeCurvesByID <- function(out, plim=NA, isotope.method=c("other", "incorporation")) {
    isotope_method <- match.arg(isotope.method)
    fitsG <- split(seq_along(out), sapply(out, function(x) x$seq$ID))
    lapply(fitsG, function(idx) {
        x <- out[idx]
        seqs <- data.frame(Sequence=sapply(x, function(x) x$seq$Sequence),
                           Error=sapply(x, function(x) {
                               if(inherits(x$fit, "error")) x$fit$message else ""}),
                           stringsAsFactors=FALSE)
        rownames(seqs) <- NULL
        fits <- do.call(rbind, lapply(which(seqs$Error==""), function(i) {
            xi <- x[[i]]
            data.frame(Sequence=xi$seq$Sequence, getFit(xi$fit), group=i, stringsAsFactors=FALSE)
        }))
        rownames(fits) <- NULL
        ## set intercepts properly if doing incorporation
        if(isotope_method == "incorporation") { 
            int <- c(pi=0, alpha=1, r=NA)
        } else {
            int <- c(pi=NA, alpha=NA, r=NA)
        }
        ## get single curve for all sequences
        curve <- tryCatch(pepcurve(fits, nab=plim, fillNA=TRUE, intercept=int), error=identity)        
        list(ID=x[[1]]$seq$ID, seqs=seqs, fits=fits, curve=curve)
    })
}

