lapply.progress <- function(x, FUN, ...) {
    t1 <- proc.time()[3]
    numjobs <- length(x)
    message("Processing ", numjobs, " items, one at a time")
    message("Started ", date())
    pb <- txtProgressBar(0, numjobs)
    result <- lapply(seq_along(x), function(i, ...) {
        tmp <- FUN(x[[i]], ...)
        setTxtProgressBar(pb, i)
        tmp
    })
    close(pb)
    t2 <- proc.time()[3]
    message("Ended   ", date(), "; elapsed time: ", round(t2-t1), " secs")
    result
}
    

## mclapply.progress <- function(x, FUN, ...,
##                               mc.preschedule = TRUE, mc.set.seed = TRUE, 
##                               mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L), 
##                               mc.cleanup = TRUE, mc.allow.recursive = TRUE) {
##     if(mc.cores==1) return(lapply.progress(x, FUN, ...))
##     t1 <- proc.time()[3]
##     numjobs <- length(x)
##     message("Processing ", numjobs, " items, in parallel on ", mc.cores, " cores")
##     message("Started ", date())
##     f <- fifo(tempfile(), open="w+b", blocking=TRUE)
##     if (inherits(parallel:::mcfork(), "masterProcess")) {
##         ## Child
##         pb <- txtProgressBar(0, numjobs)
##         progress <- 0
##         while (progress < numjobs && !isIncomplete(f)) {
##             msg <- readBin(f, "double")
##             progress <- progress + as.numeric(msg)
##             setTxtProgressBar(pb, progress)
##         }
##         close(pb)
##         parallel:::mcexit()
##     }
##     result <- mclapply(x, function(xi, ...) {
##         tmp <- FUN(xi,...)
##         writeBin(1, f)
##         tmp
##     }, mc.preschedule=mc.preschedule, mc.set.seed=mc.set.seed, mc.silent=mc.silent,
##        mc.cores=mc.cores, mc.cleanup=mc.cleanup, mc.allow.recursive=mc.allow.recursive)
##     writeBin(numjobs, f)
##     close(f)
##     t2 <- proc.time()[3]
##     message("Ended   ", date(), "; elapsed time: ", round(t2-t1), " secs")
##     result
## }

prepare_dir <- function(x) {
  dir <- dirname(x)
  if(!file.exists(dir)) { dir.create(dir, recursive=TRUE) }   
}

save.table <- function(..., file) {
  prepare_dir(file)
  write.table(..., file = file)
}

# to help in saving lattice graphics
save.graphic <- function(x, file, dir, color = TRUE, device=c("png", "pdf", "svg"),
                         panel.size=1, panel.width=panel.size, panel.height=panel.size,
                         border.width=1, border.height=1,
                         res=300,
                         width=NULL, height=NULL, ...) {
    out <- tryCatch({
        device <- match.arg(device)
        if(device=="pdf") res <- 1
        if(inherits(x, "error")) stop(x$message)
        file <- paste(file, device, sep=".")
        if(!missing(dir)) f <- file.path(dir, file)
        else f <- file
        prepare_dir(f)
        x <- plot(x, ...)                        
        d <- dim(x)
        dw <- d[1]
        dh <- if(length(d)==1) { 1 } else { d[2] }
        width <- (dw*panel.width + 2*border.width)
        height <- (dh*panel.height + 2*border.height)
        if(device %in% c("pdf", "svg")) {
            lattice::trellis.device(device, color = color, file = f, width=width, height=height)
        } else if(device=="png") {
            lattice::trellis.device(device, color = color, file = f, width=width, height=height, res=res, units="in")
        }
        plot(x, layout=c(dh, dw),
             panel.width=list(x=panel.width, units="in"),
             panel.height=list(x=panel.height, units="in"))
        dev.off()
        list(file=file, width=width, height=height, device=device)

    }, error=function(e) list(error=paste("Unable to plot:", e$message)))
    class(out) <- c("img", class(out))
    out
}
