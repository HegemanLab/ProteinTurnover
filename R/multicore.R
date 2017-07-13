#' lapply over multiple cores with progress bar
#'
#' lapply over multiple cores with progress bar
#'
#' @param x the elements to loop over
#' @param FUN the function to apply
#' @param \dots parameters passed to FUN
#' @param mc.preschedule passed to mclapply
#' @param mc.set.seed passed to mclapply
#' @param mc.silent passed to mclapply
#' @param mc.cores passed to mclapply
#' @param mc.cleanup passed to mclapply
#' @param mc.allow.recursive passed to mclapply
#' @return the result from mclapply(x, FUN, ...)
#' @export
mclapply.progress <- function(x, FUN, ...,
                              mc.preschedule = TRUE, mc.set.seed = TRUE, 
                              mc.silent = FALSE, mc.cores = getOption("mc.cores", 1L), 
                              mc.cleanup = TRUE, mc.allow.recursive = TRUE) {
    if(mc.cores==1) return(lapply.progress(x, FUN, ...))
    t1 <- proc.time()[3]
    numjobs <- length(x)
    message("Processing ", numjobs, " items, in parallel on ", mc.cores, " cores")
    message("Started ", date())
    f <- fifo(tempfile(), open="w+b", blocking=TRUE)
    .fff <- function(i) {
        if (i==1) {
            ## Child
            pb <- txtProgressBar(0, numjobs)
            progress <- 0
            while (progress < numjobs && !isIncomplete(f)) {
                msg <- readBin(f, "double")
                progress <- progress + as.numeric(msg)
                setTxtProgressBar(pb, progress)
            }
            close(pb)
        } else {
            parallel::mclapply(x, function(xi, ...) {
                tmp <- FUN(xi, ...)
                writeBin(1, f)
                tmp
            }, mc.preschedule=mc.preschedule, mc.set.seed=mc.set.seed, mc.silent=mc.silent,
                               mc.cores=mc.cores, mc.cleanup=mc.cleanup,
                               mc.allow.recursive=mc.allow.recursive)
        }
    }
    result <- parallel::mclapply(1:2, .fff)[[2]]
    writeBin(numjobs, f)
    close(f)
    t2 <- proc.time()[3]
    message("Ended   ", date(), "; elapsed time: ", round(t2-t1), " secs")
    result
}

## library(parallel)
## fff <- function(i, n=10000, N=500) {for(k in 1:n) rnorm(N); i}
## tmp <- mclapply.progress(1:100, fff, mc.cores=4)
## tmp <- mclapply.progress2(1:100, fff, mc.cores=4)
## system.time(tmp <- mclapply(1:100, fff, mc.cores=4))
