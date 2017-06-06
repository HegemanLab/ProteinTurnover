

## Read ProteinPilot files
## Returns data frame with: ID, Sequence, RT, MW, Z, TimePoint, Conf
readProteinPilotFiles <- function(files, times, dir, verbose=TRUE) {
    stopifnot(length(files) == length(times))
    paths <- files
    if(!missing(dir)) {paths <- file.path(dir, files)}    
    if(verbose) message("reading ProteinPilot files ...")
    r <- do.call(rbind, lapply(seq_along(files), function(i) {
        foo <- readProteinPilotFile(file=paths[i])
        foo$TimePoint <- times[i]
        if(verbose) message("read ", files[i])
        foo
    }))
    rownames(r) <- NULL
    r
}

## Read single ProteinPilot file
## Returns data frame with: ID, Sequence, RT, MW, Z, Conf
readProteinPilotFile <- function(file) {
    r <- utils::read.table(file, header=TRUE, sep="\t", quote="", as.is=TRUE, check.names=FALSE)
    ## restrict to rows above confidence.level and desired columns
    r <-    r[, c("Accessions", "Sequence", "Time", "Prec MW", "Theor z", "Conf")]
    names(r) <- c("ID",         "Sequence", "RT",   "MW",      "Z",       "Confidence")
    ## convert time to seconds
    r$RT <- as.integer(round(r$RT*60))
    ## we need to deal with multiple IDs (accessions) and choose just the first one.    
    ## also, remove the last _ARATH
    ##       fix names so they can be used as filenames  
    ##       and remove the initial part of the name  
    r$ID <- sapply(strsplit(r$ID, split=";"), `[`, 1)
    r$ID <- gsub("_ARATH$", "", r$ID)
    r$ID <- gsub("[|,]", "_",   r$ID)    
    r$ID <- sub("^[^_]*_", "",  r$ID)
    r
}

#' Read Scaffold file
#' 
#' Read Scaffold file
#'
#' @param file name of Scaffold file
#' @param times data collection times, as a vector
#' @param dir directory that file is in
#' @param verbose if verbose output desired
#' @return data frame with: ID, Sequence, RT, MW, Z, TimePoint, Conf
#' @export
readScaffoldFile <- function(file, times, dir, verbose=TRUE) {
    if(verbose) message("reading Scaffold file ...")
    path <- file
    if(!missing(dir)) path <- file.path(dir, file)
    r <- utils::read.csv(path, as.is=TRUE, check.names=FALSE)
    r <- r[,c(
              "Protein accession numbers",
              "Peptide sequence",
              "Spectrum name",            
              "Actual peptide mass (AMU)",
              "Spectrum charge",
              "Biological sample category",
        "Peptide identification probability"
              )]
    names(r) <- c("ID", "Sequence", "RT", "MW", "Z", "TimePoint", "Confidence")
    r <- r[r$ID!="", ]
    r$Confidence <- as.numeric(sub("%", "", r$Confidence))
    r$MW <- as.numeric(sub(",", "", r$MW))
    r$RT <- as.numeric(sub(".*Elution: ([0-9\\.]+) min.*", "\\1", r$RT))*60
    r$RT <- as.integer(round(r$RT))
    r$TimePoint <- as.numeric(sub("h", "", r$TimePoint))
    r$ID <- gsub('[|,]', '_', r$ID)
    stopifnot(all(times %in% unique(r$TimePoint)))
    if(verbose) message("read ", file)        
    r
}

#' Extract peptide sequence information from a data file
#'
#' Extract peptide sequence information from a data file
#'
#' @param d the data file, from reading the Scaffold or ProteinPilot files
#' @param isotope the desired isotope
#' @param confidence.level the desired peptide identification confidence value
#' @param time.zero the time to treat as time "zero"
#' @param time.unit the name of the time unit, for display
#' @param mz.tol the m/z window
#' @param rt.tol the retention time window
#' @param extra.channels how many extra channels to include
#' @param aa.table the table to take element information from
#' @param verbose if verbose output desired
#' @return a list of peptide sequences, as from makeSequence
#' @export
makeSequences <- function(d, isotope, confidence.level, time.zero, time.unit,
                          mz.tol, rt.tol, extra.channels, 
                          aa.table=aaTable, verbose=TRUE) {
    if(verbose) message("finding Sequences ...")       
    ds <- split(d, paste(d$ID, d$Sequence, sep="--"))
    pars <- list(isotope=isotope, confidence.level=confidence.level,
                 time.zero=time.zero, time.unit=time.unit, 
                 mz.tol=mz.tol, rt.tol=rt.tol, extra.channels=extra.channels)    
    lapply(ds, makeSequence, pars, aa.table=aa.table)
}

#' Extract peptide sequence information for a single sequence from a data file
#'
#' Extract peptide sequence information for a single sequence from a data file
#' 
#' counts how many elements are present in the amino acid sequence
#' then builds out the channels based on the isotope under study
#' 
#' @return list with ID, Sequence, N, C, Channels, a list of data frames
#'   (one for each TimePoint) with dataframe par (RT, MZ) with
#'   possible multiple rows, if multiple found, and matriz MZ
#'   with MZ for each channel (as columns) for each row of par
#' 
#' if error, has "error" term with description of error,
#'   then doesn't have Times
#' @param x the data to get the information from
#' @param pars a list with parameter information, as passed from makeSequences
#' @param aa.table the table to take element information from
#' @export
makeSequence <- function(x, pars, aa.table=aaTable) {
    out <- list(ID=x$ID[1], Sequence=x$Sequence[1])
    x <- x[,c("TimePoint", "RT", "MW", "Z", "Confidence")]
    x <- x[order(x$TimePoint, x$RT),]
    x$data <- NA
    x$Total <- NA
    x$Error <- NA
    x$Error[x$Confidence < pars$confidence.level] <- "Below confidence threshold."
    rownames(x) <- NULL
    out$Times <- x
    out$pars <- pars
    ## get elements in sequence
    el <- tryCatch(sequenceElementCount(out$Sequence, isotope=pars$isotope,
                                        aa.table=aa.table),
                   error=identity)
    if(inherits(el, "error")) {
        out$error <- el$message
        return(out)
    } else {
        out$Elements <- el
    }
    checkTimes(out)
}

checkTimes <- function(x) {
    ## skip if already had error
    if(!is.null(x$error)) return(x)
    ts <- x$Times$TimePoint[is.na(x$Times$Error)]
    ## only move forward if there are zero time entries
    ## also need data at least one other time point
    if(length(ts)==0) {
        x$error <- "No valid entries found."
    } else if (!x$pars$time.zero %in% ts) {
        x$error <- "No valid entries at time zero."
    } else if (length(unique(ts)) < 2) {
        x$error <- "No valid entries other than time zero."
    }
    x
}

#' Convert raw mzXML files to Rdata files
#'
#' Convert raw mzXML files to Rdata files, and save in the same directory.
#'
#' @param dir directory that files are in
#' @param files names of mzXML files
#' @param Rfiles name of Rdata files
#' @return names of Rdata files
#' @export
rawToRDS <- function(dir, files, Rfiles) {
    if (!requireNamespace("xcms", quietly = TRUE)) {
      stop("xcms package needed. Please install it.",
        call. = FALSE)
    }
    stopifnot(file.exists(file.path(dir, files)))
    stopifnot(length(files)==length(Rfiles))
    for(i in seq_along(files)) {
        message("reading ", files[i])
        foo <- xcms::xcmsRaw(file.path(dir, files[i]))
        message("writing ", Rfiles[i])
        saveRDS(foo, file.path(dir, Rfiles[i]))
    }
    Rfiles
}

#' Get EIC data on specified sequences
#'
#' Get EIC data on specified sequences. This extracts the raw eic values, not the profiled values in class.
#' 
#' @return this returns a list of lists, one for each protein, each with
#' a list for each sequence of that protein, each with
#'   \item{ID}{(a single character string)}
#'   \item{Sequence}{(a single character string)}
#'   \item{Times}{}
#'   \item{Elements}{}
#'   \item{data}{data frame with:
#'      MZ: the mz,
#'      RT: the retention times,
#'      Count: matrix with columns for channels and rows for retention times.}
#' @param dat peptide sequences to extract, from makeSequences
#' @param files data files to read, or optionally, the objects themselves
#' @param type the type of data file, or R object
#' @param times the data collection times
#' @param dir the directory the files are in
#' @param verbose if verbose output desired
#' @export
addEICs <- function(dat, files, type=c("mzXML", "rds", "obj"),
                    times, dir, verbose=TRUE) {
    if (!requireNamespace("xcms", quietly = TRUE)) {
      stop("xcms package needed. Please install it.",
        call. = FALSE)
    }
    dat <- lapply(dat, checkTimes)
    stopifnot(length(files)==length(times))
    type <- match.arg(type)   
    for(i in seq_along(files)) {
        time <- times[i]
        if(type=="mzXML") {
            if(verbose) message("reading ", files[i])
            paths <- files
            if(!missing(dir)) paths <- file.path(dir, files)
            robj <- xcms::xcmsRaw(paths[i])
        } else if(type=="rds") {
            if(verbose) message("reading ", files[i])
            paths <- files
            if(!missing(dir)) paths <- file.path(dir, files)            
            robj <- readRDS(paths[i])
        } else { robj <- files[[i]] }
        if(verbose) message("adding EICs for Time ", time)
        if(verbose) pb <- utils::txtProgressBar(min = 0, max = length(dat), style = 1)
        for(j in seq_along(dat)) {
            dat[[j]] <- addEIC(dat[[j]], robj=robj, time=time)
            if(verbose) utils::setTxtProgressBar(pb, j) 
        }
        if(verbose) close(pb)
    }
    ## run checkTimes again as maybe some that looked okay couldn't be read
    dat <- lapply(dat, checkTimes)
    dat
}

## N14 <- 14.0030740052
## N15 <- 15.0001088984
## ISOTOPE.DELTA <- N15 - N14
## C12 <- 12.0
## C13 <- 13.0033548378
## ISOTOPE.DELTA <- C13 - C12
getCount <- function(robj, RT, MW, Z, nchannels, isotope, rt.tol, mz.tol, 
                     delta=list(N=0.9970348932, C=1.0033548378)) {
    if (!requireNamespace("xcms", quietly = TRUE)) {
      stop("xcms package needed. Please install it.",
        call. = FALSE)
    }
    stopifnot(isotope %in% names(delta))
    H_WEIGHT <- 1.007277
    mz <- (MW + (0:(nchannels-1))*delta[[isotope]]) / Z + H_WEIGHT
    rtr <- RT + rt.tol * c(-1, 1)
    count <- tryCatch({
        ## get first MZ to get RT and to get matrix size
        mzr <- mz[1] + mz.tol * c(-1, 1)
        tmp <- xcms::rawEIC(robj, mzrange=mzr, rtrange=rtr)
        RT <- robj@scantime[tmp$scan]
        out <- matrix(0, ncol=length(mz), nrow=length(tmp$intensity))
        out[,1] <- tmp$intensity
        ## now get the rest of the MZs
        for(i in 2:length(mz)) {
            mzr <- mz[i] + mz.tol * c(-1, 1)
            out[,i] <- xcms::rawEIC(robj, mzrange=mzr, rtrange=rtr)$intensity
        }
        out
    }, error=identity)
    if(inherits(count, "error")) {
        list(MZ=mz, RT=RT, Count=NA, error=count$message)        
    } else {
        list(MZ=mz, RT=RT, Count=count)
    }
}

addEIC <- function(x, robj, time, mz.tol=x$pars$mz.tol, rt.tol=x$pars$rt.tol,
                   extra.channels=x$pars$extra.channels) {
    if(!is.null(x$error)) return(x)
    k <- which(x$Times$TimePoint==time & is.na(x$Times$Error))
    if(length(k)==0) return(x)
    ## may be more than one row at this timepoint
    ## get the data for each row
    nel <- names(x$Elements)
    nc <- sum(x$Elements[nel==nel[1]]) + extra.channels
    counts <- vector("list", length(k))
    for(idx in seq_along(k)) {
        i <- k[idx]
        tmp <- getCount(robj, RT=x$Times$RT[i], MW=x$Times$MW[i], Z=x$Times$Z[i],
                        rt.tol=rt.tol, mz.tol=mz.tol, 
                        nchannels=nc, isotope=nel[1])
        if(!is.null(tmp$error)) {
            x$Times$Error[i] <- sub(" *\\n$", ".", tmp$error)
        }
         counts[[idx]] <- tmp
    }
    tot <- sapply(counts, function(x) sum(x$Count))
    x$Times$Total[k] <- tot
    x$Times$Error[which(x$Times$Total <= 0)] <- "No data found."
    ## get the max total
    m <- which.max(tot)
    ## make sure there is one and it's positive
    if(length(m) == 0) { return(x) }
    if(tot[m] <= 0) { return(x) }
    ## if so, choose the max and store in "data" list
    ## also keep track in Times which it was
    m0 <- k[m]
    nn <- as.character(time)
    x$Times$data[m0] <- nn
    if(is.null(x$data)) x$data <- list()    
    x$data[[nn]] <- counts[[m]]
    x
}

sequenceElementCount <- function(sequence, isotope, aa.table=aaTable) {
    ions <- c("C","H","N","O","S")
    ions <- ions[order(ions!=isotope)]
    sequence <- strsplit(sequence, '')[[1]]
    mm <- match(sequence, aa.table$Code)
    if(any(is.na(mm))){
        stop("contains unknown elements")
    }    
    out <- aa.table[mm, ions]
    out <- colSums(out)
    ## check for any that should be in natural abundance only
    nn <- paste0("n", isotope)
    if(nn %in% names(aa.table)) {
        reagent.count <- sum(aa.table[[nn]][mm])
        if (reagent.count > 0) {
            out[1] <- out[1] - reagent.count
            names(reagent.count) <- isotope
            out <- c(out, reagent.count)
        }
    }
    out
}

