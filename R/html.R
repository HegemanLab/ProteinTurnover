###############################################################################
# R2HTML Functions
################################################################################


HTML.error <- function(x, file="") {
    cat("<p>", x$message, "</p>\n", sep="", file=file, append=TRUE)
}

HTML.img <- function(img, scale=60, file="") {
    if(!is.null(img$error)) {
        cat("<p>", img$error, "</p>\n", sep="", file=file, append=TRUE)
    } else {
        cat('<img src="', img$file, '" width=', img$width*60,',  height=', img$height*60,
            '>\n', file=file, append=TRUE, sep="")
    }
}

myHTMLtable <- function(x, digits = getOption("R2HTML.format.digits"), 
                        nsmall = getOption("R2HTML.format.nsmall"),
                        big.mark = getOption("R2HTML.format.big.mark"), 
                        big.interval = getOption("R2HTML.format.big.interval"),
                        decimal.mark = getOption("R2HTML.format.decimal.mark"),
                        align='left', Border=0, row.names=FALSE, file="", ...) {
    if(inherits(x, "error")) {
        R2HTML::HTML(x, file=file)
    } else {
        R2HTML::HTML(format(x, digits=digits, nsmall=nsmall,big.mark=big.mark,
                    big.interval=big.interval,decimal.mark=decimal.mark, trim=TRUE),
             align=align, Border=Border, row.names=row.names, file=file, ...)
    }
}
  
prepareHTML <- function(x, dir, subdir, name) {
    if(missing(name)) {
        seq <- gsub("([a-z])","\\1!", x$seq$Sequence)
        name <- paste0(x$seq$ID, "--", seq)
    }
    dir0 <- dir
    dir <- file.path(dir, subdir)
    if(!file.exists(dir)) { dir.create(dir, recursive=TRUE) }   
    out <- list(ID=x$seq$ID,
                Sequence=x$seq$Sequence,
                pars = x$seq$pars,
                fit  = getFit(x$fit),
                par  = if(inherits(x$fit, "error")) x$fit else rbind(est=x$fit$par),
                score  = if(inherits(x$fit, "error")) x$fit else rbind(unlist(score.visual(x$fit))),
                data = getData(x$seq),
                dir  = dir0,
                subdir = subdir, 
                file = paste0(name, ".html"),
                eic.plot = save.graphic(x$relAb, dir=dir,
                    file = paste0(name, "_eic")),
                reg.plot = save.graphic(x$relAb, type="regression", dir=dir,
                    file = paste0(name, "_reg")),
                fit.plot = save.graphic(x$fit, main=NULL, sub=NULL, dir=dir, panel.size=2,
                    file = paste0(name, "_fit")) )
    if(!is.null(x$fit$vcov)) {
        out$par.se <- sqrt(diag(x$fit$vcov))
    }
    if(!is.null(x$curve)) {
        out$curve <- getCurve(x$curve)
        out$curve.plot = save.graphic(x$curve, main=NULL, dir=dir, panel.size=3,
            file = paste0(x$seq$Sequence, "_curve"))
    }
    out
}
outputHTML <- function(x) {
    f <- file.path(x$dir, x$subdir, x$file)
    SummaryReportBegin(f)
    R2HTML::HTML.title(paste0(x$ID, ": ", x$Seq),file=f,HR=2)
    if(!is.null(x$curve)) {
        R2HTML::HTML.title("Fitted Curves", file=f, HR=2)
        myHTMLtable(x$curve, file=f)
        HTML.img(x$curve.plot, file=f)
    }
    R2HTML::HTML.title("Model Fits", file=f, HR=2)
    if(!is.null(x$par.se)) {
        myHTMLtable(rbind(x$par, se=x$par.se), file=f)
    } else {
        myHTMLtable(x$par, file=f)
    }
    R2HTML::HTML.title("Visual Score", file=f, HR=2)
    myHTMLtable(x$score, file=f)
    R2HTML::HTML.title("Model Plot", file=f, HR=2)
    HTML.img(x$fit.plot, file=f)
    R2HTML::HTML.title("Available Data", file=f, HR=2)    
    myHTMLtable(x$data, file=f)
    R2HTML::HTML.title("EIC Data", file=f, HR=2)
    HTML.img(x$eic.plot, file=f)
    R2HTML::HTML.title("Relative Abundance Fits", file=f, HR=2)
    HTML.img(x$reg.plot, file=f)
    R2HTML::HTML.title("Parameters Used", file=f, HR=2)
    myHTMLtable(data.frame(cbind(value=x$pars)), row.names=TRUE, file=f, col.names=FALSE)
    cat("</body></html>", file=f, append=TRUE)
    invisible(x)
}

## make individual pages
sequencesHTML <- function(seqs, dir, subdir="sequences") {
    ok <- !sapply(seqs, function(x) inherits(x$fit, "error"))
    seqs1 <- seqs[!ok]
    seqs2 <- seqs[ok]
    message("Creating output for sequences with errors (these are faster).")
    out1 <- mclapply.progress(seqs1, function(seq) {
        outputHTML(prepareHTML(seq, dir=dir, subdir=subdir))
    }, mc.preschedule=FALSE)
    message("Creating output for sequences without errors.")
    out2 <- mclapply.progress(seq_along(seqs2), function(j) {
        outputHTML(prepareHTML(seqs2[[j]], dir=dir, subdir=subdir))
    }, mc.preschedule=FALSE)
    out <- vector(mode="list", length=length(seqs))
    out[!ok] <- out1
    out[ok] <- out2
    out
}

#' Make HTML output for a list of fits
#' 
#' Make HTML output for a list of fits
#'
#' @param seqs A list of completed fits. Each element should be a list with three elements: \itemize{
#'   \item{seq}{the original sequence information}
#'   \item{relAb}{the relative abundance}
#'   \item{fit}{the fit}
#' }
#' @param dir the directory the output should be in
#' @param file the main file name (without the .html added)
#' @param subdir the subdirectory the peptide sequence results should go in
#' @return information on the individual pages, invisibly
#' @export
makeSequenceHTML <- function(seqs, dir, file, subdir=paste0(file,"-seqs")) {
    out <- sequencesHTML(seqs, dir=dir, subdir=subdir)
    summaryHTML(seqs, out, dir=dir, file=paste0(file, ".html"))
    invisible(out)
}


summaryHTML <- function(seqs, out, dir, file) {
    f <- file.path(dir, file)
    pp <- getPars(seqs)
    ff <- sapply(out, function(x) file.path(x$subdir, x$file))
    pp <- cbind(Num=1:length(seqs), Link=paste0('<a href="', ff, '">link</a>'), pp)
    SummaryReportBegin(f)
    R2HTML::HTML.title("All Sequences", file=f, HR=2)
    myHTMLtable(pp, file=f)
    cat("</body></html>", file=f, append=TRUE)    
}

SummaryReportBegin <- function(file="summary_report.html",title="Half-Life Summary Report") {
  cat(paste('<html><head><title>', title, '</title>
    <style type="text/css"> 
        h2 {  background-color:#CCCCCC;}
        h3 {  background-color:#EBEBEB;}
tr    
{
	font-family: "Lucida Sans Unicode", "Lucida Grande", Sans-Serif;
	font-size: 12px;
	background: #fff;
	border-collapse: collapse;
	text-align: left;
        white-space: nowrap;
}

th
{
	font-size: 14px;
	font-weight: normal;
	color: #039;
	padding: 2px 2px;
	border-bottom: 1px solid #6678b1;
}
td
{
	padding: 2px 5px 0px 0px;
}

        
        footer {text-align:right; font-style:italics;}
        .protein_seq {
            display: inline-block;
        }
        .character { width: inherit; display: inherit}
        span { 
            background-color: #D8BFD8;
            color: #800000;
            border: 1px solid black;
        }
    </style>
    </head>',sep=""),
    file=file,append=FALSE)
}

