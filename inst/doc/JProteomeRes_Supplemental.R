## ---- setup, include=FALSE-----------------------------------------------
 knitr::opts_chunk$set(eval=FALSE)

## ------------------------------------------------------------------------
#  system.file(package="ProteinTurnover", "extscripts", "JProteomeRes_Supplemental.R")

## ------------------------------------------------------------------------
#  f <- system.file(package="ProteinTurnover", "extscripts", "JProteomeRes_Supplemental.R")
#  file.copy(f, ".")

## ------------------------------------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("xcms")

## ------------------------------------------------------------------------
#  devtools::install_github("HegemanLab/ProteinTurnover")

## ------------------------------------------------------------------------
#  R.version.string
#  # [1] "R version 3.3.2 (2016-10-31)"
#  packageVersion("ProteinTurnover")
#  # [1] '1.1.1'

## ---- eval=TRUE----------------------------------------------------------
library(ProteinTurnover)

## ------------------------------------------------------------------------
#  p <- list()

## ------------------------------------------------------------------------
#  p$num.cores <- 1

## ------------------------------------------------------------------------
#  p$dir.data <- "mzs"

## ------------------------------------------------------------------------
#  p$dir.results <- "Turnover"

## ------------------------------------------------------------------------
#  p$scaffoldfile <- "soluble.csv"

## ------------------------------------------------------------------------
#  p$times <- c(0,4,8,16,24,32,40,48)
#  p$rawfiles <- paste0("hr_", p$times, ".mzXML")
#  p$rdsfiles <- paste0("hr_", p$times, ".rds")

## ------------------------------------------------------------------------
#  p$isotope <- "N"
#  p$confidence.level <- 80
#  p$time.zero <- 0
#  p$time.unit <- "Hour"
#  p$mz.tol <- 0.005
#  p$rt.tol <- 30
#  p$extra.channels <- 5

## ------------------------------------------------------------------------
#  p$regression.model <- "lm"   # "lm" as linear regression

## ------------------------------------------------------------------------
#  p$isotope.method <- "incorporation"

## ------------------------------------------------------------------------
#  p$M.model <- "one"

## ------------------------------------------------------------------------
#  p$alpha.model <- "log2k"     #make log2 value of turnover rate k

## ------------------------------------------------------------------------
#  p$plim <- NA

## ------------------------------------------------------------------------
#  dir.create(p$dir.results)

## ------------------------------------------------------------------------
#  rawToRDS(p$dir.data, p$rawfiles, p$rdsfiles)
#  # reading hr_0.mzXML (time: 1m 38.3s)
#  # writing hr_0.rds (time: 53.2s)
#  # reading hr_4.mzXML (time: 1m 39.7s)
#  # writing hr_4.rds (time: 54.2s)
#  # reading hr_8.mzXML (time: 1m 49.0s)
#  # writing hr_8.rds (time: 56.9s)
#  # reading hr_16.mzXML (time: 1m 19.1s)
#  # writing hr_16.rds (time: 41.9s)
#  # reading hr_24.mzXML (time: 1m 42.0s)
#  # writing hr_24.rds (time: 56.3s)
#  # reading hr_32.mzXML (time: 1m 41.6s)
#  # writing hr_32.rds (time: 56.4s)
#  # reading hr_40.mzXML (time: 1m 42.0s)
#  # writing hr_40.rds (time: 56.6s)
#  # reading hr_48.mzXML (time: 1m 46.8s)
#  # writing hr_48.rds (time: 53.6s)
#  # Warning message:
#  # 'loadRcppModules' is deprecated.
#  # Use 'loadModule' instead.
#  # See help("Deprecated")

## ------------------------------------------------------------------------
#  dat <- readScaffoldFile(p$scaffoldfile, times=p$times, dir=p$dir.data)
#  # reading Scaffold file ...
#  # read soluble.csv

## ------------------------------------------------------------------------
#  seqs <- makeSequences(dat,
#                        isotope=p$isotope,
#                        confidence=p$confidence.level,
#                        time.zero=p$time.zero,
#                        time.unit=p$time.unit,
#                        mz.tol=p$mz.tol,
#                        rt.tol=p$rt.tol,
#                        extra.channels=p$extra.channels)
#  # finding Sequences ...

## ------------------------------------------------------------------------
#  seqs <- addEICs(seqs, files=p$rdsfiles, type="rds", times=p$times, dir=p$dir.data)
#  # reading hr_0.rds (time: 27.0s)
#  # adding EICs for Time 0 (time: 4.5s)
#  # reading hr_4.rds (time: 27.4s)
#  # adding EICs for Time 4 (time: 3.7s)
#  # reading hr_8.rds (time: 28.9s)
#  # adding EICs for Time 8 (time: 2.9s)
#  # reading hr_16.rds (time: 21.0s)
#  # adding EICs for Time 16 (time: 2.2s)
#  # reading hr_24.rds (time: 28.3s)
#  # adding EICs for Time 24 (time: 2.5s)
#  # reading hr_32.rds (time: 3.2s)
#  # adding EICs for Time 32 (time: 2.0s)
#  # reading hr_40.rds (time: 3.2s)
#  # adding EICs for Time 40 (time: 1.7s)
#  # reading hr_48.rds (time: 28.6s)
#  # adding EICs for Time 48 (time: 1.8s)

## ------------------------------------------------------------------------
#  save(dat, seqs, p, file=file.path(p$dir.results, "eicdata.Rdata"))

## ------------------------------------------------------------------------
#  dir.results <- "Turnover"
#  load(file.path(dir.results, "eicdata.Rdata"))

## ------------------------------------------------------------------------
#  if(!is.null(p$num.cores)) { options(mc.cores=p$num.cores) }

## ------------------------------------------------------------------------
#  seqs <- mclapply.progress(seqs, function(seq) {
#    relAb <- makeRelAb(seq, regression.model=p$regression.model)
#    fit <- makeFit(seq=seq, relAb=relAb, M=p$M.model, alpha=p$alpha.model, se=TRUE, isotope.method=p$isotope.method)
#    list(seq=seq, relAb=relAb, fit=fit)
#  })
#  # Processing 1854 items, in parallel on 8 cores
#  # Started Sat Jul 15 16:03:51 2017
#  # ================================================================================
#  # Ended   Sat Jul 15 16:06:08 2017; elapsed time: 2m 17.0s

## ------------------------------------------------------------------------
#  save(dat, seqs, p, file=file.path(p$dir.results, "fit-log2k.Rdata"))

## ------------------------------------------------------------------------
#  dir.results <- "Turnover"
#  load(file.path(dir.results, "fit-log2k.Rdata"))

## ------------------------------------------------------------------------
#  fits <- getPars(seqs)
#  # There were 50 or more warnings (use warnings() to see the first 50)
#  write.csv(fits, file=file.path(p$dir.results, "fits-log2k.csv"), na="")

## ------------------------------------------------------------------------
#  dir.results <- "Turnover"
#  load(file.path(dir.results, "fit-log2k.Rdata"))

## ------------------------------------------------------------------------
#  if(!is.null(p$num.cores)) { options(mc.cores=p$num.cores) }

## ------------------------------------------------------------------------
#  out <- makeSequenceHTML(seqs, dir=p$dir.results, file="results-log2k")
#  # Creating output for sequences with errors (these are faster).
#  # Processing 866 items, in parallel on 8 cores
#  # Started Sat Jul 15 16:06:20 2017
#  # ================================================================================
#  # Ended   Sat Jul 15 16:06:55 2017; elapsed time: 35.0s
#  # Creating output for sequences without errors.
#  # Processing 988 items, in parallel on 8 cores
#  # Started Sat Jul 15 16:06:55 2017
#  # ================================================================================
#  # Ended   Sat Jul 15 16:19:49 2017; elapsed time: 12m 54.0s
#  # There were 50 or more warnings (use warnings() to see the first 50)

## ------------------------------------------------------------------------
#  dir.results <- "Turnover"
#  load(file.path(dir.results, "fit-log2k.Rdata"))

## ---- eval=TRUE, include=FALSE-------------------------------------------
s <- ProteinTurnover:::example2

## ------------------------------------------------------------------------
#  s <- seqs[[2]]

## ---- eval=TRUE----------------------------------------------------------
s$fit
score.visual(s$fit)

## ---- fig.width=15, fig.height=5, eval=TRUE------------------------------
regressionPlot(s$relAb)

## ---- fig.width=15, fig.height=5, eval=TRUE------------------------------
plot(s$relAb)

## ---- eval=TRUE, fig.width=7, fig.height=3-------------------------------
plot(s$fit)

## ------------------------------------------------------------------------
#  dir.results <- "Turnover"
#  load(file.path(dir.results, "eicdata.Rdata"))

## ------------------------------------------------------------------------
#  p$alpha.model <- "many"

## ------------------------------------------------------------------------
#  if(!is.null(p$num.cores)) { options(mc.cores=p$num.cores) }

## ------------------------------------------------------------------------
#  seqs <- mclapply.progress(seqs, function(seq) {
#    relAb <- makeRelAb(seq, regression.model=p$regression.model)
#    fit <- makeFit(seq=seq, relAb=relAb, M=p$M.model, alpha=p$alpha.model, se=TRUE, isotope.method=p$isotope.method)
#    list(seq=seq, relAb=relAb, fit=fit)
#    curve <- makeCurve(fit, plim=p$plim, isotope.method=p$isotope.method)
#    list(seq=seq, relAb=relAb, fit=fit, curve=curve)
#  })
#  # Processing 1854 items, in parallel on 8 cores
#  # Started Sat Jul 15 16:20:01 2017
#  # ================================================================================
#  # Ended   Sat Jul 15 16:23:44 2017; elapsed time: 3m 43.0s

## ------------------------------------------------------------------------
#  save(dat, seqs, p, file=file.path(p$dir.results, "fit-many.Rdata"))

## ------------------------------------------------------------------------
#  fits <- getPars(seqs)
#  # There were 50 or more warnings (use warnings() to see the first 50)
#  write.csv(fits, file=file.path(p$dir.results, "fits-many.csv"), na="")

## ------------------------------------------------------------------------
#  out <- makeSequenceHTML(seqs, dir=p$dir.results, file="results-many")
#  # Creating output for sequences with errors (these are faster).
#  # Processing 853 items, in parallel on 8 cores
#  # Started Sat Jul 15 16:23:55 2017
#  # ================================================================================
#  # Ended   Sat Jul 15 16:24:13 2017; elapsed time: 18.0s
#  # Creating output for sequences without errors.
#  # Processing 1001 items, in parallel on 8 cores
#  # Started Sat Jul 15 16:24:13 2017
#  # ================================================================================
#  # Ended   Sat Jul 15 16:38:23 2017; elapsed time: 14m 10.0s
#  # There were 50 or more warnings (use warnings() to see the first 50)

## ---- eval=TRUE, include=FALSE-------------------------------------------
s <- ProteinTurnover:::example2

## ------------------------------------------------------------------------
#  s <- seqs[[2]]

## ---- eval=TRUE-----------------------------------------------------------------------------------
options(width=100)
s$fit
score.visual(s$fit)

## ---- eval=TRUE, fig.width=7, fig.height=3--------------------------------------------------------
plot(s$fit)

## ---- eval=TRUE, fig.width=6, fig.height=3--------------------------------------------------------
plot(s$curve)

## ---- include=FALSE, eval=FALSE-------------------------------------------------------------------
#  download.file("http://users.stat.umn.edu/~rend0020/ProteinTurnoverScript/Turnover/fit-log2k.Rdata", "fit-log2k.Rdata")
#  load("fit-log2k.Rdata")
#  example1 <- seqs[[2]]
#  
#  download.file("http://users.stat.umn.edu/~rend0020/ProteinTurnoverScript/Turnover/fit-many.Rdata", "fit-many.Rdata")
#  load("fit-many.Rdata")
#  example2 <- seqs[[2]]
#  
#  aaTable <- ProteinTurnover:::aaTable
#  LocusNumbers <- ProteinTurnover:::LocusNumbers
#  devtools::use_data(aaTable, LocusNumbers, example1, example2, internal = TRUE)
#  
#  file.remove("fit-log2k.Rdata")
#  file.remove("fit-many.Rdata")

