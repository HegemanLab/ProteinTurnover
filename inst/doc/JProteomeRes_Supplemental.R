## ---- setup, include=FALSE-----------------------------------------------
 knitr::opts_chunk$set(eval=FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  system.file(package="ProteinTurnover", "extscripts", "JProteomeRes_Supplemental.R")

## ---- eval=FALSE---------------------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("xcms")

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github(HegemanLab/ProteinTurnover)

## ------------------------------------------------------------------------
#  library(ProteinTurnover)

## ------------------------------------------------------------------------
#  p <- list()

## ------------------------------------------------------------------------
#  p$num.cores <- 1

## ------------------------------------------------------------------------
#  p$dir.data <- "C:\\Users\\USERNAME\\Documents\\mzs"

## ------------------------------------------------------------------------
#  p$dir.results <- "C:\\Users\\USERNAME\\Documents\\Turnover"

## ------------------------------------------------------------------------
#  p$scaffoldfile <- "soluble-FDR.csv"

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

## ------------------------------------------------------------------------
#  dat <- readScaffoldFile(p$scaffoldfile, times=p$times, dir=p$dir.data)

## ------------------------------------------------------------------------
#  seqs <- makeSequences(dat,
#                        isotope=p$isotope,
#                        confidence=p$confidence.level,
#                        time.zero=p$time.zero,
#                        time.unit=p$time.unit,
#                        mz.tol=p$mz.tol,
#                        rt.tol=p$rt.tol,
#                        extra.channels=p$extra.channels)

## ------------------------------------------------------------------------
#  seqs <- addEICs(seqs, files=p$rdsfiles, type="rds", times=p$times, dir=p$dir.data)

## ------------------------------------------------------------------------
#  save(dat, seqs, p, file=file.path(p$dir.results, "eicdata.Rdata"))

## ------------------------------------------------------------------------
#  dir.results <- "C:\\Users\\USERNAME\\Documents\\Turnover"
#  load(file.path(dir.results, "eicdata.Rdata"))

## ------------------------------------------------------------------------
#  if(!is.null(p$num.cores)) { options(mc.cores=p$num.cores) }

## ------------------------------------------------------------------------
#  seqs <- mclapply.progress(seqs, function(seq) {
#    relAb <- makeRelAb(seq, regression.model=p$regression.model)
#    fit <- makeFit(seq=seq, relAb=relAb, M=p$M.model, alpha=p$alpha.model, se=TRUE, isotope.method=p$isotope.method)
#    list(seq=seq, relAb=relAb, fit=fit)
#  })

## ------------------------------------------------------------------------
#  save(dat, seqs, p, file=file.path(p$dir.results, "fit-log2k.Rdata"))

## ------------------------------------------------------------------------
#  dir.results <- "C:\\Users\\USERNAME\\Documents\\Turnover"
#  load(file.path(dir.results, "fit-log2k.Rdata"))

## ------------------------------------------------------------------------
#  fits <- getPars(seqs)
#  write.csv(fits, file=file.path(p$dir.results, "fits-log2k.csv"), na="")

## ------------------------------------------------------------------------
#  dir.results <- "C:\\Users\\USERNAME\\Documents\\Turnover"
#  load(file.path(dir.results, "fit-log2k.Rdata"))

## ------------------------------------------------------------------------
#  if(!is.null(p$num.cores)) { options(mc.cores=p$num.cores) }

## ------------------------------------------------------------------------
#  out <- makeSequenceHTML(seqs, dir=p$dir.results, file="results-log2k")

## ------------------------------------------------------------------------
#  dir.results <- "C:\\Users\\USERNAME\\Documents\\Turnover"
#  load(file.path(dir.results, "fit-log2k.Rdata"))

## ------------------------------------------------------------------------
#  s <- seqs[[2]]

## ------------------------------------------------------------------------
#  s$fit
#  score.visual(s$fit)

## ------------------------------------------------------------------------
#  regressionPlot(s$relAb)

## ------------------------------------------------------------------------
#  plot(s$relAb)

## ------------------------------------------------------------------------
#  plot(s$fit)

## ------------------------------------------------------------------------
#  plot(s$curve)

## ------------------------------------------------------------------------
#  dir.results <- "C:\\Users\\USERNAME\\Documents\\Turnover"
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

## ------------------------------------------------------------------------
#  save(dat, seqs, p, file=file.path(p$dir.results, "fit-many.Rdata"))

## ------------------------------------------------------------------------
#  fits <- getPars(seqs)
#  write.csv(fits, file=file.path(p$dir.results, "fits-many.csv"), na="")

## ------------------------------------------------------------------------
#  out <- makeSequenceHTML(seqs, dir=p$dir.results, file="results-many")

