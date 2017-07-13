####################################################################
## SECTION 1: Install R packages
####################################################################
## Install package "ProteinTurnover" and other required packages
## R-3.2.0 or higher is required

## To access mzXML files, the package "xcms" from Bioconductor is also required. For installation:
source("http://bioconductor.org/biocLite.R")
biocLite("xcms")  

## The ProteinTurnover package is available on GitHub, https://github.com/HegemanLab/ProteinTurnover

#####################################################################
## SECTION 2: Prepare data files
#####################################################################
## All MS data files were converted from .RAW to .mzXML files 
## Peptide identification file was exported from Scaffold spectrum report and saved as .csv file
# This .csv file needs to contain peptide retention time in column "Spectrum name" following this format "Elution: xx min"
# Please also refer Supporting Information Table S-3 for the required column headers in the .csv file
## All data files should be placed in the same data directory

#####################################################################
## SECTION 3: Run ProteinTurnover
#####################################################################
## load library
library(ProteinTurnover)

## set up list to hold parameters for this run
p <- list()

## Set the number of cores to use
## On Windows OS, this should be 1; on Linux, it can be set as 8
p$num.cores <- 1

## Set the data file directory 
## Set the directory that result files should be stored in
## Set scaffold file, which should be store in the data directory
## Note: For Windows OS users, replace all the \ with \\ from the copy of directory path

p$dir.data <- "C:\\Users\\USERNAME\\Documents\\mzs"  
# This is an example data directory 

p$dir.results <- "C:\\Users\\USERNAME\\Documents\\Turnover"  
# This is an example result directory

p$scaffoldfile <- "soluble-FDR.csv"   
# This is an example peptide file

## Set times as the labeling time points 
## set names of mzXML files that correspond to the above time points
## set names of rds (R data files) that correspond to the above times
# In this example, mzXML files in the data directory are named as hr_0.mzXML, hr_8.mzXML, etc. by users
# Note: R is case sensitive, including the file extension name

p$times <- c(0,4,8,16,24,32,40,48)
p$rawfiles <- paste0("hr_", p$times, ".mzXML") 
p$rdsfiles <- paste0("hr_", p$times, ".rds")

## Set parameters for use in the fit
p$isotope <- "N"
p$confidence.level <- 80  # Based on the peptide identification confidence value 

p$time.zero <- 0
p$time.unit <- "Hour" 
p$mz.tol <- 0.005      # m/z windows set as plus or minus 0.005
p$rt.tol <- 30         
# Retention time winodw set as plus or minus 30 sec

p$extra.channels <- 5  # Add extra 5 isotopic channels

## Set paprameter to get the relative abundance from EIC (relAbForTimes)
## available options are ("lm", "rlm", "lqs", "rq", "sum","log")
p$regression.model <- "lm"   # "lm" as linear regression

## Set up parameters for fitting
## options are "incorporation" for stable isotope label incorporation and "dilution" for label dilution

p$isotope.method <- "incorporation"

## Set up betabinomial options for the "M" parameter in Eq.(5)
## Available options are ("none", "one", "many")
# if users choose "none", then the binomial model will be used

p$M.model <- "one"

## Set the alpha parameters, Eq.(6)
## For incorporation, available options are ("many", "k", "kplateau", or "log2k")
## For dilution, available options are ("one", "many")
## "many" always means the algorithm will fit a different value of that parameter for each time point

p$alpha.model <- "log2k"     #make log2 value of turnover rate k

## passed to the "nab" parameter of pepcurve
p$plim <- NA

#######################################################
## create results directory
## will give warning if already exists
dir.create(p$dir.results)

#######################################################
## convert raw mzXML files to Rdata files
## This only need be done once
## This speeds up the steps below; it can be skipped if the raw files are used in addEICs
rawToRDS(p$dir.data, p$rawfiles, p$rdsfiles)

#######################################################
## read in all the data
## read the scaffold file
dat <- readScaffoldFile(p$scaffoldfile, times=p$times, dir=p$dir.data)

## Prepare the data structure for each sequence
seqs <- makeSequences(dat,
                      isotope=p$isotope,
                      confidence=p$confidence.level,
                      time.zero=p$time.zero,
                      time.unit=p$time.unit,                      
                      mz.tol=p$mz.tol,
                      rt.tol=p$rt.tol,
                      extra.channels=p$extra.channels)

## Get the EIC data for each sequence
seqs <- addEICs(seqs, files=p$rdsfiles, type="rds", times=p$times, dir=p$dir.data)

## Save what we've done so far
save(dat, seqs, p, file=file.path(p$dir.results, "eicdata.Rdata"))

########################################################
## Get fits for each sequence
## Reload what we've done so far
# Skip these two command lines if users just run the above code
dir.results <- "C:\\Users\\USERNAME\\Documents\\Turnover"  
load(file.path(dir.results, "eicdata.Rdata"))       

## Set up multiple cores
if(!is.null(p$num.cores)) { options(mc.cores=p$num.cores) } 

## Get the distribution abundance proportion and fit for each sequence
## mclapply.progress is just like lapply except works in parallel and shows progress
## makeRelAb is a wrapper for relAbForTimes, to handle any error
## makeFit is a wrapper for pepfit, to handle any error
# The algorithm chooses parameters that make the smallest difference between the observed proportion for each isotopic channel to the predicted proportion.
# In makeFit, if users prefer using relative abundance instead of proportion (the current default setting) for any reason, they can add "useRelAb=TRUE" as a new parameter 
# Across time points, the proportionality constant differs, so using relative abundance gives more weight to the time points with higher total values of relative abundance.
# Using proportions will treat each time point equally

seqs <- mclapply.progress(seqs, function(seq) {
  relAb <- makeRelAb(seq, regression.model=p$regression.model)
  fit <- makeFit(seq=seq, relAb=relAb, M=p$M.model, alpha=p$alpha.model, se=TRUE, isotope.method=p$isotope.method)
  list(seq=seq, relAb=relAb, fit=fit)    
})


## Save results so far
save(dat, seqs, p, file=file.path(p$dir.results, "fit-log2k.Rdata"))

########################################################
## Get parameters and save to a file
## Load results
# Skip these two command lines if users just run the above code
dir.results <- "C:\\Users\\USERNAME\\Documents\\Turnover"  
load(file.path(dir.results, "fit-log2k.Rdata"))     

## Get the parameter fits for those succeeded and write to file
fits <- getPars(seqs)
write.csv(fits, file=file.path(p$dir.results, "fits-log2k.csv"), na="")

########################################################
## Make HTML output for each seq
## Load what we've done so far
# Skip these two command lines if users just run the above code
dir.results <- "C:\\Users\\USERNAME\\Documents\\Turnover"  
load(file.path(dir.results, "fit-log2k.Rdata"))     

## set up multiple cores
if(!is.null(p$num.cores)) { options(mc.cores=p$num.cores) } 

## make the output
out <- makeSequenceHTML(seqs, dir=p$dir.results, file="results-log2k")

#####################################################################
## SECTION 4: Explore each peptide's fit as desired
#####################################################################
## Load saved result
dir.results <- "C:\\Users\\USERNAME\\Documents\\Turnover"  
load(file.path(dir.results, "fit-log2k.Rdata"))

## Pick a fit to explore
## Put the index number of whichever peptide you want to explore, eg. No.2 peptide in the output file
s <- seqs[[2]]

## Get usual output for the fit
s$fit
score.visual(s$fit)

## Make desired pictures for this No.2 peptide
# Make regression plots
regressionPlot(s$relAb)

# Make eic plots
plot(s$relAb)

# Make MLE plots
plot(s$fit)

# Make the plot of non-linear regression curve
plot(s$curve)

#####################################################################
## SECTION 5: Use alpha model as "many" instead
#####################################################################
## Load saved EIC data
dir.results <- "C:\\Users\\USERNAME\\Documents\\Turnover"  
load(file.path(dir.results, "eicdata.Rdata"))

## set alpha model differently
p$alpha.model <- "many"

## set up multiple cores
if(!is.null(p$num.cores)) { options(mc.cores=p$num.cores) } 

## get the distribution abundance proportion, fit, and curve for each sequence
seqs <- mclapply.progress(seqs, function(seq) {
  relAb <- makeRelAb(seq, regression.model=p$regression.model)
  fit <- makeFit(seq=seq, relAb=relAb, M=p$M.model, alpha=p$alpha.model, se=TRUE, isotope.method=p$isotope.method)
  list(seq=seq, relAb=relAb, fit=fit)    
  curve <- makeCurve(fit, plim=p$plim, isotope.method=p$isotope.method)
  list(seq=seq, relAb=relAb, fit=fit, curve=curve)
})

## save results so far
save(dat, seqs, p, file=file.path(p$dir.results, "fit-many.Rdata"))

## get the parameter fits for those succeeded and write to file
fits <- getPars(seqs)
write.csv(fits, file=file.path(p$dir.results, "fits-many.csv"), na="")

## output to HTML
out <- makeSequenceHTML(seqs, dir=p$dir.results, file="results-many")


