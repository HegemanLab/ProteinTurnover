# 4/1, isotopes7d
# go back to convolving using outer as the inexactness in the
# builtin convolve function is apparently causing issues

# 4/24, isotopes7e
# add option to choose convolve function

# 7/24, isotopes 7f
# take out above option; put in fast exact function
# (though not quite as fast as FFT)
# make option work by setting conv to the desired function

# 10/13, isotopes 8
# add C code for convolution

# Exact Convolution Function

#' @useDynLib ProteinTurnover
convC <- function(x, y, scale=TRUE) {
    if(scale) {
        x<-x/sum(x)
        y<-y/sum(y)
    }
       .C("convolve",
          as.double(x),
          as.integer(length(x)),
          as.double(y),
          as.integer(length(y)),
          xy = double(length(x) + length(y) - 1))$xy
     }

getabundance<-function(newab) {
  if(missing(newab)) newab<-NULL
  abundance<-list( C=c(98.930, 1.070),
                   H=c(99.985, 0.015),
                   O=c(99.757, 0.038, 0.205),
                   N=c(99.632, 0.368),
                   S=c(95.020, 0.750, 4.210, 0.000, 0.02),
                   P=100,
                   Si=c(92.2101, 4.6999, 3.090) )
  abundance <- lapply(abundance, function(x) x/sum(x))

  # put new abundance in abundance list
  if(!is.null(newab)) {
    stop.dup <- duplicated(names(newab))
    stopif("Abundance for element repeated", stop.dup, names(newab))
    for(i in seq_along(newab)) {
      abundance[[names(newab)[i]]]<-newab[[i]]
    }
  }

  stop.neg <- sapply(abundance, function(x) any(x<0))
  stopif("Element has negative abundance: ", stop.neg)

  warn.one <- sapply(abundance, function(x) {
    s <- sum(x)
    !isTRUE(all.equal(s, 1)) & !isTRUE(all.equal(s, 100))
  })
  warnif("Abundance of elements corrected to sum to 100%", warn.one)
  abundance <- lapply(abundance, function(x) x/sum(x))

  abundance
}

getelementprob<-function(n, p) {
  if(length(p)==2) {
    elementprob<-dbinom(0:n,n,p[2])
  } else {
    elementprob<-1
    for(i in 1:n) elementprob<-convC(elementprob, p)
  }
  elementprob
}

mytext<-function(text, x, n=names(x)) { paste(text, ": ", paste(unique(n[x]), collapse=" "), sep="") }
warnif<-function(text, x, n=names(x)) { if(any(x)) warning(mytext(text, x, n), call.=FALSE) }
stopif<-function(text, x, n=names(x)) { if(any(x)) stop(mytext(text, x, n), call.=FALSE) }                               
#' Calculate natural abundance of isotopes over several elements.
#'
#' For a set of elements, compute the total natural abundance.
#'
#' It computes the natural abundance for each element, and convolves
#' them together to get the total natural abundance.
#' 
#' @param \dots The number of each element of interest.
#' @param abundance   A list of the natural abundance for each element.
#' By default, it knows about C, H, O, N, and S.
#' @return  A numeric vector with the natural abundance of each unit mass extra,
#' starting with zero extra.
#' @examples
#' # natural abundance for 1 nitrogen element
#' calcEnv(N=1)
#'  
#' # combined natural abundance for 2 nitrogen elements
#' calcEnv(N=2)
#'  
#' # combined natural abundance for a new element Z
#' calcEnv(Z=1, abundance=list(Z=c(0.7, 0.2, 0.1)))
#'  
#' # combined natural abundance for Z and N
#' calcEnv(N=1, Z=1, abundance=list(Z=c(0.7, 0.2, 0.1)))
#' @export

calcEnv <- function(..., abundance=NULL) {
  elements <- list(...)
  if( class(elements[[1]])=="list" ) elements <- elements[[1]]
  abundance<-getabundance(abundance)
  
  # check input
  warn.dup <- duplicated(names(elements))
  warnif("Duplicated elements found", warn.dup, names(elements))

  warn.zero <- elements<=0
  elements <- elements[!warn.zero]
  warnif("Elements with counts of zero or less were removed", warn.zero)

  stop.noab <- (! names(elements) %in% names(abundance))
  stopif("No abundance given for elements", stop.noab, names(elements))
       
  # now do the calculation
  prob <- 1    
  for (j in seq_along(elements)) {
    elementprob <- getelementprob(elements[[j]], abundance[[names(elements)[j]]])
    prob <- convC(prob, elementprob)
  }
  as.numeric(prob)
}
