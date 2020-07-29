library(spectral)
library(progress)
library(zoo)
library(plyr)

lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}


normalizeProfileMatrix <- function(intMatrix)
  t(aaply(intMatrix, 2, function(clmn) {
  clmn <- na.locf(clmn)
  clmn / max(clmn)
}, .progress = "text"))


lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}


# We decide that this is OK and recalculate profilesFreq with "only" this 
# reduced set of frequencies
lombPeriodograms <- function(intMatrix, sampleList, 
                             baseTime, f) {
  sampleList$time <- sampleTimes(sampleList$mtime)
  
  sampleListFiltrate <- sampleList[sampleList$sampleType == "filtrate",,drop=FALSE]
  intMatrix <- intMatrix[as.character(sampleListFiltrate$sampleIDs),,drop=FALSE]
  
  timeSeries <- sampleListFiltrate$time - baseTime
  alply(
    intMatrix, 2, function(row) spec.lomb(y = row, x = timeSeries, f = f),
    .progress = "text")
}


lombMatrix <- function(lombPeriodograms) {
  laply(lombPeriodograms, function(spec)
    c(spec$A, cos(spec$phi), sin(spec$phi)))
}


sampleTimes <- function(rawTimes) {
  as.numeric(as.POSIXct(rawTimes))/60/24/60
}
