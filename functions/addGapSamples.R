#' Fill into equally-spaced time matrix
#'
#' The `intMatrix` for filtrate samples is filled into a 
#' matrix with equally-spaced time distances, i.e. adding
#' `NA` rows where samples are missing in between (from blanks,
#' spikes, calibration, maintenance, other problems.)
#' 
#' 
#' @param sampleList 
#' @param intMatrix 
#'
#' @return
#' @export
#'
#' @examples
addGapSamples <- function(sampleList, intMatrix)
{
  # calculate the number of skip samples to add after every row
  sampleListFiltrate <- sampleList[sampleList$sampleType == "filtrate",,drop=FALSE]
  sampleTime <- as.POSIXct(sampleListFiltrate$mtime)
  
  diffTime <- as.numeric(diff(sampleTime))
  sampleDuration <- median(diffTime)
  skipSamples <- c(0, round((diffTime - sampleDuration) / sampleDuration))
  targetRows <- seq_len(nrow(intMatrix)) + cumsum(skipSamples)
  
  # Make sure the intensity matrix is ordered correctly
  intMatrix <- intMatrix[as.character(sampleListFiltrate$sampleIDs),,drop=FALSE]
  filledRows <- nrow(intMatrix) + sum(skipSamples)
  filledMatrix <- matrix(NA, nrow=filledRows, ncol=ncol(intMatrix))
  # Vector of target row per original row
  filledMatrix[targetRows,] <- intMatrix
  return(filledMatrix)
}


# Checking time correction

.timeCorrection <- function(sampleTime)
{
  diffTime <- as.numeric(diff(sampleTime))
  sampleDuration <- median(diffTime)
  skipSamples <- c(0, round((diffTime - sampleDuration) / sampleDuration))
  targetRows <- seq_len(nrow(intMatrix)) + cumsum(skipSamples)
  
  correctedTime <- (targetRows-1) * sampleDuration
  originalTime <- cumsum(c(0, as.numeric(diffTime)))
  
  df <-data.frame(originalTime=originalTime, diffTime= c(0, diffTime), skipSamples, 
                  targetRows, correctedTime)
  df$timeShift <- df$correctedTime - df$originalTime
}
