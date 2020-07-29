library(plyr)
library(dplyr)
library(magrittr)
library(tidyr)

#' Filter profiles
#'
#' Ratio factors must be >1 summed, i.e. 
#' $blankRatio/blankRatioFactor + isRatio/isRatioFactor > 1$.
#' 
#' Note that `blankRatio` and `isRatio` are the *maximum* peak intensity
#' in filtrate samples versus the *mean* intensity in blanks or non-IS filtrate, respectively,
#' in order to not lose intense peaks that appear in few samples.
#' 
#' Note: non-IS filtrate is used as a comparison in order to not lose IS peaks, which are 
#' expected to be highly abundant in samples, but should be missing in non-IS filtrate.
#' 
#' `minCountSamples` guards against the inclusion of profiles with only spurious hits.
#' 
#' Note that we disregard spikes and calibration, therefore if a peak is not in filtrates, it should
#' get lost even though it is in calibration.
#'
#' @param peaksTable A `peaksTable` (profile peaks table with joined sample information)
#' @param minCountSamples Minimal count (number of detects) in filtrate samples
#' @param blankRatioFactor Scaling of blank ratio
#' @param isRatioFactor Scaling of IS ratio
#' @param sampleIntensityCutoff Minimum intensity of the top peak in filtrate
#'
#' @return
#' @export
#'
#' @examples
filterProfiles <- function(peaksTable,
                           settings = list(
                             minCountSamples = 8,
                             blankRatioFactor = 4,
                             isRatioFactor = 2,
                             sampleIntensityCutoff = 4e5))
{
  minCountSamples <- settings$minCountSamples
  blankRatioFactor <- settings$blankRatioFactor
  isRatioFactor <- settings$isRatioFactor
  sampleIntensityCutoff <- settings$sampleIntensityCutoff
  
  profilesFiltrate <- 
    peaksTable %>%
    group_by(sampleType) %>%
    mutate(samples = n_distinct(sampleIDs)) %>%
    filter(sampleType == "filtrate") %>%
    group_by(profileIDs) %>%
    dplyr::summarize(
      countSamples = n(),
      maxIntensitySamples = max(intensity),
      meanDetectedIntensitySamples = mean(intensity),
      meanIntensitySamples = sum(intensity/samples)
    )
  
  profilesBlank <-
    peaksTable %>%
    group_by(sampleType) %>%
    mutate(samples = n_distinct(sampleIDs)) %>%
    filter(sampleType == "blank") %>%
    group_by(profileIDs) %>%
    dplyr::summarize(
      countBlank = n(),
      meanIntensityBlank = sum(intensity/samples),
      meanIntensityBlankDetected = mean(intensity)
    )
  
  profilesNoIs <-
    peaksTable %>%
    group_by(sampleType) %>%
    mutate(samples = n_distinct(sampleIDs)) %>%
    filter(sampleType == "filtrate_noIS") %>%
    group_by(profileIDs) %>%
    dplyr::summarize(
      countIs = n(),
      meanIntensityIs = sum(intensity/samples),
      meanIntensityIsDetected = mean(intensity)
    )  
  
  
  profilesFiltering <-
    peaksTable %>%
    tidyr::expand(profileIDs) %>%
    dplyr::left_join(profilesFiltrate, by = "profileIDs") %>%
    dplyr::left_join(profilesBlank, by = "profileIDs") %>%
    dplyr::left_join(profilesNoIs, by = "profileIDs")
  
  profilesFiltered <-
    profilesFiltering %>%
    filter(countSamples >= minCountSamples) %>%
    mutate(meanIntensityBlank = coalesce(meanIntensityBlank, 1),
           meanIntensityIs = coalesce(meanIntensityIs, 1)) %>%
    mutate(blankRatio = maxIntensitySamples / meanIntensityBlank,
           isRatio = maxIntensitySamples / meanIntensityIs) %>%
    filter(blankRatio/blankRatioFactor + isRatio/isRatioFactor > 1) %>%
    filter(maxIntensitySamples > sampleIntensityCutoff)
  
  profilesFiltered 
}

#' Produce visual plots for verifying the filtering results by hand
#' 
#' From a `hitIntensities` (list of time series matrices for all hits by suspect),
#' filter out the ones which contribute at least a bit to the total profile,
#' plot them, and return a table which shows if they passed filtering or not.
#'
#' @param hitIntensities 
#' @param profilesFiltered List of filtered profiles.
#'
#' @return
#' @export
#'
#' @examples
checkFiltering <- function(hitIntensities, profilesFiltered)
{
  hits <- do.call(rbind, lapply(hitIntensities, function(x) attr(x, "hit")))
  hits %>% group_by(ID) %>% filter(int > 0.3*max(int)) ->hitsFilter
  hitsFilter$flag <- hitsFilter$profileID %in% profilesFiltered$profileIDs
  hitsFilter$RT <- hitsFilter$RT / 60
  
  hitsIntensities <- lapply(
    hitIntensities, function(h) {
      hi <- h[,colnames(h) %in% hitsFilter$profileID,drop=FALSE]
      attr(hi, "hit") <- attr(h, "hit")[colnames(h) %in% hitsFilter$profileID,,drop=FALSE]
      hi
    })
  
  
  l_ply(hitsIntensities, function(intMatrix)
    {
    intSum <- rowSums(intMatrix, na.rm=TRUE)
    intRange <- range(intSum, 0, na.rm = TRUE)
    par(mar=c(3,2,1,1)+0.1)
    plot.new()
    plot.window(xlim=c(0, nrow(intMatrix)), ylim=intRange)
    axis(1)
    axis(2)
    columnSeq <- seq_len(ncol(intMatrix))
    l_ply(columnSeq, function(iCol)
    {
      lines(intMatrix[,iCol], col=iCol)
    })
    lines(intSum, col="black", lty=3)
    hits <- attr(intMatrix, "hit")
    hitsLegend <- paste(hits$Name.y, "(", round(hits$RT / 60, 1), "min)", format(hits$int, scientific = TRUE, digits=1))
    legend("topleft", bty="n", fill=columnSeq, legend=hitsLegend)
  })
  hitsFilter[,c("Name.y", "int", "RT", "mz", "flag")] -> hitsFilterShow
  return(hitsFilterShow)
}

#' Runlength filtering of intensity matrix
#' 
#' Corresponding to Sabine Anlikers workflow, filter out time series with
#' less than n (default n=5) consecutive non-zero datapoints.
#'
#' @param intMatrix Matrix with features in the columns, timepoints in the rows
#' @param minRunlength minimum runlength required to retain the series
#'
#' @return The intensity matrix with only matching columns retained. Find the feature
#'   IDs by `colnames(result)`
#' @export
#'
#' @examples
filterRunlength <- function(intMatrix, minRunlength = 5) {
  profilesRle <- apply(intMatrix > 0, 2, rle)
  profilesMaxLen <- map_dbl(profilesRle, function(profile) max(0, profile$lengths[profile$values == TRUE]))
  intMatrix[,profilesMaxLen > minRunlength]
}



#' Quantile spread filtering
#' 
#' Corresponding to Sabine Anlikers workflow, filter out time series with a Qa to Qb spread
#' of less than r-fold, with default `a=5%, b=95% r=10` i.e. find time series with `Q_{95} / Q_{5} > 10`
#' 
#' Note that this is not working well here when time series have long zero stretches, because this will
#' give `Q_5 = 0` for a very large number of profiles, essentially any episodic profile.
#' Note also that data imputation (e.g. with the cheap MOZER I do) would very directly affect the
#' performance of this.
#' 
#' To try and improve this, we determine the quantiles ignoring the zero values. This is working better,
#' but is also harder to interpret.
#'
#' @param intMatrix 
#' @param q1 
#' @param q2 
#' @param spread 
#'
#' @return
#' @export
#'
#' @examples
filterQuantiles <- function(intMatrix, q1 = 0.05, q2 = 0.95, spread = 10, ignoreZero = TRUE) {
  profilesQQ <- apply(intMatrix, 2, function(x) quantile(x[x>0 | !ignoreZero], probs = c(q1, q2)))
  profilesSpread <- profilesQQ[2,] / (profilesQQ[1,] + 1)
  intMatrix[,profilesSpread > spread]
}