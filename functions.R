library(RMassScreening)
library(rlang)
library(plyr)

source("functions/parseTargets.R")
source("functions/addGapSamples.R")
source("functions/blankFilter.R")
source("functions/frequencies.R")
source("functions/processMS2DDA.R")
source("functions/consolidateProfiles.R")
source("functions/processMS2DIA.R")

ihs <- function(x) log(x + sqrt(x^2 + 1))


readSampleList <- function(sampleList, files, sampleType) {
  read.csv(sampleList, stringsAsFactors = FALSE) %>% 
    mutate(method = basename(method)) %>%
    left_join(sampleTypes) %>%
    #mutate(sampleType = sampleTypes[method]) %>%
    assignSamples(files, .)
}



batchPickClean <- function(files, out, ...)
{
  filesInDir <- list.files(out, full.names = TRUE)
  file.remove(filesInDir)
  return(batchPick(files, out, ...))
}

getFilesDIA <- function(sampleList, dir_in = "") {
  if(dir_in != "")
    dir_in <- paste(dir_in, "/")
  sampleList %>%
    filter(processingType == "DIA") %>%
    pull(filename) %>%
    paste0(dir_in, ., ".mzXML")
}

getSampleIdDIA <- function(sampleList) {
  sampleList %>%
    filter(processingType == "DIA") %>%
    select(filename, sampleIDs) %>%
    mutate(newID = row_number())
}

# For better performance, we should avoid re-picking the MS1 scans,
# which would make this quite fast since MS2 pick much faster.
batchPickDIAClean <- function(sampleList, dir_in, dir_out, ...)
{
  filesInDir <- list.files(dir_out, full.names = TRUE)
  file.remove(filesInDir)
  files <- sampleList %>%
    filter(processingType == "DIA") %>%
    pull(filename) %>%
    paste0(dir_in, "/", ., ".mzXML")
  #return(files)
  return(batchPickDIA(files, dir_out, ...))
}


batchPickDiff <- function(files, out, ...)
{
  filesInDir <- list.files(out)
  filesBase <- unique(gsub("(.*).MSlist.(.*)", "\\1", filesInDir))
  filesToPick <- which(!(basename(files) %in% filesBase))
  #return(files[filesToPick])
  return(batchPick(files[filesToPick], out, ...))
}
  
screenProfilesWithId <- function(profiles, suspects, polarity, 
                                 ppmLimit = getOption("RMassScreening")$screenProfiles$ppmLimit, 
                                 rtLimit = getOption("RMassScreening")$screenProfiles$rtLimit,
                                 ...)
{
  df <- screenProfiles(profiles, suspects, polarity, ppmLimit, rtLimit)
  df$polarity <- rep(polarity, nrow(df))
  cols <-dots_list(...)
  #print(str(cols))
  for(col in names(cols))
  {
    if(length(cols[[col]]) != 1)
      stop("Dot arguments must be of length 1 for id columns")
    df[,col] <- rep(as_string(cols[[col]]), nrow(df))
  }
  return(df)
}

getHitIntensities <- function(screenHits, intMatrix)
{
  # find only those hits which have data in the selected samples
  # i.e. in filtrate samples, rather than blanks or calibration samples
  cols <- as.character(screenHits$profileID)
  screenHits$inSamples <- cols %in% colnames(intMatrix)
  screenHits <- screenHits[screenHits$inSamples,,drop=FALSE]
  # return a list of submatrices, each submatrix containing
  # all profiles for all hits on a compound
  cpdSplit <- split(screenHits, screenHits$ID)
  cpdMatrix <- llply(cpdSplit, function(cpdHit)
  {
    mat <- intMatrix[,as.character(cpdHit$profileID),drop=FALSE]
    attr(mat, "hit") <- cpdHit
    mat
  })
  # remove elements with zero hits
  return(cpdMatrix)
}

profilesQC <- function(name, profiles, settings)
{
  profTotal <- nrow(profiles$index_prof)
  profMulti <- sum(profiles$index_prof[,"number_peaks_total"] > 1)
  profFive <- sum(profiles$index_prof[,"number_peaks_total"] > 5)
  minpeak <- settings$enviPick$minpeak
  maxint <-  settings$enviPick$maxint
  SN <- settings$enviPick$SN
  SB <- settings$enviPick$SB
  return(c(
    name=name, 
    minpeak=minpeak, maxint=maxint, SN=SN, SB=SB,
    profTotal=profTotal, profMulti=profMulti, profFive = profFive))
}