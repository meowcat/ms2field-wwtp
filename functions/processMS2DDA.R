library(tidyr)
library(purrr)
library(dplyr)
library(mzR)

extractHeadersDDA <- function(sampleList, dirInput, polarity, verbose = FALSE) {
  sampleListDda <- sampleList %>% 
    filter(processingType == "DDA") %>%
    filter(sampleType == "filtrate")
  #return(paste0(dirInput, "/", sampleListDda$file_full))
  sampleListDda$fileDda <- paste0(dirInput, "/", sampleListDda$filename, ".mzXML")
  filesOK <- file.exists(sampleListDda$fileDda)
  if(any(!filesOK))
    stop("extractHeadersDDA: Some raw files are not in the expected location")
  polarityHeader <- ifelse(polarity == "+", 1, 0)
  # Extract all MS2 headers from the specified files - here these are all the files
  # marked DDA and coming from filtrate (e.g. no calibration or blanks)
  ms2 <- sampleListDda$fileDda %>% map(
    function(msfile) {
      d <- openMSfile(msfile)
      h <- header(d)
      h <- h %>% filter(msLevel == 2, polarity == polarityHeader)
      close(d)
      if(verbose) 
        message(msfile)
      h #scans(h$)
    })
  names(ms2) <- sampleListDda$sampleIDs
  ms2 %>% bind_rows(.id = "sampleIDs")
  # note: map_dfr does almost the same
}


getSpectraOneFile <- function(file, headersDDA, scan="seqNum", limit = 0) {
  d <- openMSfile(file)
  s <- peaks(d, headersDDA[[scan]])
  names(s) <- headersDDA$specID
  close(d)
  #return(s)
  s %>% 
    map(as_tibble) %>% 
    map2(headersDDA$specID, ~ .x %>% mutate(specID = .y)) %>%
    bind_rows() %>%
    dplyr::filter(V2 > limit)
}


getSpectraDDA <- function(headersDDA, files, multicore = 5, outfile = "ms2.log", scan="seqNum", limit = 0) {
  filesTab <- tibble(sampleIDs = as.character(seq_along(files)), file_ = files)
  headersDDA <- headersDDA %>% left_join(filesTab, by="sampleIDs")
  
  headersPerFile <- headersDDA %>% group_by(file_) %>% group_split() %>%
    set_names(map(., pull, "file_") %>% map(unique))
  
  cl <- makeCluster(multicore, outfile = outfile)
  clusterEvalQ(cl, library(mzR))
  clusterEvalQ(cl, library(tidyverse))
  clusterExport(cl, "getSpectraOneFile")
  
  spectra <- parLapply(cl, headersPerFile, function(h) {
    file <- unique(h$file_)
    message(file)
    getSpectraOneFile(file, h, scan = scan, limit = limit)
  })
  stopCluster(cl)
  
  spectra <- bind_rows(spectra)
  spectra
}

processSpectraDDA <- function(headersDDA, files, multicore = 5, outfile = "ms2.log",
                              dmass = 10, dret = 20, limit_ms1 = 1e6) {
  headersDDA <- headersDDA %>% dplyr::mutate(specID = row_number())
  spectra <- getSpectraDDA(headersDDA, files, multicore, outfile)
  precursors <- getSpectraDDA(headersDDA, files, multicore, outfile, scan = "precursorScanNum", limit = limit_ms1)
  profilesDDA <- profileSpectraDDA(headersDDA, dmass, dret)
  profilesDDA$spectraData <- spectra %>%
    left_join(profilesDDA$spectraHeaders %>%
                select(specID, specGroupID, sampleIDs),
              by="specID")
  profilesDDA$precursorData <- precursors %>%
    left_join(profilesDDA$spectraHeaders %>%
                select(specID, specGroupID, sampleIDs),
              by="specID")
  profilesDDA
}
  
profileSpectraDDA <- function(headersDDA, dmass = 10, dret = 20) {
  headersDDAlist <- map(c(1,0), ~ headersDDA %>% filter(polarity == .x))
  headersDDAlist <- headersDDAlist %>% purrr::keep(~ nrow(.x) > 0)
  # Fill a profiles-like object with spectra (per scan, i.e. 1 pos 1 neg, as in the standard workflow),
  # then process i.e. compute profiles as usual
  profilesDDA <- headersDDAlist %>% map(fillSpectra)
  profilesDDA <- RMassScreening::computeProfiles(profilesDDA, dmass = dmass, dret = dret)
  # Join the two scans together into a single profiles-like object
  # with "index_prof" and "peaks" as the original object.
  # Note that in this notation 
  # * the profile_ID is the group of spectra that belong together, more aptly renamed to specGroupID
  # * the sampleIDs is a single spectrum! and gets the more appropriate second name specID
  profilesDDA <- list(
    index_prof = profilesDDA %>% 
      map("index_prof") %>% 
      map(as_tibble) %>% 
      bind_rows(.id="scan") %>% 
      dplyr::mutate(specGroupID = paste(scan, profile_ID, sep="-")),
    
    # This "peaks" is included for compatibility! It's quite confusing...
    peaks = profilesDDA %>% 
      map("peaks") %>% 
      map(as_tibble) %>% 
      bind_rows(.id="scan") %>%
      dplyr::mutate(specGroupID = paste(scan, profileIDs, sep="-"),
                    specID = sampleIDs)
  )
  # The two entries of the profiles object are copied to more useful versions
  # labeled according to what is really in there, and without a lot of useless columns.
  # MS2 spectra header information is also added.
  profilesDDA$spectraHeaders <- profilesDDA$peaks %>%
    select(-partitionIDs, -profileIDs, -peakIDs, -sampleIDs) %>%
    left_join(headersDDAlist %>% bind_rows(), by="specID")
  stopifnot(nrow(profilesDDA$spectraHeaders %>% filter(scan.x != scan.y)) == 0)
  profilesDDA$spectraHeaders <- profilesDDA$spectraHeaders %>% dplyr::select(-scan.x) %>% dplyr::rename(scan = scan.y)
  profilesDDA
  
}


screenSpectraDDA <- function(profileIndex, headersDDA, screenDDASettings) {
  
  container <- list(index_prof = profileIndex)
  suspects <- headersDDA %>% mutate(mass = precursorMZ, ret = retentionTime)
  screen <- screenProfiles(container, suspects, "",
                           screenDDASettings$dppm, 
                           screenDDASettings$drt)
  # save(screen, file="results/screenMS2.RData")
  # #screen %>% select(profileID) %>% n_distinct()
  # 
  # ramclustPeaksMS2 <- screen %>% 
  #   group_by(profileID) %>% 
  #   dplyr::summarize(
  #     ms2count = n(),
  #     dppmMin = min(dppm),
  #     drtMin = min(dRT),
  #     precursorIntMax = max(precursorIntensity),
  #     intMax = max(totIonCurrent)
  #   ) %>% 
  #   dplyr::rename(profile_ID = profileID) %>%
  #   right_join(ramclustPeaks)
  # 
  # ramclustPeaksMS2 <- ramclustPeaksMS2 %>%
  #   group_by(component) %>%
  #   dplyr::mutate(compMs2count = sum(ms2count, na.rm = TRUE))
  # 
  # loadd(sampleList)
  
  # screen <- screen %>% 
  #   mutate(filename = sub("(.*)\\.(.*)", "\\1", filename)) %>%
  #   left_join(sampleList %>% select(filename, sampleIDs))
}


fillSpectra <- function (headersDDA) 
{
  
  profiles <- list(0)
  profiles[[1]] <- data.frame(TRUE, FALSE, FALSE, FALSE)
  colnames(profiles[[1]]) <- c("peaks?", "agglom?", "profiled?", 
                               "trends?")
  profiles[[2]] <- 0
  profiles[[3]] <- 0
  profiles[[4]] <- 0
  profiles[[5]] <- 0
  profiles[[6]] <- 0
  profiles[[7]] <- 0
  profiles[[8]] <- 0
  profiles[[9]] <- 0
  names(profiles) <- c("state", "peaks", "datetime", "sampleID", 
                       "place", "index_agglom", "index_prof", "parameters", 
                       "type")
  
  peaks <- matrix(0, nrow=nrow(headersDDA), ncol=8)
  colnames(peaks) <- c("m/z", "intensity", "RT", "peakIDs", 
                       "componentIDs", "sampleIDs", "partitionIDs", "profileIDs")
  peaks[,1] <- headersDDA$precursorMZ
  peaks[,2] <- headersDDA$precursorIntensity
  peaks[,3] <- headersDDA$retentionTime
  peaks[,4] <- 1
  peaks[,6] <- headersDDA$specID
  peaks <- peaks[order(peaks[, 1], decreasing = FALSE), ]
  profiles[[2]] <- peaks
  datetime <- as.POSIXct.numeric(seq_len(nrow(headersDDA)), origin = "1960-01-01")
  profiles[[3]] <- datetime
  
  return(profiles)
}

# Convenience functions to extract a spectrum or spectrum group from profilesDDA

getSpectrum <- function(profilesDDA, specID_, level = 2) {
  if(level == 2)
    df <- profilesDDA$spectraData %>% dplyr::filter(specID == specID_)
  else
    df <- profilesDDA$precursorData %>% dplyr::filter(specID == specID_)
  h <- profilesDDA$spectraHeaders %>% dplyr::filter(specID == specID_)
  attr(df, "header") <- h
  df
}

getSpecGroup <- function(profilesDDA, specGroup) {
  df <- profilesDDA$spectraData %>% dplyr::filter(specGroupID == specGroup)
  h <- profilesDDA$spectraHeaders %>% dplyr::filter(specGroupID == specGroup)
  attr(df, "header") <- h
  df
}

