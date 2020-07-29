# From the raw DDA files, extract up to 5 MS2 per profile.
# Additionally, for every fragment in a DDA MS2, extract the chromatogram
# of that ion from the neighboring DIA acquisition, ideally giving 
# chromatograms for every single fragment.
# Correspondingly, check the EIC correlation for each fragment
# (Note: low-intensity fragments are lost here, also collision energy
# plays a role. So the filtered MS2 are used only in particular cases)

library(drake)
source("functions.R")
library(reshape2)
library(furrr)
loadd(sampleList)
source("parameters.R")
library(MSMSsim)
loadd(files)
library(RMassBank)

load("results/MS2-20200109.RData")


# find associated DIA files for DDA files
sampleNextDIA <- sampleList %>% 
  dplyr::filter(sampleType == "filtrate") %>%
  mutate(time = as.POSIXct(mtime)) %>%
  arrange(time) %>%
  mutate(DIA_ = if_else(processingType == "DIA", sampleIDs, NA_integer_),
         DIA_preceding = DIA_,
         DIA_following = DIA_) %>%
  fill(DIA_preceding, .direction = "down") %>%
  fill(DIA_following, .direction = "up") %>%
  fill(DIA_following, .direction = "down") %>%
  fill(DIA_preceding, .direction = "up")

# find EICs for a given MS2 spectrum in a DIA file
findEICDIA <- function(ms, h, scansDIA, f, d, p, hDIA_,
                       rtWindow = 15, # seconds per side
                       ppmEic = 5,
                       level = 2) {
  
  precursorEps <- 1 
  # this is a safety thing
  # if for whatever reason the precursor in the headers are rounded differently
  # than the precursors from the scansDIA table
  
  h_ <- h
  precursor <- h_$precursorMZ
  polarity_ <- h_$polarity
  polarity__ <- if_else(polarity_  == 1, "pos", "neg")
  sampleID_ <- h_$sampleIDs
  retentionTime_ <- h_$retentionTime
  
  # find the attached file (i.e. as requested either the preceding
  # or the following DIA acquisition relative to a (DDA) spectrum)
  
  # find the correct scan type to extract: 
  # right polarity, within range, closest center mass
  precursorMZ_ <- scansDIA %>% dplyr::filter(polarity == polarity__,
                                      precursor >= min,
                                      precursor <= max) %>%
    arrange(abs(precursor - center)) %>%
    slice(1) %>%
    pull("center")
  
  hSub <- hDIA_ %>% 
    dplyr::filter(
      (level == 1) | (abs(precursorMZ - precursorMZ_) < precursorEps),
      msLevel == level, 
      polarity == polarity_,
      abs(retentionTime - retentionTime_) < rtWindow) %>%
    mutate(msLevel = 1)
  pSub <- p[hSub$acquisitionNum]
  eics <- lapply(ms$V1, function(mz) {
    eic <- findEIC(d, mz, RMassBank::ppm(mz, ppmEic, p = TRUE), headerCache = hSub, peaksCache = pSub)
    eic <- eic %>% dplyr::mutate(precursorScan = if_else(rep(level == 2, nrow(hSub)), hSub$precursorScanNum, hSub$acquisitionNum))
    eic
  })
  # find all precursorScans that match to an MS2 (DIA) scan from this dataset
  precursorScans <- eics %>% bind_rows() %>% pull(precursorScan) %>% unique()
  hMS1 <- hDIA_ %>% dplyr::filter(
    msLevel == 1,
    acquisitionNum %in% precursorScans)
  pMS1 <- p[hMS1$acquisitionNum]
  eicMS1 <- findEIC(d, precursor, RMassBank::ppm(precursor, ppmEic, p = TRUE),
                    headerCache = hMS1, peaksCache = pMS1)
  eicMS1$scan <- hMS1$acquisitionNum
  eicMS1$precursorScan <- eicMS1$scan
  specList <- c(list(eicMS1), eics)
  #names(specList) <- c("precursor", seq_along(eics))
  #attr(specList, "sampleIdDIA") <- DIAsampleID
  specList
}



# fetchTopSpectra <- function(specGroupID, profilesDDA, n=5) {
#   ## get all spectra for this specGroup
#   specGroup <- getSpecGroup(profilesDDA, specGroupID)
#   # Find the specID for the top-5 intensity spectra from this group
#   h <- attr(specGroup, "header") %>% arrange(desc(precursorIntensity)) %>% slice(seq_len(n))
#   specIDs_ <- (h %>% pull("specID")) #[seq_len(min(n, nrow(h)))]
#   # Pull those 5 spectra from profilesDDA
#   spectraMS2 <- map(specIDs_, ~ profilesDDA %>% getSpectrum(.x))
#   spectraMS1 <- map(specIDs_, ~ profilesDDA %>% getSpectrum(.x, level = 1))
#   return(list(h=h, spectraMS2=spectraMS2, spectraMS1=spectraMS1))
# }


checkTopSpectra <- function(topSpectra,   
                            isolationWindow = 1, precursorWindow = 5) {

  # Find the specID for the top-5 intensity spectra from this group
  h <- topSpectra$h %>% bind_rows()
  spectraMS2 <- topSpectra$MS2
  spectraMS1 <- topSpectra$MS1
  # Make sure those are all similar to each other
  similarity <- outer(
    seq_along(spectraMS2),
    seq_along(spectraMS2),
    Vectorize(function(sp1, sp2) {
      OrgMSSim(spectraMS2[[sp1]] %>% dplyr::select(V1, V2),
               spectraMS2[[sp2]] %>% dplyr::select(V1, V2))$score
    })
  )
  
  mz_ <- h$precursorMZ[[1]]
  dmz_ <- RMassBank::ppm(mz_, precursorWindow, p = TRUE)/2
  
  purity <- spectraMS1 %>% 
    map(function(sp) {
      spIso <- sp %>% dplyr::filter(between(V1, mz_ - isolationWindow/2, mz_ + isolationWindow/2 ))
      spIon <- sp %>% dplyr::filter(between(V1, mz_ - dmz_, mz_ + dmz_))
      sum(spIon$V2) / sum(spIso$V2)
    })
  
  # Tack on the header and a peak index onto each spectrum
  spectraMS2 <- spectraMS2 %>% map( ~ .x %>% dplyr::mutate(peakIndex = row_number()))
  spectraMS1 <- spectraMS1 %>% map( ~ .x %>% dplyr::mutate(peakIndex = row_number()))
  
  spectraMS2 <- map2(spectraMS2, seq_len(nrow(h)), function(sp, n) {
    attr(sp, "header") <- h[n,]
    sp
  })
  spectraMS1 <- map2(spectraMS1, seq_len(nrow(h)), function(sp, n) {
    attr(sp, "header") <- h[n,]
    sp
  })
  
  
  # Build a matrix out of the EICs
  eics2 <- topSpectra$eics2_ %>%
    map( ~ .x %>% bind_rows(.id = "mzIndex") %>%
    acast(precursorScan ~ as.numeric(as.character(mzIndex)), value.var = "intensity",
          drop = FALSE))
  eics1 <- topSpectra$eics1_  %>%
    map( ~ .x %>% bind_rows(.id = "mzIndex") %>%
           acast(precursorScan ~ as.numeric(as.character(mzIndex)), value.var = "intensity",
                 drop = FALSE))
  # Find correlations for fragment EICs to precursor, and attach to data frame
  spectraMS2 <- map2(spectraMS2, eics2, function(spectrum, eics) {
    eicCor_ <- cor(eics[,1], eics[,-1]) %>% replace_na(0)
    spectrum %>% mutate(eicCor = eicCor_ %>% as.vector())
  })
  spectraMS1 <- map2(spectraMS1, eics1, function(spectrum, eics) {
    eicCor_ <- cor(eics[,1], eics[,-1]) %>% replace_na(0)
    spectrum %>% mutate(eicCor = eicCor_ %>% as.vector())
  })
  
  spectrum <- map2(spectraMS2, spectraMS1, 
                  ~ bind_rows(.x %>% mutate(level = 2),
                              .y %>% mutate(level = 1)))
  attr(spectrum, "eics2") <- eics2
  attr(spectrum, "eics1") <- eics1
  attr(spectrum, "header") <- h
  attr(spectrum, "source") <- h$sampleIDs
  attr(spectrum, "purity") <- purity
  attr(spectrum, "similarity") <- similarity
  spectrum
}

# Fast way to fetch spectra circumventing one-by-one pulling
# They end up nested in the columns MS1 and MS2
# (Note: To some degree, I am just replicating the Spectra package here! rformassspectrometry.com)
n <- 5
spectraToFetch <- spectraDDA$spectraHeaders %>% 
  right_join(componentsDDA %>% select(specGroupID, matrixID), by="specGroupID") %>% 
  group_by(matrixID, specGroupID) %>%
  arrange(desc(precursorIntensity)) %>%
  slice(seq_len(n)) %>%
  group_split() 
spectraTF <- spectraToFetch %>% bind_rows()
spectraMS2 <- spectraTF %>% select(-specGroupID, -sampleIDs) %>% #select(matrixID, specID, precursorIntensity) %>% 
  left_join(spectraDDA$spectraData, by="specID") %>%
  group_by(matrixID, specGroupID, specID) %>%
  nest(V1, V2, .key="MS2")
spectraMS1 <- spectraTF %>% select(-specGroupID, -sampleIDs) %>% #select(matrixID, specID, precursorIntensity) %>% 
  left_join(spectraDDA$precursorData, by="specID") %>%
  group_by(matrixID, specGroupID, specID) %>%
  nest(V1, V2, .key="MS1")
spectra <- spectraTF %>% left_join(spectraMS2) %>% left_join(spectraMS1) %>%
  mutate(nMS1 = map_lgl(MS1, is.null), nMS2 = map_lgl(MS2, is.null)) %>%
  filter(!nMS1, !nMS2) %>%
  filter(between(precursorMZ, min(scansDIA$min), max(scansDIA$max)))

#future::plan("sequential")
future::plan("multiprocess", workers = 7)

spectra <- spectra %>% group_by(sampleIDs) %>% #filter(sampleIDs %in% c("810", "811", "813")) %>%
  group_split() %>%
  future_map(function(df) {
    
    fileSel <- "DIA_preceding"
    DIAsampleID <- sampleNextDIA %>% 
      dplyr::filter(sampleIDs == unique(df$sampleIDs)) %>%
      pull(fileSel) # This must be a single row...
    message(DIAsampleID)
    f <- files[[DIAsampleID]]
    
    # open the attached file
    message(f)
    d <- openMSfile(f)
    hDIA_ <- header(d)
    p <- makePeaksCache(d, hDIA_)
    
    
    
    h_ <- df %>% dplyr::mutate(i = row_number()) %>% 
      dplyr::select(-MS1, -MS2) %>% nest(-i, .key = "h")
    eics1_ <- df$MS1 %>% map2(h_$h,  ~ findEICDIA(.x, .y, scansDIA, f, d, p, hDIA_, level=1) )
    eics2_ <- df$MS2 %>% map2(h_$h,  ~ findEICDIA(.x, .y, scansDIA, f, d, p, hDIA_) )
    
    df$h <- h_$h
    df$eics1_ <- eics1_
    df$eics2_ <- eics2_
    df
  }) %>% bind_rows()


spectra <- spectra %>% 
  group_by(specGroupID, matrixID) %>% 
  group_split() %>%
  future_map(checkTopSpectra)


save(spectra, file="results/MS2-extracted-20200110-02.RData")

