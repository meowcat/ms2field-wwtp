# Screen feature list for in-house target compounds and internal standards
# Produce a result datafile for the shinyViewer
# (and the basis for internal standard-based matrix effect correction)

library(drake)
library(MASS)
library(RMassScreening)
library(tidyverse)
source("functions.R")
loadd(intMatrixAllMS1_SB1MP3e20_profileSettings_set1)
loadd(componentsMS1_SB1MP3e20_profileSettings_set1)


# Feature table: rename columns such that RMassScreening can work with it,
# split into positive and negative mode data

screenProfiles <- componentsMS1_SB1MP3e20_profileSettings_set1$profiles %>%
  mutate(profile_ID_orig = profile_ID, 
         profile_ID = matrixID)

screenProfilesPos <- componentsMS1_SB1MP3e20_profileSettings_set1$profiles %>%
  filter(scan == 1) %>%
  mutate(profile_ID = matrixID)

screenProfilesNeg <- componentsMS1_SB1MP3e20_profileSettings_set1$profiles %>%
  filter(scan == 2) %>%
  mutate(profile_ID = matrixID)


# Load compound database, which is split into compounds and ions, and separate accordingly,
# also split into positive and negative

traceFinderRaw <- readLines("input/20190409_EawagMasterDB_RTAtlantis.csv")
traceFinderCpdType <- read.csv(stringsAsFactors = FALSE,
                               text = traceFinderRaw[6261:7722]) %>%
  dplyr::select(name = Compound.Name, category = Category, compoundType = Compound.Type, mix = Compound.Group,
         formula = Compound.Formula)

traceFinderDb <- read.csv(stringsAsFactors=FALSE,
                          text = traceFinderRaw[6:6256]) %>%
  dplyr::filter(Workflow == "TargetPeak", !is.na(Retention.Time)) %>%
  dplyr::mutate(rtRankDb = rank(Retention.Time)) %>%
  dplyr::select(mass = m.z, 
                name = Compound.Name,
                adduct = Adduct, chargeState = Charge.State,
                retAtlantis = Retention.Time, 
                rtRankDb,
                PeakPolarity) %>%
  inner_join(traceFinderCpdType) %>% filter(!is.na(mass)) %>%
  dplyr::mutate(id = row_number())

dmass <- outer(traceFinderDb$mass, traceFinderDb$mass, `-`) %>% abs
isobaric <- which(dmass < 0.001, arr.ind = TRUE) %>% 
  as_tibble %>% 
  filter(row > col) %>% 
  as.matrix() %>% as.vector() %>% 
  unique()

traceFinderDb <- traceFinderDb %>% dplyr::mutate(isIsobaric = row_number() %in% isobaric)

traceFinderDbPos <- traceFinderDb %>% filter(PeakPolarity == "Positive")
traceFinderDbNeg <- traceFinderDb %>% filter(PeakPolarity == "Negative")

# Screen for internal standards, plot retention time rank fit

screenIstd <- screenProfiles(list(index_prof=screenProfilesPos),
                             traceFinderDbPos %>% dplyr::filter(compoundType == "eInternalStandard"),
                             polarity = "",
                             ppmLimit = 3.5) %>%
  dplyr::group_by(id) %>% dplyr::mutate(relint = int / max(int)) %>%
  left_join(screenProfiles %>% as.data.frame() %>%
              dplyr::mutate(rtRank = rank(mean_RT)) %>%
              dplyr::select(profileID = profile_ID, rtRank))

rtFit <- ggplot(screenIstd %>% filter(relint > 0.99), aes(x=rtRankDb, y=rtRank)) + 
  geom_point(aes(col = abs(dppm), shape = factor(isIsobaric))) +
  scale_color_gradient(low = "blue", high = "red")

rtRankFit <- screenIstd %>% filter(rtRankDb < 750, relint > 0.99) %>%
  rlm(rtRank ~ rtRankDb,
      data = .,
      w = exp(-abs(dppm - median(dppm)))
  )
predictRange <- 0:800
rtRankLine <- data.frame(rtRankDb = predictRange) %>%
  mutate(rtRank = predict(rtRankFit, newdata=.))
rtFit + geom_line(data = rtRankLine)


screenIstd <- screenIstd %>%
  ungroup() %>%
  mutate(rtRankPred = predict(rtRankFit, newdata=.)) %>%
  mutate(dRtRank = rtRank - rtRankPred )

screenIstd$ID <- screenIstd$id
screenIstd$Name.y <- screenIstd$name
hitIntensities <- getHitIntensities(screenIstd, intMatrixAllMS1_SB1MP3e20_profileSettings_set1)
write.csv(screenIstd, file="results/screenIstd_pos.csv")
rmarkdown::render("reportTemplates/targetReports.Rmd",
                  output_file = paste0("results/targets-istd-pos.html"),
                  output_dir = "results",
                  params = list(hitIntensities = hitIntensities))



screenIstdNeg <- screenProfiles(list(index_prof=screenProfilesNeg),
                                traceFinderDbNeg %>% dplyr::filter(compoundType == "eInternalStandard"),
                                polarity = "",
                                ppmLimit = 4.5) %>%
  dplyr::group_by(id) %>% dplyr::mutate(relint = int / max(int)) %>%
  left_join(screenProfiles %>% as.data.frame() %>%
              dplyr::mutate(rtRank = rank(mean_RT)) %>%
              dplyr::select(profileID = profile_ID, rtRank))

screenIstdNeg$ID <- screenIstdNeg$id
screenIstdNeg$Name.y <- screenIstdNeg$name
write.csv(screenIstdNeg, file="results/screenIstd_neg.csv")
hitIntensities <- getHitIntensities(screenIstdNeg, intMatrixAllMS1_SB1MP3e20_profileSettings_set1)
rmarkdown::render("reportTemplates/targetReports.Rmd",
                  output_file = paste0("results/targets-istd-neg.html"),
                  output_dir = "results",
                  params = list(hitIntensities = hitIntensities))


ggplot(screenIstd %>% filter(relint > 0.99), aes(x=rtRankDb, y=rtRank)) + 
  geom_point(aes(color = factor(abs(dRtRank) < 5e4),
                 shape = factor(isIsobaric)))



# The internal standards were screened on the entire set of profiles,
# the suspects will only be screened on profiles that are in samples (as opposed to
# spikes, calibration etc.)
# Note that the rt rank is still coming from the entire profile set (so it is comparable to the ISTD data)

screenParams <- tibble(
  data = list(screenProfilesPos, screenProfilesNeg),
  db = list(traceFinderDbPos, traceFinderDbNeg)
)

screenSuspects <- screenParams %>% pmap_dfr(function(data, db) {
  screenProfiles(list(index_prof = data),
                 db %>% dplyr::filter(compoundType == "eTargetCompound"),
                 polarity = "",
                 ppmLimit = 3.5) %>%
    dplyr::group_by(id) %>% dplyr::mutate(relint = int / max(int)) %>%
    left_join(data %>% as.data.frame() %>%
                dplyr::mutate(rtRank = rank(mean_RT)) %>%
                dplyr::select(profileID = matrixID, rtRank)) %>%
    rename(matrixID = profileID)
})

suspectPlot <- ggplot(screenSuspects %>% filter(relint > 0.99), aes(x=rtRankDb, y=rtRank)) + 
  geom_point(aes(col = abs(dppm), shape = factor(isIsobaric))) +
  scale_color_gradient(low = "blue", high = "red")
suspectPlot
screenSuspects <- screenSuspects %>%
  ungroup() %>%
  mutate(rtRankPred = predict(rtRankFit, newdata=.)) %>%
  mutate(dRtRank = rtRank - rtRankPred )


ggplot(screenSuspects %>% filter(relint > 0.99, !isIsobaric), aes(x=rtRankPred, y=rtRank)) + 
  geom_point(aes(color = factor(abs(dRtRank) < 6e4),
                 shape = factor(isIsobaric),
                 size = int
  ))

suspectsNarrow <- screenSuspects %>% 
  filter(relint > 0.99, !isIsobaric, abs(dRtRank) < 6e4) %>%
  dplyr::mutate(dppmRecal = dppm - median(dppm)) %>%
  filter(abs(dppmRecal) < 1) %>%
  dplyr::select(name, rtRank, rtRankPred, dRtRank, dppm, 
                matrixID, formula, category, dppmRecal)



ramclustPeaksMS2 <- screenProfiles %>% 
  mutate(matrixID = profile_ID,
         profile_ID = profile_ID_orig) %>%
  left_join(suspectsNarrow)

save(ramclustPeaksMS2, file="results/clustering/viewerInput-global-screening-20200109.RData")


