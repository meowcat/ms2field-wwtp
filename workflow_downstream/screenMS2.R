# Find matches between the MS1 profiles (components$profiles) and the MS2 groups (spectraDDA).
# MS1 profiles indexed by components$profiles$matrixID
# MS2 groups indexed by spectraDDA$index_prof$specGroupID

library(drake)
library(MASS)
library(RMassScreening)
library(tidyverse)
source("functions.R")
loadd(componentsMS1_SB1MP3e20_profileSettings_set1)
loadd(spectraDDA)



#load("results/clustering/viewerInput-global-screening.RData")
componentsDDA <- map(c(1,2),
                     ~ screenProfiles(
                       list(index_prof = componentsMS1_SB1MP3e20_profileSettings_set1$profiles %>%
                              filter(scan==.x) %>%
                              select(-profile_ID) %>%
                              rename(profile_ID = matrixID)),
                       
                       spectraDDA$index_prof %>%
                         filter(scan==as.character(.x)) %>%
                         select(mean_mz, mean_RT, specGroupID,
                                scan,mean_int) %>%
                         mutate(mass = mean_mz,
                                ret = mean_RT),
                       polarity = "",
                       ppmLimit = 4, rtLimit = 15)
) %>% bind_rows() %>% rename(matrixID = profileID) %>% select(-mass, -ret)

save(componentsDDA, file="results/MS2-20200109.RData")
