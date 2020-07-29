# Extract component intensities,
# and scale both feature and component intensities to a maximum of 1
# (not cutting the minimum to 0).

loadd(intMatrixAllMS1_SB1MP3e20_profileSettings_set1)
loadd(componentsMS1_SB1MP3e20_profileSettings_set1)
source("functions.R")

getComponentsMatrix <- function(components, intMatrix) {
  componentsMatrix_ <- components$profiles %>%
    filter(!is.na(component)) %>%
    group_by(component) %>%
    group_split() %>%
    set_names(map(., pull, "component") %>% map(unique)) %>%
    map(pull, "matrixID") %>%
    map( ~ intMatrix[,.x]) %>%
    map(rowSums)
  return(do.call(cbind, componentsMatrix_))
}

componentsMatrix <- getComponentsMatrix(componentsMS1_SB1MP3e20_profileSettings_set1,
                                        intMatrixAllMS1_SB1MP3e20_profileSettings_set1)

rcMatrixNorm <- apply(componentsMatrix, 2, function(col) col/max(col))
intMatrixNorm <- apply(intMatrixAllMS1_SB1MP3e20_profileSettings_set1, 2, function(col) col/max(col))
save(rcMatrixNorm, intMatrixNorm, file = "results/clustering/viewerInput-global-20200109.RData")



######### Below here is QC that isn't needed for the viewer


loadd(targetCpdList)
screenHits  <- screenProfiles(list(index_prof = components$profiles %>% 
                                     mutate(profile_ID = matrixID) %>% 
                                     filter(scan == 1)),
                              targetCpdList, polarity = "+", ppmLimit=4)
hits <- components$profiles %>% filter(matrixID %in% screenHits$profileID)




hitIntensities <- getHitIntensities(screenHits, intMatrixAllDIA_SB1MP3e20_profileSettings_set1)

addComponent <- function(hit) {
  attr(hit, "hit") <- attr(hit, "hit") %>% 
    left_join(components$profiles %>% select(matrixID, component), by = c("profileID" = "matrixID")) %>%
    mutate(profileID = paste(profileID, component))
  hit
}

hitIntensities <- hitIntensities %>% map(addComponent)


  rmarkdown::render(knitr_in("reportTemplates/targetReports.Rmd"),
                    output_file = paste0("reports/componentHits.html"),
                    output_dir = "reports",
                    params = list(hitIntensities = hitIntensities))

  


