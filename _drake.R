source("functions.R")
source("environment.R")
source("parameters.R")
source("targets.R")
library(drake)
library(dplyr)
library(purrr)
library(tibble)
library(yaml)



plan <- drake_plan(
  # peak picking
  files = list.files(file_in(dirInput), full.names = TRUE),
  settings = target(yaml.load_file(file_in(settingsFile)),
                    transform=map(
                      .data=!!enviPickSettings, .id=pickSettings
                    )),
  
  pickPeaks = target(batchPickClean(files, 
                                    file_out(pickedDir), 
                                    multicore = 10, 
                                    settings = settings$enviPick),
                     transform=map(pickSettings, .id=pickSettings)),
  
  # sample list assignment
  sampleList = readSampleList(file_in("input/sampleList.csv"), files),
  
  # scan types list
  scans = target(scans,
                 transform = map(.data=!!polarities, .id = polarity)),
  
  # profiling, for pos and neg separately
  profilesRaw = target(
    fillProfiles(file_in(pickedDir), files, polarity=scans),
    #transform = map(.data=!!crossfactors, .id=c(polarity, pickSettings))
    transform = cross(pickPeaks, scans, .id=c(polarity, pickSettings))
  ),
  
  
  profileSettings = target(profileSettings, 
                           transform = map(.data=!!profileSettings, .id = profileSetting)),
  
  profilesComputed = target(
    computeProfiles(profilesRaw, dret=dret, dmass=dmass),
    transform=cross(profilesRaw, profileSettings,
                    .id=c(polarity, pickSettings, profileSetting))),

  
  # make a peak table and intensity matrix
  peaksTable = target(
    profilesComputed$peaks %>% as_tibble() %>% left_join(sampleList, by="sampleIDs"),
    transform=map(profilesComputed, .id=c(polarity, pickSettings, profileSetting))),

  intMatrix = target(
    acast(sampleIDs ~ profileIDs, data=peaksTable, 
          subset=.((sampleType == "filtrate")),
          value.var="intensity", fun.aggregate = mean, fill=0 ,drop=FALSE),
    transform = map(peaksTable, .id=c(polarity, pickSettings, profileSetting))
  ),
  
  profilesConsolidated = target(
    profileGroup.upwards(profilesComputed) %>% # merge highly similar profiles
      profilesMerge(peaksTable, sampleList) %>% # recalculate the $peaks table
      profilesRename(), # recalculate the $index_prof table
    transform = map(profilesComputed, peaksTable,
                    .id=c(polarity, pickSettings, profileSetting))
  ),
  
  intMatrixConsolidated = target(
    acast(sampleIDs ~ profileIDs, data=profilesConsolidated$peaks, 
          subset=.((sampleType == "filtrate")),
          value.var="intensity", fun.aggregate = mean, fill=0 ,drop=FALSE),
    transform = map(profilesConsolidated, .id=c(polarity, pickSettings, profileSetting))
  ),

  profileIndexFiltered = target(
    filterProfiles(profilesConsolidated$peaks,
                   filterProfilesSettings) %>%
      left_join(profilesConsolidated$index_prof %>% as_tibble(),                 
                by=c("profileIDs" = "profile_ID")) %>%
      rename(profile_ID = profileIDs),
    transform = map(profilesConsolidated,
                    .id=c(polarity, pickSettings, profileSetting))
  ),
  
  # generate target list and screen
  targetCpdList = target(parseTargets(file_in(targetList), file_in(cosomiDb))),
  screenHits = target(
    screenProfilesWithId(profilesComputed, targetCpdList, polarity = scans, ppmLimit=4,
                         pickSettings=as_string(quote(pickSettings))),
    transform = map(profilesComputed, .id=c(polarity, pickSettings, profileSetting)) 
  ),
  screenHitsCombined = target(
    bind_rows(screenHits, .id=c("id")),
    transform = combine(screenHits)
  ),
  # Pull out all intensity matrix rows that correspond to a match
  hitIntensities = target(
    getHitIntensities(screenHits, intMatrix),
    transform=map(screenHits, intMatrix, .id=c(polarity, pickSettings, profileSetting))
  ),
  
  targetReports = target(
    rmarkdown::render(knitr_in("reportTemplates/targetReports.Rmd"),
                      output_file = paste0("reports/", .id_chr, ".html"),
                      output_dir = "reports",
                      params = list(hitIntensities = hitIntensities)
    ),
    transform = map(hitIntensities, .id=c(polarity, pickSettings, profileSetting))
  ),
  
  # Same target reports also with all the consolidated data
  screenHitsConsolidated = target(
    screenProfilesWithId(profilesConsolidated,
                         targetCpdList, polarity = scans, ppmLimit=4,
                         pickSettings=as_string(quote(pickSettings))),
    transform = map(profilesConsolidated, .id=c(polarity, pickSettings, profileSetting)) 
  ),
  screenHitsConsolidatedCombined = target(
    bind_rows(screenHitsConsolidated, .id=c("id")),
    transform = combine(screenHitsConsolidated)
  ),
  
  # Pull out all intensity matrix rows that correspond to a match
  hitIntensitiesConsolidated = target(
    getHitIntensities(screenHitsConsolidated, intMatrixConsolidated),
    transform=map(screenHitsConsolidated, intMatrixConsolidated, .id=c(polarity, pickSettings, profileSetting))
    #transform=map(screenHits, intMatrix, .id=c(polarity, pickSettings))
  ),
  targetReportsConsolidated = target(
    rmarkdown::render(knitr_in("reportTemplates/targetReports.Rmd"),
                      output_file = paste0("reports/", .id_chr, "-consolidated.html"),
                      output_dir = "reports",
                      params = list(hitIntensities = hitIntensitiesConsolidated)
    ),
    transform = map(hitIntensitiesConsolidated, .id=c(polarity, pickSettings, profileSetting))
  ),
  
  
  # Quick assessment of information about profiles: count of total, multihit and 5+-hit profiles
  profQC = target(
    profilesQC(.id_chr, profilesComputed, settings),
    transform=map(profilesComputed, settings, .id=c(polarity, pickSettings, profileSetting))
  ),
  profSummary = target(
    bind_rows(profQC),
    transform=combine(profQC)
  ),
  
  # Handling MS2: DDA
  headersDDA = target(
    extractHeadersDDA(sampleList, dirInput, polarity=scans, verbose = FALSE),
    transform = map(scans, .id=polarity)),
  
  spectraDDA = target(
    processSpectraDDA(list(headersDDA) %>% 
      set_names(as.character(seq_along(.))) %>%
      bind_rows(.id = "scan"),
      files
      # additional arguments left to default for now
      ),
    transform = combine(headersDDA)),
    
  # Handling MS2: DIA
  pickPeaksDIA = target(batchPickDIAClean(sampleList,
                                    file_in(dirInput),
                                    paste0("dia-", pickedDir), 
                                    multicore = 10, 
                                    settings = settings$enviPick),
                     transform=map(pickSettings, .id=pickSettings)),
  
  scansDIA = target(scanDIA,
                 transform = map(.data=!!scansDIA, .id = scanDIA)),
  
  profilesRawDIA = target(
    getFilesDIA(sampleList) %>%
      fillProfiles(paste0("dia-", pickedDir)
                   , . , pattern=scanDIA),
    #transform = map(.data=!!crossfactors, .id=c(polarity, pickSettings))
    transform = cross(pickPeaksDIA, scansDIA, .id=c(scanDIA, pickSettings))
  ),
  
  
  profilesComputedDIA = target(
    computeProfiles(profilesRawDIA, dret=dret, dmass=dmass) %>%
      replaceSampleIDs(sampleList %>% getSampleIdDIA()),
    transform=cross(profilesRawDIA, profileSettings,
                    .id=c(scanDIA, pickSettings, profileSetting))),
            
  
  # make a peak table and intensity matrix
  peaksTableDIA = target(
    profilesComputedDIA$peaks %>% as_tibble() %>% 
      right_join(sampleList %>% 
                  filter(processingType == "DIA"), 
                by="sampleIDs"),
    transform=map(profilesComputedDIA, .id=c(scanDIA, pickSettings, profileSetting))),
  
  intMatrixDIA = target(
    acast(sampleIDs ~ profileIDs, data=peaksTableDIA, 
          subset=.((sampleType == "filtrate")),
          value.var="intensity", fun.aggregate = mean, fill=0 ,drop=FALSE),
    transform = map(peaksTableDIA, .id=c(scanDIA, pickSettings, profileSetting))
  ),
  
  
  profilesConsolidatedDIA = target(
    profileGroup.upwards(profilesComputedDIA)  %>% # merge highly similar profiles
      profilesMerge(peaksTableDIA, 
                    sampleList %>% filter(processingType == "DIA")) %>% # recalculate the $peaks table
      profilesRename(), # recalculate the $index_prof table,
    transform = map(profilesComputedDIA, peaksTableDIA,
                    .id=c(scanDIA, pickSettings, profileSetting))
  ),
  
  
  intMatrixConsolidatedDIA = target(
    acast(sampleIDs ~ profileIDs, data=profilesConsolidatedDIA$peaks, 
          subset=.((sampleType == "filtrate")),
          value.var="intensity", fun.aggregate = sum, fill=0 ,drop=FALSE),
    transform = map(profilesConsolidatedDIA, .id=c(scanDIA, pickSettings, profileSetting))
  ),
  
  
  profileIndexFilteredDIA = target(
    filterProfiles(profilesConsolidatedDIA$peaks,
                   filterProfilesDIASettings) %>%
      left_join(profilesConsolidatedDIA$index_prof %>%
                  as_tibble(), 
                by=c("profileIDs" = "profile_ID")) %>%
      rename(profile_ID = profileIDs),
    transform = map(peaksTableConsolidatedDIA, profilesConsolidatedDIA,
                    .id=c(scanDIA, pickSettings, profileSetting))
  ),
  
  intMatrixCombinedDIA = target(
    list(intMatrixConsolidatedDIA),
    transform = combine(intMatrixConsolidatedDIA,
                        .by = c(pickSettings, profileSettings),
                        .id = c(pickSettings, profileSettings))),
  
  intMatrixCombinedMS1 = target(
    list(intMatrixConsolidated),
    transform = combine(intMatrixConsolidated,
                        .by = c(pickSettings, profileSettings),
                        .id = c(pickSettings, profileSettings))),
  
  intMatrixAllMS1 = target(
    combineIntMatrix(sampleList %>% 
                       filter(sampleType == "filtrate") %>%
                       filter(!(filename %in% problematicSamples)),
                     c(intMatrixCombinedMS1)),
    transform = map(intMatrixCombinedMS1,
                    .id=c(pickSettings, profileSettings))
  ),
  
  intMatrixAllDIA = target(
    combineIntMatrix(sampleList %>% 
                       filter(processingType == "DIA", 
                              sampleType == "filtrate") %>%
                       filter(!(filename %in% problematicSamples)),
                     c(intMatrixCombinedMS1, intMatrixCombinedDIA)),
    transform = map(intMatrixCombinedMS1, intMatrixCombinedDIA,
                    .id=c(pickSettings, profileSettings))
  ),
  
  profilesListDIA = target(
    list(profilesConsolidatedDIA) %>% map("index_prof"),
    transform = combine(profilesConsolidatedDIA,
                        .by = c(pickSettings, profileSettings),
                        .id = c(pickSettings, profileSettings))),

  profilesListMS1 = target(
    list(profilesConsolidated) %>% map("index_prof"),
    transform = combine(profilesConsolidated ,
                        .by = c(pickSettings, profileSettings),
                        .id = c(pickSettings, profileSettings))),
  
  # Profile index for all MS1 and all DIA scans
  profilesDIA = target(
    bind_rows(c(profilesListMS1, profilesListDIA), .id="scan") %>%
      dplyr::mutate(matrixID = paste(scan, profile_ID, sep="-")),
    transform = map(profilesListMS1, profilesListDIA,
                    .id=c(pickSettings, profileSettings))
  ),
  
  # Profile index for all MS1 scans only
  profilesMS1 = target(
    bind_rows(c(profilesListMS1), .id="scan") %>%
      dplyr::mutate(matrixID = paste(scan, profile_ID, sep="-")),
    transform = map(profilesListMS1,
                    .id=c(pickSettings, profileSettings))
  ),
  
  # Nearest neighbors for MS1 scans only
  profilesNNMS1 = target(
    processDIA(intMatrixAllMS1, profilesMS1, 
               minDataCount = 6,
               sr = 0.5, st = 10,
               drtmax = 15, neighbors = 25,
               minC11 = 8, maxNeqRatio = 0.5,
               multicore = 10),
    transform = map(profilesMS1, intMatrixAllMS1,
                    .id=c(pickSettings, profileSettings))
  ),
  
  componentsMS1 = target(
    postprocessDIA(profilesNNMS1, profilesMS1, maxWeight = 0.5),
    transform = map(profilesNNMS1, profilesMS1,
                    .id=c(pickSettings, profileSettings))
  ),
  
  trace = TRUE
)


if(length(targets) > 0)
{
  config <- drake_config(plan, targets = targets)
} else {
  config <- drake_config(plan)
}

