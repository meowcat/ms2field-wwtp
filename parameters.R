library(tibble)
library(tidyr)
library(plyr)

enviPickSettingFiles <- list.files("input" ,"settings-enviPick-(.*).ini$")
enviPickSettings <- suppressWarnings(as_tibble(do.call(rbind,llply(
  enviPickSettingFiles, function(settingsFile)
  {
    re <- "settings-enviPick-(.*).ini"
    settingsName <- sub(re, "\\1", settingsFile)
    if(!is.na(as.numeric(settingsName)))
      settingsName <- paste0("pick", settingsName)
    list(
      pickSettings = settingsName,
      settingsFile = paste0("input/", settingsFile),
      pickedDir = paste0("picked/", settingsName)
    )
  }
))))

enviPickSettings$pickSettings <- rlang::syms(enviPickSettings$pickSettings)

methodsFiltrate <- c("20190204_PAL_filtrate_DDA.meth", "20190204_PAL_filtrate_DIA.meth",
                     "20190204_PAL_filtrate_DDA_IL.meth")


sampleTypes <-  read.csv("input/sampleTypes.csv", stringsAsFactors = FALSE, sep=";")

# # 
# polarities <- tibble(polarity=c("pos"),
#                    scans=c("+")
#                 )
polarities <- tibble(polarity=c("pos", "neg"),
                     scans=c("+", "-")
)
polarities$polarity <- rlang::syms(polarities$polarity)


# this is more clear than data.frame(dret=..., dmass=...)
profileSettings <- as_tibble(rbind(
  list(profileSetting="set1", dret=60, dmass=4)
))


filterProfilesSettings <- list(
  minCountSamples = 8,
  blankRatioFactor = 4,
  isRatioFactor = 2,
  sampleIntensityCutoff = 4e5)

filterProfilesDIASettings <- list(
  minCountSamples = 4,
  blankRatioFactor = 4,
  isRatioFactor = 2,
  sampleIntensityCutoff = 1e5)


screenDDASettings <- list(
  dppm = 5,
  drt = 40
)
#crossfactors <- crossing(enviPickSettings, polarities)

# Create directories for picked files, this is not yet automated
l_ply(enviPickSettings$pickedDir, function(pickedDir)
{
  dir.create(pickedDir, recursive = TRUE, showWarnings = FALSE)
})


# Create directories for dia-picked files, this is not yet automated
# and needs to be solved better
l_ply(enviPickSettings$pickedDir, function(pickedDir)
{
  dir.create(paste0("dia-", pickedDir), recursive = TRUE, showWarnings = FALSE)
})


# DIA scantypes to process in profiling
# c("1-pos--", "2-pos-145-110", "2-pos-355-110", "1-neg--", "2-neg-145-110",  "2-neg-355-110", "2-pos-250-110", "2-pos-705-600", "2-neg-250-110",  "2-neg-705-600")
# we remove the MS1 scans here, since we only want to use the MS2
# and for the time being, only positive mode
.scansDIA <- c(
  # "1-pos--",
  "2-pos-145-110", "2-pos-355-110",
  # "1-neg--", 
  "2-neg-145-110",  "2-neg-355-110", 
  "2-pos-250-110", "2-pos-705-600", 
  "2-neg-250-110",  "2-neg-705-600"
)

scansDIA <- readScantype(.scansDIA)
scansDIA$scanDIA <- .scansDIA

# Midnight of the first day of the measurement, in days after 1970/01/01
baseDay <- 17946
# Frequency set to use for Lomb periodogram
lombFrequencies <- c(0, sort(1/lseq(1/21, 12, 200)))


problematicSamples <- c("20190219_010", "20190219_011", "20190219_046",
                        "20190219_047", "20190219_642",
                        "20190219_931", "20190219_932",
                        "20190219_933", "20190219_934",
                        "20190219_935", "20190219_936",
                        "20190219_937", "20190219_938",
                        "20190219_939", "20190219_940",
                        "20190219_127",   "20190219_253",
                        "20190219_585")

