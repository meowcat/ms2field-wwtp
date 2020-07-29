

copySpectrum <- function(spec, file="clipboard", cor = 0.5, type = c("mona", "massbank", "metlin") ) {
  #text_ <- ""
  text <- textConnection("text_", open = "w")
  writeLines(attr(spec, "header")$precursorMZ %>% as.character(), con = text)
  data <- spec %>% dplyr::filter(eicCor > cor, level == 2) %>% select(V1, V2)
  if(type == "mona") {
    writeLines(paste(data$V1, data$V2), con=text)
  }
  else if(type == "massbank") {
    writeLines(paste(data$V1, data$V2 / max(data$V2) * 999), con=text)
  }
  else if(type == "metlin") {
    data <- data %>% 
      arrange(desc(V2)) %>%
      slice(1:30) %>% 
      arrange(V1)
    writeLines(paste(data$V1, data$V2, sep=","), con=text)
  }
  close(text)
  writeLines(text_, con=file)
}



exportSirius <- function(h, cutoffRelative = 0.05, cutoffMs1 = 0.01, cor = 0.5, ms1range = c(-1, 5), file = NULL)
{
  sp <- h$spectrum[[1]]
  buffer <- c()
  
  buffer <- c(buffer, paste0(">compound ", h$name))
  buffer <- c(buffer, paste0(">parentmass ", h$precursorMZ))
  buffer <- c(buffer, paste0(">ionization ", ifelse(h$polarity == 1, "[M+H]+", "[M-H]-")))
  buffer <- c(buffer, "")
  buffer <- c(buffer, paste0(">collision ", h$collisionEnergy))
  
  df <- sp %>% dplyr::filter(level == 2, eicCor > cor)
  table <- paste(df$V1, df$V2, sep=" ")
  buffer <- c(buffer, table)
  buffer <- c(buffer, "")
  buffer <- c(buffer, paste0(">ms1peaks"))
  df <- sp %>% dplyr::filter(level == 1, eicCor > cor, between(V1, h$precursorMZ + ms1range[[1]], h$precursorMZ + ms1range[[2]]))
  table <- paste(df$V1, df$V2, sep=" ")
  buffer <- c(buffer, table)
  if(!is.null(file))
    writeLines(buffer, file)
  return(buffer)
}


exportMsp <- function(h, cutoffRelative = 0.05, cutoffMs1 = 0.01, cor = 0.5, ms1range = c(-1, 5), file = NULL)
{
  sp <- h$spectrum[[1]]
  buffer <- c()
  buffer <- c(buffer, paste0("Name: ", h$name))
  buffer <- c(buffer, paste0("Precursor_type: ", ifelse(h$polarity == 1, "[M+H]+", "[M-H]-")))
  buffer <- c(buffer, paste0("Spectrum_type: MS2"))
  buffer <- c(buffer, paste0("PrecursorMZ: ", h$precursorMZ))
  df <- sp %>% dplyr::filter(level == 2, eicCor > cor)
  buffer <- c(buffer, paste0("Num peaks: ", nrow(df)))
  table <- paste(df$V1, df$V2, sep=" ")
  buffer <- c(buffer, table)
  if(!is.null(file))
    writeLines(buffer, file)
  return(buffer)
}


exportMgf <- function(h, cutoffRelative = 0.05, cutoffMs1 = 0.01, cor = 0.5, ms1range = c(-1, 5), file = NULL) {
  sp <- h$spectrum[[1]]
  buffer <- c()
  buffer <- c(buffer, "BEGIN IONS")
  buffer <- c(buffer, paste0("TITLE=", h$name))
  buffer <- c(buffer, paste0("RTINSECONDS=", h$retentionTime))
  buffer <- c(buffer, paste0("PEPMASS=", h$precursorMZ))
  buffer <- c(buffer, paste0("FEATURE_ID=", h$name))
  df <- sp %>% dplyr::filter(level == 2, eicCor > cor)
  table <- paste(df$V1, round(df$V2, 0), sep=" ")
  buffer <- c(buffer, table)
  buffer <- c(buffer, "END IONS")
  if(!is.null(file))
    writeLines(buffer, file)
  return(buffer)
}

