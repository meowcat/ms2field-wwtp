parseTargets <- function(targets, cosomiDb)
{
  targetCpds <- read.csv(targets, sep="_", stringsAsFactors = FALSE, header=FALSE)
  colnames(targetCpds) <- c("Name", "ID", "polarity", "_pol")
  targetCpds <- targetCpds[,c("Name", "ID", "polarity"),drop=FALSE]
  cosomi <- read.csv(cosomiDb, stringsAsFactors = FALSE)
  targetCpds <- merge(targetCpds, cosomi, by.x="ID", by.y="UCHEM.ID", all.x = TRUE)
  
  targetIsIds <- targetCpds$Own.Labelled.Compound.In.Use.Id
  targetIs <- data.frame(
    "Name" = "---",
    "ID" = targetIsIds,
    "polarity" = targetCpds$polarity)
  targetIs <- merge(targetIs, cosomi, by.x="ID", by.y="UCHEM.ID", all.x = TRUE)
  
  targetAll <- rbind(targetCpds, targetIs)
  # targetsFail <- targets[is.na(targets[,"Name.y"]),,drop=FALSE]
  # # only empty rows
  targetAll$name <- targetAll$name.y
  targetAll$mass <- targetAll$Monoisotopic.Mass
  return(targetAll)
}

