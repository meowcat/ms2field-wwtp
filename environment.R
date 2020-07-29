# Folder in which the workflow scripts are placed
dirBase <- ""
# Folder in which the raw files are placed 
dirInput <- ""

# Target compounds to be tracked (name in underlying DB)
targetList <- "input/singer-cpds.txt"
# Compound database
cosomiDb <- "input/20190410_Cosomi.csv"

setwd(dirBase)

library(RMassScreening)


