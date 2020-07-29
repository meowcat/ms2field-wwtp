# Matrix effect correction by running mean with internal standard
# Note: This produces input for clustering/helpers-corr.R.
# The final data for clustering is composed from the corrected AND uncorrected data,
# using the lower-intensity frequency component each (see clustering/lomb-hclust-filtered.R)
# However, in the viewer itself and in the plots, only the uncorrected data is shown
# since there is no obvious way of "mixing" the data in the time (not frequency) domain.



library(drake)
source("functions.R")
loadd(intMatrixAllMS1_SB1MP3e20_profileSettings_set1)
loadd(componentsMS1_SB1MP3e20_profileSettings_set1)
loadd(profilesMS1_SB1MP3e20_profileSettings_set1)
intMatrix <- intMatrixAllMS1_SB1MP3e20_profileSettings_set1

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

# These features were identified to be internal standards with screening and manual revision:

istdPos <- c("1-2787", "1-3017", "1-3344", "1-5523", "1-12733", "1-4639", "1-2759", "1-2483", "1-3934",
             "1-13284", "1-9533", "1-2540", "1-11024", "1-44682", "1-8006", "1-2857", "1-5195", "1-6199",
             "1-1506", "1-3677", "1-4109", "1-5610", "1-868", "1-4788", "1-8901", "1-5349", "1-7770",
             "1-2236", "1-2186", "1-3401", "1-14175", "1-3313", "1-7357", "1-4844", "1-2271", "1-10935",
             "1-931", "1-5590", "1-1349", "1-9485", "1-8236", "1-1930", "1-2173", "1-2596", "1-2101", "1-6627",
             "1-2912", "1-7682", "1-2265", "1-2903", "1-117", "1-12082")
istdNeg <- c("2-1193", "2-4998", "2-9254", "2-12815", "2-4858", "2-3370", "2-6314", "2-4283", "2-2246")


intensityCorrectionMatrix_pos <- intMatrix[, istdPos]
intensityCorrectionMatrix_neg <- intMatrix[, istdNeg]

intensityCorrectionMatrix_pos  <- log10(intensityCorrectionMatrix_pos)
intensityCorrectionMatrix_neg <- log10(intensityCorrectionMatrix_neg)

intensityCorrectionMatrix_pos <- apply(intensityCorrectionMatrix_pos, 2,
                                       function(x) x - median(x, na.rm = TRUE))
intensityCorrectionMatrix_neg <- apply(intensityCorrectionMatrix_neg, 2,
                                       function(x) x - median(x, na.rm = TRUE))

intensityCorrection_pos <- apply(intensityCorrectionMatrix_pos, 1, 
                                 function(x) median(x, na.rm = TRUE))
intensityCorrection_neg <- apply(intensityCorrectionMatrix_neg, 1, 
                                 function(x) median(x, na.rm = TRUE))


plot.new()
plot.window(xlim=c(1, 902), ylim=range(intensityCorrectionMatrix_pos, finite=TRUE))
axis(1); axis(2)
apply(intensityCorrectionMatrix_pos, 2, function(x) lines(x))
lines(intensityCorrection_pos, col="red")

plot.new()
plot.window(xlim=c(1, 902), ylim=range(intensityCorrectionMatrix_neg, finite=TRUE))
axis(1); axis(2)
apply(intensityCorrectionMatrix_neg, 2, function(x) lines(x))
lines(intensityCorrection_neg, col="red")



intensityCorrection_pos <- 1/(10 ^ intensityCorrection_pos)
intensityCorrection_neg <- 1/(10 ^ intensityCorrection_neg)

plot.new()
plot.window(xlim=c(1, 902), ylim=range(intensityCorrection_pos,
                                       intensityCorrection_neg,
                                       finite=TRUE))
axis(1); axis(2)
lines(intensityCorrection_pos, col="blue")
lines(intensityCorrection_neg, col="red")
# Should we smooth the negative?
intensityCorrection_neg_smooth <- stats::smooth(intensityCorrection_neg)
lines(intensityCorrection_neg_smooth, col="orange")
intensityCorrection_pos_smooth <- stats::smooth(intensityCorrection_pos)
lines(intensityCorrection_pos_smooth, col="green")

# We know that scan 1 comes first and scan 2 later, but let's be strict
intMatrixColumnScan <- 
  profilesMS1_SB1MP3e20_profileSettings_set1[match(
    colnames(intMatrix), profilesMS1_SB1MP3e20_profileSettings_set1$matrixID
  ), "scan"]
intMatrix[, intMatrixColumnScan == 1] <- 
  intMatrix[, intMatrixColumnScan == 1] * intensityCorrection_pos_smooth
intMatrix[, intMatrixColumnScan == 2] <- 
  intMatrix[, intMatrixColumnScan == 2] * intensityCorrection_neg_smooth



componentsMatrix <- getComponentsMatrix(componentsMS1_SB1MP3e20_profileSettings_set1,
                                        intMatrix)

rcMatrixNorm <- apply(componentsMatrix, 2, function(col) col/max(col))
intMatrixNorm <- apply(intMatrix, 2, function(col) col/max(col))
save(intensityCorrection_pos, intensityCorrection_pos_smooth, 
     intensityCorrection_neg, intensityCorrection_neg_smooth, 
     file="results/intensityCorrection.RData")
save(rcMatrixNorm, intMatrixNorm, intMatrix, file = "results/clustering/viewerInput-corr-20200110.RData")



