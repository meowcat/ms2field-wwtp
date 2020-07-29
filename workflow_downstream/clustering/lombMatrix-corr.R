library(magrittr)
library(dplyr)
library(purrr)
library(tidyr)
source("workflowActive/clustering/helpers-corr.R")



# Midnight of the first day of the measurement, in days after 1970/01/01
baseDay <- 17946

# sampleList %>% filter(sampleType == "filtrate") %>% pull(mtime)  %>% parse_datetime() %>% plot()
# This shows that there are only mini-interruptions here

# Frequency set to use for Lomb periodogram
lombFrequencies <- c(0, sort(1/lseq(1/21, 12, 200)))

# Build periodograms for every RAMClust component
# (the trace used is the highest-intensity profile)
rcPeriodograms <- lombPeriodograms(rcMatrixNorm, 
                                   sampleListOK, 
                                   baseTime = baseDay,
                                   f = lombFrequencies)
rcLombMatrix <- lombMatrix(rcPeriodograms)
rcLombMatrixScale <- scale(rcLombMatrix)
rcLombMatrixScale <- rcLombMatrixScale[,-c(202,403)]
save(rcLombMatrix, rcLombMatrixScale, file="results/rcLombMatrix-corr-20200110.RData")


# 
# pca_vis <- princomp(rcLombMatrixScale)
# plot(pca_vis$scores[,c(1,2)], col=cluster_col(km$cluster), pch = cluster_pch(km$cluster))
# plot(pca_vis$scores[,c(1,3)], col=cluster_col(km$cluster), pch = cluster_pch(km$cluster))
# plot(pca_vis$scores[,c(2,3)], col=cluster_col(km$cluster), pch = cluster_pch(km$cluster))
# pairs(pca_vis$scores[,1:6],  col=cluster_col(km$cluster), pch = cluster_pch(km$cluster))
# 
# 
# ica_vis <- ica::icafast(rcLombMatrixScale, nc=6)
# pairs(ica_vis$S,  col=cluster_col(km$cluster), pch = cluster_pch(km$cluster))
# 
# plot(c(1:50) %% 10, col=cluster_col(c(1:50)), pch = cluster_pch(c(1:50)))
# 
# 
# save(tsne_vis, umap_cos_vis, umap_n2_vis, lv_vis, pca_vis, ica_vis, vis, file="results/clustering/visualizations.RData")
