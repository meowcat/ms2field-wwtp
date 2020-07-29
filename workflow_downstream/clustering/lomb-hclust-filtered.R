library(magrittr)
library(dplyr)
library(purrr)
library(tidyr)
source("workflowActive/clustering/helpers.R")
load("results/rcLombMatrix-corr-20200110.RData")
rcLombMatrixCorr <- rcLombMatrix
load("results/rcLombMatrix-20200109.RData")
rcLombMatrix[,1:201] <- pmin(
  rcLombMatrix[,1:201],
  rcLombMatrixCorr[,1:201]
)
rcLombMatrixScale <- scale(rcLombMatrix)
rcLombMatrixScale <- rcLombMatrixScale[,-c(202,403)]

seed <- 869853L
set.seed(seed)
rcDist <- dist(rcLombMatrixScale)
hc <- fastcluster::hclust(rcDist, "average")
save(hc, rcDist, file="results/clustering/data-hc-filtered-20200111-01.RData")

range(hc$height)
h <- c(10, 20, 30, 40)
# clusters <- map(h, ~ dynamicTreeCut::cutreeDynamicTree(
#   hc, maxTreeHeight = .x, deepSplit = FALSE, minModuleSize = 2))


clusters <- cutreeHybrid(hc, as.matrix(rcDist))
save(hc, rcDist, clusters, file="results/clustering/data-hc-filtered-20200111-01b.RData")
# This is identical to the "01" save file,  as you can check with some all.equal stuff

hc$cluster <- clusters$labels
#printClusters(clusters, name = "kmeans-filtered", n=10)
exportClusters(hc, NA, "hc-filtered-20200111-01")

# re-export the data with dampened intensities to account for matrix effect
source("workflowActive/clustering/helpers-corr.R")
exportClusters(hc, NA, "hc-filtered-corr-20200113-01")

