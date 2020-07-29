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

# Create a set of visualizations

# UMAP visualization: cosine
library(uwot)
min_dist = c(0.001,0.005,0.01,0.05,0.1,0.5)

set.seed(seeds[[2]])
umap_cos_vis <- map(
  min_dist,
  ~ uwot::umap(rcLombMatrixScale, metric="cosine", scale="none", min_dist = .x, nn_method="annoy")
) %>% set_names(paste0("umap_cos_", seq_along(min_dist)))

set.seed(seeds[[3]])
umap_n2_vis <- map(
  min_dist,
  ~ uwot::umap(rcLombMatrixScale, metric="euclidean", scale="none", min_dist = .x, nn_method="annoy")
) %>% set_names(paste0("umap_n2_", seq_along(min_dist)))

library(Rtsne)
perplexity = c(5,10,15,20,25,30, 35, 40, 45)
set.seed(seeds[[4]])
tsne_vis <- map(
  perplexity,
  ~ Rtsne(rcLombMatrixScale, perplexity = .x, num_threads = 8)
) %>% map(`$`, "Y") %>% set_names(paste0("tsne_", perplexity))

library(largeVis)
perplexity = c(10, 25, 50, 75, 100)
set.seed(seeds[[5]])
lv_vis <- map(perplexity, ~ largeVis(t(rcLombMatrixScale), perplexity = .x)) %>%
  map(`$`, "coords") %>% map(t) %>% set_names(paste0("largeVis_", perplexity))

vis <- c(umap_cos_vis, umap_n2_vis, tsne_vis, lv_vis) 

save(vis, file="results/clustering/visualizations-filtered-20200110-01.RData")
