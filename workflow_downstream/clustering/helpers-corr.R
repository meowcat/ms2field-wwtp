library(drake)
library(tidyverse)
load("results/clustering/viewerInput-corr-20200110.RData")
loadd(sampleList)
source("functions.R")
source("parameters.R")

# some seeds for random generators to use before non-deterministic functions
# > cat(deparse(sample.int(1e6, 10)))
#seeds <- c(982765L, 346680L, 932367L, 57356L, 100155L, 898623L, 811843L,  372530L, 610931L, 659153L)
seeds <- c(600723L, 155268L, 914927L, 417585L, 920786L, 163209L, 521476L,  113475L, 250119L, 869853L)

sampleListOK <- sampleList %>% filter(!(filename %in% problematicSamples))


printClusters <- function(clusters, name, n) {
  
  sampleProfiles <- clusters$cluster %>% 
    map(~ split(clusters$labels, .x)) %>% 
    map_depth(2, sample, size=n) %>%
    map_depth(2, ~ rcMatrixNorm[,.x])
  
  clustersCount <- map(clusters$cluster, table)
  
  
  for(s in seq_along(clusters$cluster)) {
    pdf(paste0("results/sampleProfiles/", name, "-", names(sampleProfiles)[[s]], ".pdf"))
    for(j in seq_along(sampleProfiles[[s]])) {
      subsample <- sampleProfiles[[s]][[j]]
      plot.new()
      plot.window(xlim = range(0,902), ylim=range(0,1))
      axis(1)
      axis(2)
      for(i in seq_len(ncol(subsample)))
        lines(subsample[,i], col=i)
      title(main=paste0("hmax = ", names(sampleProfiles)[[s]], ", cluster ", j-1))
      legend("topleft", bty='n', legend=clustersCount[[s]][[j]])
    }
    dev.off()
  }
  
}



exportClusters <- function(clusters, index = NA, name) {
  
  if(!is.na(index))
    clustCuts <- clusters$cluster[[index]]
  else
    clustCuts <- clusters$cluster
  
  rcClustCuts <- rcMatrixNorm %>% t() %>% as.data.frame() %>% 
    split(clustCuts) %>% map(as.matrix) %>% 
    map(function(profileCluster) {
      clusterMean <- apply(profileCluster, 2, mean, na.rm=TRUE)
      clusterSd <- apply(profileCluster, 2, sd, na.rm=TRUE)
      cbind(clusterMean, clusterSd, clusterMean + clusterSd, clusterMean - clusterSd)
    })
  
  save(clustCuts, rcClustCuts, file=paste0("results/clustering/viewerInput-", name, ".RData"))
}



