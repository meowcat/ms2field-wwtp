# Export collected MS2 spectra (up to 5 per feature) for each frequency cluster
# as MSP, SIRIUS MS format, and MGF for GNPS
# Also export a CSV table per cluster with all features in the cluster.

source("functions/exportMS2.R")
load("results/MS2-extracted-20200110-02.RData")
load("results/clustering/viewerInput-global-screening-20200109.RData")
load("results/clustering/viewerInput-hc-filtered-20200111-01.RData")
load("results/clustering/viewerInput-global-20200109.RData")



library(tidyverse)

spectra2 <- spectra %>% map(function(s) {
  h <- attr(s, "header")
  pu <- attr(s, "purity")
  h$purity <- as.numeric(unlist(pu))
  h$name <- paste0("feature-", h$matrixID, "-spec-", h$specGroupID, "-", 
                   h$sampleIDs, "-", h$acquisitionNum, "@", round(h$retentionTime/60, 2),
                   "-", round(h$purity, 2))
  h$spectrum <- s
  return(h)
})

spectra3 <- spectra2 %>% bind_rows()

correspondence <- ramclustPeaksMS2 %>%
  select(matrixID, component, name, mean_int, mean_RT, mean_mz) %>%
  mutate(cluster = clustCuts[component])

spectra4 <- spectra3 %>% left_join(correspondence, by="matrixID") %>%
  mutate(name = paste0("cl", cluster, "-", name.x, "=", name.y))

# all(colnames(rcMatrixNorm) == as.character(seq_along(colnames(rcMatrixNorm))))
# is TRUE, so we don't need rcMatrixNorm, components are column names

  
spectraMsp <- spectra4 %>% group_by(cluster) %>%
  dplyr::mutate(suspect = (name.y != "")) %>%
  arrange(suspect, desc(mean_int)) %>% 
  group_split() %>%
  map(
    ~ .x %>% rowwise() %>% group_split() %>%
      map(exportMsp)
  )

spectraMsp %>% iwalk(
  ~ writeLines(unlist(.x), con=paste0("results/msp/spectra-cluster", .y, ".msp"))
)


spectraSirius <- spectra4 %>% group_by(cluster) %>%
  dplyr::mutate(suspect = (name.y != "")) %>%
  arrange(suspect, desc(mean_int)) %>% 
  group_split() %>%
  map(
    ~ .x %>% rowwise() %>% group_split() %>%
      map(exportSirius)
  )

spectraSirius %>% iwalk(
  ~ writeLines(unlist(.x), con=paste0("results/sirius/spectra-cluster", .y, ".ms"))
)



spectraSiriusFullMS1 <- spectra4 %>% group_by(cluster) %>%
  dplyr::mutate(suspect = (name.y != "")) %>%
  arrange(suspect, desc(mean_int)) %>% 
  group_split() %>%
  map(
    ~ .x %>% rowwise() %>% group_split() %>%
      map(exportSirius, ms1range = c(-1000,1000))
  )

spectraSiriusFullMS1 %>% iwalk(
  ~ writeLines(unlist(.x), con=paste0("results/sirius/spectra-cluster-fullMS1-", .y, ".ms"))
)


spectraMgf <- spectra4 %>% group_by(cluster) %>%
  dplyr::mutate(suspect = (name.y != "")) %>%
  arrange(suspect, desc(mean_int)) %>% 
  group_split() %>%
  imap(
    ~ c(list(paste0("COM=spectra-cluster-", .y )),
      .x %>% rowwise() %>% group_split() %>%
        map(exportMgf)
    )
  )

spectraMgf %>% iwalk(
  ~ writeLines(unlist(.x), con=paste0("results/mgf/spectra-cluster", .y, ".mgf"))
)



spectraMspFull <- spectra4 %>% group_by(cluster) %>%
  dplyr::mutate(suspect = (name.y != "")) %>%
  arrange(suspect, desc(mean_int)) %>% 
  group_split() %>%
  map(
    ~ .x %>% rowwise() %>% group_split() %>%
      map(exportMsp, cor=-1)
  )

spectraMspFull %>% iwalk(
  ~ writeLines(unlist(.x), con=paste0("results/msp/spectra-cluster-full", .y, ".msp"))
)



ramclustPeaksMS2 %>%
  mutate(cluster = clustCuts[component]) %>%
  group_by(cluster) %>% 
  arrange(mean_mz) %>%
  group_split() %>%
  imap(~write.csv(.x, file=paste0("results/clusterPeaks/cluster",.y,".csv")))
