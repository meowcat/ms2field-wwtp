# consolidate "similar" profiles

# from Jen, 19.9.2019
## grouping profiles which should belong together
## applied after normal profiling to combat "peak splitting" and poor intensity repeatability

profileGroup <- function(profiles, dppm = 2, drt = 15){
  
  l <- nrow(profiles$index_prof)
  
  #peaklist <- as.data.frame(profiles$index_prof[,c("mean_mz", "mean_RT", "max_int", "number_peaks_total")])
  #peaklist$grouped <- NA
  peaklist <- profiles$index_prof %>% 
    as_tibble() %>%
    mutate(grouped = NA)
  
  peaklist <- peaklist %>% arrange(desc(max_int))
  # # Compute in matrix form and hope this is faster :)
  # peaklist <- profiles$index_prof[,c("profile_ID", "mean_mz", "mean_RT",
  #                                    "max_int", "number_peaks_total")]
  # peaklist <- cbind(peaklist, rep(NA, l))
  #                     
  # rownames(peaklist) <- profiles$index_prof[,"profile_ID"]
  # colnames(peaklist) <- c("profile_ID", "mean_mz", "mean_RT",
  #                         "max_int", "number_peaks_total", "grouped")
  # peaklist <- peaklist[order(peaklist[,"max_int"], decreasing = TRUE),]
  #n.peaks <- peaklist[,4]
  # peaklist[,1] <- as.numeric(as.character(peaklist[,1]))
  # peaklist[,2] <- as.numeric(as.character(peaklist[,2]))
  # peaklist[,3] <- as.numeric(as.character(peaklist[,3]))
  # peaklist[,4] <- as.numeric(as.character(peaklist[,4]))
  
  groupID <- list()
  
  pb <- txtProgressBar(min = 0, max = l)
  
  for(i in seq_len(l)){
    
    #print(i)
    setTxtProgressBar(pb,i)
    
    ## checks if profile already in a group
    if(!is.na(peaklist$grouped[i])){next} 
    
    ## no grouping if profile only has 1 peak in it
    if(peaklist$number_peaks_total[i] == 1){
      
      groupID[[as.character(i)]] <- peaklist$profile_ID[i]
      # todo: yes or no?
      # peaklist[c(i, take), "grouped"] <- i
      next
    } 
    
    mz.min <- peaklist$mean_mz[i] - 
      RMassScreening::ppm(peaklist$mean_mz[i], dppm)
    mz.max <- peaklist$mean_mz[i] + 
      RMassScreening::ppm(peaklist$mean_mz[i], dppm)
    
    rt.min <- peaklist$mean_RT[i] - drt
    rt.max <- peaklist$mean_RT[i] + drt
    
    take <- peaklist %>% 
      filter(mean_mz < mz.max,
             mean_mz > mz.min) %>%
      filter(mean_RT < rt.max,
             mean_RT > rt.min) %>%
      filter(is.na(grouped)) %>%
      pull(profile_ID)
    
    #   
    #   which(peaklist[,"mean_mz"] < mz.max & 
    #                  peaklist[,"mean_mz"] > mz.min)
    # 
    # 
    # take2 <- which(peaklist[,"mean_RT"] < rt.max & 
    #                  peaklist[,"mean_RT"] > rt.min)
    # 
    # take <- intersect(take1, take2)
    # 
    
    # Sanity check: the profile should be in its own match set
    stopifnot(peaklist$profile_ID[[i]] %in% take)
    groupID[[as.character(i)]] <- take #peaklist[take,"profile_ID"]
    
    peaklist <- peaklist %>% mutate(grouped = ifelse(profile_ID %in% take, i, grouped))
    # peaklist$grouped[c(take)] <- i
    # if(length(take) == 1){next}else{
    #   
    #   #peaklist[i,5] <- paste(rownames(peaklist[take,]), collapse = ';')
    #   ## indicates profileIDs to be grouped
    #   
    #   peaklist[take,5] <- "x" ## marks profiles which are already in a group
    #   
    # }
    
  }
  
  close(pb)
  return(groupID)
  
}



profileGroup.upwards <- function(profiles, dppm = 2, drt = 25){
  
  l <- nrow(profiles$index_prof)
  
  #peaklist <- as.data.frame(profiles$index_prof[,c("mean_mz", "mean_RT", "max_int", "number_peaks_total")])
  #peaklist$grouped <- NA
  peaklist <- profiles$index_prof %>% 
    as_tibble() %>%
    mutate(grouped = NA)
  
  peaklist <- peaklist %>% arrange(desc(max_int))
  # # Compute in matrix form and hope this is faster :)
  # peaklist <- profiles$index_prof[,c("profile_ID", "mean_mz", "mean_RT",
  #                                    "max_int", "number_peaks_total")]
  # peaklist <- cbind(peaklist, rep(NA, l))
  #                     
  # rownames(peaklist) <- profiles$index_prof[,"profile_ID"]
  # colnames(peaklist) <- c("profile_ID", "mean_mz", "mean_RT",
  #                         "max_int", "number_peaks_total", "grouped")
  # peaklist <- peaklist[order(peaklist[,"max_int"], decreasing = TRUE),]
  #n.peaks <- peaklist[,4]
  # peaklist[,1] <- as.numeric(as.character(peaklist[,1]))
  # peaklist[,2] <- as.numeric(as.character(peaklist[,2]))
  # peaklist[,3] <- as.numeric(as.character(peaklist[,3]))
  # peaklist[,4] <- as.numeric(as.character(peaklist[,4]))
  
  groupID <- list()
  
  pb <- txtProgressBar(min = 0, max = l)
  
  j <- 1
  ambiguousMerges <- 0
  
  for(i in seq_len(l)){
    
    # Determine boundaries for merging profiles
    mz.min <- peaklist$mean_mz[i] - 
      RMassScreening::ppm(peaklist$mean_mz[i], dppm)
    mz.max <- peaklist$mean_mz[i] + 
      RMassScreening::ppm(peaklist$mean_mz[i], dppm)
    
    rt.min <- peaklist$mean_RT[i] - drt
    rt.max <- peaklist$mean_RT[i] + drt
    
    # Try to find if there is a group this peak can join
    hitGroups <- peaklist %>% 
      filter(!is.na(grouped)) %>%
      filter(mean_mz < mz.max,
             mean_mz > mz.min) %>%
      filter(mean_RT < rt.max,
             mean_RT > rt.min)

    nGroupHits <- (hitGroups %>% pull(grouped) %>% n_distinct()) 
    # Multiple groups: join the closest one in RT
    if(nGroupHits > 1) {
      hitGroup <- hitGroups %>%
        mutate(drt = abs(mean_RT - peaklist$mean_RT[i])) %>%
        arrange(drt) %>%
        slice(1) %>%
        pull(grouped)
      ambiguousMerges <- ambiguousMerges + 1  
    }
    else if(nGroupHits == 1)
      hitGroup <- hitGroups %>% slice(1) %>% pull(grouped)
    else {
      hitGroup <- j
      j <- j + 1
    }
    peaklist$grouped[i] <- hitGroup

    #print(i)
    setTxtProgressBar(pb,i)
    
  }
  
  close(pb)
  profiles$grouped <- peaklist
  profiles$ambiguousMerges <- ambiguousMerges
  
  return(profiles)
  
}



#' Consolidate profiles: apply grouping
#' 
#' 
#'
#' @param profiles 
#'
#' @return
#' @export
#'
#' @examples
profilesMerge <- function(profiles, peaksTable, sampleList) {
  peaksTable <- peaksTable %>%
    left_join(profiles$grouped %>% 
                select(profile_ID, grouped), 
              by=c("profileIDs" = "profile_ID"))
  peaksSummary <- peaksTable %>% 
    dplyr::group_by(sampleIDs, grouped) %>%
    dplyr::summarize(RT = weighted.mean(RT, intensity),
              `m/z` = weighted.mean(`m/z`, intensity),
              intensity = sum(intensity)) %>%
    left_join(sampleList, by="sampleIDs") %>%
    dplyr::rename(profileIDs = grouped) %>%
    dplyr::ungroup()
  profiles$peaks <- peaksSummary
  return(profiles)
}

profilesRename <- function(profiles) {
  profiles$index_prof <- profiles$grouped %>% 
    dplyr::select(-profile_ID) %>%
    dplyr::rename(profile_ID = grouped) %>%
    dplyr::group_by(profile_ID) %>% 
    dplyr::arrange(desc(mean_int)) %>%
    dplyr::mutate(mean_int = sum(mean_int)) %>%
    slice(1)
  return(profiles)
}


consolidationResults <- function(profiles) {
  cnt <- profiles$grouped %>%
    dplyr::group_by(grouped) %>% dplyr::summarize(count = n()) %>% pull(count)
  return(list(
    fivenum = fivenum(cnt),
    profilesCount = c(length(cnt), nrow(profiles$grouped))
  ))
}


replaceSampleIDs <- function(profiles, correction) {
  # Note: for compatibility reasons, we write into the original data format
  # and do not already convert into tibbles...
  stopifnot(all(correction$newID == seq_along(correction$newID)))
  profiles$peaks[,"sampleIDs"] <- 
    correction$sampleIDs[profiles$peaks[,"sampleIDs"]]
  return(profiles)
}