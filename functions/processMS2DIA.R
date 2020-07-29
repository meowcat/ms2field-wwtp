library(readr)
#library(fastliclust)
library(tidyverse)
library(weights)
library(igraph)

combineIntMatrix <- function(sampleList, ...) {
  # get all DIA sample IDs in the sample list as characters for row selection
  rows <- sampleList %>%
    arrange(parse_datetime(mtime)) %>%
    pull(sampleIDs) %>%
    as.character()
  # From each matrix, select the "DIA filtrate" sample rows as rows,
  # set the colnames to "scan-profileID" format,
  # and cbind them to a single matrix
  allCols <- map(as.list(...), ~ .[rows,]) %>%
    map2(., seq_along(.),  ~ set_colnames(.x, paste(.y, colnames(.x), sep="-"))) %>%
    do.call(cbind, .)
  allCols
}

#' Quick MOZER missing value replacement
#' 
#' Values below eps (typically zeroes) are replaced by a value calculated from
#' FUN and the number of missing values in a column. The value of `FUN(column)`
#' is computed from non-"zeroes". For a fraction of `f` zeroes in the column,
#' every zero is replaced by `(1-f) * FUN(column)`.
#' 
#' In this way, zero-heavy columns are filled with a low value that has a relation
#' to the data in the column, and the range of *imputed* values in each column is 
#' limited based on how rare non-zeroes are. I.e. for 90% zero values, the range of
#' minimum to median forced by imputation is limited to a factor of 10, whereas 
#' for 90% zero values, the range is limited to a factor of 100. 
#' 
#' Note that this applies for imputed values only, and correspondingly, real values
#' lower than the imputed zero may still occur. This is relevant when doing things
#' such as consolidating profiles, which can lead to very low non-missing outliers 
#' in otherwise high-intensity profiles when there is a missed integration for the
#' main peak. 
#' 
#' A solution would be to set eps (per column) such that it will impute at least 
#' the lowest 5 percentile of datapoints.
#' 
#' Columns with few missing values are filled with a value close to the median 
#' (when using the default `median` for `FUN`), representing "missing at random" values. 
#' This is crude and could use some real imputation.
#' 
#' The use of `median` over `min` favors a replacement with non-outliers to keep
#' the statistical parameters of the column relatively close; 
#' 
#'
#' @param intMatrix A matrix with variables in the columns, samples in the rows
#' @param eps A scalar or function describing the imputation level. Zero by default, such that
#'   only zeroes will be imputed
#' @param FUN Function to determine the base value for replacement. May return a scalar or a
#'   vector of the length of its input. By default, the median.
#'
#' @return
#' @export
#'
#' @examples
quickMozer <- function(intMatrix, eps = 0, FUN = median) {
  .eps <- eps
  if(!is_function(eps))
    .eps <- function(col) eps
  apply(intMatrix, 2, function(col) {
    col[col <= .eps(col)] <- NA
    col %>% replace_na((1 - sum(is.na(col)) / length(col)) * FUN(col[!is.na(col)]))
  })
}

filterAutocorrelations <- function(intMatrix, lag = 3, demean = TRUE) {
  #intMatrix[intMatrix <= eps] <- NA
  autocorrelations <- apply(intMatrix, 2, function(col)
    acf(col, lag.max=lag, plot=FALSE, demean = demean)$acf)
  autocorrelationSums <- colSums(autocorrelations[-1,]) %>% replace_na(0)
  autocorrelations[1,] <- autocorrelationSums
  return(t(autocorrelations))
  #return(intMatrix[,autocorrelationSums > cutoff])
}

imputeZeroes <- function(intMatrix, eps = 0) {
  # Remove all-zero columns since they only cause trouble. They should ideally be removed before this!
  zerocols <- apply(intMatrix, 2, function(col) sum(col > 0) == 0)
  # Impute zeroes with MOZER, see in docs
  intMatrix <- intMatrix[,!zerocols] %>% quickMozer(eps = eps) 
  intMatrix
}




#' Calculate correlations for edgelist
#' 
#' Given the edgelist in `profilesLinks`, computes the correlations for
#' the sample columns in `intMatrix`, whereby the `intMatrix` columns
#' are additionally indexed through the `labels` table (i.e. edge `i` references
#' column `labels[[i]]`, not column `i`).
#' 
#' Note: the edgelist is assumed to be strictly lower triangular. Otherwise, behaviour
#' is undefined.
#' 
#' @param profilesLinks Two-column (`c1, c2`) lower triangular (`c2 > c1`) edgelist.
#' @param friends The number of target edges per source edge
#' @param labels Reindexing for addressing subsets of intMatrix. By default, the identity map
#' @param intMatrix Matrix with features(rows) x samples(columns)
#' @param progress 
#'
#' @return
#' @export
#'
#' @examples
computeCor <- function(profilesLinks, friends, intMatrix, 
                       labels = seq_len(ncol(intMatrix)),
                       progress = TRUE) {
  # Calculate correlation in i-blocks
  if(progress)
    corProgress <- progress::progress_bar$new(total = length(friends))
  
  score <- numeric(nrow(profilesLinks))
  
  # Precompute matrix column indices to avoid hashtable lookup
  matrixCols <- match(labels, colnames(intMatrix))
  
  
  pos <- 1
  for(i in seq_along(friends)) {
    if(friends[i] > 0) {
      targets <- seq(from = pos, to = pos + friends[i] - 1)
      score[targets] <- 
        as.vector(cor(intMatrix[,matrixCols[[i]]], 
                      intMatrix[,matrixCols[profilesLinks[targets, 2]]]))
      pos <- pos + friends[[i]]
      if(progress)
        corProgress$tick()
    }
  }
  
  score
}


#' Compute RAMClust score
#'
#' @param cor Correlations vector
#' @param drt delta-RT vector
#' @param sr Weight for correlations, lower is stricter
#' @param st Weight for RT difference, lower is stricter
#'
#' @return
#' @export
#'
#' @examples
computeScore <- function(cor, drt,  sr = 0.5, st = 20) {
  score <- 1 - (exp(-(1-cor)^2 / sr^2)) *
    (exp(-drt^2 / st^2))
}


ramclustBlocks <- function(profiles, intMatrix,
                           sr = 0.5, st = 10,
                           drtmax = 15, neighbors = 25,
                           minC11 = 8, maxNeqRatio = 0.5,
                           multicore = 10,
                           outfile = "ramclustBlocks.log") {
  
  # Divide RT range into overlapping blocks of drtmax length
  rtlim <- range(profiles$mean_RT)
  blocks <- (rtlim[[2]] %/% (drtmax / 2)) + 1
  rtStart <- rtlim[[1]]
  rtEnd <- rtStart + drtmax
  
  estimated_blocks <- 2  * ( rtlim[[2]] - rtlim[[1]] ) / drtmax 
  message("estimated blocks:", estimated_blocks)
  block <- 0
  
  nnlists <- list()
  
  cl <- makeCluster(multicore, outfile = outfile)
  clusterEvalQ(cl, library(weights))
  
  while((rtStart < rtlim[[2]])) {
    block <- block + 1
    message(paste("block", block))
    profilesBlock <- profiles %>% dplyr::filter(
      dplyr::between(mean_RT, rtStart, rtEnd))
    profileNames <- profilesBlock$matrixID
    # extract intensity matrix and presence-absence matrix fopr
    # the profiles in discussion
    intSubmatrix <- intMatrix[,profileNames]
    notSubmatrix <- (intSubmatrix == 0)
    # Filter criteria for common and noncommon features:
    # c11 is where both profiles are nonzero,
    # c00 is where both profiles are zero,
    c11 <- t(intSubmatrix > 0) %*% (intSubmatrix > 0)
    # c00 <- t(intSubmatrix == 0) %*% (intSubmatrix == 0)
    c10 <- t(intSubmatrix > 0) %*% (intSubmatrix == 0)
    c01 <- t(intSubmatrix == 0) %*% (intSubmatrix > 0)

    dt <- outer(profilesBlock$mean_RT, profilesBlock$mean_RT, `-`)
    
    # Calculate scores only when there are sufficient
    # common non-zero points to be relevant,
    # then fill into the correlation matrix
    corSubmatrix <- matrix(0, nrow = ncol(intSubmatrix), ncol = ncol(intSubmatrix))
    validPoints <- (c11 > minC11) & 
      (((c10*c01) / (c11*c11)) < maxNeqRatio) # &
    #abs(dt) < drtmax
    calcIndices <- which( validPoints , arr.ind = TRUE)
    corPts <- which( validPoints)
    
    pb <- progress_bar$new(total = nrow(calcIndices))
    
    message(paste("Calculating", length(corPts), "of", length(corSubmatrix), "points (",
                  length(corPts) / length(corSubmatrix),")"))
    
    clusterExport(cl, "intSubmatrix", envir = environment())
    clusterExport(cl, "notSubmatrix", envir = environment())
    corData <- parApply(cl, calcIndices, 1, function(row) {
      #pb$tick()
      wtd.cors( c(intSubmatrix[,row[[1]]], 0), 
                c(intSubmatrix[,row[[2]]], 0),
                c(1-(notSubmatrix[,row[[1]]] * notSubmatrix[,row[[2]]]), 1)
      )
    })
    
    corSubmatrix[corPts] <- corData
    
    dimnames(corSubmatrix) <- list(
      colnames(intSubmatrix),
      colnames(intSubmatrix))
    
    dt <- outer(profilesBlock$mean_RT, profilesBlock$mean_RT, `-`)
    score <- computeScore(corSubmatrix, dt, sr = sr, st = st)
    # Filter out matches that shouldn't count:
    # the self-loop (diagonal)
    diag(score) <- 1

    # find (maximally) the specified number of nearest neighbors
    nn <- apply(score, 2, function(col) 
      as.list(sort(col)[seq_len(min(neighbors, nrow(score)))]))
    nn <- map(nn, unlist)
    nnlists[[as.character(block)]] <- nn
    #pb$tick()
    rtStart <- rtStart + drtmax/2
    rtEnd <- rtEnd + drtmax/2
  }
  stopCluster(cl)
  
  # Join the individual blocks together and select top-n neighbors
  # again, since we might have doubled the number here
  # Note, this is not using unique() at any time, so when there are
  # overlapping blocks, will this give the duplicate neighbor??
  nnlist <- nnlists %>% 
    map_dfr( ~ imap_dfr(
      .,~ tibble(e1 = .y, e2 = names(.x), weight = .x)),
      .id = "block")
  nnlist <- nnlist %>% 
    select(-block) %>% 
    group_by(e1) %>% 
    arrange(weight) %>%
    slice(seq_len(neighbors)) %>% 
    #mutate(w = 1-w) %>%
    ungroup()
  
  return(nnlists)
}




#' Title
#'
#' @param intMatrix 
#' @param profiles 

#'
#' @return
#' @export
#'
#' @examples
processDIA <- function(intMatrix, profiles, 
                       minDataCount = 6,
                       sr = 0.5, st = 10,
                       drtmax = 15, neighbors = 25,
                       minC11 = 8, maxNeqRatio = 0.5,
                       multicore = 10
                       # eps = 0,
                       # lag = 3, cutoff = 1, useAcf = 1,
                       # drtmax = 20, minCor = 0.25,
                       # sr = 0.5, st = 20,
                       # hmax = 0.6, deepSplit = TRUE, minModuleSize = 2
                       ) {

  # Remove data without enough data points  
  intMatrixCount <- colSums(intMatrix > 0)
  intMatrix <- intMatrix[, intMatrixCount >= minDataCount]
  
  # Remove profile IDs that are not in any filtrate samples, since they are not in the matrix
  # Notably not because of the combineIntMatrix step, but because of the acast step when 
  # intMatrixConsolidated is built
  
  cols <- tibble(matrixID = colnames(intMatrix)) %>% dplyr::mutate(ii = row_number())
  profiles_ <- profiles %>% 
    left_join(cols, by="matrixID") %>%
    ungroup() %>%
    dplyr::filter(!is.na(ii)) %>%
    dplyr::select(-ii) %>%
    dplyr::arrange(mean_RT) %>%
    dplyr::mutate(i = row_number())
  
  # reorder matrix in ascending RT order to allow direct addressing
  # Note: this is currently not used, since we use $matrixID to access the data
  intMatrix <- intMatrix[,profiles_$matrixID]
  
  profiles_ <- profiles_ %>%
    filter(matrixID %in% colnames(intMatrix))
  # make sure all profiles are present
  stopifnot(nrow(profiles_) == ncol(intMatrix))
  
  nnlist <- ramclustBlocks(profiles_, intMatrix,
                            sr, st,
                            drtmax, neighbors,
                            minC11, maxNeqRatio,
                            multicore)
  
  return(nnlist)
}

postprocessDIA <- function(nnlist, profiles, maxWeight) {
  # Crop to top-1 neighbor because it works so well :P
  nnlist <- nnlist %>% map_dfr( ~ imap_dfr(., ~ tibble(e1 = .y, e2 = names(.x), weight = .x)), .id = "block")
  
  nntop <- nnlist %>% 
    group_by(e1) %>% 
    arrange(weight) %>%
    slice(1) %>% 
    ungroup()

  
  # Make the graph consisting of all edges between top
  # neighbor profiles. 
  # Note: we computed the entire matrix, not just the upper
  # or lower triangle, so links will exist in both
  # directions, and every profile will have its nearest neighbor
  # linked to it.
  allProfiles <- union(nntop$e1, nntop$e2)
  profilesGraph <- profiles %>% 
    dplyr::filter(matrixID %in% allProfiles) %>%
    left_join(nntop %>% dplyr::rename(matrixID = e1), by = "matrixID") %>%
    dplyr::arrange(desc(mean_int)) %>%
    ungroup() %>%
    dplyr::mutate(i = row_number())
  
  profilesGraphIds <- profilesGraph$matrixID
  profilesGraph$e2i <- match(profilesGraph$e2, profilesGraphIds)
  stopifnot(!any(is.na(profilesGraph$e2i)))
  
  edges <- profilesGraph %>% select(e1 = i, e2 = e2i, weight = weight)

  edges_ <- edges %>% filter(weight < maxWeight)
  edgeGraph <- graph_from_data_frame(edges_, directed = FALSE)
  components <- tibble(i = as.numeric(names(V(edgeGraph))), component = membership(components(edgeGraph)))
  profilesGraph_ <- profilesGraph %>% select(-component) %>% left_join(components, by="i")
  return(list(profiles = profilesGraph_, maxWeight = maxWeight, n = count_components(edgeGraph)))
# 
#   componentsSummary <- componentsDf %>% 
#     group_by(run) %>%
#     mutate(hasComponent = !is.na(component)) %>% 
#     group_by(hasComponent, add = TRUE) %>% 
#     dplyr::summarize(totalInt = sum(mean_int), n=n(),
#                      profiles = n_distinct(component)) %>%
#     mutate(maxWeight = map_dbl(componentsData, "maxWeight")[run])
#   
#   componentsData[[3]]$profiles -> components
  
  
}


