# Viewer for results


# Load basic profiles and intensity matrix
load("results/clustering/viewerInput-global-20200109.RData")
# # Replace ramclustPeaksMS2 with the suspect screening results:
load("results/clustering/viewerInput-global-screening-20200109.RData")
# # Load clustering to overlay on the data
#load("results/clustering/viewerInput-km-filtered-20200109-03.RData")
load("results/clustering/viewerInput-hc-filtered-20200111-01.RData")
# Load visualization of clusters
load("results/clustering/visualizations-filtered-20200110-01.RData")

library(drake)
library(tidyverse)
library(shiny)
library(DT)
library(MSnbase)
library(circlize)
library(gplots)
library(ggplot2)
library(reshape2)
library(cowplot)


loadd(sampleList)

# this requires:
# intMatrixNorm, a matrix with features in columns normalized to [0..1].
# rcClustCuts
# clustCuts
# rcMatrixNorm

colorFun <- colorRamp2(c(0,1), c("darkgrey", "red"))
mapColor <- function(x, fun = colorFun, range = c(5,9)) {
  colorFun((log10(x) - range[[1]]) / (range[[2]]-range[[1]]))
}

# Colors for cluster IDs
cluster_col <- function(x) 2+ x %% 6
cluster_pch <- function(x) x %/% 6

profilesSelectedColumns <- c("component", "profile_ID", "matrixID", "mean_mz", "mean_int", "mean_RT",
  "number_peaks_total", "compMs2count",
  # suspect screening hits
  "name", "rtRank", "rtRankPred", "dppm", "formula", "category"
) 
profilesSelectedColumns <- c(
  intersect(colnames(ramclustPeaksMS2), profilesSelectedColumns), "RT_min")


ui <- fixedPage(
  tabsetPanel(
    tabPanel("cluster",
             fixedRow(
               column(6, numericInput("cluster", "Cluster",1,  1, length(rcClustCuts), step=1)),
               column(6, selectInput("vis", "Visualization", choices = names(vis), selected = names(vis)[[1]]))
             ),
             fixedRow(
               # clusterPlot: plots "averaged cluster" +- sd
               column( 6, 
                       plotOutput("clusterPlot",
                                  width = "500px",
                                  height = "400px"
                       )),
               # profilesPlot: plots RT-mz distribution of components in cluster
               # (the most intense profile in the component is used)
               column( 6, tabsetPanel(
                 tabPanel("mz-RT",
                          plotOutput("profilesPlot",
                                     width = "500px",
                                     height = "400px",
                                     click = "profilesToggle"
                                  )), # tabPanel "mz-RT"
                 tabPanel("cluster visualization",
                          plotOutput("clusterVisPlot",
                                     width = "500px",
                                     height = "400px",
                                     click = "clusterVisToggle"
                                     )) # tabPanel "cluster-vis"
                          )) # tabsetPanel and column
             ),
             # profilesTable: shows components in a cluster
             fixedRow(
               DT::dataTableOutput("profilesTable")
             )
    ),
    tabPanel("profile",
             # profilePlot: shows time series of selected profile(s)
             fixedRow(
               plotOutput("profilePlot",
                          width = "700px",
                          height = "400px")),
             # showMS2: checkbox to toggle display of MS2 shot positions, if present
             fixedRow(
               checkboxInput("showMS2", "show MS2 positions")
             ),
             fixedRow(
               DT::dataTableOutput("componentsTable")
             )
    ),
    tabPanel("MS2",
             fixedRow(
               DT::dataTableOutput("ms2Table")
             ),
             fixedRow(
               actionButton("extractMS2", "Extract MS2")
             ),
             fixedRow(
               column(4, selectInput("selectMS2", "Select MS2", choices = c())),
               column(4, downloadButton("mspDownload", label = "MSP")),
               column(4, downloadButton("siriusDownload", "SIRIUS")),
               column(4, checkboxInput("ms2RemovePrecursor", "Remove precursor"))
             ),
             fixedRow(
               plotOutput("ms2Plot",
                          width = "700px",
                          height = "400px")
             ),
             fixedRow(
               DT::dataTableOutput("ms2Data")
             )
             
    )
  ))




# 
# largeVisPlot <- neighborsLargeVis$coords %>%
#   t() %>% as.data.frame() %>%
#   mutate(cluster = largeVisClusters$cluster,
#          profile = colnames(rcLombMatrixScale))

server <- function(input, output, session) {
  
  store <- reactiveValues(ms2 = list())
  
  clusterName <- reactive({
    i_ <- input$cluster
    i <- names(rcClustCuts)[[i_]]
    return(i)
  })
  
  output$clusterPlot <- renderPlot({
    i <- clusterName()
    i_ <- input$cluster
    df <- tibble(
      val = rcClustCuts[[i]][,1],
      sampleIDs = as.integer(rownames(rcClustCuts[[i]]))
    )
    df <- df %>% left_join(sampleList, by="sampleIDs")
    df <- df %>% mutate(time = parse_datetime(mtime)) %>%
      mutate(wday=time  %>% as.POSIXlt() %>% `$`("wday")) %>%
      mutate(weekend = if_else(wday %in% c(0,6),1,0))
    
    
    
    # find the 5 most intense components
    pt <- profilesTable()
    ids <- pt %>% group_by(component) %>% arrange(desc(mean_int)) %>% slice(1) %>%
      ungroup() %>% arrange(desc(mean_int)) %>% slice(1:5) %>% pull(matrixID)
    
    profs <- intMatrixNorm[,ids] %>% melt(varnames = c("sampleIDs", "matrixID"))
    
    plot1 <- df %>% ggplot() + aes(x=time, y=val, col=weekend) + geom_line()
    plot2 <- profs %>% ggplot + aes(x=sampleIDs, y=value, col=matrixID) + geom_line()
    cowplot::plot_grid(plot1, plot2, nrow=2)
    
    
    
    #plot(rcClustCuts[[i]][,1], type="l", lwd=2, col="red", ylim=range(rcClustCuts[[i]][,c(1,3,4)]))
    #lines(rcClustCuts[[i]][,3], lwd=1)
    #lines(rcClustCuts[[i]][,4], lwd=1)
    #clusterIntensities <- rcMatrixNorm[,clustCuts == i_]
    #hist2d(rcMatri)
    # imageX <- seq_len(nrow(clusterIntensities))
    # imageY <- seq(from=0, to=1, length.out = 20)
    # clusterHistogram <- apply(clusterIntensities, 1,
    #                           function(col) {
    #                             hist(col[col>0], breaks=imageY, plot=FALSE)$counts
    #                           })
    # image(imageX, imageY, t(clusterHistogram))
    #title(paste0(i, " (", sum(clustCuts == i_) ," features)"))
  })
  
  
  profilesTable <- reactive({
    i <- clusterName()
    profiles <- colnames(rcMatrixNorm)[clustCuts == as.integer(i)]
    ramclustPeaksMS2 %>%
      dplyr::filter(component %in% profiles) %>%
      dplyr::mutate(RT_min = mean_RT / 60) %>%
      dplyr::mutate(suspect = (name != "")) %>%
      arrange(suspect, desc(mean_int))
  })
  
  output$profilesPlot <- renderPlot({
    par(mar=c(5,4,1,2) + 0.1)
    dataInCluster <- profilesTable() %>%
      arrange(desc(mean_int))
    plot(dataInCluster$mean_RT, dataInCluster$mean_mz, col = mapColor(dataInCluster$mean_int),
         pch = 16, xlim=c(200,950))
    rows <- input$profilesTable_rows_selected
    profiles <- profilesTable()[rows,,drop=FALSE]
    if(nrow(profiles) > 0)
      points(profiles$mean_RT, profiles$mean_mz, col = "black",
             pch = 8)
  })
  
  vis_ <- reactive({
    vis[[input$vis]]
  })
  
  visSelected <- reactive({
    vis_()[clustCuts == as.integer(clusterName()),,drop=FALSE]
  })
  
  output$clusterVisPlot <- renderPlot({
    par(mar=c(5,4,1,2) + 0.1)
    plot(vis_(), col = cluster_col(clustCuts), pch = cluster_pch(clustCuts))
    points(visSelected(), pch=21, col="black", bg="orange")
    
    rows <- input$profilesTable_rows_selected
    components <- profilesTable()[rows,,drop=FALSE] %>% pull(component)
    if(length(rows) > 0) {
      message(components)
      points(vis_()[components[!is.na(components)],,drop=FALSE], pch = 8, col="black")
    }
      
  })
  
  proxy <- dataTableProxy('profilesTable')
  
  
  observeEvent(input$profilesToggle, {
    toggleRow <- which(
      nearPoints(profilesSelected(), input$profilesToggle, xvar = "mean_RT", yvar = "mean_mz",
                 allRows = TRUE)$selected_)
    selected <- input$profilesTable_rows_selected
    
    selectNew <- setdiff(union(toggleRow,selected), intersect(toggleRow, selected))
    selectRows(proxy, selectNew)
  })
  
  
  observeEvent(input$clusterVisToggle, {
    toggleRow <- which(
      nearPoints(visSelected() %>% as_tibble(), input$clusterVisToggle, xvar = "c1", yvar = "c2",
                 allRows = TRUE)$selected_)
    selected <- input$profilesTable_rows_selected
    
    selectNew <- setdiff(union(toggleRow,selected), intersect(toggleRow, selected))
    selectRows(proxy, selectNew)
  })
  
  
  profilesSelected <- reactive({
    profilesTable() %>%
      dplyr::select(profilesSelectedColumns
      ) 
  })
  
  output$profilesTable <- DT::renderDataTable(
    { profilesSelected() %>%
        datatable(selection = "multiple") %>%
        formatRound("mean_mz", 4) %>%
        formatRound("mean_RT", 0) %>%
        formatRound("RT_min", 2) %>%
        formatSignif("mean_int", 2)
    })
  componentsSelected <- reactive({
    rows <- input$profilesTable_rows_selected
    if(length(rows) > 0)
      components <- as.data.frame(profilesSelected())[rows, "component"]
    else
      components <- as.data.frame(profilesSelected())
    ramclustPeaksMS2 %>%
      dplyr::filter(component %in% components) %>%
      arrange(component, mean_mz) %>%
      dplyr::select(component, #sampleID,
                    matrixID,
                    profile_ID, mean_mz, 
                    mean_RT, mean_int, number_peaks_total, #ms2count
                    ) %>%
      dplyr::mutate(RT_min = mean_RT / 60)
  })
  
  output$componentsTable <- DT::renderDataTable({
    componentsSelected()  %>%
      datatable(selection = "multiple") %>%
      formatRound("mean_mz", 4) %>%
      formatRound("mean_RT", 0) %>%
      formatRound("RT_min", 2) %>%
      formatSignif("mean_int", 2)
  })
  
  output$profilePlot <- renderPlot({
    # extract intensity profiles for all selected profile_IDs
    # (or all, if no profile_ID is selected in the table)
    rows <- input$componentsTable_rows_selected
    if(length(rows) > 0)
      profiles <- as.data.frame(componentsSelected())[rows, "matrixID"]
    else
      profiles <- as.data.frame(componentsSelected())[,"matrixID"]
    profilesData <- intMatrixNorm[
      ,as.character(profiles),drop=FALSE]
    # Plot colored, with profile_ID color legend
    plot.new()
    plot.window(xlim=c(0, nrow(profilesData)), ylim=range(0,1))
    axis(1)
    axis(2)
    box(bty="l")
    for(i in seq_len(ncol(profilesData)))
      lines(profilesData[,i], col=i)
    legend("topleft", legend = profiles, fill = seq_along(profiles))
    # If MS2 overplot is on: find matching MS2 and plot positions
    if(input$showMS2) {
      for (i in seq_along(profiles)) {
        ms2table() %>% 
          dplyr::mutate(relPrecursorInt = profilesData[index,i]) %>%
          dplyr::mutate(relint = totIonCurrent / precursorIntensity * relPrecursorInt) %>% 
          lines(relint ~ index, data = ., type='h', col=i)
      }
    }
  })
  
  #output$ms2table <- 
  ms2table <- reactive({
    rows <- input$componentsTable_rows_selected
    profiles <- as.data.frame(componentsSelected())[rows, "matrixID", drop=FALSE]
    ms2table <- componentsDDA %>% 
      dplyr::filter(matrixID %in% profiles) %>%
      #dplyr::mutate(index = match(sampleIDs, as.integer(rownames(intMatrixNorm)))) %>%
      dplyr::mutate(RT_min = RT / 60) %>%
      dplyr::select(specGroupID, matrixID,  RT_min, mean_int,
             mz, dppm, dRT, RT) 
    ms2table
  })
  
  output$ms2Table <-  DT::renderDataTable({
    ms2table()  %>%
      datatable(selection = "multiple") %>%
      formatRound("RT", 0) %>%
      formatRound("mz", 4) %>%
      formatRound("RT_min", 2) %>%
      formatSignif("mean_int", 2) %>%
      formatRound("dRT", 1) %>%
      formatRound("dppm", 1)
  })
  
  
  observeEvent(input$extractMS2, {
    rows <- input$ms2Table_rows_selected
    ms2data <- ms2table()[rows,,drop=FALSE]
    specGroupID <- ms2data$specGroupID
    store$ms2 <- map(specGroupID, ~ getBestSpectrum(.x, spectraDDA, sampleNextDIA, scansDIA, files, n = 5,
                                file = c("DIA_preceding"), isolationWindow = 1, precursorWindow = 5))
    updateSelectInput(session, "selectMS2", choices = seq_along(store$ms2))
  })
  
  
  selectedMS2 <- reactive({
    store$ms2[[as.numeric(input$selectMS2)]]
  })
  
  output$ms2Plot <- renderPlot({
    ms2 <- selectedMS2()
    ms2 %>% ggplot() + 
      aes(x=V1, xend= V1, yend=V2, col=eicCor) +
      geom_segment(y=0) +
      facet_grid(rows = vars(level), scales = "free") +
      scale_color_gradient( low = "red", high = "blue", limits = c(0, 0.5))
    #plot(ms2$V1, ms2$V2, type='h')
  })
  
  output$ms2Data <- DT::renderDataTable({
    as.data.frame(selectedMS2())
  })
  
  output$mspDownload <- downloadHandler(
    function() paste0(attr(selectedMS2(), "header")$specGroupID, ".msp"),
    function(f) writeLines(exportMsp(selectedMS2()), con=f)
  )
  
  output$siriusDownload <- downloadHandler(
    function() paste0(attr(selectedMS2(), "header")$specGroupID, ".ms"),
    function(f) writeLines(exportSirius(selectedMS2()), con=f)
  )
  
    
}

shinyApp(ui, server, options = list(width="1000px"))
