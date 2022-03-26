#####################
#Shiny App Server
#####################
server <- function(
  input = input, 
  output = output, 
  session = session
){
  
  ### What would need to be changed in server.R
  # tileSize <- 
  # ArchRProj <- 
  # groupBy <- 
    
  output$Metadata <- renderRHandsontable({
    groups <- gtools::mixedsort(unique(ccd[,input$grouping]))
    mdata <- data.frame(
      groupBy = input$grouping,
      include = rep(TRUE,length(groups)), 
      group = groups, 
      color = paletteDiscrete(values = groups)[groups], 
      nCells = as.vector(table(ccd[,input$grouping])[groups]),
      medianTSS = getGroupSummary(ArchRProj = ArchRProj, select = "TSSEnrichment", summary = "median", groupBy = input$grouping)[groups],
      medianFragments = getGroupSummary(ArchRProj = ArchRProj, select = "nFrags", summary = "median", groupBy = input$grouping)[groups],
      stringsAsFactors = FALSE
    )
    rownames(mdata) <- NULL
    rhandsontable(mdata)
  })
  
  #Update Sliders
  observeEvent(input$range_min, {
    updateSliderInput(session, "range",
                      value = c(input$range_min,max(input$range)))
  })
  
  observeEvent(input$range_max, {
    updateSliderInput(session, "range",
                      value = c(input$range_min,input$range_max))
  })
  
  observeEvent(input$range , {
    
    updateNumericInput(session, "range_min", value = min(input$range))
    updateNumericInput(session, "range_max", value = max(input$range))
    
  }, priority = 200)
  
  output$checkbox <- renderUI({
    choice <- gtools::mixedsort(unique(ccd[,input$grouping,drop=TRUE]))
    checkboxGroupInput("checkbox","Select Groups", choices = choice, selected = choice)
  })    
  
  #################################
  # Inputs that cause re-plotting
  #################################
  #toListen <- reactive({
  #  list(input$restartButton, input$name)
  #})
  
  restartFN <- observeEvent(input$restartButton, {
    
    if (input$name == ""){
      
      output$ATAC <- renderPlot({
        p <- ggplot() +
          xlim(c(-5,5)) + ylim(c(-5,5)) +
          geom_text(size=20, aes(x = 0, y = 0, label = "Please Supply\nA Valid Gene!")) + theme_void()
        print(p)
      })
      
    }else{
      output$ATAC <- renderPlot({
        
        withProgress(message = 'Plotting', style = "notification", value = 0, {
          
          #Get Region if Gene Symbol
          region <- geneAnnotation$genes
          
          if(tolower(input$name) %ni% tolower(mcols(region)$symbol)){
            p <- ggplot() +
              xlim(c(-5,5)) + ylim(c(-5,5)) +
              geom_text(size=20, aes(x = 0, y = 0, label = "Please Supply\nA Valid Gene!")) + theme_void()
            return(print(p))
          }
          
          region <- region[which(tolower(mcols(region)$symbol) %in% tolower(input$name))]
          region <- region[order(match(tolower(mcols(region)$symbol), tolower(input$name)))]
          region1 <- GenomicRanges::resize(region, 1, "start")
          strand(region1) <- "*"
          
          #Extend Region
          #region <- extendGR(region, upstream = -min(input$range) * 1000, downstream = max(input$range) * 1000)
          #Pre-Load full window for even faster plotting
          region <- extendGR2(region1, upstream = 250000, downstream = 250000)
          tmpArchRRegion <<- extendGR2(region1, 
                                       upstream = -min(isolate(input$range)) * 1000, 
                                       downstream = max(isolate(input$range)) * 1000
          )
          region <- tmpArchRRegion
          
          setProgress(0.1)
          
          #User Inputs
        #  groupBy <- isolate(input$grouping)
          
          groupDF <- tryCatch({
            isolate(hot_to_r(input$Metadata))
          },error=function(x){
            groups <- gtools::mixedsort(unique(ccd[,isolate(input$grouping)]))
            mdata <- data.frame(
              groupBy = input$grouping,
              include = rep(TRUE,length(groups)), 
              group = groups, 
              color = paletteDiscrete(values = groups)[groups], 
              nCells = as.vector(table(ccd[,input$grouping])[groups]),
              medianTSS = getGroupSummary(ArchRProj = ArchRProj, select = "TSSEnrichment", summary = "median", groupBy = input$grouping)[groups],
              medianFragments = getGroupSummary(ArchRProj = ArchRProj, select = "nFrags", summary = "median", groupBy = input$grouping)[groups],
              stringsAsFactors = FALSE
            )
            rownames(mdata) <- NULL
            mdata
          })
          
          if(groupDF$groupBy[1] != groupBy){
            groups <- gtools::mixedsort(unique(ccd[,isolate(input$grouping)]))
            groupDF <- data.frame(
              groupBy = groupBy,
              include = rep(TRUE,length(groups)), 
              group = groups, 
              color = paletteDiscrete(values = groups)[groups], 
              stringsAsFactors = FALSE
            )
            rownames(groupDF) <- NULL
          }
          
          useGroups <- groupDF[groupDF[,"include"],"group"]
          
          
          if(!all(.isColor(groupDF[groupDF[,"include"], "color"]))){
            p <- ggplot() +
              xlim(c(-5,5)) + ylim(c(-5,5)) +
              geom_text(size=20, aes(x = 0, y = 0, label = "Error Colors from Metadata\n not real R colors!")) + theme_void()
            return(print(p))               
          }
          
          if(any(useGroups %ni% ccd[, groupBy])){
            p <- ggplot() +
              xlim(c(-5,5)) + ylim(c(-5,5)) +
              geom_text(size=20, aes(x = 0, y = 0, label = "Error Groups from Metadata\n not present in groupBy!")) + theme_void()
            return(print(p)) 
          }
          
          pal <- groupDF[groupDF[,"include"], "color"]
          names(pal) <- groupDF[groupDF[,"include"], "group"]
          
          ylim <- c(0, isolate(input$ymax))
          normMethod <- isolate(input$normATAC)
          tileSize <- isolate(input$tile_size)
          
          p <- .bulkTracks(
            ArchRProj = ArchRProj, 
            region = region, 
            tileSize = tileSize, 
            useGroups = useGroups,
            groupBy = groupBy,
            threads = threads, 
            minCells = minCells,
            ylim = ylim,
            baseSize = baseSize,
            borderWidth = borderWidth,
            tickWidth = tickWidth,
            facetbaseSize = facetbaseSize,
            normMethod = normMethod,
            geneAnnotation = geneAnnotation,
            title = "",
            pal = pal, 
            tstart = NULL,
            logFile = logFile
          ) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
          
          #p <- p + .suppressAll(scale_x_continuous(limits = c(start(tmpArchRRegion), end(tmpArchRRegion)), expand = c(0,0)))
          
          tmpArchRP <<- p
          
          setProgress(0.5)
          
          if(!is.null(features)){
            
            f <- .featureTracks(
              features = features, 
              region = tmpArchRRegion,
              facetbaseSize = facetbaseSize, 
              hideX = TRUE, 
              title = "Peaks",
              logFile = logFile
            ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
            
            #f <- f + .suppressAll(scale_x_continuous(limits = c(start(tmpArchRRegion), end(tmpArchRRegion)), expand = c(0,0)))
            
          }
          
          if(!is.null(loops)){
            
            l <- .loopTracks(
              loops = loops, 
              region = tmpArchRRegion, 
              facetbaseSize = facetbaseSize,
              hideX = TRUE, 
              hideY = TRUE,
              title = "Loops",
              logFile = logFile
            ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
          }
          
          setProgress(0.6)
          
          g <- .geneTracks(
            geneAnnotation = geneAnnotation, 
            region = tmpArchRRegion, 
            facetbaseSize = facetbaseSize,
            labelSize = 3,
            title = "Genes",
            logFile = logFile
          ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
          
          #g <- .suppressAll(g + scale_x_continuous(limits = c(start(tmpArchRRegion), end(tmpArchRRegion)), expand = c(0,0)))
          
          setProgress(0.8)
          
          if(!is.null(loops)){
            if(!is.null(features)){
              suppressWarnings(print(ggAlignPlots(p, f, l, g, sizes = c(10, 1.5, 3, 4),type = "v", draw = TRUE)))
            }else{
              suppressWarnings(print(ggAlignPlots(p, l, g, sizes = c(10, 3, 4),type = "v", draw = TRUE)))
            }
          }else{
            if(!is.null(features)){
              suppressWarnings(print(ggAlignPlots(p, f, g, sizes = c(10, 2, 4),type = "v", draw = TRUE)))
            }else{
              suppressWarnings(print(ggAlignPlots(p, g, sizes = c(10, 4),type = "v", draw = TRUE)))
            }
          }
          
          setProgress(1)
          
        })
        
      })
    }    
    
  })
  
  #######################################
  # When Download Is Initiated
  #######################################
  
  # downloadHandler contains 2 arguments as functions, namely filename, content
  output$down <- downloadHandler(
    
    filename <- function(){
      paste0("ArchRBrowser-",input$name,"-",seqnames(tmpArchRRegion)[1],":",start(tmpArchRRegion)[1],"-",end(tmpArchRRegion)[1],".pdf")
    },
    
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      
      withProgress(message = 'Creating PDF', style = "notification", value = 0, {
        
        if(!exists("tmpArchRP")){
          
          #User Inputs
          groupBy <- isolate(input$grouping)
          groupDF <- isolate(hot_to_r(input$Metadata))
          useGroups <- groupDF[groupDF[,"include"],"group"]
          
          .isColor <- function(x = NULL) {
            unlist(lapply(x, function(y) tryCatch(is.matrix(col2rgb(y)), 
                                                  error = function(e) FALSE)))
          }
          
          if(!all(.isColor(groupDF[groupDF[,"include"], "color"]))){
            p <- ggplot() +
              xlim(c(-5,5)) + ylim(c(-5,5)) +
              geom_text(size=20, aes(x = 0, y = 0, label = "Error Colors from Metadata\n not real R colors!")) + theme_void()
            return(print(p))               
          }
          
          if(any(useGroups %ni% ccd[, groupBy])){
            p <- ggplot() +
              xlim(c(-5,5)) + ylim(c(-5,5)) +
              geom_text(size=20, aes(x = 0, y = 0, label = "Error Groups from Metadata\n not present in groupBy!")) + theme_void()
            return(print(p)) 
          }
          
          pal <- groupDF[groupDF[,"include"], "color"]
          names(pal) <- groupDF[groupDF[,"include"], "group"]
          
          ylim <- c(0, isolate(input$ymax))
          normMethod <- isolate(input$normATAC)
          tileSize <- isolate(input$tile_size)
          
          p <- .bulkTracks(
            ArchRProj = ArchRProj, 
            region = tmpArchRRegion, 
            tileSize = tileSize, 
            useGroups = useGroups,
            groupBy = groupBy,
            threads = threads, 
            minCells = minCells,
            ylim = ylim,
            baseSize = baseSize,
            borderWidth = borderWidth,
            tickWidth = tickWidth,
            facetbaseSize = facetbaseSize,
            normMethod = normMethod,
            geneAnnotation = geneAnnotation,
            title = "",
            pal = pal, 
            tstart = NULL,
            logFile = logFile
          ) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
          
        }else{
          
          print("Using previous ggplot")
          
          p <- tmpArchRP
          
        }
        
        setProgress(0.5)
        
        if(!is.null(features)){
          
          f <- .featureTracks(
            features = features, 
            region = tmpArchRRegion,
            facetbaseSize = facetbaseSize, 
            hideX = TRUE, 
            title = "Peaks",
            logFile = logFile
          ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
          
        }
        
        setProgress(0.6)
        
        if(!is.null(loops)){
          
          l <- .loopTracks(
            loops = loops, 
            region = tmpArchRRegion, 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            hideY = TRUE,
            title = "Loops",
            logFile = logFile
          ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
        }
        
        setProgress(0.7)
        
        g <- .geneTracks(
          geneAnnotation = geneAnnotation, 
          region = tmpArchRRegion, 
          facetbaseSize = facetbaseSize,
          labelSize = 3,
          title = "Genes",
          logFile = logFile
        ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
        
        setProgress(0.8)
        
        pdf(file = file, width = input$plot_width, height = input$plot_height)
        
        if(!is.null(loops)){
          if(!is.null(features)){
            a <- suppressWarnings(ggAlignPlots(p, f, l, g, sizes = c(10, 1.5, 3, 4),type = "v", draw = FALSE))
          }else{
            a <- suppressWarnings(ggAlignPlots(p, l, g, sizes = c(10, 3, 4),type = "v", draw = FALSE))
          }
        }else{
          if(!is.null(features)){
            a <- suppressWarnings(ggAlignPlots(p, f, g, sizes = c(10, 2, 4),type = "v", draw = FALSE))
          }else{
            a <- suppressWarnings(ggAlignPlots(p, g, sizes = c(10, 4),type = "v", draw = FALSE))
          }
        }
        
        suppressWarnings(grid::grid.draw(a))
        dev.off()
        
        setProgress(1)
        
      })
      
    })
  
  
  exitFN <- observeEvent(input$exitButton, {
    if(exists("tmpArchRRegion")){
      .suppressAll(rm(tmpArchRRegion))
    }
    if(exists("tmpArchRP")){
      .suppressAll(rm(tmpArchRP))
    }
    shiny::stopApp()
  })
  
}


