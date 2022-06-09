####################################################################
# Signal Track Plotting Methods
####################################################################

#' Launch ArchR Genome Browser
#' 
#' This function will open an interactive shiny session in style of a browser track. It allows for normalization of the signal which
#' enables direct comparison across samples. Note that the genes displayed in this browser are derived from your `geneAnnotation`
#' (i.e. the `BSgenome` object you used) so they may not match other online genome browsers that use different gene annotations.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param features A `GRanges` object containing the "features" to be plotted via the "featureTrack". This should be thought of as a
#' bed track. i.e. the set of peaks obtained using `getPeakSet(ArchRProj))`. 
#' @param loops A `GRanges` object containing the "loops" to be plotted via the "loopTrack".
#' This `GRanges` object start represents the center position of one loop anchor and the end represents the center position of another loop anchor. 
#' A "loopTrack" draws an arc between two genomic regions that show some type of interaction. This type of track can be used 
#' to display chromosome conformation capture data or co-accessibility links obtained using `getCoAccessibility()`. 
#' @param minCells The minimum number of cells contained within a cell group to allow for this cell group to be plotted. This argument
#' can be used to exclude pseudo-bulk replicates generated from low numbers of cells.
#' @param baseSize The numeric font size to be used in the plot. This applies to all plot labels.
#' @param borderWidth The numeric line width to be used for plot borders.
#' @param tickWidth The numeric line width to be used for axis tick marks.
#' @param facetbaseSize The numeric font size to be used in the facets (gray boxes used to provide track labels) of the plot.
#' @param geneAnnotation The `geneAnnotation` object to be used for plotting the "geneTrack" object. See `createGeneAnnotation()` for more info.
#' @param browserTheme A `shinytheme` from shinythemes for viewing the ArchR Browser. If not installed this will be NULL.
#' To install try devtools::install_github("rstudio/shinythemes").
#' @param threads The number of threads to use for parallel execution.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
ArchRBrowser <- function(
  ArchRProj = NULL,
  features = getPeakSet(ArchRProj),
  loops = getCoAccessibility(ArchRProj),
  minCells = 25,
  baseSize = 10,
  borderWidth = 0.5,
  tickWidth = 0.5,
  facetbaseSize = 12,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  browserTheme = "cosmo",
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("ArchRBrowser")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = features, name = "features", valid = c("granges", "grangeslist", "null"))
  .validInput(input = loops, name = "loops", valid = c("granges", "grangeslist", "null"))
  .validInput(input = minCells, name = "minCells", valid = c("integer"))
  .validInput(input = baseSize, name = "baseSize", valid = c("integer"))
  .validInput(input = borderWidth, name = "borderWidth", valid = c("numeric"))
  .validInput(input = tickWidth, name = "tickWidth", valid = c("numeric"))
  .validInput(input = facetbaseSize, name = "facetbaseSize", valid = c("numeric"))
  geneAnnotation <- .validGeneAnnotation(geneAnnotation)
  .validInput(input = browserTheme, name = "browserTheme", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  .startLogging(logFile=logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "ArchRBrowser Input-Parameters", logFile = logFile)

  .requirePackage("shiny", installInfo = 'install.packages("shiny")')
  .requirePackage("rhandsontable", installInfo = 'install.packages("rhandsontable")')

  #Determine Grouping Methods
  ccd <- getCellColData(ArchRProj)
  discreteCols <- lapply(seq_len(ncol(ccd)), function(x){
    .isDiscrete(ccd[, x])
  }) %>% unlist %>% {colnames(ccd)[.]}
  if("Clusters" %in% discreteCols){
    selectCols <- "Clusters"
  }else{
    selectCols <- "Sample"
  }

  #Extend where upstream can be negative for browser
  extendGR2 <-  function(gr = NULL, upstream = NULL, downstream = NULL){
    .validInput(input = gr, name = "gr", valid = c("GRanges"))
    .validInput(input = upstream, name = "upstream", valid = c("integer"))
    .validInput(input = downstream, name = "downstream", valid = c("integer"))
    #Get Info From gr
    st <- start(gr)
    ed <- end(gr)
    #https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
    isMinus <- BiocGenerics::which(strand(gr) == "-")
    isOther <- BiocGenerics::which(strand(gr) != "-")
    #Forward
    st[isOther] <- st[isOther] - upstream
    ed[isOther] <- ed[isOther] + downstream
    #Reverse
    ed[isMinus] <- ed[isMinus] + upstream
    st[isMinus] <- st[isMinus] - downstream
    #If Any extensions now need to be flipped.
    end(gr) <- pmax(st, ed)
    start(gr) <- pmin(st, ed)
    return(gr)
  }


  #####################
  #Shiny App UI
  #####################
  if(!requireNamespace("shinythemes", quietly = TRUE)){
    .logMessage("shinythemes not found! To see a nice theme use :\n\tinstall.packages('shinythemes')\nContinuing wihtout shinythemes!", verbose = verbose, logFile = logFile)
    theme <- NULL
  }else{
    theme <- shinythemes::shinytheme(browserTheme)
  }

  ui <- fluidPage(
    theme = theme,
    titlePanel(
        h1(div(HTML(paste0("<b>ArchR Browser v1 : nCells = ", formatC(nCells(ArchRProj), format="f", big.mark = ",", digits=0), "</b>"))), align = "left")
    ),
    sidebarLayout(
      sidebarPanel(
        tags$head(
          tags$style(HTML("hr {border-top: 1px solid #000000;}"))
        ),
        tags$head(
          tags$style(HTML('#exitButton{background-color:#D60000}'))
        ), 
        tags$head(
          tags$style(HTML('#restartButton{background-color:#02A302}'))
        ),
        tags$head(
          tags$style(HTML('#plot_height{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#plot_width{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#ymax{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#tile_size{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#range_min{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#range_max{height: 35px}'))
        ),
        actionButton(inputId = "exitButton", label = "Exit Browser", icon = icon("times-circle")),
        br(),
        br(),
        actionButton(inputId = "restartButton", label = "Plot Track!", icon = icon("play-circle")),
        #div(style="display:inline-block;width:50%;text-align: center;",actionButton("exitButton", label = "Exit Browser", icon = icon("paper-plane"))),
        #br(),
        #div(style="display:inline-block;width:50%;text-align: center;",actionButton("restartButton", label = "Plot Track!", icon = icon("paper-plane"))),
        br(),
        br(),
        selectizeInput("name",
                       label = "Gene Symbol",
                       choices = as.vector(geneAnnotation$genes$symbol)[!is.na(as.vector(geneAnnotation$genes$symbol))],
                       multiple = FALSE,
                       options = list(placeholder = 'Select a Center'),
                       selected = "CD4"
                      ),
        selectizeInput("grouping",
           label = "groupBy",
           choices = discreteCols,
           multiple = FALSE,
           options = list(placeholder = 'Select Grouping'),
           selected = selectCols
        ),
        sliderInput("range", "Distance From Center (kb):", min = -250, max = 250, value = c(-50,50)),
        splitLayout(cellWidths = c("50%","50%"),
          numericInput("range_min", "Distance (-kb):", min = -250, max = 250, value = -50),
          numericInput("range_max", "Distance (+kb):", min = -250, max = 250, value = 50)
        ),
        splitLayout(cellWidths = c("50%","50%"),
          numericInput("tile_size", "TileSize:", min = 10, max = 5000, value = 250),
          numericInput("ymax", "Y-Max (0,1):", min = 0, max = 1, value = 0.99)
        ),
        hr(),
        downloadButton(outputId = "down", label = "Download the Track!"),
        br(),
        br(),
        splitLayout(cellWidths = c("50%","50%"),
          numericInput("plot_width", "Width", min = 0, max = 250, value = 8),
          numericInput("plot_height", "Height", min = 0, max = 250, value = 12)
        ),
        width = 2, height = "750px", position = "left"
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Plot",
            plotOutput(outputId = "ATAC", width= "800px", height = "725px")
          ),
          tabPanel("Additional Params",
            br(),
            br(),
            selectizeInput("normATAC",
               label = "normMethod",
               choices = c("ReadsInTSS", "ReadsInPromoter", "nFrags", "None"),
               multiple = FALSE,
               options = list(placeholder = 'Select NormMethod'),
               selected = "ReadsInTSS"
            ),
            br(),
            div(HTML("<b>Group Metadata</b>")),
            rHandsontableOutput("Metadata", width= "1085px",height = "800px")
          )
        )
      )
    )
  )

  #####################
  #Shiny App Server
  #####################
  server <- function(
        input = input, 
        output = output, 
        session = session
    ){

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
            groupBy <- isolate(input$grouping)

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

  shiny::runGadget(ui, server)

}

#' @export
ArchRBrowserTrack <- function(...){
    .Deprecated("plotBrowserTrack")
    plotBrowserTrack(...)
}

#' Plot an ArchR Region Track
#' 
#' This function will plot the coverage at an input region in the style of a browser track. It allows for normalization of the signal
#' which enables direct comparison across samples. Note that the genes displayed in these plots are derived from your `geneAnnotation`
#' (i.e. the `BSgenome` object you used) so they may not match other online genome browsers that use different gene annotations.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param region A `GRanges` region that indicates the region to be plotted. If more than one region exists in the `GRanges` object,
#' all will be plotted. If no region is supplied, then the `geneSymbol` argument can be used to center the plot window at the
#' transcription start site of the supplied gene.
#' @param groupBy A string that indicates how cells should be grouped. This string corresponds to one of the standard or
#' user-supplied `cellColData` metadata columns (for example, "Clusters"). Cells with the same value annotated in this metadata
#' column will be grouped together and the average signal will be plotted.
#' @param useGroups A character vector that is used to select a subset of groups by name from the designated `groupBy` column in
#' `cellColData`. This limits the groups to be plotted.
#' @param plotSummary A character vector containing the features to be potted. Possible values include "bulkTrack" (the ATAC-seq signal),
#' "scTrack" (scATAC-seq signal), "featureTrack" (i.e. the peak regions), "geneTrack" (line diagrams of genes with introns and exons shown. 
#' Blue-colored genes are on the minus strand and red-colored genes are on the plus strand), and "loopTrack" (links between a peak and a gene).
#' @param sizes A numeric vector containing up to 3 values that indicate the sizes of the individual components passed to `plotSummary`.
#' The order must be the same as `plotSummary`.
#' @param features A `GRanges` (for a single feature track) or `GRangesList` (for multiple feature tracks) object containing the "features" to
#' be plotted via the "featureTrack". This should be thought of as a bed track. i.e. the set of peaks obtained using `getPeakSet(ArchRProj))`.
#' If you provide a `GRangesList`, then each element of that object must be named and this name will be used on the plot.
#' For example - `GRangesList("peaks" = peak_gr, "other" = other_gr)`.
#' @param loops A `GRanges` object containing the "loops" to be plotted via the "loopTrack".
#' This `GRanges` object start represents the center position of one loop anchor and the end represents the center position of another loop anchor. 
#' A "loopTrack" draws an arc between two genomic regions that show some type of interaction. This type of track can be used 
#' to display chromosome conformation capture data or co-accessibility links obtained using `getCoAccessibility()`. 
#' @param geneSymbol If `region` is not supplied, plotting can be centered at the transcription start site corresponding to the gene symbol(s) passed here.
#' @param useMatrix If supplied geneSymbol, one can plot the corresponding GeneScores/GeneExpression within this matrix. I.E. "GeneScoreMatrix"
#' @param log2Norm If supplied geneSymbol, Log2 normalize the corresponding GeneScores/GeneExpression matrix before plotting.
#' @param upstream The number of basepairs upstream of the transcription start site of `geneSymbol` to extend the plotting window.
#' If `region` is supplied, this argument is ignored.
#' @param downstream The number of basepairs downstream of the transcription start site of `geneSymbol` to extend the plotting window.
#' If `region` is supplied, this argument is ignored.
#' @param tileSize The numeric width of the tile/bin in basepairs for plotting ATAC-seq signal tracks. All insertions in a single bin will be summed.
#' @param minCells The minimum number of cells contained within a cell group to allow for this cell group to be plotted. This argument can be
#' used to exclude pseudo-bulk replicates generated from low numbers of cells.
#' @param normMethod The name of the column in `cellColData` by which normalization should be performed. The recommended and default value
#' is "ReadsInTSS" which simultaneously normalizes tracks based on sequencing depth and sample data quality.
#' @param threads The number of threads to use for parallel execution.
#' @param ylim The numeric quantile y-axis limit to be used for for "bulkTrack" plotting. This should be expressed as `c(lower limit, upper limit)` such as `c(0,0.99)`. If not provided, the y-axis limit will be c(0, 0.999).
#' @param pal A custom palette (see `paletteDiscrete` or `ArchRPalettes`) used to override coloring for groups.
#' @param baseSize The numeric font size to be used in the plot. This applies to all plot labels.
#' @param scTileSize The width of the tiles in scTracks. Larger numbers may make cells overlap more. Default is 0.5 for about 100 cells.
#' @param scCellsMax The maximum number of cells for scTracks.
#' @param borderWidth The numeric line width to be used for plot borders.
#' @param tickWidth The numeric line width to be used for axis tick marks.
#' @param facetbaseSize The numeric font size to be used in the facets (gray boxes used to provide track labels) of the plot.
#' @param geneAnnotation The `geneAnnotation` object to be used for plotting the "geneTrack" object. See `createGeneAnnotation()` for more info.
#' @param title The title to add at the top of the plot next to the plot's genomic coordinates.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
plotBrowserTrack <- function(
  ArchRProj = NULL, 
  region = NULL, 
  groupBy = "Clusters",
  useGroups = NULL, 
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(10, 1.5, 3, 4),
  features = getPeakSet(ArchRProj),
  loops = getCoAccessibility(ArchRProj),
  geneSymbol = NULL,
  useMatrix = NULL,
  log2Norm = TRUE,
  upstream = 50000,
  downstream = 50000,
  tileSize = 250, 
  minCells = 25,
  normMethod = "ReadsInTSS",
  threads = getArchRThreads(), 
  ylim = NULL,
  pal = NULL,
  baseSize = 7,
  scTileSize = 0.5,
  scCellsMax = 100,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  verbose = TRUE,
  logFile = createLogFile("plotBrowserTrack")
  ){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  .validInput(input = region, name = "region", valid = c("granges","null"))
  .validInput(input = groupBy, name = "groupBy", valid = "character")
  .validInput(input = useGroups, name = "useGroups", valid = c("character", "null"))
  .validInput(input = plotSummary, name = "plotSummary", valid = "character")
  .validInput(input = sizes, name = "sizes", valid = "numeric")
  .validInput(input = features, name = "features", valid = c("granges", "grangeslist", "null"))
  .validInput(input = loops, name = "loops", valid = c("granges", "grangeslist", "null"))
  .validInput(input = geneSymbol, name = "geneSymbol", valid = c("character", "null"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character", "null"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = upstream, name = "upstream", valid = c("integer"))
  .validInput(input = downstream, name = "downstream", valid = c("integer"))
  .validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  .validInput(input = minCells, name = "minCells", valid = c("integer"))
  .validInput(input = normMethod, name = "normMethod", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = ylim, name = "ylim", valid = c("numeric", "null"))
  .validInput(input = pal, name = "pal", valid = c("palette", "null"))
  .validInput(input = baseSize, name = "baseSize", valid = "numeric")
  .validInput(input = scTileSize, name = "scTileSize", valid = "numeric")
  .validInput(input = scCellsMax, name = "scCellsMax", valid = "integer")
  .validInput(input = borderWidth, name = "borderWidth", valid = "numeric")
  .validInput(input = tickWidth, name = "tickWidth", valid = "numeric")
  .validInput(input = facetbaseSize, name = "facetbaseSize", valid = "numeric")
  geneAnnotation <- .validGeneAnnotation(geneAnnotation)
  .validInput(input = title, name = "title", valid = "character")

  tstart <- Sys.time()
  .startLogging(logFile=logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "plotBrowserTrack Input-Parameters", logFile = logFile)

  ##########################################################
  # Get Region Where Plot Will Occur (GenomicRanges)
  ##########################################################
  .logDiffTime("Validating Region", t1=tstart, verbose=verbose, logFile=logFile)
  if(is.null(region)){
    if(!is.null(geneSymbol)){
      region <- geneAnnotation$genes
      region <- region[which(tolower(mcols(region)$symbol) %in% tolower(geneSymbol))]
      region <- region[order(match(tolower(mcols(region)$symbol), tolower(geneSymbol)))]
      print(region)
      region <- GenomicRanges::resize(region, 1, "start")
      strand(region) <- "*"
      region <- extendGR(region, upstream = upstream, downstream = downstream)
    }
  }
  region <- .validGRanges(region)
  .logThis(region, "region", logFile = logFile)

  if(is.null(geneSymbol)){
    useMatrix <- NULL
  }

  if(!is.null(useMatrix)){
    featureMat <- .getMatrixValues(
      ArchRProj = ArchRProj,
      matrixName = useMatrix,
      name = mcols(region)$symbol
    )
    if(log2Norm){
      featureMat <- log2(featureMat + 1) 
    }
    featureMat <- data.frame(t(featureMat))
    featureMat$Group <- getCellColData(ArchRProj, groupBy, drop = FALSE)[rownames(featureMat), 1]
  }

  ggList <- lapply(seq_along(region), function(x){

    plotList <- list()

    ##########################################################
    # Bulk Tracks
    ##########################################################
    if("bulktrack" %in% tolower(plotSummary)){
      .logDiffTime(sprintf("Adding Bulk Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$bulktrack <- .bulkTracks(
        ArchRProj = ArchRProj, 
        region = region[x], 
        tileSize = tileSize, 
        groupBy = groupBy,
        threads = threads, 
        minCells = minCells,
        pal = pal,
        ylim = ylim,
        baseSize = baseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth,
        facetbaseSize = facetbaseSize,
        normMethod = normMethod,
        geneAnnotation = geneAnnotation,
        title = title,
        useGroups = useGroups,
        tstart = tstart,
        logFile = logFile) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
    }
    
    ##########################################################
    # Single-cell Tracks
    ##########################################################
    if("sctrack" %in% tolower(plotSummary)){
      .logDiffTime(sprintf("Adding SC Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$sctrack <- .scTracks(
        ArchRProj = ArchRProj, 
        region = region[x], 
        tileSize = tileSize, 
        groupBy = groupBy,
        threads = threads, 
        minCells = 5,
        maxCells = scCellsMax,
        pal = pal,
        baseSize = baseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth,
        scTileSize = scTileSize,
        facetbaseSize = facetbaseSize,
        geneAnnotation = geneAnnotation,
        title = title,
        useGroups = useGroups,
        tstart = tstart,
        logFile = logFile) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
    }

    ##########################################################
    # Feature Tracks
    ##########################################################
    if("featuretrack" %in% tolower(plotSummary)){
      if(!is.null(features)){
        .logDiffTime(sprintf("Adding Feature Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
        plotList$featuretrack <- .featureTracks(
            features = features, 
            region = region[x], 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            title = "Peaks",
            logFile = logFile) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
      }
    }

    ##########################################################
    # Loop Tracks
    ##########################################################
    if("looptrack" %in% tolower(plotSummary)){
      if(!is.null(loops)){
        .logDiffTime(sprintf("Adding Loop Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
        plotList$looptrack <- .loopTracks(
            loops = loops, 
            region = region[x], 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            hideY = TRUE,
            title = "Loops",
            logFile = logFile) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
      }
    }

    ##########################################################
    # Gene Tracks
    ##########################################################
    if("genetrack" %in% tolower(plotSummary)){
      .logDiffTime(sprintf("Adding Gene Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$genetrack <- .geneTracks(
        geneAnnotation = geneAnnotation, 
        region = region[x], 
        facetbaseSize = facetbaseSize,
        title = "Genes",
        logFile = logFile) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
    }

    ##########################################################
    # Time to plot
    ##########################################################
    plotSummary <- tolower(plotSummary)
    names(sizes) <- plotSummary
    sizes <- sizes[order(plotSummary)]
    plotSummary <- plotSummary[order(plotSummary)]

    # nullSummary <- unlist(lapply(seq_along(plotSummary), function(x) is.null(eval(parse(text=paste0("plotList$", plotSummary[x]))))))
    # if(any(nullSummary)){
    #   sizes <- sizes[-which(nullSummary)]
    # }
    sizes <- sizes[tolower(names(plotList))]

    if(!is.null(useMatrix)){

      suppressWarnings(.combinedFeaturePlot(
        plotList = plotList,
        log2Norm = log2Norm,
        featureMat = featureMat,
        feature = region[x]$symbol[[1]],
        useMatrix = useMatrix,
        pal = pal,
        sizes = sizes,
        baseSize = baseSize,
        facetbaseSize = facetbaseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth
      ))

    }else{

      .logThis(names(plotList), sprintf("(%s of %s) names(plotList)",x,length(region)), logFile=logFile)
      .logThis(sizes, sprintf("(%s of %s) sizes",x,length(region)), logFile=logFile)
      #.logThis(nullSummary, sprintf("(%s of %s) nullSummary",x,length(region)), logFile=logFile)
      .logDiffTime("Plotting", t1=tstart, verbose=verbose, logFile=logFile)
      
      tryCatch({
        suppressWarnings(ggAlignPlots(plotList = plotList, sizes=sizes, draw = FALSE))
      }, error = function(e){
        .logMessage("Error with plotting, diagnosing each element", verbose = TRUE, logFile = logFile)
        for(i in seq_along(plotList)){
          tryCatch({
            print(plotList[[i]])
          }, error = function(f){
            .logError(f, fn = names(plotList)[i], info = "", errorList = NULL, logFile = logFile)
          })
        }
        .logError(e, fn = "ggAlignPlots", info = "", errorList = NULL, logFile = logFile)
      })

    }

  })

  if(!is.null(mcols(region)$symbol)){
    names(ggList) <- mcols(region)$symbol
  }else{
    if(length(ggList) == 1){
      ggList <- ggList[[1]]
    }
  }

  .endLogging(logFile=logFile)

  ggList

}
    
#######################################################
# Bulk Aggregated ATAC Track Methods
#######################################################
.bulkTracks <- function(
  ArchRProj = NULL, 
  region = NULL, 
  tileSize = 100, 
  minCells = 25,
  groupBy = "Clusters",
  useGroups = NULL,
  normMethod = "ReadsInTSS",
  threads = 1, 
  ylim = NULL,
  baseSize = 7,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  pal = NULL,
  tstart = NULL,
  verbose = FALSE,
  logFile = NULL
  ){

  .requirePackage("ggplot2", source = "cran")

  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  df <- .groupRegionSumArrows(
    ArchRProj = ArchRProj, 
    groupBy = groupBy, 
    normMethod = normMethod,
    useGroups = useGroups,
    minCells = minCells,
    region = region, 
    tileSize = tileSize, 
    threads = threads,
    verbose = verbose,
    logFile = logFile
  )
  .logThis(split(df, df[,3]), ".bulkTracks df", logFile = logFile)

  ######################################################
  # Plot Track
  ######################################################
  if(!is.null(ylim)){
    ylim <- quantile(df$y, ylim)
    df$y[df$y < ylim[1]] <- ylim[1]
    df$y[df$y > ylim[2]] <- ylim[2]
  }else{
    ylim <- c(0,quantile(df$y, probs=c(0.999)))
    df$y[df$y < ylim[1]] <- ylim[1]
    df$y[df$y > ylim[2]] <- ylim[2]
  }
  uniqueGroups <- gtools::mixedsort(unique(paste0(df$group)))
  if(!is.null(useGroups)){
    uniqueGroups <- unique(useGroups)
  }
  df$group <- factor(df$group, levels = uniqueGroups)
  title <- paste0(as.character(seqnames(region)),":", start(region)-1, "-", end(region), " ", title)

  allGroups <- gtools::mixedsort(unique(getCellColData(ArchRProj = ArchRProj, select = groupBy, drop = TRUE)))

  if(is.null(pal)){
    pal <- suppressWarnings(paletteDiscrete(values = allGroups))
  }
  
  #Plot Track
  p <- ggplot(df, aes_string("x","y", color = "group", fill = "group")) + 
    geom_area(stat = "identity") + 
    facet_wrap(facets = ~group, strip.position = 'right', ncol = 1) +
    ylab(sprintf("Coverage\n(Norm. ATAC Signal Range (%s-%s) by %s)", round(min(ylim),2), round(max(ylim),2), normMethod)) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
    scale_y_continuous(limits = ylim, expand = c(0,0)) +
    theme_ArchR(baseSize = baseSize,
                baseRectSize = borderWidth,
                baseLineSize = tickWidth,
                legendPosition = "right",
                axisTickCm = 0.1) +
    theme(panel.spacing= unit(0, "lines"),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text = element_text(
            size = facetbaseSize, 
            color = "black", 
            margin = margin(0,0.35,0,0.35, "cm")),
            strip.text.y = element_text(angle = 0),
          strip.background = element_rect(color="black")) +
    guides(fill = "none", colour = "none") + ggtitle(title)

  p

}

##############################################################################
# Create Average Tracks from Arrows
##############################################################################
.groupRegionSumArrows <- function(
  ArchRProj = NULL,
  useGroups = NULL,
  groupBy = NULL,
  region = NULL,
  tileSize = NULL,
  normMethod = NULL,
  verbose = FALSE,
  minCells = 25,
  maxCells = 500,
  threads = NULL,
  logFile = NULL
  ){

  #Group Info
  cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  if(!is.null(minCells)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% names(table(cellGroups)[table(cellGroups) >= minCells]),,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  if(!is.null(useGroups)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% useGroups,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  tabGroups <- table(cellGroups)
  
  if(any(tabGroups > maxCells)){
    cellGroups2 <- getCellColData(ArchRProj, groupBy, drop = FALSE)
    splitGroups <- split(rownames(cellGroups2), cellGroups2[,1])
    useCells <- lapply(seq_along(splitGroups), function(x){
      if(length(splitGroups[[x]]) > maxCells){
        sample(splitGroups[[x]], maxCells)
      }else{
        splitGroups[[x]]
      }
    }) %>% unlist
    ArchRProj@cellColData <- ArchRProj@cellColData[useCells,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
    tabGroups <- table(cellGroups)
  }

  cellsBySample <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj, "Sample", drop = TRUE))
  groupsBySample <- split(cellGroups, getCellColData(ArchRProj, "Sample", drop = TRUE))
  uniqueGroups <- gtools::mixedsort(unique(cellGroups))
  
  #Tile Region
  regionTiles <- seq(trunc(start(region) / tileSize), trunc(end(region) / tileSize) + 1) * tileSize
  ArrowFiles <- getArrowFiles(ArchRProj)
  ArrowFiles <- ArrowFiles[names(cellsBySample)]

  groupMat <- .safelapply(seq_along(ArrowFiles), function(i){
    .logMessage(sprintf("Getting Region From Arrow Files %s of %s", i, length(ArrowFiles)), logFile = logFile)
    tryCatch({
      .regionSumArrows(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
    }, error = function(e){
      errorList <- list(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
      .logError(e, fn = ".groupRegionSumArrows", info = .sampleName(ArrowFiles[i]), errorList = errorList, logFile = logFile)
    })
  }, threads = threads) %>% Reduce("+" , .)

  #Plot DF
  df <- data.frame(which(groupMat > 0, arr.ind=TRUE))
  df$y <- groupMat[cbind(df[,1], df[,2])]

  #Minus 1 Tile Size
  dfm1 <- df
  dfm1$row <- dfm1$row - 1
  dfm1$y <- 0

  #Plus 1 Size
  dfp1 <- df
  dfp1$row <- dfp1$row + 1
  dfp1$y <- 0

  #Create plot DF
  df <- rbind(df, dfm1, dfp1)
  df <- df[!duplicated(df[,1:2]),]
  df <- df[df$row > 0,]
  df$x <- regionTiles[df$row]
  df$group <- uniqueGroups[df$col]

  #Add In Ends
  dfs <- data.frame(
    col = seq_along(uniqueGroups), 
    row = 1, 
    y = 0,
    x = start(region),
    group = uniqueGroups
  )

  dfe <- data.frame(
    col = seq_along(uniqueGroups),
    row = length(regionTiles),
    y = 0,
    x = end(region),
    group = uniqueGroups
  )
  
  #Final output
  plotDF <- rbind(df,dfs,dfe)
  plotDF <- df[order(df$group,df$x),]
  plotDF <- df[,c("x", "y", "group")]
  
  #Normalization 
  g <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  if(tolower(normMethod) == "readsintss"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "readsinpromoter"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "nfrags"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "ncells"){
      groupNormFactors <- table(g)
  }else if(tolower(normMethod) == "none"){
      groupNormFactors <- rep(10^4, length(g))
      names(groupNormFactors) <- g
  }else{
    stop("Norm Method Not Recognized : ", normMethod)
  }

  #Scale with Norm Factors
  scaleFactors <- 10^4 / groupNormFactors
  matchGroup <- match(paste0(plotDF$group), names(scaleFactors))
  plotDF$y <- plotDF$y * as.vector(scaleFactors[matchGroup])

  return(plotDF)

}

.regionSumArrows <- function(
  ArrowFile = NULL,
  region = NULL,
  regionTiles = NULL,
  tileSize = NULL,
  cellNames = NULL,
  cellGroups = NULL,
  uniqueGroups = NULL,
  logFile = NULL
  ){
  
  cellFragsRegion <- .getFragsFromArrow(
      ArrowFile = ArrowFile, 
      chr = paste0(seqnames(region)), 
      cellNames = cellNames, 
      out = "GRanges"
    ) %>% subsetByOverlaps(., region, ignore.strand = TRUE)
  
  #Starts
  ts <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ids <- which(ts > 0)
  
  #Ends
  te <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ide <- which(te > 0)
  
  #Match
  matchID <- S4Vectors::match(mcols(cellFragsRegion)$RG, cellNames)
  
  #Sparse Matrix
  mat <- Matrix::sparseMatrix(
    i = c(ts[ids], te[ide]),
    j = c(matchID[ids], matchID[ide]),
    x = rep(1,  length(ids) + length(ide)),
    dims = c(length(regionTiles), length(cellNames))
  )
  colnames(mat) <- cellNames
  
  mat@x[mat@x > 1] <- 1

  #Create Group Matrix
  groupMat <- matrix(0, nrow = length(regionTiles), ncol = length(uniqueGroups))
  colnames(groupMat) <- uniqueGroups
  uniqueGroups <- uniqueGroups[uniqueGroups %in% unique(cellGroups)]
  for(i in seq_along(uniqueGroups)){
    groupMat[,uniqueGroups[i]] <- Matrix::rowSums(mat[,which(cellGroups == uniqueGroups[i]),drop=FALSE])
  }

  return(groupMat)

}

#######################################################
# Gene Tracks
#######################################################
.geneTracks <- function(
  geneAnnotation = NULL, 
  region = NULL, 
  baseSize = 9, 
  borderWidth = 0.4, 
  title = "Genes",
  geneWidth = 2, 
  exonWidth = 4, 
  labelSize = 2,
  facetbaseSize,
  colorMinus = "dodgerblue2",
  colorPlus = "red",
  logFile = NULL
  ){

  .requirePackage("ggplot2", source = "cran")
  .requirePackage("ggrepel", source = "cran")

  #only take first region
  region <- .validGRanges(region)
  region <- .subsetSeqnamesGR(region[1], as.character(seqnames(region[1])))

  genes <- sort(sortSeqlevels(geneAnnotation$genes), ignore.strand = TRUE)
  exons <- sort(sortSeqlevels(geneAnnotation$exons), ignore.strand = TRUE)
  genesO <- data.frame(subsetByOverlaps(genes, region, ignore.strand = TRUE))

  if(nrow(genesO) > 0){

    #Identify Info for Exons and Genes
    exonsO <- data.frame(subsetByOverlaps(exons, region, ignore.strand = TRUE))
    exonsO <- exonsO[which(exonsO$symbol %in% genesO$symbol),]
    genesO$facet = title
    genesO$start <- matrixStats::rowMaxs(cbind(genesO$start, start(region)))
    genesO$end <- matrixStats::rowMins(cbind(genesO$end, end(region)))

    #Collapse Iteratively
    #backwards iteration so that the last value chosen is the lowest cluster possible to fit in.
    genesO$cluster <- 0
    for(i in seq_len(nrow(genesO))){
      if(i==1){
        genesO$cluster[i] <- 1
      }else{
        for(j in seq_len(max(genesO$cluster))){
          jEnd <- rev(genesO$end)[match(rev(seq_len(max(genesO$cluster)))[j], rev(genesO$cluster))]
          if(genesO$start[i] > jEnd + median(genesO$width)){
            genesO$cluster[i] <- rev(genesO$cluster)[match(rev(seq_len(max(genesO$cluster)))[j],rev(genesO$cluster))]
          }
        }
        if(genesO$cluster[i]==0){
          genesO$cluster[i] <- genesO$cluster[i-1] + 1
        }
      }
    }
    exonsO$cluster <- genesO$cluster[match(exonsO$symbol, genesO$symbol)]
    pal <- c("-"=colorMinus,"+"=colorPlus,"*"=colorPlus)
    
    p <- ggplot(data = genesO, aes(color = strand, fill = strand)) +
      facet_grid(facet~.) +
      #################################################
      #Limits
      #################################################
      ylim(c(0.5, max(genesO$cluster) + 0.5)) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) + 
      #################################################
      #Segment for Not Minus Stranded
      #################################################
      geom_segment(data = genesO[which(as.character(genesO$strand)!="-"),], 
        aes(x = start, xend = end, y = cluster, yend = cluster, color = strand),size=geneWidth) +
      #################################################
      #Segment for Minus Stranded
      #################################################
      geom_segment(data = genesO[which(as.character(genesO$strand)=="-"),], 
        aes(x = end, xend = start, y = cluster, yend = cluster, color = strand),size=geneWidth) +
      #################################################
      #Segement for Exons
      #################################################
      geom_segment(data = exonsO, aes(x = start, xend = end, y = cluster, 
        yend = cluster, color = strand),size=exonWidth) +
      #################################################
      #Colors
      #################################################
      scale_color_manual(values = pal, guide = FALSE) + 
      scale_fill_manual(values = pal) +
      #################################################
      #Theme
      #################################################
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
      theme(legend.text = element_text(size = baseSize), strip.text.y = element_text(size = facetbaseSize, angle = 0)) +
      guides(fill = guide_legend(override.aes = list(colour = NA, shape = "c", size=3)), color = "none") + 
      theme(legend.position="bottom") +
      theme(legend.title=element_text(size=5), legend.text=element_text(size=7),
        legend.key.size = unit(0.75,"line"), legend.background = element_rect(color =NA), strip.background = element_blank())

    #Add Labels if There are Genes with this orientation!
    if(length(which(genesO$strand!="-")) > 0){
      p <- p + ggrepel::geom_label_repel(data=genesO[which(genesO$strand!="-"),], 
        aes(x = start, y = cluster, label = symbol, color = strand), 
          segment.color = "grey", nudge_x = -0.01*(end(region) - start(region)), nudge_y = -0.25, 
          size = labelSize, direction = "x", inherit.aes=FALSE)
    }

    #Add Labels if There are Genes with this orientation!
    if(length(which(genesO$strand=="-")) > 0){
      p <- p + ggrepel::geom_label_repel(data=genesO[which(genesO$strand=="-"),], 
        aes(x = end, y = cluster, label = symbol, color = strand), 
          segment.color = "grey", nudge_x = +0.01*(end(region) - start(region)), nudge_y = 0.25, 
          size = labelSize, direction = "x", inherit.aes=FALSE)
    }

    p <- p + theme(legend.justification = c(0, 1), 
      legend.background = element_rect(colour = NA, fill = NA), legend.position="none")

  }else{

    #create empty plot
    df <- data.frame(facet = "GeneTrack", start = 0, end = 0, strand = "*", symbol = "none")
    pal <- c("*"=colorPlus)
    p <- ggplot(data = df, aes(start, end, fill = strand)) + geom_point() +
      facet_grid(facet~.) +
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      scale_color_manual(values = pal) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  }

  if(!is.ggplot(p)){
    .logError("geneTrack is not a ggplot!", fn = ".geneTracks", info = "", errorList = NULL, logFile = logFile)
  }

  return(p)

}

#######################################################
# Feature Tracks
#######################################################
.featureTracks <- function(
  features = NULL, 
  region = NULL, 
  title = "FeatureTrack", 
  pal = NULL,
  baseSize = 9, 
  facetbaseSize = NULL,
  featureWidth = 2, 
  borderWidth = 0.4, 
  hideX = FALSE, 
  hideY = FALSE,
  logFile = NULL
  ){

  .requirePackage("ggplot2", source = "cran")

  #only take first region
  region <- .validGRanges(region)
  region <- .subsetSeqnamesGR(region[1], as.character(seqnames(region[1])))

  if(!is.null(features)){

    if(!.isGRList(features)){
      features <- .validGRanges(features)
      featureList <- SimpleList(FeatureTrack = features)
      hideY <- TRUE
    }else{
      featureList <- features
      hideY <- FALSE
    }

    #make sure all elements in featureList have a name for plot display
    for(i in seq_along(featureList)){
      if(is.null(names(featureList)[i]) || is.na(names(featureList)[i]) || nchar(names(featureList)[i]) == 0) {
        message("Warning! Object ",i," in your GRangesList (features) is not named. Generic numbering will be used.")
        names(featureList)[i] <- as.character(i)
      }
    }

    featureList <- featureList[rev(seq_along(featureList))]

    featureO <- lapply(seq_along(featureList), function(x){
      featurex <- featureList[[x]]
      namex <- names(featureList)[x]
      if(is.null(namex) || namex == "") {
        message("Warning! Object ",x," in your GRangesList (features) is not named. Generic numbering will be used.")
        namex <- as.character(x)
      }
      mcols(featurex) <- NULL
      sub <- subsetByOverlaps(featurex, region, ignore.strand = TRUE)
      if(length(sub) > 0){
        data.frame(sub, name = namex)
      }else{
        empty <- GRanges(as.character(seqnames(region[1])), ranges = IRanges(0,0))
        data.frame(empty, name = namex)
      }

    })

    featureO <- Reduce("rbind", featureO)
    
    .logThis(featureO, "featureO", logFile = logFile)

    featureO$facet <- title

    if(is.null(pal)){
      pal <- paletteDiscrete(set = "stallion", values = rev(unique(paste0(featureO$name))))
    }
    
    featureO$name <- factor(paste0(featureO$name), levels=names(featureList))

    p <- ggplot(data = featureO, aes(color = name)) +
      facet_grid(facet~.) +
      geom_segment(data = featureO, aes(x = start, xend = end, y = name, yend = name, color = name), size=featureWidth) +
      ylab("") + xlab("") + 
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      scale_color_manual(values = pal) +
      theme(legend.text = element_text(size = baseSize)) + 
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      guides(color = "none", fill = "none") + theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank())

  }else{

    #create empty plot
    df <- data.frame(facet = "FeatureTrack", start = 0, end = 0, strand = "*", symbol = "none")
    p <- ggplot(data = df, aes(start, end)) + 
      geom_point() +
      facet_grid(facet~.) +
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  }

  if(hideX){
    p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }

  if(hideY){
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  if(!is.ggplot(p)){
    .logError("featureTrack is not a ggplot!", fn = ".featureTracks", info = "", errorList = NULL, logFile = logFile)
  }

  return(p)

}

#######################################################
# Loop Tracks
#######################################################
.loopTracks <- function(
  loops = NULL, 
  region = NULL, 
  title = "LoopTrack", 
  pal = NULL,
  baseSize = 9, 
  facetbaseSize = 9,
  featureWidth = 2, 
  borderWidth = 0.4, 
  hideX = FALSE, 
  hideY = FALSE,
  logFile = NULL
  ){

  getArchDF <- function(lp, r = 100){
    angles <- seq(pi, 2*pi,length.out=100)
    rx <- (end(lp)-start(lp))/2
    rscale <- r * (rx/max(rx))
    cx <- start(lp) + rx
    if(is.null(mcols(lp)$value)){
      mcols(lp)$value <- 1
    }
    df <- lapply(seq_along(cx), function(z){
      xz <- rx[z]*cos(angles)+cx[z]
      dfz <- DataFrame(x=xz, y=rscale[z]*sin(angles), id=Rle(paste0("l",z)), value = mcols(lp)$value[z])
    }) %>% Reduce("rbind",.)
    return(df)
  }

  if(!is.null(loops)){

    if(is(loops, "GRanges")){
      loops <- SimpleList(Loops = loops)
    }else if(.isGRList(loops)){
    }else{
      stop("Loops is not a GRanges or a list of GRanges! Please supply valid input!")
    }

    valueMin <- min(unlist(lapply(loops, function(x) min(x$value))))
    valueMax <- max(unlist(lapply(loops, function(x) max(x$value))))

    loopO <- lapply(seq_along(loops), function(x){
       subLoops <- subsetByOverlaps(loops[[x]], region, ignore.strand = TRUE, type = "within") 
       if(length(subLoops)>0){
         dfx <- getArchDF(subLoops)
         dfx$name <- Rle(paste0(names(loops)[x]))
         dfx
       }else{
         NULL
       }
    }) %>% Reduce("rbind",.)
    .logThis(loopO, "loopO", logFile = logFile)

    testDim <- tryCatch({
      if(is.null(loopO)){
        FALSE
      }
      if(nrow(loopO) > 0){
        TRUE
      }else{
        FALSE
      }
    }, error = function(x){
      FALSE
    })

    if(testDim){

      loopO$facet <- title
      if(is.null(pal)){
        pal <- colorRampPalette(c("#E6E7E8","#3A97FF","#8816A7","black"))(100)
      }

      p <- ggplot(data = data.frame(loopO), aes(x = x, y = y, group = id, color = value)) + 
        geom_line() +
        facet_grid(name ~ .) +
        ylab("") + 
        coord_cartesian(ylim = c(-100,0)) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        scale_color_gradientn(colors = pal, limits = c(valueMin, valueMax)) +
        theme(legend.text = element_text(size = baseSize)) +
        theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth, legendPosition = "right") +
        theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank(),
          legend.box.background = element_rect(color = NA)) +
        guides(color= guide_colorbar(barwidth = 0.75, barheight = 3))

    }else{

      #create empty plot
      df <- data.frame(facet = "LoopTrack", start = 0, end = 0, strand = "*", symbol = "none")
      p <- ggplot(data = df, aes(start, end)) + 
        geom_point() +
        facet_grid(facet~.) +
        theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
        theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

    }

  }else{

    #create empty plot
    df <- data.frame(facet = "LoopTrack", start = 0, end = 0, strand = "*", symbol = "none")
    p <- ggplot(data = df, aes(start, end)) + 
      geom_point() +
      facet_grid(facet~.) +
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  }

  if(hideX){
    p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }

  if(hideY){
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  if(!is.ggplot(p)){
    .logError("loopTracks is not a ggplot!", fn = ".loopTracks", info = "", errorList = NULL, logFile = logFile)
  }

  return(p)

}

.subsetSeqnamesGR <- function(gr = NULL, names = NULL){
  .validInput(input = gr, name = "gr", valid = c("GRanges"))
  .validInput(input = names, name = "names", valid = c("character"))
  gr <- gr[which(as.character(seqnames(gr)) %in% names),]
  seqlevels(gr) <- as.character(unique(seqnames(gr)))
  return(gr)
}

#######################################################
# scATAC Track Methods
#######################################################

.scTracks <- function(
  ArchRProj = NULL,
  region = NULL,
  tileSize = 100,
  minCells = 5,
  maxCells = 100,
  groupBy = "Clusters",
  useGroups = NULL,
  threads = 1,
  baseSize = 7,
  scTileSize = 0.5,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  pal = NULL,
  tstart = NULL,
  verbose = FALSE,
  logFile = NULL
  ){

  .requirePackage("ggplot2", source = "cran")

  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  #Group Info
  cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  if(!is.null(minCells)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% names(table(cellGroups)[table(cellGroups) >= minCells]),,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  if(!is.null(useGroups)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% useGroups,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  tabGroups <- table(cellGroups)
  
  if(any(tabGroups > maxCells)){
    cellGroups2 <- getCellColData(ArchRProj, groupBy, drop = FALSE)
    splitGroups <- split(rownames(cellGroups2), cellGroups2[,1])
    useCells <- lapply(seq_along(splitGroups), function(x){
      if(length(splitGroups[[x]]) > maxCells){
        sample(splitGroups[[x]], maxCells)
      }else{
        splitGroups[[x]]
      }
    }) %>% unlist
    ArchRProj@cellColData <- ArchRProj@cellColData[useCells,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
    tabGroups <- table(cellGroups)
  }

  cellsBySample <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj, "Sample", drop = TRUE))
  groupsBySample <- split(cellGroups, getCellColData(ArchRProj, "Sample", drop = TRUE))
  uniqueGroups <- gtools::mixedsort(unique(cellGroups))
  
  #Tile Region
  regionTiles <- seq(trunc(start(region) / tileSize), trunc(end(region) / tileSize) + 1) * tileSize
  ArrowFiles <- getArrowFiles(ArchRProj)
  ArrowFiles <- ArrowFiles[names(cellsBySample)]

  groupMat <- .safelapply(seq_along(ArrowFiles), function(i){
    .logMessage(sprintf("Getting Region From Arrow Files %s of %s", i, length(ArrowFiles)), logFile = logFile)
    tryCatch({
      .regionSCArrows(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
    }, error = function(e){
      errorList <- list(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
      .logError(e, fn = ".groupRegionSCArrows", info = .sampleName(ArrowFiles[i]), errorList = errorList, logFile = logFile)
    })
  }, threads = threads) %>% Reduce("cbind" , .)

  groupDF <- data.frame(Matrix::summary(groupMat))
  groupDF$group <- getCellColData(ArchRProj, groupBy, drop = FALSE)[colnames(groupMat)[groupDF$j], 1]
  groupDF <- lapply(split(groupDF, groupDF$group), function(z){
    nz <- tabGroups[z$group[1]]
    nc <- length(unique(z$j))
    idx <- sort(sample(seq_len(nz), nc))
    idx[1] <- 1
    idx[length(idx)] <- nz
    z$y <- idx[match(z$j, unique(z$j))]
    z
  }) %>% Reduce("rbind", .)
  groupDF$bp <- regionTiles[groupDF$i]
  
  if(is.null(pal)){
    pal <- suppressWarnings(paletteDiscrete(values = names(tabGroups)))
  }

  nn <- paste0(names(tabGroups), ":", tabGroups)
  names(nn) <- names(tabGroups)
  groupDF$group2 <- nn[groupDF$group]
  names(pal) <- nn[names(pal)]

  title <- paste0(as.character(seqnames(region)),":", start(region)-1, "-", end(region), " ", title)
  
  #Re-Order
  groupDF$group2 <- factor(
    paste0(groupDF$group2), 
    levels = gtools::mixedsort(unique(paste0(groupDF$group2)))
  )

  p <- ggplot(groupDF, aes(x=bp, y=y, width = tileSize, fill = group2, color = group2)) + 
      geom_tile(size = scTileSize) + 
      facet_grid(group2 ~ ., scales="free_y") + 
      theme_ArchR() + 
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      ylab("Binarized SC Coverage") + 
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme_ArchR(baseSize = baseSize,
                  baseRectSize = borderWidth,
                  baseLineSize = tickWidth,
                  legendPosition = "right",
                  axisTickCm = 0.1) +
      theme(panel.spacing= unit(0, "lines"),
            axis.title.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.text = element_text(
              size = facetbaseSize, 
              color = "black", 
              margin = margin(0,0.35,0,0.35, "cm")),
              strip.text.y = element_text(angle = 0),
            strip.background = element_rect(color="black")) +
      guides(fill = "none", colour = "none") + ggtitle(title)

    p

}

.regionSCArrows <- function(
  ArrowFile = NULL,
  region = NULL,
  regionTiles = NULL,
  tileSize = NULL,
  cellNames = NULL,
  cellGroups = NULL,
  uniqueGroups = NULL,
  logFile = NULL
  ){
  
  cellFragsRegion <- .getFragsFromArrow(
      ArrowFile = ArrowFile, 
      chr = paste0(seqnames(region)), 
      cellNames = cellNames, 
      out = "GRanges"
    ) %>% subsetByOverlaps(., region, ignore.strand = TRUE)
  
  #Starts
  ts <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ids <- which(ts > 0)
  
  #Ends
  te <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ide <- which(te > 0)
  
  #Match
  matchID <- S4Vectors::match(mcols(cellFragsRegion)$RG, cellNames)
  
  #Sparse Matrix
  mat <- Matrix::sparseMatrix(
    i = c(ts[ids], te[ide]),
    j = as.vector(c(matchID[ids], matchID[ide])),
    x = rep(1,  length(ids) + length(ide)),
    dims = c(length(regionTiles), length(cellNames))
  )
  colnames(mat) <- cellNames
  
  mat@x[mat@x > 1] <- 1

  return(mat)

}



####################################
# Combined Feature Plot
####################################

.combinedFeaturePlot <- function(
  plotList = NULL,
  useMatrix = NULL,
  featureMat = NULL,
  log2Norm = FALSE,
  feature = NULL,
  pal = NULL,
  sizes = NULL,
  baseSize = NULL,
  facetbaseSize = NULL,
  borderWidth = NULL,
  tickWidth = NULL
  ){

  .requirePackage("patchwork", installInfo = "devtools::install_github('thomasp85/patchwork')")

  if(is.null(pal)){
    pal <- paletteDiscrete(values=featureMat$Group, set = "stallion")
  }

  if(log2Norm){
    title <- paste0("Log2 ", useMatrix, " : ", feature)
  }else{
    title <- paste0("Raw ", useMatrix, " : ", feature) 
  }

  featurePlot <- ggGroup(
      x = featureMat$Group,
      y = featureMat[,feature],
      groupOrder = gtools::mixedsort(paste0(unique(featureMat$Group))),
      pal = pal
    ) + 
    facet_wrap(x~., ncol=1,scales="free_y",strip.position="right") +
    guides(fill = "none", colour = "none") +
    theme_ArchR(baseSize = baseSize,
              baseRectSize = borderWidth,
              baseLineSize = tickWidth,
              legendPosition = "right",
              axisTickCm = 0.1) +
    theme(panel.spacing= unit(0, "lines"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_text(
          size = facetbaseSize, 
          color = "black", 
          margin = margin(0,0.35,0,0.35, "cm")),
          strip.text.y = element_text(angle = 0),
        strip.background = element_rect(color="black")) +
    theme(plot.margin = unit(c(0.35, 0.15, 0.35, 0.15), "cm")) +
    ggtitle(title)

  if(any(tolower(names(plotList)) %in% "bulktrack")){

    idx <- which(tolower(names(plotList)) == "bulktrack")
    
    p <- plotList[[idx]] + featurePlot + plot_spacer()
    
    plotList[idx] <- NULL
    
    for(i in seq_along(plotList)){
      p <- p + plotList[[i]] + plot_spacer() + plot_spacer()
    }
    
    p <- p + plot_layout(
      ncol = 3,
      widths = c(3, 1, 0.2), 
      heights = sizes
    )

  }else{


    idx <- which(tolower(names(plotList)) == "sctrack")
    
    p <- plotList[[idx]] + featurePlot + plot_spacer()
    
    plotList[idx] <- NULL
    
    for(i in seq_along(plotList)){
      p <- p + plotList[[i]] + plot_spacer() + plot_spacer()
    }
    
    p <- p + plot_layout(
      ncol = 3,
      widths = c(3, 1, 0.2), 
      heights = sizes
    )

  }

  p

}







