####################################################################
# Signal Track Plotting Methods
####################################################################

#' Launch ArchR Genome Browser
#' 
#' This function will open an interactive shiny session in style of a browser track. It allows for normalization of the signal which
#' enables direct comparison across samples.
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
  threads = getArchRThreads()
  ){

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
    message("shinythemes not found! To see a nice theme use :\n\tinstall.packages('shinythemes')\nContinuing wihtout shinythemes!")
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
               choices = c("ReadsInTSS", "nFrags"),
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
            region1 <- resize(region, 1, "start")
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
                tstart = NULL
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
                  title = "Peaks"
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
                title = "Loops") + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
            }

            setProgress(0.6)

            g <- .geneTracks(
              geneAnnotation = geneAnnotation, 
              region = tmpArchRRegion, 
              facetbaseSize = facetbaseSize,
              labelSize = 3,
              title = "Genes"
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
                tstart = NULL
              ) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))

          }else{

            print("Using previous ggplot")

            p <- tmpArchRP

          }

          #p <- p + .suppressAll(scale_x_continuous(limits = c(start(tmpArchRRegion), end(tmpArchRRegion)), expand = c(0,0)))

          setProgress(0.5)

          if(!is.null(features)){

            f <- .featureTracks(
                features = features, 
                region = tmpArchRRegion,
                facetbaseSize = facetbaseSize, 
                hideX = TRUE, 
                title = "Peaks"
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
              title = "Loops") + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
          }

          setProgress(0.7)

          g <- .geneTracks(
            geneAnnotation = geneAnnotation, 
            region = tmpArchRRegion, 
            facetbaseSize = facetbaseSize,
            labelSize = 3,
            title = "Genes"
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

#' Plot an ArchR Region Track
#' 
#' This function will plot the coverage at an input region in the style of a browser track. It allows for normalization of the signal
#' which enables direct comparison across samples.
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
#' "featureTrack" (i.e. the peak regions), and "geneTrack" (line diagrams of genes with introns and exons shown. Blue-colored genes
#' are on the minus strand and red-colored genes are on the plus strand).
#' @param sizes A numeric vector containing up to 3 values that indicate the sizes of the individual components passed to `plotSummary`.
#' The order must be the same as `plotSummary`.
#' @param features A `GRanges` object containing the "features" to be plotted via the "featureTrack". This should be thought of as a
#' bed track. i.e. the set of peaks obtained using `getPeakSet(ArchRProj))`. 
#' @param loops A `GRanges` object containing the "loops" to be plotted via the "loopTrack".
#' This `GRanges` object start represents the center position of one loop anchor and the end represents the center position of another loop anchor. 
#' A "loopTrack" draws an arc between two genomic regions that show some type of interaction. This type of track can be used 
#' to display chromosome conformation capture data or co-accessibility links obtained using `getCoAccessibility()`. 
#' @param geneSymbol If `region` is not supplied, plotting can be centered at the transcription start site corresponding to the gene symbol(s) passed here.
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
#' @param ylim The numeric quantile y-axis limit to be used for for "bulkTrack" plotting. If not provided, the y-axis limit will be c(0, 0.999).
#' @param baseSize The numeric font size to be used in the plot. This applies to all plot labels.
#' @param borderWidth The numeric line width to be used for plot borders.
#' @param tickWidth The numeric line width to be used for axis tick marks.
#' @param facetbaseSize The numeric font size to be used in the facets (gray boxes used to provide track labels) of the plot.
#' @param geneAnnotation The `geneAnnotation` object to be used for plotting the "geneTrack" object. See `createGeneAnnotation()` for more info.
#' @param title The title to add at the top of the plot next to the plot's genomic coordinates.
#' @export
ArchRBrowserTrack <- function(
  ArchRProj = NULL, 
  region = NULL, 
  groupBy = "Clusters",
  useGroups = NULL, 
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(10, 1.5, 3, 4),
  features = getPeakSet(ArchRProj),
  loops = getCoAccessibility(ArchRProj),
  geneSymbol = NULL,
  upstream = 50000,
  downstream = 50000,
  tileSize = 250, 
  minCells = 25,
  normMethod = "ReadsInTSS",
  threads = getArchRThreads(), 
  ylim = NULL,
  baseSize = 7,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = ""
  ){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  .validInput(input = region, name = "region", valid = c("granges","null"))
  .validInput(input = groupBy, name = "groupBy", valid = "character")
  .validInput(input = useGroups, name = "useGroups", valid = c("character", "null"))
  .validInput(input = plotSummary, name = "plotSummary", valid = "character")
  .validInput(input = sizes, name = "sizes", valid = "numeric")
  #.validInput(input = features, name = "features", valid = c("granges", "grangeslist", "null")) #JJJ 
  #.validInput(input = loops, name = "loops", valid = c("granges", "grangeslist", "null")) #JJJ 
  .validInput(input = geneSymbol, name = "geneSymbol", valid = c("character", "null"))
  .validInput(input = upstream, name = "upstream", valid = c("integer"))
  .validInput(input = downstream, name = "downstream", valid = c("integer"))
  .validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  .validInput(input = minCells, name = "minCells", valid = c("integer"))
  .validInput(input = normMethod, name = "normMethod", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = ylim, name = "ylim", valid = c("numeric", "null"))
  .validInput(input = baseSize, name = "baseSize", valid = "numeric")
  .validInput(input = borderWidth, name = "borderWidth", valid = "numeric")
  .validInput(input = tickWidth, name = "tickWidth", valid = "numeric")
  .validInput(input = facetbaseSize, name = "facetbaseSize", valid = "numeric")
  geneAnnotation <- .validGeneAnnotation(geneAnnotation)
  .validInput(input = title, name = "title", valid = "character")

  tstart <- Sys.time()

  ##########################################################
  # Get Region Where Plot Will Occur (GenomicRanges)
  ##########################################################
  .messageDiffTime("Validating Region", tstart)
  if(is.null(region)){
    if(!is.null(geneSymbol)){
      region <- geneAnnotation$genes
      region <- region[which(tolower(mcols(region)$symbol) %in% tolower(geneSymbol))]
      region <- region[order(match(tolower(mcols(region)$symbol), tolower(geneSymbol)))]
      print(region)
      region <- resize(region, 1, "start")
      strand(region) <- "*"
      region <- extendGR(region, upstream = upstream, downstream = downstream)
    }
  }
  region <- .validGRanges(region)

  ggList <- lapply(seq_along(region), function(x){

    plotList <- list()

    ##########################################################
    # Bulk Tracks
    ##########################################################
    if("bulktrack" %in% tolower(plotSummary)){
      .messageDiffTime(sprintf("Adding Bulk Tracks (%s of %s)",x,length(region)), tstart)
      plotList$bulktrack <- .bulkTracks(
        ArchRProj = ArchRProj, 
        region = region[x], 
        tileSize = tileSize, 
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
        title = title,
        useGroups = useGroups,
        tstart = tstart) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
    }
    
    ##########################################################
    # Feature Tracks
    ##########################################################
    if("featuretrack" %in% tolower(plotSummary)){
      if(!is.null(features)){
        .messageDiffTime(sprintf("Adding Feature Tracks (%s of %s)",x,length(region)), tstart)
        plotList$featuretrack <- .featureTracks(
            features = features, 
            region = region[x], 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            title = "Peaks") + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
      }
    }

    ##########################################################
    # Feature Tracks
    ##########################################################
    if("looptrack" %in% tolower(plotSummary)){
      if(!is.null(loops)){
        .messageDiffTime(sprintf("Adding Loop Tracks (%s of %s)",x,length(region)), tstart)
        plotList$looptrack <- .loopTracks(
            loops = loops, 
            region = region[x], 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            hideY = TRUE,
            title = "Loops") + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
      }
    }

    ##########################################################
    # Gene Tracks
    ##########################################################
    if("genetrack" %in% tolower(plotSummary)){
      .messageDiffTime(sprintf("Adding Gene Tracks (%s of %s)",x,length(region)), tstart)
      plotList$genetrack <- .geneTracks(
        geneAnnotation = geneAnnotation, 
        region = region[x], 
        facetbaseSize = facetbaseSize,
        title = "Genes") + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
    }

    ##########################################################
    # Time to plot
    ##########################################################
    plotSummary <- tolower(plotSummary)
    sizes <- sizes[order(plotSummary)]
    plotSummary <- plotSummary[order(plotSummary)]

    nullSummary <- unlist(lapply(seq_along(plotSummary), function(x) is.null(eval(parse(text=paste0("plotList$", plotSummary[x]))))))
    if(any(nullSummary)){
      sizes <- sizes[-which(nullSummary)]
    }

    .messageDiffTime("Plotting", tstart)
    suppressWarnings(ggAlignPlots(plotList = plotList, sizes=sizes, draw = FALSE))

  })

  if(!is.null(mcols(region)$symbol)){
    names(ggList) <- mcols(region)$symbol
  }else{
    if(length(ggList) == 1){
      ggList <- ggList[[1]]
    }
  }

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
  verbose = FALSE
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
    verbose = verbose
  )

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
    guides(fill = FALSE, colour = FALSE) + ggtitle(title)

  p

}

##############################################################################
# Create Average Tracks from Coverages
##############################################################################
.groupRegionSumCoverages <- function(
  ArchRProj = NULL,
  groupBy = NULL,
  useGroups = NULL,
  minCells = 25,
  region = NULL,
  tileSize = NULL,
  normMethod = NULL,
  verbose = NULL,
  threads = NULL
  ){

  coverageMetadata <- .getCoverageMetadata(
    ArchRProj = ArchRProj, 
    groupBy = groupBy, 
    useGroups = useGroups,
    minCells = minCells
  )

  cellGroups <- .getCoverageParams(
    ArchRProj = ArchRProj, 
    groupBy = groupBy
  )[["cellGroups"]] %>% unlist
  
  groupRegionRle <- .groupRegionCoverages(
    coverageMetadata = coverageMetadata, 
    region = region, 
    tileSize = tileSize, 
    buffer = tileSize * 5, 
    threads = threads
  )
  groupNames <- names(groupRegionRle)
  
  #Normalization 
  g <- names(unlist(cellGroups, use.names = TRUE))
  if(tolower(normMethod) == "readsintss"){
      v <- getCellColData(ArchRProj, normMethod, drop = FALSE)[unlist(cellGroups),]
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "nfrags"){
      v <- getCellColData(ArchRProj, normMethod, drop = FALSE)[unlist(cellGroups),]
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "ncells"){
      groupNormFactors <- table(g)
  }else{
    stop("Norm Method Not Recognized : ", normMethod)
  }
  
  #Scale with Norm Factors
  scaleFactors <- 10^4 / groupNormFactors

  #Normalize
  groupRegionRle <- lapply(seq_along(groupRegionRle), function(x){
      groupRegionRle[[x]] * scaleFactors[names(groupRegionRle)[x]]
  })

  #Group And Average
  groupRegionRle <- split(groupRegionRle, coverageMetadata$Group)
  groupRegionList <- lapply(seq_along(groupRegionRle), function(x){
    Reduce("+", groupRegionRle[[x]]) / length(groupRegionRle[[x]])
  })
  names(groupRegionList) <- names(groupRegionRle)

  #Tile Region
  tileSize <- floor(tileSize / 2)
  regionTiles <- seq(trunc(start(region) / tileSize) - 1, trunc(end(region) / tileSize) + 1) * tileSize

  plotDF <- lapply(seq_along(groupRegionList), function(x){
    data.frame(x = regionTiles, y = as.vector(groupRegionList[[x]][regionTiles]), group = names(groupRegionList)[x])
  }) %>% Reduce("rbind", .)

  plotDF

}

.groupRegionCoverages <- function(
  coverageMetadata = NULL,
  region = NULL,
  tileSize = 100,
  buffer = 1000,
  threads = 1
  ){
  
  region <- .validGRanges(region[1])
  coverageFiles <- coverageMetadata$File
  names(coverageFiles) <- coverageMetadata$Name
  
  covList <- .safelapply(seq_along(coverageFiles), function(x){
    .getCoverageFromRegion(coverageFiles[x], region, tileSize, buffer)
  }, threads = threads) %>% {as(.,"RleList")}
  names(covList) <- names(coverageFiles)
  
  covList

}

.getCoverageFromRegion <- function(
  coverageFile = NULL,
  region = NULL,
  tileSize = NULL,
  buffer = NULL
  ){
  chr <- as.character(seqnames(region))
  cov <- Rle(
    lengths = h5read(coverageFile, paste0("Coverage/",chr,"/Lengths")), 
    values = h5read(coverageFile, paste0("Coverage/",chr,"/Values"))
  )
  w <- sum(runLength(cov))
  idx <- cumsum(runLength(cov)) %>% {which(.  >= start(region) - buffer & . <= end(region) + buffer)}
  runValue(cov)[-idx] <- 0
  covRanges <- ranges(cov)
  mcols(covRanges)$values <- runValue(cov)
  covRanges <- covRanges[mcols(covRanges)$values > 0]
  start(covRanges) <- trunc(start(covRanges) / tileSize) * tileSize
  end(covRanges) <- (trunc(end(covRanges) / tileSize) + 1) * tileSize - 1
  o <- coverage(covRanges, weight = mcols(covRanges)$values, width = w)
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
  verbose = NULL,
  minCells = 25,
  threads = NULL
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
  cellsBySample <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj, "Sample", drop = TRUE))
  groupsBySample <- split(cellGroups, getCellColData(ArchRProj, "Sample", drop = TRUE))
  uniqueGroups <- gtools::mixedsort(unique(cellGroups))
  
  #Tile Region
  regionTiles <- seq(trunc(start(region) / tileSize), trunc(end(region) / tileSize) + 1) * tileSize
  ArrowFiles <- getArrowFiles(ArchRProj)
  ArrowFiles <- ArrowFiles[names(cellsBySample)]

  groupMat <- .safelapply(seq_along(ArrowFiles), function(i){
    gmi <- .regionSumArrows(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
    )
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
  }else if(tolower(normMethod) == "nfrags"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "ncells"){
      groupNormFactors <- table(g)
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
  uniqueGroups = NULL
  ){
  
  cellFragsRegion <- .getFragsFromArrow(
      ArrowFile = ArrowFile, 
      chr = paste0(seqnames(region)), 
      cellNames = cellNames, 
      out = "GRanges"
    ) %>% subsetByOverlaps(., region, ignore.strand = FALSE)
  
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
  colorPlus = "red"
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
      guides(fill = guide_legend(override.aes = list(colour = NA, shape = "c", size=3)), color = FALSE) + 
      theme(legend.position="bottom") +
      theme(legend.title=element_text(size=5), legend.text=element_text(size=7),
        legend.key.size = unit(0.75,"line"), legend.background = element_rect(color =NA), strip.background = element_blank())

    #Add Labels if There are Genes with this orientation!
    if(length(which(genesO$strand!="-")) > 0){
      p <- p + ggrepel::geom_label_repel(data=genesO[which(genesO$strand!="-"),], 
        aes(x = start, y = cluster, label = symbol, color = strand, fill = NA), 
          segment.color = "grey", nudge_x = -0.01*(end(region) - start(region)), nudge_y = -0.25, 
          size = labelSize, direction = "x")
    }

    #Add Labels if There are Genes with this orientation!
    if(length(which(genesO$strand=="-")) > 0){
      p <- p + ggrepel::geom_label_repel(data=genesO[which(genesO$strand=="-"),], 
        aes(x = end, y = cluster, label = symbol, color = strand, fill = NA), 
          segment.color = "grey", nudge_x = +0.01*(end(region) - start(region)), nudge_y = 0.25, 
          size = labelSize, direction = "x")
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
  hideY = FALSE
  ){

  .requirePackage("ggplot2", source = "cran")

  #only take first region
  region <- .validGRanges(region)
  region <- .subsetSeqnamesGR(region[1], as.character(seqnames(region[1])))

  if(!is.null(features)){

    if(!inherits(features,"GRangesList") & !inherits(features,"GenomicRangesList")){
      features <- .validGRanges(features)
      featureList <- GenomicRanges::GenomicRangesList(features)
      names(featureList) <- "FeatureTrack"
      hideY <- TRUE
    }else{
      featureList <- features
      hideY <- FALSE
    }
    featureList <- featureList[rev(seq_along(featureList))]

    featureO <- lapply(seq_along(featureList), function(x){
      featurex <- featureList[[x]]
      namex <- names(featureList)[x]
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
    featureO$facet <- title

    if(is.null(pal)){
      pal <- paletteDiscrete(set = "stallion", values = rev(unique(paste0(featureO$name))))
    }

    p <- ggplot(data = featureO, aes(color = name)) +
      facet_grid(facet~.) +
      geom_segment(data = featureO, aes(x = start, xend = end, y = name, yend = name, color = name), size=featureWidth) +
      ylab("") + xlab("") + 
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      scale_color_manual(values = pal) +
      theme(legend.text = element_text(size = baseSize)) + 
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      guides(color = FALSE, fill = FALSE) + theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank())

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
  hideY = FALSE
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

    if(inherits(loops, "GRanges")){
      loops <- GenomicRanges::GenomicRangesList(loops)
      names(loops) <- "Loops" 
    }else if(all(unlist(lapply(loops, function(x) inherits(x, "GRanges"))))){
    }else{
      stop("Loops is not a GRanges or a list of GRanges! Please supply valid input!")
    }

    valueMin <- min(unlist(lapply(loops, function(x) min(x$value))))
    valueMax <- max(unlist(lapply(loops, function(x) max(x$value))))

    loopO <- lapply(seq_along(loops), function(x){
       subLoops <- subsetByOverlaps(loops[[x]], region, ignore.strand = TRUE, type = "within")    
       dfx <- getArchDF(subLoops)
       dfx$name <- Rle(paste0(names(loops)[x]))
       return(dfx)
    }) %>% Reduce("rbind",.)

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

  if(hideX){
    p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }

  if(hideY){
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
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

