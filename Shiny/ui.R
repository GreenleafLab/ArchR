library(shinybusy)

# This file contains UI widgets. 

# Umap plotting ----------------------------------------------------------------------
umap_panel <- tabPanel(id="umap_panel",
                       
                       titlePanel(h5("scClusters")),
                       sidebarPanel(
                         titlePanel(h3('UMAP 1', align = 'center')),
                         width = 3,
                         h4(''),
                         hr(style = "border-color: grey"),
                         
                         selectizeInput(
                           'matrix_UMAP1_forComparison',
                           label = 'UMAP type',
                           choices =  c("Clusters","Constrained","Constrained remap","Sample","Unconstrained","Gene Score Matrix","Gene Integration Matrix","Motif Matrix"),
                           selected ="Clusters"
                         ),
                         
                         conditionalPanel(condition = "input.matrix_UMAP1_forComparison=='Gene Score Matrix' ||input.matrix_UMAP1_forComparison=='Gene Integration Matrix' ||  input.matrix_UMAP1_forComparison=='Motif Matrix'",
                                          selectizeInput(
                                            'UMAP1_forComparison',
                                            label = 'UMAP 1',
                                            choices = "",
                                            selected = NULL
                                          )),
                         
                         splitLayout(cellWidths = c("30%","30%","40%"),
                                     numericInput("UMAP1_plot_width", "Width", min = 0, max = 250, value = 8),
                                     numericInput("UMAP1_plot_height", "Height", min = 0, max = 250, value = 12),
                                     selectizeInput(
                                       'plot_choice_download_UMAP1',
                                       label = "Format",
                                       choices = c(".pdf",".png",".tiff"),
                                       selected = ".pdf"),
                                     tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
                         ),
                         
                         downloadButton(outputId = "download_UMAP1", label = "Download UMAP 1"),
                         
                         titlePanel(h3('UMAP 2', align = 'center')),
                         hr(style = "border-color: grey"),
                         selectizeInput(
                           'matrix_UMAP2_forComparison',
                           label = 'UMAP type',
                           choices = c("Clusters","Constrained","Constrained remap","Sample","Unconstrained","Gene Score Matrix","Gene Integration Matrix","Motif Matrix"),
                           selected ="Clusters"
                         ),
                         
                         conditionalPanel(condition = "input.matrix_UMAP2_forComparison=='Gene Score Matrix' ||input.matrix_UMAP2_forComparison=='Gene Integration Matrix' ||  input.matrix_UMAP2_forComparison=='Motif Matrix'",
                                          selectizeInput(
                                            'UMAP2_forComparison',
                                            label = 'UMAP 2',
                                            choices ="",
                                            selected = NULL
                                          )),
                         
                         
                         splitLayout(cellWidths = c("30%","30%","40%"),
                                     numericInput("UMAP2_plot_width", "Width", min = 0, max = 250, value = 8),
                                     numericInput("UMAP2_plot_height", "Height", min = 0, max = 250, value = 12),
                                     selectizeInput(
                                       'plot_choice_download_UMAP2',
                                       label = "Format",
                                       choices = c(".pdf",".png",".tiff"),
                                       selected = ".pdf"),
                                     tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
                         ),
                         downloadButton(outputId = "download_UMAP2", label = "Download UMAP 2"),
                         
                       ),
                       
                       mainPanel(
                         verbatimTextOutput("text"),
                         fluidRow(h5("Dimension Reduction scClusters UMAPs"
                         )),
                         fluidRow(helpText("Users can view and compare side-by-side UMAPs' representing identified scATAC-seq clusters,
                            origin of sample, unconstrained and constrained integration with scRNA-seq datasets, and integrated remapped clusters.", style = "font-family: 'Helvetica Now Display Bold'; font-si20pt"),
                         ),
                         fluidRow(
                           column(6,plotOutput("UMAP_plot_1")),  ##%>% withSpinner(color="#0dc5c1")
                           column(6,plotOutput("UMAP_plot_2"))
                         )
                       )
)

# Plot Browser:scATAC Clusters --------------------------------------------------------

scATACbrowser_panel <- tabPanel(
  
  titlePanel(h5("scATAC-seq peak browser")),
  
  sidebarPanel(
    titlePanel(h5('Gene Name', align = 'center')),
    width = 3,
    h4(''),
    hr(style = "border-color: grey"),
    
    actionButton(inputId = "restartButton", label = "Plot Track", icon = icon("play-circle")),
    
    
    checkboxGroupInput(inputId = "selectPlotSummary", label = "Select track plots",
                       choices = c("Feature" = "featureTrack", "Loop" = "loopTrack", "Gene" = "geneTrack"),
                       selected = c("featureTrack", "loopTrack", "geneTrack"),
                       inline = TRUE),
    
    selectizeInput(
      'browserContent',
      label = 'Type',
      choices = c("Unconstrained","Constrained"),
      selected = "Unconstrained"
    ),
    
    selectizeInput(
      'gene_name',
      label = 'Gene Name',
      choices = sort(gene_names),
      selected = sort(sort(gene_names))[1]
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
    
    hr(style = "border-color: grey"),
    
    splitLayout(cellWidths = c("30%","30%","40%"),
                numericInput("plot_width", "Width", min = 0, max = 250, value = 8),
                numericInput("plot_height", "Height", min = 0, max = 250, value = 12),
                selectizeInput(
                  'plot_choice_download_peakBrowser',
                  label = "Format",
                  choices = c(".pdf",".png",".tiff"),
                  selected = ".pdf"),
                tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
    ),
    downloadButton(outputId = "down", label = "Download"),
    
  ),
  
  mainPanel(fluidRow(h5("Peak browser of scATAC-seq clusters"
  )),
  plotOutput("browser_atacClusters")
  )
)

ui <- shinyUI(fluidPage(
  add_busy_spinner(spin = "radar", color = "#CCCCCC", onstart = TRUE, height = "55px", width = "55px"),
  
  navbarPage( 
    umap_panel,
    scATACbrowser_panel,
    title ="ShinyArchR Export",
    tags$head(tags$style(".shiny-output-error{color: grey;}"))
  ),
  
  tags$footer(HTML("<p><i>This webpage was made using</i> <a href='https://github.com/GreenleafLab/ArchR' target=\"_blank\">ArchR Browser</a>.</p>"),
              align = "left", style = "
              position:relative;
              bottom:0;
              color: black;
              padding: 10px;
              z-index: 1000;")
)
)
