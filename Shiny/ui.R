library(shinybusy)

# This file contains UI widgets. 

# EMBEDING plotting ----------------------------------------------------------------------
EMBED_panel <- shiny::tabPanel(id="EMBED_panel",
                        
                        shiny::titlePanel(htmltools::h5("scClusters")),
                        shiny::sidebarPanel(
                          shiny::titlePanel(htmltools::h3('EMBEDDING 1', align = 'center')),
                          width = 3,
                          htmltools::h4(''),
                          htmltools::hr(style = "border-color: grey"),
                          
                          shiny::selectizeInput(
                            'matrix_EMBED1_forComparison',
                            label = 'EMBEDDING type',
                            choices =  base::c(EMBEDs_dropdown, matrices_dropdown),
                            selected = NULL
                          ),
                          
                          shiny::conditionalPanel(
                            condition = '!(input.matrix_EMBED1_forComparison %in% EMBEDs_dropdown)',
                            shiny::selectizeInput(
                              'EMBED1_forComparison',
                              label = 'EMBEDDING 1',
                              choices = "",
                              selected = NULL
                            )),
                          
                          shiny::splitLayout(cellWidths = c("30%","30%","40%"),
                                      shiny::numericInput("EMBED1_plot_width", "Width", min = 0, max = 250, value = 8),
                                      shiny::numericInput("EMBED1_plot_height", "Height", min = 0, max = 250, value = 12),
                                      shiny::selectizeInput(
                                        'plot_choice_download_EMBED1',
                                        label = "Format",
                                        choices = c(".pdf",".png",".tiff"),
                                        selected = ".pdf"),
                                      tags$head(tags$style(htmltools::HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
                          ),
                          
                          shiny::downloadButton(outputId = "download_EMBED1", label = "Download EMBEDDING 1"),
                          
                          shiny::titlePanel(htmltools::h3('EMBEDDING 2', align = 'center')),
                          htmltools::hr(style = "border-color: grey"),
                          shiny::selectizeInput(
                            'matrix_EMBED2_forComparison',
                            label = 'EMBEDDING type',
                            choices = base::c(EMBEDs_dropdown, matrices_dropdown),
                            selected =NULL
                          ),
                          
                          shiny::conditionalPanel(condition = '!(input.matrix_EMBED2_forComparison %in% EMBEDs_dropdown)',
                                           shiny::selectizeInput(
                                             'EMBED2_forComparison',
                                             label = 'EMBEDDING 2',
                                             choices ="",
                                             selected = NULL
                                           )),
                          
                          
                          shiny::splitLayout(cellWidths = c("30%","30%","40%"),
                                      shiny::numericInput("EMBED2_plot_width", "Width", min = 0, max = 250, value = 8),
                                      shiny::numericInput("EMBED2_plot_height", "Height", min = 0, max = 250, value = 12),
                                      shiny::selectizeInput(
                                        'plot_choice_download_EMBED2',
                                        label = "Format",
                                        choices = base::c(".pdf",".png",".tiff"),
                                        selected = ".pdf"),
                                      tags$head(tags$style(htmltools::HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
                          ),
                          shiny::downloadButton(outputId = "download_EMBED2", label = "Download EMBEDDING 2"),
                          
                        ),
                        
                        shiny::mainPanel(
                          shiny::verbatimTextOutput("feat"),
                          shiny::verbatimTextOutput("text"),
                          shiny::fluidRow(htmltools::h5("Dimension Reduction scClusters EMBEDs"
                          )),
                          shiny::fluidRow(shiny::helpText("Users can view and compare side-by-side EMBEDs' representing identified scATAC-seq clusters,
                            origin of sample, unconstrained and constrained integration with scRNA-seq datasets, and integrated remapped clusters.", style = "font-family: 'Helvetica Now Display Bold'; font-si20pt"),
                          ),
                          shiny::fluidRow(
                            shiny::column(6,shiny::plotOutput("EMBED_plot_1")),  ##%>% withSpinner(color="#0dc5c1")
                            shiny::column(6,shiny::plotOutput("EMBED_plot_2"))
                          )
                        )
)

# Plot Browser:scATAC Clusters --------------------------------------------------------

scATACbrowser_panel <- shiny::tabPanel(
  
  shiny::titlePanel(htmltools::h5("scATAC-seq peak browser")),
  
  shiny::sidebarPanel(
    shiny::titlePanel(htmltools::h5('Gene Name', align = 'center')),
    width = 3,
    htmltools::h4(''),
    htmltools::hr(style = "border-color: grey"),
    
    shiny::actionButton(inputId = "restartButton", label = "Plot Track", icon = shiny::icon("play-circle")),
    
    
    shiny::checkboxGroupInput(inputId = "selectPlotSummary", label = "Select track plots",
                       choices = c("Feature" = "featureTrack", "Loop" = "loopTrack", "Gene" = "geneTrack"),
                       selected = c("featureTrack", "loopTrack", "geneTrack"),
                       inline = TRUE),
    
    shiny::selectizeInput(
      'browserContent',
      label = 'Type',
      choices = EMBEDs_dropdown,
      selected = EMBEDs_dropdown[1]
    ),
    
    shiny::selectizeInput(
      'gene_name',
      label = 'Gene Name',
      choices = base::sort(gene_names),
      selected = base::sort(base::sort(gene_names))[1]
    ),
    
    shiny::sliderInput("range", "Distance From Center (kb):", min = -250, max = 250, value = c(-50,50)),
    shiny::splitLayout(cellWidths = c("50%","50%"),
                shiny::numericInput("range_min", "Distance (-kb):", min = -250, max = 250, value = -50),
                shiny::numericInput("range_max", "Distance (+kb):", min = -250, max = 250, value = 50)
    ),
    shiny::splitLayout(cellWidths = c("50%","50%"),
                shiny::numericInput("tile_size", "TileSize:", min = 10, max = 5000, value = 250),
                shiny::numericInput("ymax", "Y-Max (0,1):", min = 0, max = 1, value = 0.99)
    ),
    
    htmltools::hr(style = "border-color: grey"),
    
    shiny::splitLayout(cellWidths = c("30%","30%","40%"),
                shiny::numericInput("plot_width", "Width", min = 0, max = 250, value = 8),
                shiny::numericInput("plot_height", "Height", min = 0, max = 250, value = 12),
                shiny::selectizeInput(
                  'plot_choice_download_peakBrowser',
                  label = "Format",
                  choices = c(".pdf",".png",".tiff"),
                  selected = ".pdf"),
                tags$head(tags$style(htmltools::HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
    ),
    shiny::downloadButton(outputId = "down", label = "Download"),
    
  ),
  
  shiny::mainPanel(shiny::fluidRow(htmltools::h5("Peak browser of scATAC-seq clusters"
  )),
  shiny::plotOutput("browser_atacClusters")
  )
)

ui <- shiny::shinyUI(shiny::fluidPage(
  shinybusy::add_busy_spinner(spin = "radar", color = "#CCCCCC", onstart = TRUE, height = "55px", width = "55px"),
  
  shiny::navbarPage( 
    EMBED_panel,
    scATACbrowser_panel,
    title ="ShinyArchR Export",
    tags$head(tags$style(".shiny-output-error{color: grey;}"))
  ),
  
  tags$footer(htmltools::HTML("<p><i>This webpage was made using</i> <a href='https://github.com/GreenleafLab/ArchR' target=\"_blank\">ArchR Browser</a>.</p>"),
              align = "left", style = "
              position:relative;
              bottom:0;
              color: black;
              padding: 10px;
              z-index: 1000;")
)
)
