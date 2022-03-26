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
      # tags$head(
      #   tags$style(HTML('#tile_size{height: 35px}'))
      # ),
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
      # selectizeInput("grouping",
      #                label = "groupBy",
      #                choices = discreteCols,
      #                multiple = FALSE,
      #                options = list(placeholder = 'Select Grouping'),
      #                selected = selectCols
      # ),
      sliderInput("range", "Distance From Center (kb):", min = -250, max = 250, value = c(-50,50)),
      splitLayout(cellWidths = c("50%","50%"),
                  numericInput("range_min", "Distance (-kb):", min = -250, max = 250, value = -50),
                  numericInput("range_max", "Distance (+kb):", min = -250, max = 250, value = 50)
      ),
      splitLayout(cellWidths = c("50%","50%"),
                #  numericInput("tile_size", "TileSize:", min = 10, max = 5000, value = 250),
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