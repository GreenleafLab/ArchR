# Load libraries so they are available
# Run the app through this file.
base::source("ui.R")
base::source("server.R")
shiny::shinyApp(ui:ui, server:shinyServer)
# http://127.0.0.1:6747