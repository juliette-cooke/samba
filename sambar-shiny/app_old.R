library(shiny)
library(sambar)
options(shiny.maxRequestSize=50*1024^2)

ui <- fluidPage(
  fileInput(inputId = "sdatain", label = "Upload sampling data here", multiple = T,
            accept = c(".csv"))
)

server <- function(input, output) {
  observeEvent(input$sdatain, print(input$sdatain))
}

#shinyApp(ui = ui, server = server)