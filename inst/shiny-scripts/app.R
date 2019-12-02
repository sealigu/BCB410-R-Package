library(shiny)

ui <- fluidPage(
  titlePanel("UbiquitinAnalysis: R package for Plotting Protein Relationship"),

  sidebarLayout(
    sidebarPanel(
          textInput(inputId = "UP1",
                  label = h3("Input an UniProt ID"),
                  placeholder = NULL),
          textInput(inputId = "UP2",
                  label = h3("Input an UniProt ID"),
                  placeholder = NULL),
          submitButton("Submit")
      ),
      mainPanel(
        imageOutput("plot1", width = "100%", height = "400px"),
        imageOutput("plot2", width = "100%", height = "400px"),
        verbatimTextOutput("percentage")
      )
  )
)

server <- function(input, output) {

  output$plot1 <- renderPlot({
      UbiquitinAnalysis::PlotProteinInteractions()
  })

  output$plot2 <- renderPlot({
    UbiquitinAnalysis::PlotResModification()
  })

  output$percentage <- renderText({
    p3 <- UbiquitinAnalysis::GetSimilarPercentage()
    paste("Protein Sequences Similar Percentage: ",
          p3)
  })
}

shinyApp(ui = ui, server = server)
