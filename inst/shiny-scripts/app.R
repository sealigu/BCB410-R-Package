library(shiny)

ui <- fluidPage(
  titlePanel("UbiquitinAnalysis: R package for Plotting Protein Relationship"),

  sidebarLayout(
    sidebarPanel(
          textInput(inputId = "UniProt_ID1",
                  label = h3("Input an UniProt ID"),
                  placeholder = NULL),
          textInput(inputId = "UniProt_ID2",
                  label = h3("Input an UniProt ID"),
                  placeholder = NULL),
          submitButton("Submit")
      ),
      mainPanel(
#        imageOutput("plotgraph", width = 500, height = 400),
        imageOutput("plot1", width = "100%", height = "400px"),
        imageOutput("plot2", width = "100%", height = "400px"),
        verbatimTextOutput("percentage")
      )
  )
)

server <- function(input, output) {
  p1 = PlotProteinInteractions()
  p2 = PlotResModification()
  plotgraph = c(p1, p2)
  percentage = GetSimilarPercentage()

  output$plot1 <- renderPlot({ p1 })
  output$plot2 <- renderPlot({ p2 })
#  output$plotgraph <- renderPlot({ plotgraph })

  output$percentage <- renderText({
    paste("Protein Sequences Similar Percentage: ", percentage)
  })
}

shinyApp(ui = ui, server = server)
