library(shiny)

ui <- fluidPage(
  titlePanel("UbiquitinAnalysis: R package for Plotting Protein Relationship"),

  sidebarPanel(
    textInput(inputId = "UniProt_ID",
              label = h3("Input an UniProt ID"),
              placeholder = NULL),
    textInput(inputId = "UniProt_ID",
              label = h3("Input an UniProt ID"),
              placeholder = NULL)
  ),

  mainPanel(
    plotOutput("plot1"),
    plotOutput("plot2"),
    textOutput("percentage")
  )
)

server <- function(input, output) {
  plot1 = PlotProteinInteractions()
  plot2 = PlotResModification()
  percentage = GetSimilarPercentage()
  output$plot1 <- renderPlot({
    plot1
  })
  output$plot2 <- renderPlot({
    plot2
  })
  output$percentage <- renderText({
    paste("Similar Percentage: ", percentage)
  })
}

shinyApp(ui = ui, server = server)
