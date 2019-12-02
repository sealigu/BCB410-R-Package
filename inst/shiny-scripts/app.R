library(shiny)
library(d3heatmap)

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
      d3heatmapOutput("plot1", width = "100%", height = "400px"),
      d3heatmapOutput("plot2", width = "100%", height = "400px"),
      verbatimTextOutput("percentage")
    )
  )
)

server <- function(input, output) {

  output$plot1 <- renderD3heatmap({
    validate({
      need(input$UP1, 'Check there is at least one letter!')
      need(input$UP2, 'Check there is at least one letter!')
    })
    p1 <- UbiquitinAnalysis::PlotProteinInteractions(input$UP1, input$UP2)
    d3heatmap(p1, dendrogram = "none", colors = c("white", "lightblue"))
  })

  output$plot2 <- renderD3heatmap({
    validate({
      need(input$UP1, 'Check there is at least one letter!')
      need(input$UP2, 'Check there is at least one letter!')
    })
    p2 <- UbiquitinAnalysis::PlotResModification(input$UP1, input$UP2)
    d3heatmap(p2, dendrogram = "none", colors = c("white", "gray88"))
  })

  output$percentage <- renderText({
    validate({
      need(input$UP1, 'Check there is at least one letter!')
      need(input$UP2, 'Check there is at least one letter!')
    })

    p3 <- UbiquitinAnalysis::GetSimilarPercentage(input$UP1, input$UP2)
    paste("Protein Sequences Similar Percentage: ", p3)
  })
}

shinyApp(ui = ui, server = server)
