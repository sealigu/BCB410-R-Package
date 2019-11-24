#' Launch the shiny app for package UbiquitinAnalysis
#' 
#' A function that launches the shiny app for the package
#' The code has been placed in \code{./inst/shiny-scripts}.
#' 
#' @return No return value but open up a shiny page.
#' 
#' @examples
#' \dontrun{
#' runUA()
#' }
#' 
#' @importFrom shiny runApp

runUA <- function() {
  appDir <- system.file("shiny-scripts",
      package = "UbiquitinAnalysis")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}