#' Run the Shiny Application
#'
#' @importFrom shiny shinyApp runApp
#' 
#' @export
#' @return Launches RShiny app
run_app <- function() {
    app = shinyApp(
      ui = app_ui, 
      server = app_server
    )
    runApp(app)
}