#' Run the Shiny Application
#'
#' @importFrom shiny shinyApp runApp
#' @importFrom golem with_golem_options
#' 
#' @export
#' @return Launches RShiny app
run_app <- function(
  ...
) {
  with_golem_options(
      app = shinyApp(
        ui = app_ui, 
        server = app_server
      ), 
      golem_opts = list(...)
  )
}