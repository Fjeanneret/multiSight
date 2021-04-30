#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#' @importFrom shiny callModule stopApp
#' @importFrom R6 R6Class 
#' @import S4Vectors
#' @noRd
app_server <- function(input, output, session ) {
    myDataProcess <- R6Class("MyDataProcess", list(
        organismDb = NULL,
        data = list(),
        classification = list(),
        enrichment = list(),
        networkInference = list()
    ))
    obj <- NULL
    obj <<- myDataProcess$new()
    
    # First level callModules
    callModule(mod_dashboard_structure_server, "dashboard_structure_ui_1")
    # startSignal <- callModule(mod_home_server, "home_ui_1", obj)
    startSignal <- callModule(mod_user_input_server, "user_input_ui_1")
    bioDB <- callModule(mod_understand_server, "understand_ui_1", startSignal)
    callModule(mod_MLmodels_server, "MLmodels_ui_1", startSignal, bioDB)
    # callModule(mod_hypothesisGenerator_server, "hypothesisGenerator_ui_1", 
    # obj)
    session$onSessionEnded(stopApp)
}
