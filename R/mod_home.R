#' #' home UI Function
#' #'
#' #' @description A shiny Module.
#' #'
#' #' @param id,input,output,session Internal parameters for {shiny}.
#' #'
#' #' @noRd 
#' #'
#' #' @importFrom shiny NS tagList 
#' mod_home_ui <- function(id){
#'   ns <- NS(id)
#'   
#'   mod_user_input_ui(ns("user_input_ui_1"))
#'   # mod_user_output_ui("user_output_ui_1")
#' }
#'     
#' #' home Server Function
#' #'
#' #' @noRd 
#' mod_home_server <- function(input, output, session, obj){
#'   ns <- session$ns
#'  
#'   multiOmicData <- callModule(mod_user_input_server, "user_input_ui_1")
#'   obj$data$wholeData <- multiOmicData
#'   
#'   # callModule(mod_user_output_server, "user_output_ui_1")
#'   
#'   return(input)
#' }
#'     
#' ## To be copied in the UI
#' # mod_home_ui("home_ui_1")
#'     
#' ## To be copied in the server
#' # callModule(mod_home_server, "home_ui_1")
#'  
