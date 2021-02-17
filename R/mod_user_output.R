#' #' user_output UI Function
#' #'
#' #' @description A shiny Module.
#' #'
#' #' @param id,input,output,session Internal parameters for {shiny}.
#' #'
#' #' @noRd 
#' #'
#' #' @importFrom shiny NS tagList 
#' mod_user_output_ui <- function(id){
#'   ns <- NS(id)
#'   box(width = 10, title = "Analysis Results", align = "center",
#'       solidHeader = TRUE, 
#'       status = "warning",
#'       
#'       ## Download report
#'       radioButtons(ns('format'), 'Document format', c('HTML', 'Word'),
#'                    inline = TRUE),
#'       downloadButton(ns('downloadReport'), "Generate report")
#'       
#'       ## Download RData results
#'   )
#' }
#'     
#' #' user_output Server Function
#' #'
#' #' @noRd 
#' mod_user_output_server <- function(input, output, session){
#'   ns <- session$ns
#'  
#'   ## Automatic results report
#'   output$downloadReport <- downloadHandler(
#'     filename = function() {
#'       paste('MONIA_analysisReport', sep = '.', switch(
#'         input$format, HTML = 'html', Word = 'docx'
#'       ))
#'     },
#'     
#'     content = function(file) {
#'       src <- normalizePath('report.Rmd')
#'       
#'       # temporarily switch to the temp dir, in case you do not have write
#'       # permission to the current working directory
#'       owd <- setwd(tempdir())
#'       on.exit(setwd(owd))
#'       file.copy(src, 'report.Rmd', overwrite = TRUE)
#'       
#'       params <- list(obj = obj, format = input$format)
#'       out <- render('report.Rmd', switch(
#'         input$format,
#'         HTML = html_document(toc = TRUE, 
#'                              theme = "united", 
#'                              highlight = "tango", 
#'                              toc_float = TRUE),
#'         Word = word_document(),
#'       ))
#'       file.rename(out, file)
#'     }
#'   )
#' }
#'     
#' ## To be copied in the UI
#' # mod_user_output_ui("user_output_ui_1")
#'     
#' ## To be copied in the server
#' # callModule(mod_user_output_server, "user_output_ui_1")
#'  
