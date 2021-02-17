#' user_input UI Function
#'
#' @description A shiny Module to manage user data.
#' All multi omic data, with Y in last position, are included in the R6 
#' object obj
#' 
#' @param id,input,output,session Internal parameters for {shiny}.
#' @param obj R6 object to wrap all data from different analysis.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList fluidRow column downloadButton fixedRow
#' fileInput selectInput actionButton radioButtons span icon br
#' @importFrom shinydashboard box
#' @importFrom caret nearZeroVar 
mod_user_input_ui <- function(id){
    ns <- NS(id)
    
    fluidRow(
        col_12(align = "center",
            box(width = 10, title = span(icon("file-import"), "Data input"), 
                solidHeader = TRUE, 
                status = "warning",
            
            col_12(align = "center",
                column(6, 
                    box(width = 12, solidHeader = TRUE, status = "info",
                        title = "Omic 1",
                        fileInput(ns("Omic1"), "", width = "100%"),
                        selectInput(ns("sepInput_Omic1"), 
                                    label = "Column separator:",
                                    choices = c(":" = ":",
                                                ";" = ";",
                                                "," = ",",
                                                "\\t" = "\t"),
                                    selected = ",")
                    )
                ),
              col_6(
                  box(width = 12, solidHeader = TRUE, status = "info",
                      title = "Omic 2",
                      fileInput(ns("Omic2"), "", width = "100%"),
                      selectInput(ns("sepInput_Omic2"), 
                                  label = "Column separator:",
                                  choices = c(":" = ":",
                                              ";" = ";",
                                              "," = ",",
                                              "\\t" = "\t"),
                                  selected = ",")
                )
              ),
              uiOutput(ns("inputs")),
              box(width = 12, collapsed = TRUE, status="info",  
                  column(6, actionButton(ns("addInput"),"Add Input", 
                                         icon = icon("plus-square"))),
                  column(6, actionButton(ns("rmvInput"),"Remove Input", 
                                         icon = icon("minus-square")))
              ),
            ),
            column(6,  offset = 3,
                box(width = 12, solidHeader = TRUE, status = "danger",
                    title = "Samples classes vector",
                    fileInput(ns("OmicClass"), "File", width = "100%"),
                    selectInput(ns("sepInput_OmicClass"), 
                                label = "Column separator:",
                                choices = c(":" = ":",
                                            ";" = ";",
                                            "," = ",",
                                            "\\t" = "\t"),
                                selected = ",")
                  )
            ),
            fixedRow(
                column(6, offset = 3,
                    selectInput(ns("organismDb"), 
                        label = span(("Organism"), icon("paw")), 
                        choices = organismList(), ## see input_utils 
                        selected = "org.Hs.eg.db")
                )
            ),
            fixedRow(
                column(6, offset = 3, 
                    actionButton(ns("start"),"Start Analysis", 
                                 icon = icon("play-circle")))
                ),
                column(6, offset = 3, 
                    uiOutput(ns("errorData"))
                ),
                column(6, offset = 3, 
                    uiOutput(ns("warnData"))
                )
          ),
          box(width = 10, title = "Analysis Results", align = "center",
              solidHeader = TRUE, 
              status = "warning",
              
              ## Download report
              radioButtons(ns('format'), 
                           'Document format', 
                           c('HTML', 'Word'),
                           inline = TRUE),
              
              downloadButton(ns('downloadReport'), 
                             "Generate report"),
              br(),
              ## Download RData results
              downloadButton(ns('downloadRes'), 
                             'Download Results', 
                             class="dlButton")
          )
      )
)}
    
#' user_input Server Function
#' 
#' @param input,output,session Internal parameters for shiny
# #' @param obj R6 object to wrap all data from different analysis.
#' 
#' @importFrom shiny reactiveValues observeEvent renderUI span downloadHandler
#' @importFrom rmarkdown render html_document word_document
#' 
#' @return Launches multi-omic data analysis and saves results.
mod_user_input_server <- function(input, output, session){
    ns <- session$ns

    ids <- reactiveValues()
    ids$nb <- c(1, 2)
    
    ## When user wants to add more input file
    observeEvent(input$addInput,{
      if (max(ids$nb)<2){
        
      }else{
          ids$nb <- c(ids$nb, max(ids$nb)+1)
      }
      output$inputs <- renderUI({
          tagList(
              lapply(3:length(ids$nb), function(i){
                  column(6,
                      box(width = 12, solidHeader = TRUE, status = "info",
                      title = paste0("Omic ", ids$nb[i]), 
                      fileInput(ns(paste0("Omic", ids$nb[i])), 
                                label = "",
                                width = "100%",
                                accept = c("text/csv",
                                          "text/comma-separated-values",
                                          "text/plain",
                                          ".csv")),
                      selectInput(ns(paste0("sepInput_Omic", ids$nb[i])), 
                          label = "Column separator:",
                          choices = c(":" = ":",
                                      ";" = ";",
                                      "," = ",",
                                      "\\t" = "\\t"),
                          selected = ","
                            
                )
              )
            )
          })
        )
      })
    })
    
    ## When user wants to remove input file
    observeEvent(input$rmvInput,{
        if (max(ids$nb)<=3){
            output$inputs <- NULL
            ids$nb <- c(1,2)
        }else{
            ids$nb <- ids$nb[-max(ids$nb)]
            output$inputs <- renderUI({
                tagList(
                    lapply(3:length(ids$nb), function(i){
                        column(6, 
                            box(width = 12, 
                                solidHeader = TRUE, 
                                status = "info",
                                title = paste0("Omic ", ids$nb[i]), 
                                fileInput(ns(paste0("Omic", ids$nb[i])), 
                                      label = "",
                                      width = "100%",
                                      accept = c("text/csv",
                                                "text/comma-separated-values,
                                                text/plain",
                                                ".csv")),
                                selectInput(ns(paste0("sepInput_Omic", 
                                                      ids$nb[i])), 
                                            label = "Column separator:",
                                            choices = c(":" = ":",
                                                ";" = ";",
                                                "," = ",",
                                                "\t" = "\t"),
                                            selected = ","
                            )
                  ))
                })
              )
            })
        }
    })
    
    ## When analysis starts
    observeEvent(input$start,{
        if(is.null(ids$nb))
        {
        }else
        {
            out <- tryCatch(
            {
                obj$data$wholeData <- loadMultiOmicData(input, 
                                                        output, 
                                                        session,
                                                        ids$nb)
            },
              error=function(e) {
                  msg <- paste0(e$message,
                    "\nError - Please check your data file structures.")
                  message(msg)
                  output$errorData <- renderUI({
                    span(msg, style="color:red")
                    })
                  # Choose a return value in case of error
                  return(NULL)
              },
              warning=function(w) 
                {
                    message("Warning in file loading")
                    msg <- "\nWarning - Please check your data file format."
                    output$warnData <- renderUI({
                        span(msg, style="color:red")
                      })
                    # Choose a return value in case of warning
                    return(NULL)
              })
          }
    })
    
    ## Automatic results report
    output$downloadReport <- downloadHandler(
        filename = function() {
            paste('MONIA_analysisReport', sep = '.', switch(
                input$format, HTML = 'html', Word = 'docx'
            ))
      },
      
      content = function(file) {
          src <- normalizePath('report.Rmd')
          
          # temporarily switch to the temp dir, in case you do not have write
          # permission to the current working directory
          owd <- setwd(tempdir())
          on.exit(setwd(owd))
          file.copy(src, 'report.Rmd', overwrite = TRUE)
          
          params <- list(obj = obj, format = input$format)
          out <- render('report.Rmd', switch(
              input$format,
              HTML = html_document(toc = TRUE, 
                                   theme = "united", 
                                   highlight = "tango", 
                                   toc_float = TRUE),
              Word = word_document(),
          ))
          file.rename(out, file)
      }
    )
    
    output$downloadRes <- downloadHandler(
      filename <- function(){
        "Results.RData"
      },
      
      content = function(file) {
        save(obj, file = file)
      }
    )
    
    return(input)
}
 
