#' hypothesisGenerator UI Function
#'
#' @description A shiny Module to create biological networks.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList column sliderInput 
#' textInput actionButton htmlOutput 
#' @importFrom shinydashboard box
#' @importFrom DT dataTableOutput
#' @importFrom networkD3 forceNetworkOutput
mod_hypothesisGenerator_ui <- function(id){
  ns <- NS(id)
  tagList(
      # box(title ="Diablo features", width = 12, #height = "800px",
        box(title = "Correlation", solidHeader = TRUE, width = 12,
            collapsible = TRUE, collapsed = TRUE, status = "success",
              column(width = 7,
                  forceNetworkOutput(ns("correlationNetwork")),
                  sliderInput(ns("correlationNumThreshold"),
                                 label = "Correlation absolute value threshold",
                                 value = 0.2, 
                                 min = 0.1, 
                                 max = 1, 
                                 step = 0.05)
              ),
              ## table
              column(width = 5,
                  dataTableOutput(ns("correlationCouple"), width = "100%"),
                  textInput(ns("PubMedRequest_correlation"), 
                            "Selected features query:", 
                            value = "", 
                            width = "100%", 
                            placeholder = NULL)
              ),
          ),
          box(title = "Partial Correlation", solidHeader = TRUE, width = 12, 
            collapsible = TRUE, collapsed = TRUE, status = "success",
              column(width = 7,
                  forceNetworkOutput(ns("partialCorNetwork")),
                  sliderInput(ns("partialCorNumThreshold"),
                               label = "Partial correlation absolute value 
                              threshold",
                               value = 0.2, 
                               min = 0.1, 
                               max = 1, 
                               step = 0.05)
              ),
              ## table
              column(width = 5,
                     dataTableOutput(ns("partialCorCouple"), width = "100%"),
                     textInput(ns("PubMedRequest_partialCor"), 
                               "Selected features query:", 
                               value = "", 
                               width = "100%", 
                               placeholder = NULL)
              ),
          ),
          box(title = "Mutual Information", solidHeader = TRUE, width = 12, 
            collapsible = TRUE, collapsed = TRUE, status = "success",
              column(width = 7,
                  forceNetworkOutput(ns("mutualInfNetwork")),
                  sliderInput(ns("mutualInfNumThreshold"),
                               label = "Mutual information absolute value 
                              threshold",
                               value = 0.2, 
                               min = 0.1, 
                               max = 2, 
                               step = 0.1)
              ),
              ## table
              column(width = 5,
                     dataTableOutput(ns("mutualInfCouple"), width = "100%"),
                     textInput(ns("PubMedRequest_mutualInf"), 
                               "Selected features query:", 
                               value = "", 
                               width = "100%", 
                               placeholder = NULL)
              ),
          ),
          column(10, 
              textInput(ns("BuildPubMedRequest"), 
                        "PubMed query", 
                        value = "", 
                        width = "100%", 
                        placeholder = NULL),
              actionButton(ns("btn_allRequest"),
                           "Update query with selected relations"),
              actionButton(ns("btn_search"),
                           "Search in PubMed")
          ),
          htmlOutput(ns("pubmed"))
    )
}
    
#' hypothesisGenerator Server Function
#' @param input,output,session Internal parameters for {shiny}.
#' @param obj R6 object to wrap all data from different analysis.
#' @param biosignature Discriminant biological features selected. 
#' @param featMethod character string to specify feature selection method 
#' (e.g. "diablo");
#'
#' @importFrom DT renderDataTable
#' 
#' @return Network inferences and PubMed queries results. 
mod_hypothesisGenerator_server <- function(input, 
                                           output, 
                                           session, 
                                           obj,
                                           biosignature, 
                                           featMethod){
  
    ##################################################################
    ##                             Data                             ##
    ##################################################################
    ## data without classes Y vector
    seqOmicData <- seq(1, length(obj$data$dataTrain)-1)
    omicData <- obj$data$wholeData[seqOmicData]
    
    message("Network building...")
    
    concatenatedOmicData <- getDataSelectedFeatures(omicData, biosignature)
    concatMatrixName <- paste0("concatMatrixFeats_", featMethod)
    obj$data[[concatMatrixName]] <- concatenatedOmicData
    
    
    #################################################################
    ##                         Correlation                         ##
    #################################################################
    message("Correlation...")
    ## Correlation reactive outputs  
    corrTable <- buildNetworkInference(input, 
                                      output, 
                                      session, 
                                      "correlation", 
                                      featMethod,
                                      concatenatedOmicData)
    output$correlationCouple <- renderDataTable(corrTable())
    
    
    #################################################################
    ##                     Partial Correlation                     ##
    #################################################################
    message("Partial correlation...")
    ## Partial correlation reactive outputs  
    pcorTable <- buildNetworkInference(input, 
                                       output, 
                                       session, 
                                       "partialCor", 
                                       featMethod,
                                       concatenatedOmicData)                                                    
    output$partialCorCouple <- renderDataTable(pcorTable())
    
    
    
    ##################################################################
    ##                      Mutual Information                      ##
    ##################################################################
    message("Mutual information...")
    ## Mutual Information reactive outputs  
    miTable <-  buildNetworkInference(input, 
                                      output, 
                                      session, 
                                      "mutualInf", 
                                      featMethod,
                                      concatenatedOmicData)
    output$mutualInfCouple <- renderDataTable(miTable())
    
    
    
    ##################################################################
    ##                         PUBMED query                         ##
    ##################################################################
    ## Correlation query
    buildPubmedQuery(input, output, session, "correlation", corrTable())
    
    ## Partial correlation query
    buildPubmedQuery(input, output, session, "partialCor", pcorTable())
    
    ## Mutual Information query
    buildPubmedQuery(input, output, session, "mutualInf", miTable())
    
    
    ## Build final Pubmed request 
    mergePubmedQuery(input, output, session)
    
    ## Launch pubmed request
    observeEvent(input$btn_search, {
        my_query <- input$BuildPubMedRequest
        if (my_query != " ")
        {
            ## Html output of 5 pubmed titles and abstracts according to query
            message("Pubmed query...")
          print(featMethod)
            pubmedInsight(session, input, output, my_query, featMethod)
        }
    })
}