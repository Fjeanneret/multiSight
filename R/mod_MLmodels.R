#' MLmodels UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList column sliderInput uiOutput
#' textInput actionButton htmlOutput verbatimTextOutput
#' @importFrom shiny tabPanel h4 
#' @importFrom shinydashboard box
#' @importFrom DT dataTableOutput renderDataTable
mod_MLmodels_ui <- function(id){
    ns <- NS(id)
    tagList(
        box(title = "Biosigner - SVM & RF models", 
            width = 10, solidHeader = TRUE, 
            collapsible = TRUE, collapsed = TRUE, 
            status = "warning",
                tabBox(width = 12,
                    tabPanel("Features", 
                        uiOutput(ns("WarnNotTwoClass")),
                        uiOutput(ns("biosignerFeatDetails"))
                    ),
                    tabPanel("Performances",
                        uiOutput(ns("WarnNotTwoClassP")),
                        h4("SVM"),
                        verbatimTextOutput(ns("biosignerPerfsvm")),
                        h4("RF"),
                        verbatimTextOutput(ns("biosignerPerfrf"))
                    )
                ),
        ),
        box(title = "Diablo - sPLS-DA model", 
            width = 10, solidHeader = TRUE,
            collapsible = TRUE, collapsed = TRUE, 
            status = "warning",
                tabBox(width = 12,
                    tabPanel("Features", 
                        uiOutput(ns("diabloFeatDetails"))
                    ),
                    tabPanel("Performances",
                        verbatimTextOutput(ns("diabloPerfsplda"))
                    )
                )
        )
    )
}
    
#' MLmodels Server Function
#'
#' @importFrom shiny observeEvent req renderUI renderPrint span
#'
#' @noRd 
mod_MLmodels_server <- function(input, output, session, startSignal, bioDB)
{
    ns <- session$ns
    
    observeEvent(startSignal$start,{
      req(length(obj$data$wholeData) > 0)
      
      ## Split data sets in train and test sets
      dataSplitted <- splitDatatoTrainTest(obj$data$wholeData, freq = 0.8)
      
      obj$data$dataTrain <- dataTrain <- dataSplitted$data.train
      obj$data$dataTest <- dataTest <- dataSplitted$data.test
      class <- obj$data$wholeData$Y
      nclass <- length(levels(factor(class)))
      
      #################################################################
      ##                       Biosigner model                       ##
      #################################################################
      if (nclass == 2)
      {
          message("Biosigner model building...")
          biosignerModelRes <- runSVMRFmodels_Biosigner(dataTrain)
          obj$classification$biosignerResult <- biosignerModelRes
          # runSVMRFmodels_Biosigner(obj)
          
          message("Biosigner's SVM-RF models assessing...")
          modelList <- biosignerModelRes$model
          # assessPerformance_Biosigner(obj)
          biosignerPerf <- assessPerformance_Biosigner(modelList, 
                                                            dataTest)
          message("...Biosigner's models assessing ended")
          obj$classification$biosignerResult$performance <- biosignerPerf
          
          ## biosigner results
          biosignerSig <- biosignerModelRes$biosignature
          biosignerFeatTables <- 
            displayFeatDetails(featuresList = biosignerSig, 
                                modelMethod = "biosigner", 
                                obj = obj)
          ## Shiny outputs
          output$biosignerFeatDetails <- renderUI({biosignerFeatTables})
          output$biosignerPerfsvm <- renderPrint(biosignerPerf$svm)
          output$biosignerPerfrf <- renderPrint(biosignerPerf$rf)
      }
      else
      {
          msg <- "More than 2 classes, not supported by Biosigner yet. 
          Biosigner will not compute."
          message(msg)
          output$WarnNotTwoClass <- WarnNotTwoClassP <- renderUI({
              span(msg, style="color:blue")
          })
      }
  
      
      ##################################################################
      ##                     Diablo sPLS-DA model                     ##
      ##################################################################
      message("Diablo model building...")
      diabloRes <- runSPLSDA(dataTrain)
      obj$classification$diabloResult <- diabloRes
      sPLSDAmodel <- diabloRes$model
      diabloSig <- diabloRes$biosignature
      
      message("Diablo's sPLS-DA model assessing...")
      diabloPerf <- assessPerformance_Diablo(sPLSDAmodel, 
                                             dataTest)
      message("...Diablo's model assessing ended")
      
      ## diablo results
      obj$classification$diabloResult$performance <- diabloPerf
      obj$classification$diabloResult$biosignature <- diabloSig
      diabloFeatTables <- displayFeatDetails(featuresList = diabloSig, 
                                             modelMethod = "diablo", 
                                             obj = obj)
      
      ## Shiny outputs
      output$diabloPerfsplda<- renderPrint(diabloPerf)
      output$diabloFeatDetails <- renderUI({diabloFeatTables})
      })
    
    
      ##################################################################
      ##         Reload tables with DESeq2 values if computed         ##
      ##################################################################
      ## Updates feature tables in classification tab with logF and padj values
      ### Biosigner
      observeEvent(bioDB$runEnrichDeseq,{
          class <- obj$data$wholeData$Y
          nclass <- length(levels(factor(class)))
          if (nclass == 2)
          {
              message("[Classification tab] Biosigner ",
                      "features tables updating...")
              biosignerSig <- obj$classification$biosignerResult$biosignature
              DESEQtables <- obj$enrichment$deseq$DEtable
              biosignerFeatTables <- 
                displayFeatDetails(biosignerSig,
                                   modelMethod ="biosigner",
                                   DESeqTables = DESEQtables, 
                                   obj = obj)       
              output$biosignerFeatDetails <- renderUI({biosignerFeatTables})
          }
      })
      
      ### Diablo
      observeEvent(bioDB$runEnrichDeseq,{
          message("[Classification tab] Diablo features tables updating...")
          diabloSig <- obj$classification$diabloResult$biosignature
          DESEQtables <- obj$enrichment$deseq$DEtable
          diabloFeatTables <- displayFeatDetails(diabloSig, 
                                                 "diablo",
                                                 DESEQtables, 
                                                 obj = obj)
          output$diabloFeatDetails <- renderUI({diabloFeatTables})
      })
}

## To be copied in the UI
# mod_MLmodels_ui("MLmodels_ui_1")
    
## To be copied in the server
# callModule(mod_MLmodels_server, "MLmodels_ui_1")
 
