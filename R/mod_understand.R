#' understand UI Function
#'
#' @description A shiny Module.
#'
#' @param id Internal parameter for {shiny}.
#'
#' @importFrom shiny NS tagList h3 span icon checkboxGroupInput numericInput
#' 
#' @return Displays UI output for Biological Insigths tab.
mod_understand_ui <- function(id){
  ns <- NS(id)
  tagList(
          box(width = 10, solidHeader = TRUE, 
              status = "warning",
              collapsible = TRUE,
              title = "Databases & parameters",
              box(width = 6, align = "center",
                  h3("Annotation databases\n"),
                  box(collapsible = TRUE,
                      title = span(icon("database")), 
                      status = "primary",
                      width = NULL,
                      
                      ## Pathways databases choice
                      checkboxGroupInput(ns("pathways_enrichDatabase"), 
                          span(tagList(icon("book"), "Pathways")), 
                          choices = c("Kegg" = "kegg", 
                              "Reactome" = "reactome", 
                              "wikiPathways" = "wikiPathways"), 
                          selected = c("reactome")
                      ),
                      ## GO ontology databases choice
                      checkboxGroupInput(ns("GO_enrichDatabase"), 
                          span(tagList(icon("book"), "Gene Ontology (GO)")), 
                          choices = c("Molecular functions (GO)" = "MF", 
                              "Cellular components (GO)" = "CC",
                              "Biological processes (GO)" = "BP"), 
                      )
                  )
              ),
          box(width = 6, align = "center",
                     h3("Enrichment parameters\n"),
                     box(collapsible = TRUE,
                         title = span(icon("gear")), 
                         status = "primary", 
                         width = NULL,
                         
                         column(width = 6,
                                numericInput(ns("minGSSize"), 
                                    "Minimal number of genes in pathways:", 
                                    value = 5, 
                                    min = 1, 
                                    max = 5000, 
                                    step = 1)
                         ),
                         column(width = 6,
                                numericInput(ns("maxGSSize"), 
                                    "Maximal number of genes in pathways:", 
                                    value = 800, 
                                    min = 1, 
                                    max = 5000, 
                                    step = 1)
                         ),
                         column(width = 12,
                                selectInput(ns("pvAdjust"), 
                                    "P-value correction method:", 
                                    choices = c("Benjamini-Hochberg" = "BH",
                                        "Holm" = "holm",
                                        "Hochberg" = "hochberg", 
                                        "Hommel" = "hommel", 
                                        "Bonferroni" = "bonferroni", 
                                        "BY" = "BY",
                                        "FDR" = "fdr", 
                                        "none"),
                                    selected = "BH")
                         ),
                         col_12(
                           sliderInput(ns("plotStoufferThreshold"),
                                       label = "Stouffer-s value threshold for 
                                       Enrichment maps",
                                       value = 0.1,
                                       min = 0.001,
                                       max = 1,
                                       step = 0.05)
                         )
                     ),
              ),
              uiOutput(ns("sourceID")),
          ),
          
      
      # tabBox(width = 10,
        # tabPanel(title = "Deseq2", width = 12,
      box(width = 10, solidHeader = TRUE, 
          status = "primary",
          collapsible = TRUE,
          title = "Deseq2 analysis & enrichment",
            uiOutput(ns("warnDeseq")),
            uiOutput(ns("errorDeseq")),
            sliderInput(ns("DEtable_padj"),
                        label = "p-adj value threshold for Differential 
                        Expression table",
                        value = 0.05,
                        min = 0.001,
                        max = 1,
                        step = 0.05),
            tabBox(width = 12,
                tabPanel("DESeq2", 
                         uiOutput(ns("uiDeseq"))),
                tabPanel("pathways", 
                         uiOutput(ns("deseq_pathways_enrichResults"))),
                tabPanel("go",
                         uiOutput(ns("deseq_GO_enrichResults")))
            )
          ),
        # tabPanel(title = "Diablo", width = 12,
      box(width = 10, solidHeader = TRUE, 
          status = "danger",
          collapsible = TRUE,
          title = "Diablo - sPLS-DA features enrichment",
            tabBox(width = 12,
                tabPanel("pathwaysTab", 
                         uiOutput(ns("diablo_pathways_enrichResults"))),
                tabPanel("goTab", 
                         uiOutput(ns("diablo_GO_enrichResults")))
            )
          )
      )
    # )
}

#' understand Server Function
#'
#' @param input,output,session Internal parameters for {shiny}.
# #' @param obj R6 object to wrap all data from different analysis.
#' @param startSignal input$start from Start button in Home tab.
#' 
#' @importFrom shiny req renderUI span icon
#' @importFrom anyLib anyLib
#' 
#' @return Biological insights from databases chosen by user by hypergeomtric 
#' tests (ORA) on DESeq2 or diablo (sPLS-DA) selected features.
mod_understand_server <- function(input, output, session, startSignal){
    ns <- session$ns
    
    ## When application is launched after data inputs
    observeEvent(startSignal$start,{
        req(length(obj$data$wholeData) > 0)
      
        ## Org.db download according user home input 
        orgDbinput <- obj$organismDb <- startSignal$organismDb
        anyLib::anyLib(orgDbinput, autoUpdate = TRUE) 
        orgDb <- eval(parse(text = orgDbinput, keep.source=FALSE))
        
        ## database sources available for organism chosen
        featureIDSources <- keytypes(orgDb)
        featureIDSources <- c("Do not enrich", featureIDSources)
        
        ## create selection input interface to get database sources for omics 
        omicDataNumber <- length(obj$data$wholeData)-1
        omicDataNames <- names(obj$data$wholeData)
        output$sourceID <- renderUI({
            box(width = 12, 
                lapply(
                    seq(1, omicDataNumber),  function(i){
                      omicNameID <- paste0(omicDataNames[i], "_db")
                      column(6,
                          selectInput(inputId = ns(omicNameID), 
                                      label = span(icon("database"), 
                                                   omicNameID), 
                                      choices = featureIDSources, 
                                      selected = featureIDSources[1])
                      )
                }), #end of lapply for each omic data block
                column(12, align = "center", 
                       actionButton(ns("runEnrichDeseq"), 
                                    "Run DESEQ2 enrichment"),
                       actionButton(ns("runEnrichDiablo"), 
                                    "Run Diablo features enrichment"),
                       uiOutput(ns("uiNoConv")),
                       uiOutput(ns("errorBase")),
                       uiOutput(ns("errorOrg")),
                )
            )
        }) # end of renderUI
    }) 
  
    #################################################################
    ##                       Diablo features                       ##
    #################################################################
    observeEvent(input$runEnrichDiablo,{
      organismDb <- obj$organismDb
      pathwaysDatabasesChosen <- input$pathways_enrichDatabase
      goDatabasesChosen <- input$GO_enrichDatabase
      output$uiNoConv <- NULL
      
      
      anyDB <- length(c(pathwaysDatabasesChosen, goDatabasesChosen)) > 0
      req(anyDB) ## enrichment not launched if no database chosen
      ## Features IDs converting to ENTREZID for enrichment
      omicSignature <- obj$classification$diabloResult$biosignature
      
      ## Chosen databases ID by user to convert feature ids
      omicDataNames <- names(omicSignature)
      dBinputsVec <- getDbFromInput(input, session, omicDataNames)
      obj$enrichment$diablo$convertDb <- dBinputsVec
      
      req(length(dBinputsVec)>0)
      
      message("Converting features by bitr()...")
      convertedFeatures <- 
          convertToEntrezid(omicSignature, dBinputsVec, organismDb)
      obj$enrichment$diablo$featureConverted <- convertedFeatures
      checkConvFeat(output, session, dBinputsVec, convertedFeatures)
      
      ## Diablo enrichment
      runMultiOmicEnrichment(input, output, session, "diablo") ## _utils.R
    })
  
  
    ##################################################################
    ##                            DESEQ2                            ##
    ##################################################################
    observeEvent(input$runEnrichDeseq,{
        organismDb <- obj$organismDb
        pathwaysDatabasesChosen <- input$pathways_enrichDatabase
        goDatabasesChosen <- input$GO_enrichDatabase
        output$uiNoConv <- NULL
        output$uiDeseq <- NULL
        
        
        anyDB <- length(c(pathwaysDatabasesChosen, goDatabasesChosen)) > 0
        req(anyDB) ## enrichment not launched if no database chosen
        
        ## DESEQ2 differential expression analysis.
        omicData <- obj$data$wholeData # with Y
        deseqRes <- builDeseqAnalysis(omicData, input, session, output)
        
        ## displays DESeq tables in ui
        deseqTables <- deseqRes$DEtable
        output$uiDeseq <- renderUI({ displayDESeqtables(deseqTables) })
        ## Features IDs converting to ENTREZID for enrichment
        omicSignature <- deseqRes$selectedFeatures
        req(omicSignature)
        
        ## Chosen databases ID by user to convert feature ids
        omicDataNames <- names(omicSignature)
        dBinputsVec <- getDbFromInput(input, session, omicDataNames)
        
        req(length(dBinputsVec)>0)
        
        convertedFeatures <- 
            convertToEntrezid(omicSignature, dBinputsVec, organismDb)
        ## Displays error if features mapping returns null feature vectors for 
        # at least one omic block.
        message("Check selected features converted...")
        checkConvFeat(output, session, dBinputsVec, convertedFeatures)
        
        ## DESEQ enrichment
        message("Run DESeq2's features multi-omic enrichment...")
        runMultiOmicEnrichment(input, output, session, "deseq") ## _utils.R
        
        obj$enrichment$deseq$featureConverted <- convertedFeatures
        obj$enrichment$deseq$DEtable <- deseqRes$DEtable
        View(deseqRes)
        
    })
    
    return(input)
}

    
## To be copied in the UI
# mod_understand_ui("understand_ui_1")
    
## To be copied in the server
# callModule(mod_understand_server, "understand_ui_1")
 
