#' Understand utils function
#'
#' @description To replace geneID column by whole geneset 
#' not only selected features mapped to pathways 
#' to build a relevant enrichment map
#'
#' @param enrichResult enrichResult class object from enrichment analysis
#' 
#' @noRd
#' 
#' @return enrichResult
fillEnrichGeneID <- function(enrichResult)
{
    ## enrichResult@geneSets contains whole geneset for x pathways
    ## not only selected features mapped to pathways
    geneSetVector <- vapply(enrichResult@result$ID, function(x) {
        geneset <- enrichResult@geneSets[[x]] 
        g <- paste(geneset, collapse = "/")
        return(g)},
        FUN.VALUE = character(1))
    
    enrichResult@result$geneSet <- geneSetVector
    return (enrichResult)
}

#' Understand utils function
#'
#' @description To build termsim matrix to plot enrichment Map. JC's method.
#' From emaplot pairwise_termsim function.
#' 
#' @param enrichmentTable Table with stouffer values
#' @param geneSets geneSets with pathways genes list
#'
#' @noRd
#' 
#' @return termsim 
termsimJC <- function(enrichmentTable, geneSets)
{
    overlap_ratio <- function(x, y) {
        x <- unlist(x)
        y <- unlist(y)
        length(intersect(x, y))/length(unique(c(x,y)))
    }
    id <- enrichmentTable[, "ID"]
    geneSets <- geneSets[id]
    n <- nrow(enrichmentTable)
    enrichTermSim <- matrix(NA, nrow=n, ncol=n)
    colnames(enrichTermSim) <- rownames(enrichTermSim) <- 
        enrichmentTable$Description
    for (i in seq_len(n-1)) {
        for (j in (i+1):n) {
            enrichTermSim[i, j] <- overlap_ratio(geneSets[id[i]], 
                                                 geneSets[id[j]])
        }
    }
    return(enrichTermSim)
}

#' Understand utils function
#'
#' @description To build Enrichment Map from
#' multi-omic table with Stouffer p-values.
#'
#' @param enrichmentStouffer Table with stouffer values
#' @param pvStouffer Stouffer's value informed by user input.
#'
#' @importFrom enrichplot emapplot
#' @importFrom ggsci scale_color_gsea 
#' @importFrom ggplot2 ggtitle ggplot annotate theme_void
#' 
#' @noRd
#' 
#' @return emplot 
emFromStouffer <- function(enrichmentStouffer, pvStouffer)
{
    plotTitle <- paste0("Enrichment Map - Stouffer's Weighted p-value <= ", 
                      pvStouffer)
    enrichmentTable <- enrichmentStouffer@result
    geneSets <- enrichmentStouffer@geneSets
    enrichmentTable <- enrichmentTable %>%
      filter(.data$StoufferWeighted <= pvStouffer)
    enrichmentStouffer@result <- enrichmentTable
    if (nrow(enrichmentTable) > 1)
    {
        # enrichTermSim <- pairwise_termsim(enrichmentStouffer,
                                          # method = "JC")

        enrichmentStouffer@termsim <- termsimJC(enrichmentTable, 
                                                geneSets)
        enrichmentStouffer@method <- 'JC'
        emplot <- emapplot(enrichmentStouffer,
                         layout="kk", 
                         min_edge = 0.2,
                         cex_label_category = 0.75) + 
        scale_color_gsea(reverse = TRUE) +
        ggtitle(plotTitle)
  }
  else 
  {
    text = paste0("\n   No one or only one Stouffer weighted's 
                  value\n less or equal to your threshold: ", 
                  pvStouffer)
    emplot <- ggplot() + 
      annotate("text", x = 4, y = 25, size=5, label = text) + 
      theme_void()
    
  }
  return(emplot)
}

#' Understand utils function
#'
#' @description To get user's database choices from Biological Insights tab.
#'
#' @param input,session Internal parameters for shiny
#' @param omicDataNames Omic data set names to build input id to retrieve
#' databases chosen in Biological Insights tab.
#'
#' @noRd
#' 
#' @return Database names vector. One by omic block (e.g. for 2 omics 
#' c("SYMBOL", "ENSEMBL")) 
getDbFromInput <- function(input, session, omicDataNames)
{
    inputIDs <- paste0(omicDataNames, "_db")
    i <- 1
    dBinputsVec <- list()
    while (i <= length(omicDataNames))
    {
        selectInputID <- inputIDs[[i]]
        omicName <- omicDataNames[i]
        dbFrom <- input[[selectInputID]] # ui input
        if(dbFrom != "Do not enrich")
        {
            dBinputsVec[[omicName]] <- dbFrom 
        }else
        {
            dBinputsVec[[omicName]] <- NULL
        }
        i = i+1
    }
    return(dBinputsVec)
}

#' Understand utils function
#'
#' @description To convert features to entrezid to enrich.
#'
#' @param featList Feature lists from each omic data block.
#' @param fromDbList Database names vector. One by omic block (e.g. for 2 omics 
#' list(omic1 = "SYMBOL", omic2 = "ENSEMBL")). Returned by getDbFromInput() 
#' for app. 
#' @param organismDb Organism database to convert features.
#'
#' @examples
#' \donttest{
#' library(org.Mm.eg.db) # Organism's database
#' splittedData <- splitDatatoTrainTest(omic2, 0.8)
#' data.train <- splittedData$data.train
#' data.test <- splittedData$data.test
#' 
#' diabloRes <- runSPLSDA(data.train)
#' diabloModels <- diabloRes$model #sPLS-DA model using all omics.
#' diabloFeats <- diabloRes$biosignature #selected features for each omic.
#' id_db <- list(omic1 = "ENSEMBL", omic2 = "ENSEMBL")
#' convFeat <- convertToEntrezid(diabloFeats, id_db, "org.Mm.eg.db")
#' }
#' 
#' featList <- list(Omic1 = c("ENSMUSG00000039621", 
#'                            "ENSMUSG00000038733", 
#'                            "ENSMUSG00000062031"), 
#'                  Omic2 = c("ENSMUSG00000031170", 
#'                            "ENSMUSG00000077495", 
#'                            "ENSMUSG00000042992"))
#' dbList <- list(Omic1 = "ENSEMBL",
#'                  Omic2 = "ENSEMBL")
#' 
#' convFeat <- convertToEntrezid(featList, dbList, "org.Mm.eg.db")
#' 
#' @importFrom clusterProfiler bitr
#' 
#' @export
#' 
#' @return featConverted 
convertToEntrezid <- function(featList, fromDbList, organismDb)
{
  omicDataNumber <- length(featList)
  omicDataNames <- names(featList)
  
  ## orgDb to convert
  orgDb <- eval(parse(text = organismDb, keep.source=FALSE))
  
  ## For each omic data features
  featConvList <- list()
  i <- 1
  while (i <= omicDataNumber)
  {
    omicName <- omicDataNames[i]
    dbFrom <- fromDbList[[omicName]] # ui input
    if(!is.null(dbFrom))
    {
      features <- featList[[i]]
      out <- tryCatch(
        {
            if (length(features) == 0)
            {
               stop(simpleError("Empty feature list"))
            }
            
            ## bitr converts features names to entrez id
            featureConverted <- bitr(features,
                                     fromType = dbFrom,
                                     toType = "ENTREZID",
                                     OrgDb = orgDb)
            if (length(featureConverted$ENTREZID) == 0)
            {
                stop(simpleError("Empty converted feature list"))
            }
            
            ## features converted toward enrichment
            featConvList[[omicDataNames[i]]] <- featureConverted$ENTREZID
        },
        error=function(e) 
        {
            msg <- e$message
            msg1 <- ": None of the features entered are valid keys for '"
            omic <- omicDataNames[i]
            dbFrom <- fromDbList[[omicName]]
            message(paste0(omic, 
                           msg1, 
                           dbFrom,
                           "'"))
            featConvList[[omicDataNames[i]]] <- NULL
          })
    }
    else 
    {
        featConvList[[omicDataNames[i]]] <- NULL
    }
    i <- i + 1
  } ## End of while, for each omic data
  
  return(featConvList)
}

#' Understand utils function
#'
#' @description To manage several enrichment database results.
#' 
#' @param output,session Internal parameters for shiny
#' @param convertedFeatList List of converted features lists returned by
#' convertToEntrezid() using bitr().
#' @param dBinputsVec Database names vector. One by omic block (e.g. for 2 
#' omics c("SYMBOL", "ENSEMBL")). Returned by getDbFromInput(). 
#' 
#' @noRd
#' 
#' @importFrom shiny renderUI span
#' 
#' @return UI output to display errors.
checkConvFeat <- function(output, session, dBinputsVec, convertedFeatList)
{
    namesOmics <- names(convertedFeatList)
    nullFeatListOmic <- vapply(convertedFeatList, is.null, c(TRUE))
    nullFeatListOmicNames <- namesOmics[nullFeatListOmic]
  
    
    if (is.null(namesOmics))
    {
        dbNull <- dBinputsVec[nullFeatListOmic]
        omics <- paste(dbNull, collapse = "', '")
        output$uiNoConv <- renderUI({
            span(
              "At least one omic's features vector entered are not valid keys 
                for ", 
              dbNull,
              "' or: empty selected features list.",
              style="color:green")
    })
    }else if (length(nullFeatListOmicNames) != 0)
    {
        dbNull <- dBinputsVec[nullFeatListOmic]
        omics <- paste(nullFeatListOmicNames, sep = ", ")
        output$uiNoConv <- renderUI({
            span(
                "At least one omic's features vector entered are not valid keys 
                for ", 
                dbNull,
                "' or: empty selected features list.",
                style="color:green")
    })
    }else
    {
        # All omic features have been converted.
    }
}


#' #' Understand utils function
#'
#' @description To manage several enrichment database results.
#' 
#' @param input,output,session Internal parameters for shiny
#' @param featSource Deseq2 or diablo features.
#' 
#' @noRd
#' 
#' @return Runs and displays pathways and Gene Ontology enrichment
#'  for each omic data block.
runMultiOmicEnrichment<- function(input, output, session, featSource)
{
    organismDb <- obj$organismDb
    minGSSize <- input$mixGSSize
    maxGSSize <- input$maxGSSize
    pvAdjust <- input$pvAdjust
    pvStouffer <- input$plotStoufferThreshold
    pathwaysDatabasesChosen <- input$pathways_enrichDatabase
    goDatabasesChosen <- input$GO_enrichDatabase
    databasesChosen <- c(pathwaysDatabasesChosen, 
                         goDatabasesChosen)
    omicFeats <- obj$enrichment[[featSource]]$featureConverted
    
    ## Multi-omic features enrichment
    # Remove previous results 
    obj$enrichment[[featSource]]$pathways <- list()
    obj$enrichment[[featSource]]$go <- list()
    if (!is.null(databasesChosen) && 
        length(omicFeats) > 0)
    {
        # Removes previous errors
        output$errorBase <- NULL
        output$errorOrg <- NULL
        
        # Enrichment analysis
        message("Multi-Omic enrichment analysis... ")
        multiEnrichRes <- runMultiEnrichment(omicSignature = omicFeats,
                                             databasesChosen = databasesChosen,
                                             organismDb = organismDb,
                                             pvAdjust = pvAdjust,
                                             minGSSize = minGSSize,
                                             maxGSSize = maxGSSize,
                                             pvStouffer = pvStouffer)
        
        ## Results
        obj$enrichment[[featSource]] <- multiEnrichRes
        obj$enrichment[[featSource]]$featureConverted <- omicFeats
        output$errorBase <- multiEnrichRes$error$errorBase
        output$errorOrg <- multiEnrichRes$error$errorOrgBase
    }
    
    dataNames <- names(omicFeats)
    
    ## export result
    outputNamePath <- paste0(featSource, "_pathways_enrichResults")
    outputNameGO <- paste0(featSource, "_GO_enrichResults")
    output[[outputNamePath]] <- NULL
    output[[outputNameGO]] <- NULL
    output[[outputNamePath]] <- uiResultExport(input,
                                               output,
                                               session, 
                                               enrichmentType = "pathways",
                                               featSource,
                                               obj$enrichment,
                                               dataNames)
    output[[outputNameGO]] <- uiResultExport(input,
                                             output, 
                                             session,
                                             enrichmentType = "go",
                                             featSource,
                                             obj$enrichment,
                                             dataNames)
}

#' Enrichment utils ui function
#'
#' @description Displays DESEQ2 analysis tables for all omics
#'
#' @param deseqTables List of DESeq() tables for all omics.
#' 
#' @noRd 
#' 
#' @importFrom DT datatable 
#' 
#' @return ui elements to display by renderUI().
displayDESeqtables <- function(deseqTables)
{
    omicNames <- names(deseqTables)
    ## box header color
    statusLabels <- c("primary", "success", "info", "warning", "danger")
    
    ## For each omic signature builds one table
    seqDeseqTables <- seq(1, length(deseqTables))
    lapply(seqDeseqTables, function(i){
        status <- sample(statusLabels, 1)
        table <- data.frame(deseqTables[[i]])
    ## Table 
    deseqTable <- datatable(table,
                            extensions = c('Responsive', 'Buttons'),
                            options = list(
                                dom = 'BRrltpi',
                                rownames = FALSE,
                                buttons = c('csv', 'excel', 'pdf')
                              )
      )
      colToFormat <- colnames(table)
      deseqTable <- deseqTable %>% formatSignif(columns = colToFormat, 
                                                digits = 5)
      ## adds color background for p.adjust values <= 0.01 and 0.05
      deseqTable <- deseqTable %>% formatStyle(
          'padj',
          color = styleInterval(c(0.01, 0.05), 
                                c('white', 'white', 'black')),
          backgroundColor = styleInterval(c(0.01, 0.05), 
                                          c('#a2465f', '#cb5658', "white")))
      
      ## UI output
      box(width = 10, 
          title = omicNames[i], 
          status = status,
          solidHeader = TRUE,
          
          DT::renderDataTable({ deseqTable }, server = FALSE)       
      )
    })
}