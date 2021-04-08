#' UI export function for understand module
#' 
#' @description UI export in several boxes of enrichment results 
#' (e.g. for 2 omics, 4 boxes: Omic1, omic2, multi-Omic table and 
#' enrichment map).
#' 
#' @param input,output,session Internal parameters for shiny
#' @param enrichmentType Pathways or Gene Ontology enrichment.
#' @param featSource Model which selected features. 
#' @param multiEnrichRes Multi enrichment results obtained by 
#'   runMultiEnrichment(). 
#' @param dataNames Omic names. 
#'
#' @return One table with relation values as rows
#'  of all omic selected features.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList renderUI renderPlot
#' @importFrom DT dataTableOutput renderDataTable JS formatRound datatable
uiResultExport <- function(input, output, session, 
                           enrichmentType, 
                           featSource,
                           multiEnrichRes,
                           dataNames) 
{
  if (enrichmentType == "pathways")
  {
    result <- multiEnrichRes[[featSource]]$pathways
    plotMethod <- "enrichMap"
  }
  else if (enrichmentType == "go")
  {
    result <- multiEnrichRes[[featSource]]$go
    plotMethod <- "gosim"
  }
  databases <- names(result)
  if (!is.null(databases))
  {
    renderUI({
        lapply(seq(1, length(databases)), function(i){
            column(width = 12,
                align = "center", offset = 1,
                ## Boxes building
                box(width = 10, title = databases[i], 
                ## check if data
                if (!is.null(names(result[[databases[i]]]$result))) 
                {
                    ## check if data after db query
                    if (is.null(result[[databases[i]]]$result)) 
                    {
                        box(width = 12, title = "Organism not supported",
                            solidHeader = TRUE,
                            collapsible = FALSE, 
                            collapsed = TRUE)
                    }
                    else{
                        ## for each simple enrichment result table 
                        lapply(seq(1,
                            length(names(result[[databases[i]]]$result))), 
                            function(j){
                            if (j <= length(dataNames))
                            {
                                box(width = 12, title = dataNames[j], 
                                    solidHeader = TRUE,
                                    collapsible = TRUE,
                                    collapsed = TRUE,
                                    
                                      renderDataTable({
                                        formatEnrichTable(
                                          result[[databases[i]]]$result[[j]], 
                                          databases[i])
                                        }, server = FALSE)
                                )
                            }else if (length(dataNames)>1) ## Multi-omic results
                            {
                                if(j == length(dataNames)+1) ## Stouffer Table
                                {
                                  box(width = 12, title = "Multi-Omic Table", 
                                      solidHeader = TRUE,
                                      collapsible = TRUE,
                                      collapsed = TRUE,
                                      
                                      renderDataTable(
                                          formatStoufferTable(
                                          result[[databases[i]]]$result$multi,
                                          databases[i], 
                                          length(dataNames)), server = FALSE,
                                      )
                                  )
                                }else if (j == length(dataNames)+2) #EnrichMap
                                {
                                  if (plotMethod == "enrichMap") 
                                    {title <- "Enrichment Map"}
                                  else if (plotMethod == "gosim") 
                                    {title <- "GO Similarities"}
                                    box(width = 12, title = title, 
                                        solidHeader = TRUE,
                                        collapsible = TRUE, 
                                        collapsed = TRUE,
      
                                        renderPlot(
                                            result[[databases[i]]]$result[[j]]
                                        ),
                                    )
                                }
                            }
                        })
                    } # end of if organism supported by database
                } # end of if databases are empties (not chosen to enrich)
                else 
                {
                  box(width = 12, 
                      title = "You need at least 1 omic data enriched", 
                      solidHeader = TRUE,
                      collapsible = FALSE, 
                      collapsed = TRUE)
                }
            ))
          })
    })
  }
}

#' Utils understand function
#'
#' @description Launches functional enrichment of features provided
#' for every databases furnished. Run by utils function runMultiOmicEnrichment 
#'
#' @param dbResult Enrichment results for a given databases for all omics.
#' @param pvStouffer Numeric values chosen by user in ui to select 
#'   pathways under threshold to draw enrichment map.
#'
#' @noRd
#'  
#' @return Stouffer s multi-omic table and multi-omic enrichment map
buildStoufferTable <- function(dbResult, pvStouffer)
{
    stoufferRes <- list()
    enrichStouffer <- stoufferTable(dbResult)
    if (length(enrichStouffer) > 1) # Stouffer results exist
    {
        stoufferRes$multi <- enrichStouffer$table
        stoufferEnrichRes <- enrichStouffer$moEnrichRes
        message("Multi-Omic Enrichment map drawing...")
        stoufferRes$enrichMap <- 
            emFromStouffer(stoufferEnrichRes, pvStouffer)
    }else ## Stouffer table does not exist
    {
        stoufferRes$multi <- NULL
        stoufferRes$enrichMap <- NULL
    }
    
    return(stoufferRes)
}


#' Utils understand function
#'
#' @description Launches functional enrichment of features provided
#' for every databases furnished. Run by utils function 
#' runMultiOmicEnrichment(). 
#'
#' @param databaseName Database name (e.g. kegg, reactome, MF) 
#'
#' @noRd
#'  
#' @return Enrichment type (e.g pathways or go)
detectTypeEnrich <- function(databaseName)
{
    ## Each list to manage several type of functions (different arguments).
    pathDb <- c("kegg",
                "reactome",
                "wikipathways")
    
    geneOntDb <- c("MF",
                   "CC",
                   "BP")
    if (databaseName %in% pathDb)
    {
        typeEnrich <- "pathways"
    }
    else if (databaseName %in% geneOntDb)
    {
        typeEnrich <- "go"
    }
    
    return(typeEnrich)
}


#' Main understand module function
#'
#' @description Launches functional enrichment of features provided
#' for every databases furnished. Run by utils function runMultiOmicEnrichment 
#'
#' @param databasesChosen Which biological database 
#' c(reactome, kegg, wikiPathways, MF, CC, BP)
#' @param omicSignature Feature lists from each omic data block.
#' @param organismDb Organism provided by user in Home tab.
#' @param pvAdjust pv adjust method (e.g "BH" for Benjamini-Hochberg)
#' @param minGSSize,maxGSSize,pvStouffer Numeric values chosen by user 
#'   in ui.
#'
#' @examples
#' \donttest{
#' data("omic2", package = "multiSight")
#' splitData <- splitDatatoTrainTest(omic2, 0.8)
#' data.train <- splitData$data.train
#' data.test <- splitData$data.test
#' 
#' diabloRes <- runSPLSDA(data.train)
#' diabloModels <- diabloRes$model #sPLS-DA model using all omics.
#' diabloFeats <- diabloRes$biosignature #selected features for each omic.
#' id_db <- list(omic1 = "ENSEMBL", omic2 = "ENSEMBL")
#' }
#' library(org.Mm.eg.db, warn.conflicts = FALSE) #Organism's database
#' featList <- list(Omic1 = c("ENSMUSG00000039621", 
#'                            "ENSMUSG00000038733", 
#'                            "ENSMUSG00000062031"), 
#'                  Omic2 = c("ENSMUSG00000031170", 
#'                            "ENSMUSG00000077495", 
#'                            "ENSMUSG00000042992"))
#' dbList <- list(Omic1 = "ENSEMBL",
#'                  Omic2 = "ENSEMBL")
#' convFeat <- convertToEntrezid(featList, dbList, "org.Mm.eg.db")
#' 
#' ## To enrich features
#' database <- c("reactome", "MF")
#' #runMultiEnrichment_result <- runMultiEnrichment(databasesChosen = database,
#' #                                  omicSignature = convFeat,
#' #                                  organismDb = "org.Mm.eg.db")
#' data(runMultiEnrichment_result, package = "multiSight")
#'   
#' @importFrom dplyr filter
#' @importFrom clusterProfiler enrichKEGG read.gmt enricher enrichGO
#' @importFrom ReactomePA enrichPathway
#' 
#' @export
#' 
#' @return Wraps in obj enrichment results for all databases and all omics.
runMultiEnrichment <- function(omicSignature,
                               databasesChosen, 
                               organismDb,
                               pvAdjust = "BH",
                               minGSSize = 5,
                               maxGSSize = 800,
                               pvStouffer = 0.1)
{
    
    ##################################################################
    ##                To manage enrichment functions                ##
    ##################################################################
    db_fct <- list(kegg = enrichKEGG,
                   reactome = enrichPathway,
                   wikipathways = enrichWP,
                   MF = enrichGO,
                   CC = enrichGO,
                   BP = enrichGO)
    
    #################################################################
    ##                 Multi functional enrichment                 ##
    #################################################################
    multiEnrichRes <- list()
    errorBase <- c() # To send to UI bases not available.
    errorOrgBase <- c() # To send to UI organism not available in db.
    for (db in databasesChosen) # Multi databases
    {
        enrichFct <- db_fct[[db]] ## Selects enrichment fct (e.g. enrichPA)
        typeEnrich <- detectTypeEnrich(db)
        dbResult <- list()
        
        ## Get organism name according to database style
        orgForAllBioDb <- (organismTable %>% # organism table in sysdata.rda"
                             filter(.data$orgDb == organismDb))
        org <- orgForAllBioDb[[db]]
        
        ## Enrichment analysis
        if(!is.null(org) || typeEnrich == "go") # spe org only for pathways
        {
            message(db, " enrichment...")
            i <- 1
            ## for each feature vector from omic data
            while (i <= length(omicSignature)) 
            {
                if (!is.null(omicSignature[[i]]) && 
                    length(omicSignature[[i]])>0)
                {
                    out <- tryCatch(
                    {
                      ## Different arguments according enrichment type 
                      if (typeEnrich == "pathways")
                      {
                          bioDbRes <- enrichFct(
                              gene = omicSignature[[i]],
                              organism = org,
                              pAdjustMethod = pvAdjust,
                              minGSSize = minGSSize,
                              maxGSSize = maxGSSize,
                        )
                      }else if (typeEnrich == "go")
                      {
                          bioDbRes <- enrichFct(
                              gene = omicSignature[[i]],
                              OrgDb = organismDb,
                              ont = db,
                              pAdjustMethod = pvAdjust,
                              minGSSize = minGSSize, 
                              maxGSSize = maxGSSize,
                              readable = TRUE
                        )  
                      }
                      enrichRes <- fillEnrichGeneID(bioDbRes)
                      dbResult$enrichObj[[paste0("enrichRes_omic", i)]] <- 
                          enrichRes
                      dbResult$result[[paste0("enrichTable_omic", i)]] <- 
                          enrichRes@result
                    },
                    error=function(cond)
                    {
                        msg <- paste0(db, 
                                "seems down, or not available throught API.")
                        errorBase <- c(errorBase, msg)
                        multiEnrichRes$error$errorBase <- errorBase
                        dbResult$enrichObj[[paste0("enrichRes_omic", i)]] <- 
                            NULL
                        dbResult$result[[paste0("enrichTable_omic", i)]] <- 
                            NULL
                    })
                }
                else
                {
                    dbResult[[paste0("enrichRes_omic", i)]] <- NULL
                }
                    i <- i+1
                }
          
            #################################################################
            ##         BUILDS Multi-Omic TABLE - Stouffer p-values         ##
            #################################################################
            stoufferRes <- buildStoufferTable(dbResult$enrichObj, pvStouffer)
            dbResult$result$multi <- stoufferRes$multi
            dbResult$result$enrichMap <- stoufferRes$enrichMap
            
            ## All results in multiEnrichRes list
            multiEnrichRes[[typeEnrich]][[db]] <- dbResult
          
        }else # Organism not in organismTable for pathways enrichment
        {
            msg <- paste0("Organism not provided by ", db ," database")
            message(msg)
            errorOrgBase <- c(errorOrgBase, msg)
            multiEnrichRes[[typeEnrich]][[db]] <- NULL
            multiEnrichRes$error$errorOrgBase <- errorOrgBase
        }
    }
    return(multiEnrichRes)
}


#' Util function From clusterProfiler 
#' 
#' @description Downloads wikiPathways data and computes ORA analysis with
#' provided features.
#' 
#' @param gene Features to enrich.
#' @param organism Value chosen by user in Home tab..
#' @param pAdjustMethod,minGSSize,maxGSSize Numeric values chosen by user in
#'   Biological Insights tab.
#'
#' @importFrom clusterProfiler read.gmt enricher
#' @importFrom rWikiPathways downloadPathwayArchive
#' @importFrom dplyr select 
#' @importFrom tidyr separate 
#' 
#' @return enrichResult object with wikiPAthways database used.
enrichWP <- function(gene,
                     organism,
                     pAdjustMethod,
                     minGSSize,
                     maxGSSize)
{
    ## WPgmtFile according to organism 
    wpgmtfile <- downloadPathwayArchive(organism=organism, format = "gmt")
    wp2gene <- clusterProfiler::read.gmt(wpgmtfile)
    wp2gene <- wp2gene %>% 
        tidyr::separate(.data$term, c("name","version","wpid","org"), "%")
    wpid2gene <- wp2gene %>% dplyr::select(.data$wpid, gene) #TERM2GENE
    wpid2name <- wp2gene %>% dplyr::select(.data$wpid, .data$name) #TERM2NAME
    
    ewp <- enricher(gene, 
                    TERM2GENE = wpid2gene, 
                    TERM2NAME = wpid2name, 
                    pAdjustMethod = pAdjustMethod,
                    minGSSize = minGSSize,
                    maxGSSize = maxGSSize)
    
    return(ewp)
}

#' To compute omic weights for Stouffer value weighted
#' 
#' @description Unique genes involved in at least one pathways are
#' taken into account for each omic feature set.
#'
#' @param enrichmentResultTables enrichResults@result list.
#'
#' @return Numeric vector with n values for n omic enrichment
#' table results provided.
#'
#' @noRd
omicWeight <- function(enrichmentResultTables)
{
    weights <- c()
    i <- 1
    while(i <= length(enrichmentResultTables))
    {
        ## character vector splitting by "/" to get pathway's genes
        # results_geneID_list <- Reduce(c,
        #   enrichmentResultTables[[i]]$geneID)
        # AllGenes <- vapply(results_geneID_list,
        #   function(x)strsplit(x, "/"))
        ## transform into vector
        # genesVec <- Reduce(c, AllGenes)
        ## unique genes number
        # genesNumber <- length(unique(genesVec))
        genesNumber <- strsplit(enrichmentResultTables[[i]][1,3], "/")[[1]][2]
        weights[i] <- as.integer(genesNumber)
        
        i <- i + 1
    }
    ## each gene number value -> proportion for each omic
    return(weights/sum(weights))
}

#' Stouffer's p-value computing with all pvalues from all enrichment
#' tables provided and all gene sets to build enrichment maps.
#'
#' @param enrichmentResult enrichResults@result list.
#'
#' @examples 
#' data(enrichResList, package = "multiSight")
#' enrichResList # list of enrichRes objects (e.g. enrichKEGG() results)
#' multiOmicRes <- stoufferTable(enrichResList)
#' multiOmicRes$table # table with stouffer's values  
#' multiOmicRes$moEnrichRes # enrichRes object for clusterProfiler plots  
#' 
#' \donttest{
#' data("omic2", package = "multiSight")
#' splitData <- splitDatatoTrainTest(omic2, 0.8)
#' data.train <- splitData$data.train
#' data.test <- splitData$data.test
#' 
#' diabloRes <- runSPLSDA(data.train)
#' diabloModels <- diabloRes$model #sPLS-DA model using all omics.
#' diabloFeats <- diabloRes$biosignature #selected features for each omic.
#' id_db <- list(omic1 = "ENSEMBL", omic2 = "ENSEMBL")
#' library(org.Mm.eg.db) # Organism's database
#' convFeat <- convertToEntrezid(diabloFeats, id_db, "org.Mm.eg.db")
#' database <- c("reactome", "MF")
#' #enrichTables <- runMultiEnrichment(databasesChosen = database,
#' #  omicSignature = convFeat,
#' #  organismDb = "org.Mm.eg.db")
#' # enrichmentTables <- enrichTables$pathways$reactome$enrichObj
#' #enrichResList # list of enrichRes objects (e.g. enrichKEGG() results)
#' data(enrichResList, package = "multiSight")
#' multiOmicTable <- stoufferTable(enrichResList)
#' }
#' 
#' @return enrichment table results merged with Stouffer's p-value
#' non-weighted and weighted.
#'
#'
#' @importFrom dplyr select arrange rename matches mutate coalesce
#' @importFrom metap sumz
#'
#' @export
stoufferTable <- function(enrichmentResult)
{
    resultCheck <- lapply(enrichmentResult, is.null)
    nonEmptyTableNumber <- length(which(resultCheck == FALSE))
    enrichTables <- 
        lapply(enrichmentResult, function(enrichResult) enrichResult@result)
    if (nonEmptyTableNumber>1)
    {
        #################################################################
        ##                 Table with Stouffer's value                 ##
        #################################################################
        enrichmentResultsMergedbyID <- Reduce(function(...)
            merge(..., by=c("ID"), all=TRUE), enrichTables)
        
        ## matches used because of merging consequences on columns names
        ## pvalue.x, pvalue.y, etc.
        id_pv <- enrichmentResultsMergedbyID %>%
            dplyr::select(.data$ID, 
                          matches("^Description"), 
                          matches("^pvalue"))
        pvalues <- id_pv %>% dplyr::select(matches("^pvalue"))
        
        ## Non-weighted Stouffer's method
        Stouffer <- apply(pvalues, 1,
            function(x) suppressWarnings(sumz(p <- x,
                                            na.action=na.exclude)$p))
        
        ## Omic's weights 
        weights <- omicWeight(enrichTables)
        
        ## Computes weighted Stouffer's p-value
        StoufferWeighted <-  apply(pvalues, 1, function(pv) 
            suppressWarnings(sumz(p = pv,
                weights = weights,
                na.action = na.exclude)$p))
        
        ## Builds result table
        stTable <- cbind(enrichmentResultsMergedbyID,
                        Stouffer,
                        StoufferWeighted)
        
        ## Biological descriptions and genesets merging
        descMerged <- as.list(dplyr::select(stTable, matches("^Description")))
        allStoufferResult <- stTable %>% 
            rename("p-value:Omic" = 
                names(dplyr::select(stTable, matches("^pvalue")))) %>%
            mutate(Description = coalesce(!!!descMerged)) %>%
            dplyr::select(.data$ID, .data$Description, matches("^p-value"), 
                    Stouffer, StoufferWeighted)
        
        ## Final result table
        stoufferTableSorted <- allStoufferResult %>% 
            arrange(StoufferWeighted)
        
        
        ##################################################################
        ##                Multi-Omic enrichResult object                ##
        ##################################################################
        ## GeneSets
        geneSetMerged <- 
            lapply(enrichmentResult, function(x) x@geneSets)
        geneSets <- (Reduce(c, geneSetMerged))
        geneSets <- geneSets[stoufferTableSorted$ID]
        ## GeneID
        pathways_geneID = vapply(geneSets, function(path)
            paste0(path, collapse = "/"), FUN.VALUE = "foo")
        stoufferTableSorted$geneID <- pathways_geneID
        ## Gene
        enrichGene <- 
            lapply(enrichmentResult, function(enrichResult) enrichResult@gene)
        gene <- Reduce(unique, enrichGene)
        ## Count
        # pathways_Count = vapply(geneSets, function(path)
            # length(path), FUN.VALUE = 5)
        # stoufferTableSorted$Count <- pathways_Count
        ## GeneRatio
        # stoufferTableSorted$GeneRatio <- pathways_Count
        geneIDtoMerge <- 
            lapply(enrichmentResult, function(x) {
                y = x@result$geneID
                names(y) <- x@result$ID
                return(y)})
        pathNameList <- lapply(geneIDtoMerge, names)
        paths <- unique(Reduce(c, pathNameList))
        countDF <- data.frame(row.names = paths, count = rep(0, length(paths)))
        toCountGene <- function(geneID)
        {
            splitRes <- str_split(geneID, "/")[[1]]
            len <- length(splitRes)
            return(len)
        }
        for (pathList in geneIDtoMerge)
        {
            countGene <- vapply(pathList, toCountGene, 1)
            countDF[names(countGene),] <- 
                countDF[names(countGene), ] + countGene
        }
        GeneRatio <- paste0(countDF$count, "/", length(gene))
        stoufferTableSorted$GeneRatio <- GeneRatio
        stoufferTableSorted$Count <- countDF$count
        ## pvalues
        stoufferTableSorted$pvalue <- stoufferTableSorted$StoufferWeighted
        stoufferTableSorted$p.adjust <- stoufferTableSorted$StoufferWeighted
        
        
        ## Universe
        universe <- enrichmentResult[[1]]@universe # same db, same universe
        
        ## EnrichResult object
        moEnrichRes <- new("enrichResult",
                   result = stoufferTableSorted,
                   pAdjustMethod = "Stouffer",
                   gene = gene,
                   universe = universe,
                   geneSets = geneSets,
                   keytype = "ENTREZID",
                   ontology = "UNKNOWN",
                   organism = "UNKNOWN")
        
        
        return(list(table = stoufferTableSorted, 
                    geneset = geneSets, 
                    moEnrichRes = moEnrichRes))
    }
    else(return(FALSE))
}

#' Enrichment table ID url
#'
#' @description Adds relative url for enrichment table ID
#'
#' @param enrichTableID_col an enrichment table ID column obtained by an 
#' enrichment function.
#' @param enrichmentDB Relative biological database.
#'
#' @noRd 
addEnrichIDUrl <- function(enrichTableID_col, enrichmentDB)
{
    keggURL <- "https://www.kegg.jp/kegg-bin/show_pathway?"
    reactomeURL <- "https://reactome.org/PathwayBrowser/#/"
    wikipURL <- "https://www.wikipathways.org/index.php/Pathway:"
    goURL <- "http://amigo.geneontology.org/amigo/term/"
    
    urls <- list(kegg = keggURL, 
                 reactome = reactomeURL, 
                 wikiPathways = wikipURL,
                 MF = goURL,
                 CC = goURL,
                 BP = goURL)
    enrichTableID_col <- paste0("<a href=", 
                             paste0(urls[[enrichmentDB]], enrichTableID_col),
                             " target='_blank'>",
                             enrichTableID_col,"</a>")
}

#' Formatting enrichment table function
#'
#' @description Formats enrichment table to display in app
#'
#' @param enrichTable,db an enrichment table obtained by an 
#' enrichment function and relative database
#'
#' @noRd 
#'
#' @importFrom DT formatStyle styleInterval formatSignif
#' @importFrom dplyr mutate 
formatEnrichTable <- function(enrichTable, db = NULL)
{
    ## adds url for pathways or go ids
    if (!is.null(db))
    {
      enrichTable$ID <- addEnrichIDUrl(enrichTable$ID, db)
    }
    
    ## adds column visibility button and responsive options
    tbl_dt <- datatable(enrichTable, 
                        escape = FALSE, 
                        extensions = c('Responsive', 'Buttons'), 
                        options = list(
                            rownames = FALSE,
                            columnDefs = list(list(
                            targets = 10,
                            render = JS(
                            "function(data, type, row, meta) {",
                            "return type === 'display' && data.length > 50 ?",
                            "'<span title=\"' + data + '\">' + 
                            data.substr(0, 50) + '...</span>' : data;",
                            "}")
                          )),
                            dom = 'Bfrtip',
                          buttons = c('csv', 'excel', 'pdf')
                        ))
    
    ## adds color background for p.adjust values <= 0.01 and 0.05
    tbl_dt <- tbl_dt %>% formatStyle(
        'p.adjust',
        color = styleInterval(c(0.01, 0.05), 
                              c('white', 'white', 'black')),
        backgroundColor = styleInterval(c(0.01, 0.05), 
                                        c('#a2465f', '#cb5658', "white"))
    ) %>%
    ## To display numeric values columns with 3 number after last zero
    formatSignif(columns = c('pvalue', 'p.adjust', 'qvalue'), digits = 3)
    
    return(tbl_dt)
}

#' Formatting Stouffer's table function
#'
#' @description Formats stouffer table to display in app
#'
#' @param stoufferTable,db an table with Stouffer's values obtained by 
#' stoufferTable() function and relative database
#'
#' @noRd 
#'
#' @importFrom DT datatable formatStyle styleInterval formatSignif
#' @importFrom dplyr mutate 
formatStoufferTable <- function(stoufferTable, db = NULL, dataNbr)
{
    ## adds url for pathways or go ids
    if (!is.null(db))
    {
        stoufferTable$ID <- addEnrichIDUrl(stoufferTable$ID, db)
    }
    
    ## adds column visibility button and responsive options
    tbl_dt <- datatable(stoufferTable, escape = FALSE, 
                        extensions = c('Responsive', 'Buttons'), 
                        options = list(
                            dom = 'Bfrtip',
                            buttons = c('csv', 'excel', 'pdf'),
                            ## To choose which column to display (3th to 9th)
                            buttons = list(
                              list(extend = 'colvis', 
                                   columns = seq(3, 4+dataNbr)))
                        ))
    
    ## adds color background for p.adjust values <= 0.05 and 0.1
    tbl_dt <- tbl_dt %>% formatStyle(c("Stouffer", "StoufferWeighted"),
        color = styleInterval(c(0.05, 0.1),
                              c("white", "white", "black")),
        backgroundColor = styleInterval(c(0.05, 0.1), 
                                        c("#31abb2", "#7bc5c9", "white"))
    ) %>%
    ## To display numeric values columns with 3 number after last zero
    formatSignif(columns = seq(3, 3+dataNbr+1), digits = 3)
    
    return(tbl_dt)
}

#' Enrichment function
#'
#' @description Runs DESEQ2 analysis on all omic data blocks.
#' 
#' @param omicDataList List of omic data blocks with Y class vector.
#' @param padjUser Threshold for p-value adjusted to select features
#' according to padj values in DESeq2 table.
#' 
#' @examples 
#' data("omic2", package = "multiSight")
#' deseqRes <- runMultiDeseqAnalysis(omic2, 0.05)
#' 
#' @export 
#' 
#' @return List of DESeq2's Differential Expression tables for all omic 
#' datasets (e.g. BaseMean, Log2FoldChange, padj columns).
#' 
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
runMultiDeseqAnalysis <- function(omicDataList, padjUser)
{
    if (is(omicDataList$Y, "data.frame")) {
        omicDataClass <- omicDataList$Y$Y
    }else{
        omicDataClass <- omicDataList$Y
    }
    resList <- list()
    i <- 1
    while (i < length(omicDataList))
    {
      
        ## data
        omicData <- omicDataList[[i]]
        sampleNames <- rownames(omicData)
        omicDataForDESEQ <- t(omicData) # features in rows
        
        ## add feature names
        message("DESeq2: Converting values to integers")
        omicDataForDESEQ <- apply((omicDataForDESEQ), 2, as.integer)
        rownames(omicDataForDESEQ) <- colnames(omicData)
        
        ## meta design for DESEQ
        meta <- data.frame(sample = sampleNames, class = factor(omicDataClass))
        
        ## DESeq2 main function to DE table
        dds <- DESeqDataSetFromMatrix(countData = omicDataForDESEQ,
                                      colData = meta, 
                                      design= ~ class) #design is class meta col
        dds <- DESeq(dds)
        res <- results(dds)
        ## add res in result list
        resList$DEtable[[names(omicDataList)[i]]] <- res
        ## biosignature by deseq_p-adj filtering
        rowPadjFiltered <- which(res$padj <= padjUser)
        DEtable_padj <- res[rowPadjFiltered, ]
        selectedFeatures <- rownames(DEtable_padj)
        ## add selected features in result list
        resList$selectedFeatures[[names(omicDataList)[i]]] <- selectedFeatures
        i <- i+1
    }
    
    return(resList)
}
      
#' Enrichment function
#'
#' @description Makes and checks DESEQ2 analysis results.
#' 
#' @param input,output,session Internal parameters for shiny.
#' @param omicDataList List of omic data blocks with Y class vector.
#' 
#' @noRd 
#' 
#' @importFrom shiny renderUI span
#' 
#' @return Sets obj values for DESEQ tables
builDeseqAnalysis <- function(omicDataList, input, session, output)
{
    padjUser <- input$DEtable_padj
    resList <- list()
    out <- tryCatch(
      {
          resList <- runMultiDeseqAnalysis(omicDataList, padjUser)
          return(resList)
      },
      warning=function(w)
      {
          msg <- w$message
          output$warnDeseq <- renderUI({
            span(msg, style="color:blue")
          })
          return(resList)
      },
      error=function(e)
      {
          msg <- e$message
          output$errorDeseq <- renderUI({
            span(msg, style="color:blue")
          })
          return(resList)
      }) 
    
    out <- tryCatch(
      {
        i <- 1
        while (i < length(omicDataList))
        {
            omic <- names(omicDataList)[[i]]
            selectedFeatures <- resList$selectedFeatures[[omic]]
            vecNull <- c()
            if (length(selectedFeatures)==0)
            {
                vecNull <- c(vecNull, omic)
            }
            i <- i+1
        }
        
        ## If at least one omic does not have selected features by threshold
        # in DESeq2 relative table.
        if (length(vecNull)>0)
        {
          
            warning(paste0(paste(vecNull, collapse = ", "),
                       ": there is not features to enrich with p.adj <= ",
                       padjUser))
        }
      },
      warning=function(w)
      {
          msg <- w$message
          output$warnDeseq <- renderUI({
            span(msg, style="color:blue")
          })
          return(resList)
      },
      error=function(e)
      {
          msg <- e$message
          output$errorDeseq <- renderUI({
            span(msg, style="color:red")
          })
          return(resList)
      })
    return(resList)
}
