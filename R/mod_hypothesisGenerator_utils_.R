#' From matrix relation values to relation as rows in dataframe. 
#'
#' @param relationMatrix A relation matrix between features 
#' (e.g. correlation matrix, partial correlation matrix)
#' @param pvAdded Optional supplementary values 
#' (e.g. partial correlation p-values).
#' @param pvThreshold p-value threshold 
#' @param numThreshold numeric threshold 
#'
#' @noRd
#' 
#' @return One table with relation values as rows
#'  of all omic selected features.
#'
#' @importFrom dplyr arrange desc
flattenCorrMatrix <- function(relationMatrix,
                              pvThreshold, 
                              numThreshold,
                              pvAdded)
{
    feature1 <- feature2 <- value <- NULL
    ut <- upper.tri(relationMatrix)
    flattenTable <- data.frame(
        feature1 = rownames(relationMatrix)[row(relationMatrix)[ut]],
        feature2 = rownames(relationMatrix)[col(relationMatrix)[ut]],
        value = (relationMatrix)[ut])
    
    if(!missing(pvAdded) && !missing(pvThreshold))
    {
        flattenTable$value <- pvAdded[ut]
        pvIndex <- which(flattenTable$value <= pvThreshold)
        
        flattenTable <- flattenTable[pvIndex,]
    }
    
    if(!missing(numThreshold))
    {
        numIndex <- which(flattenTable$value >= numThreshold)
        numIndex <- c(numIndex, which(flattenTable$value <= -(numThreshold)))
        
        flattenTable <- flattenTable[numIndex,]
    }
    
    flattenRelevantRow <- cbind(as.vector(flattenTable$feature1), 
                                as.vector(flattenTable$feature2))
    
    flattenTable <- arrange(flattenTable, desc(value))
    
    flatRes <- list(flattenTable = flattenTable, 
                    flattenRelevantRow = flattenRelevantRow)
    return(flatRes)
}

#' Create multi-omic network.
#'
#' @param matrixValue Relation matrix with n x n values (e.g correlation matrix)
#' @param relevantFeatureCouples A 2-column matrix with all feature couples
#' relevant according to thresholds provided; get by flattenCorrMatrix function.
#' @param numThreshold A numeric threshold to plot only best relation values.
#'
#' @noRd
#' 
#' @return The correlation network of selected features according to 
#' threshold values.
#' 
#' @importFrom igraph graph.adjacency vcount ecount
#'
mat2graph <- function(matrixValue, relevantFeatureCouples, numThreshold)
{
    feat4graph <- unique(Reduce(c, relevantFeatureCouples))
    mat4graph <- abs(matrixValue[feat4graph, feat4graph])
    mat4graph[mat4graph <= numThreshold] <- 0
    diag(mat4graph) <- 0
    
    graph <- graph.adjacency(mat4graph, 
                             weighted=TRUE, 
                             mode="min", 
                             diag = FALSE)
    return(graphPlot = graph)
}

#' Sets vertex colors according to omic datasets 
#' 
#' @param graph graph object obtained by correlationNI,
#'   partialCorrelationNI or mutualInformationNI.
#' @param threshold numeric threshold to display relevant edges
#' 
#' @importFrom stringr str_split
# #' @importFrom ggsci pal_npg
#' @importFrom igraph vertex_attr
#' @importFrom utils tail
#' 
#' @noRd
#' 
#' @return New graph with networkD3 properties.
buildNetworkD3 <- function(graph, threshold)
{
    ## Retrieve omic index as last split par by "_"
    nodeNames <- attr(V(graph), "names")
    nodeOmicId <- vapply(nodeNames, function(x) 
        tail(str_split(x, "_")[[1]], n=1),
        FUN.VALUE = character(1))
    omicNum <- factor(nodeOmicId)
    # omicNumLabel <- as.numeric(levels(omicNum))
    
    ## Changes omic indexes to colors
    levels(omicNum) <- paste0("Omic ", levels(omicNum))
    # col = pal_npg()(5)[omicNumLabel]
    # vertex_attr(graph)$color <- omicNum
    
    ## convert to networkD3 network
    group <- cluster_walktrap(graph)
    grp <- membership(group)
    convertedGraph <- igraph_to_networkD3(graph, 
                                          group = omicNum)
    newGraph <- forceNetwork(Links = convertedGraph$links,
                             Nodes = convertedGraph$nodes, 
                             zoom = TRUE, 
                             NodeID = "name",
                             Group = "group",
                             charge = -10,
                             opacityNoHover = FALSE,
                             opacity = 1,
                             legend = TRUE)
    # colourScale = col)
    newGraph <- newGraph
        # htmlwidgets::prependContent(newGraph, 
        #                             htmltools::tags$i(
        #                                 paste0("Threshold value: ", 
                                                             # threshold))) 
    
    return(newGraph)
}

#' InferenceNetwork util function
#'
#' @description Selects features from relation table and update pubmed query 
#' string 
#' 
#' @param input,output,session Internal parameters for {shiny}.
#' @param method Relation method (e.g. correlation, partial correlation)
#' @param table Relation table with feature names and value
#' (e.g. correlation, partial correlation)

#' @noRd 
#' 
#' @importFrom shiny updateTextInput observeEvent
#' @importFrom stringr str_split
buildPubmedQuery <- function(input, output, session, method, table)
{
    methodTableId <- paste0(method, "Couple_rows_selected")
    observeEvent(input[[methodTableId]], {
        ## Retrieve couple features
        couples <- table[input[[methodTableId]], c(1, 2)]
        newRequest_Oid <- unique(c(couples$feature1, 
                                   couples$feature2))
        ## remove omic number id 
        newRequest <- vapply(newRequest_Oid, 
                             function(x) str_split(x, "_")[[1]][1],
                             FUN.VALUE = character(1))
        ## Build pubmed query
        newRequest <- paste0("(", newRequest, ")")
        request <- paste(newRequest, collapse = " OR ")
        ## Update user text input on client side
        queryUiId <- paste0("PubMedRequest_", method)
        updateTextInput(session, queryUiId, value = request)
    })
}

#' InferenceNetwork util function
#'
#' @description Merge PubMed UI inputs generated by buildPubmedQuery() 
#' function.
#'
#' @param input,output,session Internal parameters for {shiny}.

#' @noRd 
#' 
#' @importFrom shiny updateTextInput observeEvent
mergePubmedQuery <- function(input, output, session)
{ 
    ## Build final Pubmed request 
    observeEvent(input$btn_allRequest, {
        corrQuery <- input$PubMedRequest_correlation
        pcorQuery <- input$PubMedRequest_partialCor
        miQuery <- input$PubMedRequest_mutualInf
        allRequests <- unique(c(corrQuery, pcorQuery, miQuery))
        allMethodRelations <- paste(allRequests, collapse = " ")
        ## Updates pubmed query input on client side
        updateTextInput(session, "BuildPubMedRequest", 
                        value = allMethodRelations)
    })
}

#' InferenceNetwork util function
#'
#' @description Build network inference: graph and numeric relation table.
#' 
#' @param input,output,session Internal parameters for {shiny}.
#' @param featMethod Feature selection source (e.g. diablo, biosigner)
#' @param concatenatedOmicData A concatenated matrix of all omic selected 
#' features. Returned by getDataSelectedFeatures().
#' @noRd 
#' 
#' @importFrom shiny eventReactive
#' @importFrom networkD3 renderForceNetwork
# #' @importFrom htmlwidgets onRender
buildNetworkInference <- function(input, 
                                  output, 
                                  session, 
                                  method, 
                                  featMethod,
                                  concatenatedOmicData)
{ 
    thresholdID <- paste0(method, "NumThreshold")
    eventReactive(input[[thresholdID]], {
        threshold <- input[[thresholdID]]
        ## Network Inference function according to method
        niMethodFunctions <- list("correlation" = correlationNetworkInference,
                                  "partialCor" = partialCorrelationNI,
                                  "mutualInf" = mutualInformationNI)
        netInf <-  niMethodFunctions[[method]]
        ## Correlation computing
        network <- netInf(concatenatedOmicData, 
                          threshold)
        
        ## obj
        obj$networkInference[[featMethod]]$relationMatrix[[method]] <- 
            network$num.matrix
        obj$networkInference[[featMethod]]$graph[[method]] <- 
            network$graph
        
        ## Output UI
        plotID <- paste0(method, "Plot") ####
        # output[[plotID]] <- 
        networkID <- paste0(method, "Network") ####
        output[[networkID]] <- renderForceNetwork(network$graph)
            # onRender(network$graph, 
            #          "function(el,x) 
            #              { d3.selectAll('.node').on('mouseover', null); }"))
        ## Correlation table
        relationTable <- (network$flattenTable)
        relationTable
    })
}
