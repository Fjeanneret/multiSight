#' Select only data values from diablo features selected. 
#'
#' @param matrixDataList A matrix list according to omic type with
#' sample ID as columns and features as rows.
#' @param featureList A features list for each omic type.
#'
#' @examples
#' data("omic2", package = "multiSight")
#' splitData <- splitDatatoTrainTest(omic2, 0.8)
#' data.train <- splitData$data.train
#' data.test <- splitData$data.test
#' 
#' #diabloRes <- runSPLSDA(data.train)
#' data("diabloRes", package = "multiSight")
#' diabloModels <- diabloRes$model #sPLS-DA model using all omics.
#' diabloFeats <- diabloRes$biosignature #selected features for each omic.
#' omicMatrices <- getDataSelectedFeatures(omic2, diabloFeats)
#' 
#' @return A concatenated matrix of all omic selected features.
#'
#' @export
getDataSelectedFeatures <- function(matrixDataList, featureList = NULL)
{
    if (!is.null(featureList))
    {
          
        nData <- length(featureList)
        selectedFeatData <- 
            lapply(seq(1, nData), function(j) 
                {
                    matrixDataList[[j]][, featureList[[j]]]
                })
        
        ## To indicate from which omic come from features
        i <- 1
        while(i <= nData)
        {
            newColnames <- paste0(colnames(selectedFeatData[[i]]), "_", i)
            colnames(selectedFeatData[[i]]) <- newColnames
            i <- i+1
        }
    }else {
        selectedFeatData <- matrixDataList
    }
  
    concatenatedMatrix <- do.call(cbind, selectedFeatData)
    return(concatenatedMatrix)
}


#' Correlation network inference
#' 
#' @param concatenatedMatrix A concatenated matrix of all omic selected 
#' features. Returned by getDataSelectedFeatures().
#' @param valueThreshold Correlation absolute value threshold to select 
#' relevant values.
#'
#' @examples
#' data("omic2", package = "multiSight")
#' splitData <- splitDatatoTrainTest(omic2, 0.8)
#' data.train <- splitData$data.train
#' data.test <- splitData$data.test
#' 
#' #diabloRes <- runSPLSDA(data.train)
#' data("diabloRes", package = "multiSight")
#' diabloModels <- diabloRes$model #sPLS-DA model using all omics.
#' diabloFeats <- diabloRes$biosignature #selected features for each omic.
#' omicMatrices <- getDataSelectedFeatures(omic2, diabloFeats)
#' correlationNetworkInference(omicMatrices, 0.8)
#'
#' @return Each result for a correlation network of selected features 
#' according to threshold values.
#'
#' 
#' @import igraph
#' @importFrom stats cor na.exclude
#' @importFrom stringr str_split
#' @importFrom dplyr filter
#' 
#' @export
correlationNetworkInference <- function(concatenatedMatrix, valueThreshold)
{
    ## Correlation matrix values
    corr.matrix <- round(cor(concatenatedMatrix,
                             method="pearson"), 4)
    ## Correlation best couples
    pvalueThreshold <- 1
    corrMatrixRelation <- flattenCorrMatrix(corr.matrix, 
                                            pvThreshold = pvalueThreshold,
                                            numThreshold = valueThreshold)
    ## graph
    if (nrow(corrMatrixRelation$flattenTable) > 0)
    {
        g <- mat2graph(corr.matrix, 
                          corrMatrixRelation$flattenRelevantRow, 
                          valueThreshold)
        graph <- addEdgeWidth(g, corrMatrixRelation$flattenTable)
        graph <- buildNetworkD3(graph, valueThreshold)
    }
    else
    {
        graph <- make_empty_graph() + vertices("No relation for this threshold")
    }
    
    cRes <- list(flattenTable = corrMatrixRelation$flattenTable,
                 num.matrix = corr.matrix, 
                 graph = graph)
    return(cRes)
}

#' Partial Correlation network inference
#' 
#' @param concatenatedMatrix A concatenated matrix of all omic selected 
#'   features. Returned by getDataSelectedFeatures().
#' @param valueThreshold Partial Correlation Value threshold to select 
#' relevant values
#'
#' @examples
#' data("omic2", package = "multiSight")
#' splitData <- splitDatatoTrainTest(omic2, 0.8)
#' data.train <- splitData$data.train
#' data.test <- splitData$data.test
#' 
#' #diabloRes <- runSPLSDA(data.train)
#' data("diabloRes", package = "multiSight")
#' diabloModels <- diabloRes$model #sPLS-DA model using all omics.
#' diabloFeats <- diabloRes$biosignature #selected features for each omic.
#' omicMatrices <- getDataSelectedFeatures(omic2, diabloFeats)
#' partialCorrelationNI(omicMatrices, 0.8)
#'
#' @return Each result for a partial correlation network of selected features 
#' according to threshold values.
#'
#' @importFrom ppcor pcor
#' 
#' @export
partialCorrelationNI <- function(concatenatedMatrix, valueThreshold)
{
    pcor.res <- NULL
    out <- tryCatch(
    {
        suppressWarnings({ pcor.res <- pcor(concatenatedMatrix) })
        
    }, error=function(e)
    {
      suppressWarnings({pcor.res <- pcor(concatenatedMatrix, 
                         method = "kendall") })
        return(pcor.res)
    })
    
    if (is.null(pcor.res))
    {
        pcor.res <- out    
    }
    pcor.matrix <- pcor.res$estimate
    dimnames(pcor.matrix) <- dimnames(pcor.res$p.value) <- 
        list(colnames(concatenatedMatrix), colnames(concatenatedMatrix))
    # pvalueThreshold <- 1
    
    # When Moore-Penrose generalized matrix inverse used by pcor(), check
    # valueTotalNb <- ncol(pcor.res$p.value)**2-ncol(pcor.res$p.value)
    # if (length(which(is.na(pcor.res$p.value))) != valueTotalNb)
    # {
    #     pcorMatrixRelation <- flattenCorrMatrix(round(pcor.res$estimate, 4),
    #                                           pvAdded = pcor.res$p.value,
    #                                           pvThreshold = pvalueThreshold,
    #                                           numThreshold = valueThreshold)
    # }else
    # {
        pcorMatrixRelation <- flattenCorrMatrix(round(pcor.matrix, 4), 
                                            numThreshold = valueThreshold) 
    # }
    
    ## graph
    if (nrow(pcorMatrixRelation$flattenTable) > 0)
    {
        g <- mat2graph(pcor.matrix, 
                          pcorMatrixRelation$flattenRelevantRow, 
                          valueThreshold)
        graph <- addEdgeWidth(g, pcorMatrixRelation$flattenTable)
        graph <- buildNetworkD3(graph, valueThreshold)
    }
    else
    {
        graph <- make_empty_graph() + vertices("No relation for this threshold")
    }
    
    
    pcRes <- list(flattenTable = pcorMatrixRelation$flattenTable,
                  num.matrix = pcor.matrix,
                  graph = graph)
    return(pcRes)
}

#' Mutual Information network inference
#' 
#' @param concatenatedMatrix A concatenated matrix of all omic selected 
#'   features. Returned by getDataSelectedFeatures().
#' @param valueThreshold Mutual Information Value threshold to select relevant 
#'   values
#'
#' @examples
#' data("omic2", package = "multiSight")
#' splitData <- splitDatatoTrainTest(omic2, 0.8)
#' data.train <- splitData$data.train
#' data.test <- splitData$data.test
#' 
#' #diabloRes <- runSPLSDA(data.train)
#' data("diabloRes", package = "multiSight")
#' diabloModels <- diabloRes$model #sPLS-DA model using all omics.
#' diabloFeats <- diabloRes$biosignature #selected features for each omic.
#' omicMatrices <- getDataSelectedFeatures(omic2, diabloFeats)
#' mutualInformationNI(omicMatrices, 0.8)
#'
#' @return Each result for a mutual information network of selected features 
#' according to threshold values.
#'
#' @importFrom infotheo mutinformation discretize
#' 
#' @export
mutualInformationNI <- function(concatenatedMatrix, valueThreshold)
{
    mat2discretDF <- data.frame(discretize(concatenatedMatrix), 
                                check.names = FALSE)
    mi.matrix <- round(mutinformation(mat2discretDF), 4)
    
    miMatrixRelation <- flattenCorrMatrix(mi.matrix,
                                          numThreshold = valueThreshold)
    
    ## graph
    if (nrow(miMatrixRelation$flattenTable) > 0)
    {
        g <- mat2graph(mi.matrix, 
                          miMatrixRelation$flattenRelevantRow, 
                          valueThreshold)
        graph <- addEdgeWidth(g, miMatrixRelation$flattenTable)
        graph <- buildNetworkD3(graph, valueThreshold)
    }
    else
    {
        graph <- make_empty_graph() + vertices("No relation for this threshold")
    }
    
    ## Result list
    miRes <- list(flattenTable = miMatrixRelation$flattenTable, 
                  num.matrix = mi.matrix, 
                  graph = graph)
    return(miRes)
  
}

#' Sets edge widths with numeric relation value 
#' 
#' @param graph graph object obtained by correlationNI,
#' partialCorrelationNI or mutualInformationNI.
#' @param flattenRelation values computed by relation method 
#' (e.g. correlation values) formatted in a dataframe where 
#' two first columns are features (e.g var1 and var2) and 
#' where "value" is relation value (e.g. correlation between 
#' var1 and var2).
#' 
#' @importFrom stringr str_split
#' @importFrom dplyr filter .data
#' @importFrom networkD3 igraph_to_networkD3 forceNetwork
#' 
#' @noRd
#' 
#' @return Graph with edge width according to relation method.
addEdgeWidth <- function(graph, flattenRelation)
{
    ## To retrieve value for two features of each relation
    getVertexValue <- function(vertexIdx, edges, flattenRelation)
    {
        vertexName <- attr(edges[vertexIdx],
                           "vnames")
        feats <- str_split(vertexName, "\\|")
        feat1 <- feats[[1]][1]
        feat2 <- feats[[1]][2]
        
        value <- (flattenRelation %>% filter(.data$feature1 == feat1, 
                                             .data$feature2 == feat2))$value
        return(value)
    }
    
    edges <- E(graph)
    edgesNb <- ecount(graph)
    
    i = 1
    verticeWidth <- c()
    while (i <= edgesNb)
    {
        value <- abs(getVertexValue(i, edges, flattenRelation))
        verticeWidth <- c(verticeWidth, value)
        
        i = i + 1
    }
    
    ## Set edge widths with relation values
    # E(graph)$width  <- verticeWidth
    
    return(graph)
}
