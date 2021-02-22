#' MLmodels util function
#'
#' @description To split list of omic data in data train and data test subsets.
#'
#' @param MultiOmicData list of your data blocks
#' @param freq Split proportion of train samples
#' 
#' @examples
#' data("omic2.rda")
#' splittedData <- splitDatatoTrainTest(omic2, 0.8)
#' data.train <- splittedData$data.train
#' data.test <- splittedData$data.test
#' 
#' @export
#' 
#' @return Returns two data sets: first to train model and second to assess it.
splitDatatoTrainTest <- function(MultiOmicData, freq = 0.8)
{
    ## Random splitting indexes
    if (is(MultiOmicData$Y, "data.frame")) 
        {MultiOmicData$Y <- MultiOmicData$Y$Y}
    sampleLen <- length(MultiOmicData$Y)
    sampleSeq <- seq(1, sampleLen)
    train_index <- sample(sampleSeq, freq * sampleLen)
    
    ## Numeric matrices indexing
    listLen <- length(MultiOmicData)
    data.train <- lapply(MultiOmicData[-listLen], function(x)x[train_index,])
    data.test <- lapply(MultiOmicData[-listLen], function(x)x[-train_index,])
    
    ## Class vector indexing
    data.train$Y <- MultiOmicData$Y[train_index]
    data.test$Y <- MultiOmicData$Y[-train_index]
    
    dataSplitted <- list(data.train = data.train, data.test = data.test)
    
    return(dataSplitted)
}


#' diablo util function
#'
#' @description A covariance matrix design builder function for splsda model.
#' Builds covariance matrix n*n with n number of omic data sets.
#'
#' @param MultiOmicData list of your data blocks
#'
#' @noRd 
#'
#' @return A symmetric design matrix providing correlations to maximize 
#' between omic data sets.
#'
#' @importFrom mixOmics pls 
buildCovarianceDesign <- function(MultiOmicData)
{
  
    ## Builds 0 matrix n*n with n number of omic data sets.
    seqOmicData = seq(1, length(MultiOmicData))
    design <- matrix(0, ncol = length(MultiOmicData), 
                  nrow = length(MultiOmicData),
                  dimnames = list(names(MultiOmicData)[seqOmicData], 
                  names(MultiOmicData)[seqOmicData]))
    
    i <- 1
    while(i <= length(MultiOmicData))
    {
        j <- 1
        while (j <= length(MultiOmicData))
        {
            if (i != j)
            {
                plsModel <- pls(MultiOmicData[[i]], 
                                MultiOmicData[[j]], 
                                ncomp = 1)
                corrValue <- cor(plsModel$variates$X, 
                                 plsModel$variates$Y)
                design[i, j] <- design[j, i] <- as.numeric(corrValue)
            }
            j <- j+1
        }
        
        i <- i+1
    }
    return(design)
}

#' diablo util function
#'
#' @description A component number optimization function of splsda model.
#' 
#' Computes Balanced Error Rate (BER) using 3 types of distances and 
#' uses mean of ncomp according to these 3 models.
#'
#' @param dataTrainList List of data block training part.
#' @param YClassVector List of your sample classes vector provided 
#' @param design Covariance matrix design obtained from 
#' buildCovarianceDesign function launched by runSPLSDAmodels_Diablo function.
#'
#' @importFrom shiny NS tagList 
#' @importFrom mixOmics block.splsda perf 
#' 
#' @noRd
#' 
#' @return Returns model performances plot.
runComponentNumberTest <- function(dataTrainList, YClassVector, design) 
{
    splsdaModelToTestComponentNumber <- block.splsda(
        X = dataTrainList, 
        Y = YClassVector, 
        ncomp = 8, 
        scale = TRUE, 
        near.zero.var = FALSE,
        design = design)
    
    perf.comp <- perf(splsdaModelToTestComponentNumber, 
                      validation = "Mfold", 
                      folds = 10, 
                      nrepeat = 5)
    plot(perf.comp)
    
    ncomp_max.dist <- 
      perf.comp$choice.ncomp$WeightedVote["Overall.BER", "max.dist"] 
    ncomp_centroids.dist <- 
      perf.comp$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
    ncomp_maha.dist <- 
      perf.comp$choice.ncomp$WeightedVote["Overall.BER", "mahalanobis.dist"] 
    
    ## Chooses mean ncomp according to different distance values
    ncomp <- ceiling(mean(c(ncomp_max.dist, 
                            ncomp_centroids.dist,
                            ncomp_maha.dist)))
    
    return(list("ncomp" = ncomp, "plotPerf" = perf.comp))
}

#' diablo util function
#'
#' @description A feature number tuning function of splsda model
#' for each component of each data block.
#'
#' @param dataTrainList List of data block training part.
#' @param YClassVector List of your sample classes vector provided 
#' @param ncomp Component number in splsda model obtained from 
#' runComponentNumberTest function launched by runSPLSDAmodels_Diablo 
#' function.
#' @param design Covariance matrix design obtained from 
#'
#' @importFrom shiny NS tagList 
#' @importFrom mixOmics tune.block.splsda 
#' 
#' @return Grid of number of features to select by sPLS-DA model for each 
#' component and omic.
runFeatureNumberTuning <- function(dataTrainList, YClassVector, ncomp, design)
{
    ## To provide at least 30 features to enrich 
    min_varNumber <- round(30/ncomp, 0) 
    
    ## Builds number feature selection grid for each omic block
    test.keepX <- list()
    i <- 1
    while(i <= length(dataTrainList))
    {
        test.keepX[[names(dataTrainList)[i]]] <-  c(min_varNumber, 
                                                   min_varNumber+5, 
                                                   min_varNumber+10)
        i <- i+1
    }
    
    ## Feature number selection tuning
    tune.featuresNumberSelection <- tune.block.splsda(
        X = dataTrainList, 
        Y = YClassVector, 
        ncomp = ncomp, 
        scale = TRUE,
        test.keepX = test.keepX, 
        design = design,
        validation = 'Mfold', 
        folds = 10, 
        nrepeat = 5, 
        cpu = 2, 
        dist = "mahalanobis.dist")
    
    ## Feature number grid
    list.keepX = tune.featuresNumberSelection$choice.keepX
    
    return(list.keepX)
}

#' diablo util function
#'
#' @description A feature extraction function of sparse-PLS-DA model
#' for each component of each data block.
#'
#' @param splsdaModel sPLS-DA model obtained from block.splsda 
#' function called in runSPLSDAmodels_Diablo function. 
#' 
#' @noRd
#' 
#' @return Returns features selected by diablo sPLS-DA model for 
#' each omic data set.
getSelectedFeatures <- function(splsdaModel)
{
    allVars <- list()
    i <- 1
    ## For each omic block
    while(i < length(splsdaModel$names$blocks))
    {
        block <- splsdaModel$names$blocks[i]
        block_var <- c()
        j <- 1
        ## for each component
        while(j <= splsdaModel$ncomp[1])
        {
            vars <- selectVar(splsdaModel, 
                              block = block, 
                              comp = j)[[block]]$name
            block_var <- c(block_var, vars)
            j <- j+1
        }
        allVars[[block]] = unique(block_var)
        i <- i+1
    } ## end for each omic block
    
    return(allVars)
}

#' Models util function
#'
#' @description Builds ui outputs of detailed feature tables for selected 
#' ones by machine learning models.
#'
#' @param listFeatTables datatable of selected features.
#' 
#' @importFrom DT datatable renderDT
#' 
#' @noRd 
displayFeatDetails <- function(listFeatTables)
{   
    omicNames <- names(listFeatTables)
    statusLabels <- c("primary", 
                      "success", 
                      "info", 
                      "warning", 
                      "danger")
    
    seqFeat <- seq(1, length(listFeatTables))
    lapply(seqFeat, function(i){
      table <- listFeatTables[[i]]
      status <- sample(statusLabels, 1)
      
      ## UI output
      box(width = 10, 
          title = omicNames[i], 
          status = status,
          solidHeader = TRUE,
          
          renderDT({ table }, server = FALSE)
        
    )
  })
}

#' Models util function
#'
#' @description Builds ui outputs of detailed feature tables for selected 
#' ones by machine learning models.
#'
#' @param featuresList Selected features list, one vector by omic dataset.
#' @param modelMethod Model used to save data.
#' @param DESeqTables Deseq2 results table for an omic dataset.
#' @param obj R6 object to wrap all data from different analysis
#' 
#' @importFrom DT renderDT datatable
#' 
#' @noRd 
computeFeatDetails <- function(featuresList, 
                               modelMethod, 
                               DESeqTables = NULL,
                               obj)
{   
    omicNames <- names(featuresList)
    ## load data to compute feature information
    wholeOmicData <- obj$data$wholeData
    wholeY <- obj$data$wholeData$Y
    
    ## For each omic signature build one table
    seqFeat <- seq(1, length(featuresList))
    listFeatTables <- lapply(seqFeat, function(i){
      omic <- omicNames[[i]]
      table <- buildFeatTable(featuresList[[i]], 
                              wholeOmicData[[i]], 
                              wholeY, 
                              DESeqTables[[omic]])
      ## To save in obj
      objOmicName <- omicNames[i]
      table
    })
    names(listFeatTables) <- paste0("Omic", seqFeat)
    return(listFeatTables)
}

#' Models util function
#'
#' @description To build detailed feature tables outputs for selected ones by 
#' machine learning models. Relative mean values by label class then 
#' log2FoldChange and p.adj values if DESeq2 have been computed are indicated.
#'
#' @param featVec Selected features vector of one omic dataset.
#' @param omicBlock Omic dataset of selected features.
#' @param Y Omic sample classes.
#' @param deTable (optional) Deseq2 results table for an omic dataset.
#'
#' @examples
#' splittedData <- splitDatatoTrainTest(omic2, 0.8)
#' data.train <- splittedData$data.train
#' data.test <- splittedData$data.test
#' 
#' diabloRes <- runSPLSDA(data.train)
#' diabloModels <- diabloRes$model #sPLS-DA model using all omics.
#' diabloFeats <- diabloRes$biosignature #selected features for each omic.
#' diabloFeatTable <- buildFeatTable(diabloFeats[[1]], 
#'   omic2[[1]], 
#'   omic2$Y)
#' diabloFeatTable
#' 
#' @importFrom DT datatable formatSignif formatStyle
#' 
#' @return Returns table of features selected by classification model and
#' relative values.
#' 
#' @export
buildFeatTable <- function(featVec, omicBlock, Y, deTable = NULL)
{
    ## builds table with new columns
    df <- data.frame(features = featVec)
    
    ## Mean column
    mean_col <- vapply(df$features, function(feat) 
        mean(omicBlock[, feat]), 
        FUN.VALUE = numeric(1))
    df$mean <- mean_col
    
    ## Mean by class
    wholeY <- factor(Y)
    class <- levels(wholeY)
    for (label in class)
    {
        idx_label <- wholeY == label
        mean_label <- vapply(df$features, function(feat) 
            mean(omicBlock[idx_label, feat]), 
            FUN.VALUE = numeric(1))
        
        col_name <- paste0("mean_", label)
        df[[col_name]] <- mean_label
    }
    
    ## logF and p-adj if available after DESeq2
    if (!is.null(deTable))
    {
        logF <- vapply(df$features, function(feat)
          deTable[feat, "log2FoldChange"], 
          FUN.VALUE = numeric(1))
        
        df[["log2FoldChange"]] <- logF
        
        padj <- vapply(df$features, function(feat)
          deTable[feat, "padj"], 
          FUN.VALUE = numeric(1))
        
        df[["padj"]] <- padj
    }
    
    ## Table 
    featTable <- DT::datatable(df,
                           extensions = c('Responsive', 'Buttons'),
                           options = list(
                              dom = 'Bfrtip',
                              rownames = FALSE,
                              buttons = c('csv', 'excel', 'pdf')
                           )
    )
    colToFormat <- names(df)[-1]
    featTable <- featTable %>% formatSignif(columns = colToFormat, digits = 3)
    
    if (!is.null(df$padj))
    {
        ## adds color background for p.adjust values <= 0.01 and 0.05
        featTable <- featTable %>% formatStyle(
            'padj',
            color = styleInterval(c(0.01, 0.05), 
                                  c('white', 'white', 'black')),
            backgroundColor = styleInterval(c(0.01, 0.05), 
                                            c('#a2465f', '#cb5658', "white")))
    }
    return(featTable)
}
