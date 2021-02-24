#' MLmodels diablo function
#'
#' @description A model and selection features function from mixOmics.
#'
# #' @param obj R6 object to wrap all data from different analysis
#' @param dataTrain Data train set to build classification models.
#' Returned by splitDatatoTrainTest().
#'
#' @examples
#' data("omic2", package = "multiSight")
#' splittedData <- splitDatatoTrainTest(omic2, 0.8)
#' data.train <- splittedData$data.train
#' data.test <- splittedData$data.test
#' 
#' #diabloRes <- runSPLSDA(data.train)
#' data("diabloRes", package = "multiSight")
#' diabloModels <- diabloRes$model #sPLS-DA model using all omics.
#' diabloFeats <- diabloRes$biosignature #selected features for each omic.
#' 
#' @export
#'
#' @importFrom mixOmics block.splsda perf tune.block.splsda
#' 
#' @return Returns sPLS-DA model, performances, features selected.
runSPLSDA <- function(dataTrain)
{
    ## data to train sPLS-DA model
    seq_lenDataTrain <- seq(1, length(dataTrain)-1)
    omicTrain <- dataTrain[seq_lenDataTrain] # No Y vector
    Ytrain <- factor(dataTrain$Y)
    
    ## design matrix
    design <- buildCovarianceDesign(omicTrain)
    
    ## component number test
    compTestRes <- runComponentNumberTest(dataTrainList = omicTrain, 
                                          YClassVector = Ytrain, 
                                          design = design) 
    ncomp <- compTestRes$ncomp
  
    ## Tuning number of selected features by component and by omic
    featsNumberList <- runFeatureNumberTuning(omicTrain, 
                                                  Ytrain, 
                                                  ncomp, 
                                                  design)
    
    ## Final model
    splsdaModel <- block.splsda(X = omicTrain,
                                Y = Ytrain, 
                                ncomp = ncomp, 
                                scale = TRUE,
                                keepX = featsNumberList, 
                                design = design)
    
    ## Get selected features from each data block 
    selectedFeats <- getSelectedFeatures(splsdaModel)
    
    ## Results list
    diabloRes <- list(design = design, 
                      model = splsdaModel,
                      biosignature = selectedFeats)
    
    return(diabloRes)
}


#' MLmodels diablo function
#'
#' @description A diablo models assessing function.
#'
#' @param splsdaModel sPLS-DA model computed and returned inside
#'  runSPLSDAmodels_Diablo() results.
#' @param dataTest List of new omic data sets and samples for same omics that
#' training to test model performances. Returned by splitDatatoTrainTest().
#'
#' @examples
#' data("omic2", package = "multiSight")
#' splittedData <- splitDatatoTrainTest(omic2, 0.8)
#' data.train <- splittedData$data.train
#' data.test <- splittedData$data.test
#' 
#' #diabloRes <- runSPLSDA(data.train)
#' data("diabloRes", package = "multiSight")
#' diabloModels <- diabloRes$model #sPLS-DA model using all omics.
#' diabloFeats <- diabloRes$biosignature #selected features for each omic.
#' perfDiablo <- assessPerformance_Diablo(diabloModels, data.test)
#' perfDiablo$Omic1 # sPLS-DA's perf for omic1 data block.
#'
#' @export
#'
#' @importFrom shiny NS tagList 
#' @importFrom mixOmics auroc selectVar
#' 
#' @return Confusion matrix for sPLS-DA model
assessPerformance_Diablo <- function(splsdaModel, 
                                     dataTest)
{
    ## data
    seq_lenDataTest <- seq(1, length(dataTest)-1)
    omicTest <- dataTest[seq_lenDataTest]
    Ytest <- factor(dataTest$Y)
    
    ## model
    ncomp <- splsdaModel$ncomp[[1]]
    predict.diablo <- predict(splsdaModel, newdata = omicTest)
    
    ## Formats predictions to display by confusion matrix
    ## Anticipating predictions as NA by sPLS-DA model.
    preds <- factor(predict.diablo$MajorityVote$centroids.dist[, ncomp])
    levels(Ytest) <- c(levels(Ytest), "NA")
    levels(preds) <- c(levels(preds), "NA")
    preds[is.na(preds)] <- "NA"
    
    ## Performances
    confMat <- confusionMatrix(data = table((preds), Ytest))
    return(confMat)
    
    # aurocPlot = list()
    # i=1
    # while(i <= length(dataTest))
    # {
    #   aurocPlot <- auroc(splsdaModel, 
    #     newdata = dataTest, 
    #     outcome.test = Y.test, 
    #     roc.comp = 3,
    #     roc.block = names(dataTest)[i],print = TRUE)
    #   
    #   plotName <- paste0("graph.", names(dataTest)[i])
    #   obj$classification$diabloResult$performance[[i]] <- aurocPlot[plotName]
    #   i = i+1
    # }
    
    # (obj$classification$diabloResult$performance$Omic1$graph.Omic1$comp3)
}