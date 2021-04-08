#' MLmodels biosigner function
#'
#' @description A model and selection features function from biosigner.
#'
#' @param dataTrain Data train set to build classification models. 
#' Returned by splitDatatoTrainTest(). List of Omic blocks and Y class vector.
#'
#' @examples
#' data("omic2", package = "multiSight")
#' splitData <- splitDatatoTrainTest(omic2, 0.8)
#' data.train <- splitData$data.train
#' data.test <- splitData$data.test
#' 
#' #biosignerRes <- runSVMRFmodels_Biosigner(data.train)
#' data("biosignerRes", package = "multiSight") 
#' biosignerModels <- biosignerRes$model #list of SVM/RF models for each omic.
#' biosignerFeats <- biosignerRes$biosignature #selected features for each omic.
#' 
#' @importFrom shiny NS tagList 
#' @importFrom biosigner biosign getSignatureLs
#' 
#' @return Models and features selected for each omic block in dataTrain.
#' 
#' @export
runSVMRFmodels_Biosigner <- function(dataTrain)
{
    dataBlockName <- names(dataTrain)
    
    ## svm and rf models building to feature selection
    listBiosModels <- list()
    listBiosFeats <- list()
    i <- 1
    while(i < length(dataBlockName))
    {
        svm_rf_models <- biosign(
          dataTrain[[i]], 
          dataTrain$Y, 
          bootI = 5, # bootstrap number
          methodVc = c("svm", "randomforest"),
          fig.pdfC = c("none"),
          info.txtC = c("none"))
        
        selectedFeatures <- getSignatureLs(svm_rf_models, 
                                           tierC = c("S", "AS")[2])$complete
        
        ## R6 object obj updating with biosigner results
        listBiosModels[[dataBlockName[i]]] <- svm_rf_models
        listBiosFeats[[dataBlockName[i]]] <- selectedFeatures
        
        i <- i+1
    }
    biosignerRes <- list(model = listBiosModels, 
                         biosignature = listBiosFeats)
    return(biosignerRes)
}

#' MLmodels biosigner function
#'
#' @description A biosigner models assessing function.
#' to display by confusion matrices for SVM and RF models.
#'
#' @param modelList Models list computed and returned inside
#'  runSVMRFmodels_Biosigner() results.
#' @param dataTest List of new omic datasets and samples for same omics that
#' training to test model performances. Returned by splitDatatoTrainTest().
#' 
#' @examples
#' data("omic2", package = "multiSight")
#' splitData <- splitDatatoTrainTest(omic2, 0.8)
#' data.train <- splitData$data.train
#' data.test <- splitData$data.test
#' 
#' #biosignerRes <- runSVMRFmodels_Biosigner(data.train)
#' data("biosignerRes", package = "multiSight") 
#' biosignerModels <- biosignerRes$model #list of SVM/RF models for each omic.
#' biosignerFeats <- biosignerRes$biosignature #selected features for each omic.
#' perfBiosigner <- assessPerformance_Biosigner(biosignerModels, data.test)
#' perfBiosigner$svm$rnaRead # perf for SVM for rnaRead data block.
#' perfBiosigner$rf$rnaRead # perf for RF for rnaRead data block.
#' 
#' @importFrom shiny NS tagList 
#' @importFrom biosigner predict
#' @importFrom caret confusionMatrix
#' 
#' @return List of performances for svm and rf models (confusion matrices).
#' 
#' @export
assessPerformance_Biosigner <- function(modelList, dataTest)
{
    ## svm and rf models assessing
    i <- 1
    resList <- list()
    while(i < length(dataTest))
    {
        # models <- obj$classification$biosignerResult$model[[i]]
        name <- names(dataTest)[i]
        models <- modelList[[i]]
        
        # prediction with new model with S and A features 
        preds <- predict(models, dataTest[[i]], tierMaxC = "A")
      
        resList$svm[[name]] <- confusionMatrix(data = table(preds$svm, 
                                                                dataTest$Y))
        resList$rf[[name]] <- confusionMatrix(data = table(preds$randomforest, 
                                                               dataTest$Y))
        # obj$classification$biosignerResult$performance$svm[[i]] <- 
          # confusionMatrix(data = table(preds$svm, dataTest$Y))
        # obj$classification$biosignerResult$performance$rf[[i]] <- 
          # confusionMatrix(data = table(preds$randomforest, dataTest$Y))
        i <- i+1
    }
    return(resList)
}
