## ---- echo=FALSE, results="hide", warning=FALSE-------------------------------
suppressPackageStartupMessages({library('multiSight')})

## ---- eval = FALSE------------------------------------------------------------
#  # To install this package ensure you have BiocManager installed
#  # if (!requireNamespace("BiocManager", quietly = TRUE))
#      # install.packages("BiocManager")
#  
#  # The following initializes usage of Bioc devel
#  # BiocManager::install(version='devel')
#  
#  # BiocManager::install("multiSight")

## ----model_biosigner, echo=TRUE, message=FALSE, warning=FALSE-----------------
library(multiSight)

## omic2 is multi-omic data of 2 data sets included in package
data("omic2", package = "multiSight")

## omic2 is multi-omic data of 2 data sets included in package
splittedData <- splitDatatoTrainTest(omic2, freq = 0.8)
data.train <- splittedData$data.train
data.test <- splittedData$data.test

## Build model and one biosignature by omic data set.
biosignerRes <- runSVMRFmodels_Biosigner(data.train)

## Results
biosignerModels <- biosignerRes$model #list of SVM/RF models for each omic.
biosignerFeats <- biosignerRes$biosignature #selected features for each omic.

## Asses model classification performances
biosignerPerf <- assessPerformance_Biosigner(modelList = biosignerModels, 
                                             dataTest = data.test)
print(biosignerPerf) #confusion matrices and performance metrics

## ----model_diablo, message=FALSE, warning=FALSE-------------------------------
library(multiSight)

## omic2 is multi-omic data of 2 data sets included in package
data("omic2", package = "multiSight")
data("diabloRes", package = "multiSight")

splittedData <- splitDatatoTrainTest(omic2, freq = 0.8)
data.train <- splittedData$data.train
data.test <- splittedData$data.test

## Build model and one biosignature by omic data set.
# diabloRes <- runSPLSDA(data.train)
# diabloRes #internal object of package to save time

## Results
diabloModels <- diabloRes$model #sPLS-DA model using all omics.
diabloFeats <- diabloRes$biosignature #selected features for each omic.

## Asses model classification performances
diabloPerf <- assessPerformance_Diablo(splsdaModel = diabloModels, 
                                          dataTest = data.test)
print(diabloPerf) #confusion matrices and performance metrics

## ----message=FALSE------------------------------------------------------------
library(org.Mm.eg.db)
columns(org.Mm.eg.db)

## ----multiOmicEnrichment_deseq2, message=FALSE, warning=FALSE-----------------
library(multiSight)

## omic2 is multi-omic data of 2 data sets included in package
deseqRes <- runMultiDeseqAnalysis(omicDataList = omic2, 
                                  padjUser = 0.05)
## One Differential Expression Analysis table for each omic data set
# View(deseqRes$DEtable) 
## One feature selected list for each omic according to padjust user threshold
multiOmic_biosignature <- deseqRes$selectedFeatures
# View(multiOmic_biosignature)

## Multi-omic enrichment
### convert features
dbList <- list(Omic1 = "ENSEMBL",
               Omic2 = "ENSEMBL")
convFeat <- convertToEntrezid(multiOmic_biosignature, dbList, "org.Mm.eg.db")

### ORA enrichment analysis
library(org.Mm.eg.db, warn.conflicts = FALSE)
# database <- c("kegg", "wikiPathways", "reactome", "MF", "CC", "BP")
database <- c("reactome")
multiOmicRes <- runMultiEnrichment(omicSignature = convFeat, 
                                   databasesChosen = database, 
                                   organismDb = "org.Mm.eg.db", 
                                   pvAdjust = 0.05, #default value, optional
                                   minGSSize = 5, #default value, optional
                                   maxGSSize = 800, #default value, optional
                                   pvStouffer = 0.1) #default value, optional
reacRes <- multiOmicRes$pathways$reactome
# mfRes <- multiOmicRes$go$MF

## ----multiOmicEnrichment_diablo, message=FALSE, warning=FALSE-----------------
library(multiSight)
data("diabloRes", package = "multiSight")
## omic2 is multi-omic data of 2 data sets included in package
# splittedData <- splitDatatoTrainTest(omic2, 0.8)
# data.train <- splittedData$data.train
# data.test <- splittedData$data.test
# 
# diabloRes <- runSPLSDA(data.train)
# diabloRes #internal object of package to save time
diabloModels <- diabloRes$model #sPLS-DA model using all omics.
diabloFeats <- diabloRes$biosignature #selected features for each omic.

## Multi-omic enrichment
### convert features
names(diabloFeats) #/!\use same names for dbList and omic data sets.
dbList <- list(rnaRead = "ENSEMBL", #feature names origin
               dnaRead = "ENSEMBL")
convFeat <- convertToEntrezid(diabloFeats, 
                              dbList, 
                              "org.Mm.eg.db")

### ORA enrichment analysis for omic feature lists
library(org.Mm.eg.db, warn.conflicts = FALSE)
# database <- c("kegg", "wikiPathways", "reactome", "MF", "CC", "BP")
database <- c("reactome")
multiOmicRes <- runMultiEnrichment(omicSignature = convFeat, 
                                   databasesChosen = database, 
                                   organismDb = "org.Mm.eg.db",
                                   pvAdjust = 0.05, #default value, optional
                                   minGSSize = 5, #default value, optional
                                   maxGSSize = 800, #default value, optional
                                   pvStouffer = 0.1) #default value, optional

## Results
reacRes <- multiOmicRes$pathways$reactome
mfRes <- multiOmicRes$go$MF

## ----networkInference, message=FALSE, warning=FALSE---------------------------
library(multiSight)
data("omic2", package = "multiSight")
data("diabloRes", package = "multiSight")
## omic2 is multi-omic data of 2 data sets included in package
splittedData <- splitDatatoTrainTest(omic2, 0.8)
data.train <- splittedData$data.train
data.test <- splittedData$data.test

## Build diablo models
# diabloRes <- runSPLSDA(data.train)
diabloFeats <- diabloRes$biosignature #selected features for each omic.

## Build biosigner models
biosignerRes <- runSVMRFmodels_Biosigner(data.train)
biosignerFeats <- biosignerRes$biosignature #selected features for each omic.

## Network inference
### Diablo features
concatMat_diablo <- getDataSelectedFeatures(omic2, diabloFeats)
corrRes_diablo <- correlationNetworkInference(concatMat_diablo, 0.85)
pcorRes_diablo <- partialCorrelationNI(concatMat_diablo, 0.4)
miRes_diablo <- mutualInformationNI(concatMat_diablo, 0.2)

### Biosigner features
concatMat_biosigner <- getDataSelectedFeatures(omic2, biosignerFeats)
corrRes_bios <- correlationNetworkInference(concatMat_biosigner, 0.85)
pcorRes_bios <- partialCorrelationNI(concatMat_biosigner, 0.4)
miRes_bios <- mutualInformationNI(concatMat_biosigner, 0.2)

corrRes_diablo$graph

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

