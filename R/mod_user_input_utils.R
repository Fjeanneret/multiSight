#' Remove low variance data function
#'
#' @param omicData Numeric matrix.
#' @param omicID Numeric id to display which omic is processed.
#'  
#' @description To remove features with low variance.
#'
#' @return Numeric matrix with features having variance.
#' 
#' @noRd
#' 
#' @importFrom caret nearZeroVar
zeroVarRmv <- function(omicData, omicID)
{
    message("NOTE ", omicID, ": caret::nearZeroVar processing...")
    featIdx <- caret::nearZeroVar(omicData)
    if (length(featIdx) > 0)
    {
        omicData <- omicData[, -featIdx]
        message("caret::nearZeroVar have found and removed ",
                length(featIdx)," low variance features.")
    }else{
        message(omicID, ": caret::nearZeroVar did not find low variance 
                features.")
    }
    
    return(omicData)
}  
  
#' Check multi-omic data function
#'
#' @param listOmic Omic data sets provided by user in Home tab.
#'  
#' @description To check and select sample names overlap.
#'
#' @return list of all omic data sets with samples overlap.
#' 
#' @noRd
#' 
#' @importFrom purrr map
checkSampleOverlap <- function(listOmic)
{
    out <- tryCatch(
      {
          # checking samples ID -> only all in common
          allID <- lapply(listOmic, rownames)
          ID <- Reduce(intersect, allID) # common IDs between blocks
          sampleOverlap <- length(ID)
          if (sampleOverlap==0)
          {
            msg <- "No one sample is the same between blocks.
            You could retrieve expected format in vignette."
            stop(simpleError(msg))
          }
          matListOverlap <- map(listOmic, function(x) x[ID, ])
          message("NOTE: ", 
                  sampleOverlap, 
                  " samples in common between all blocks.")
          return(matListOverlap)
      },
      error=function(e) {
          message(e$message)
          return(NULL)
      })
}

#' Load multi-omic data function
#'
#' @param input,output,session Internal parameters for shiny
#' @param ids input ids according to user input number choice.
#'  
#' @description To load all data provided by user for further analysis.
#'
#' @importFrom utils read.csv read.csv2
#' 
#' @noRd
#' 
#' @return list of all omic data sets provided by user in Home tab.
loadMultiOmicData <- function(input, output, session, ids)
{
    multiOmicData <- list()
    # Get ids for data inputs
    data_ids <- vapply(seq(1, length(ids)), 
                      function(i){
                          paste("Omic", ids[i], sep="")
                      },
                      FUN.VALUE = character(1))

    ## Get omic data inputs by input's ids i
    for(i in seq(1, length(data_ids)))
    {
        omicID <- data_ids[i]
        sepLabel <- paste0("sepInput_", omicID)
        
        omicData <- read.csv(input[[omicID]]$datapath, 
                            header = TRUE, 
                            row.names = 1, 
                            sep = input[[sepLabel]])
        
        ## nearZeroVar: remove features with very low variance
        omicData <- zeroVarRmv(omicData, omicID)
        
        multiOmicData[[omicID]] <- omicData
    }
    
    ## Get sample classes vector input
    sepClass <- input$sepInput_OmicClass
    Y <- read.csv(input[["OmicClass"]]$datapath, 
                  header = TRUE,
                  row.names = 1,
                  # colClasses = 'factor',
                  # sep = sepClass)$x
                  sep = sepClass)
    multiOmicData$Y <- Y
    
    ## Check if sample names are overlapping 
    multiOmicData <- checkSampleOverlap(multiOmicData)
    return(multiOmicData)
}

#' Organism list function
#'
#' @description To load organism vector list for user choice in Home ui
#'
#' @noRd
organismList <- function()
{
    orgList <- c("Human (org.Hs.eg.db)" = "org.Hs.eg.db",
                 "Mouse (org.Mm.eg.db)" = "org.Mm.eg.db",
                 "Rat (org.Rn.eg.db)" = "org.Rn.eg.db",
                 "Yeast (org.Sc.sgd.db)" = "org.Sc.sgd.db",
                 "Fly (org.Dm.eg.db)" = "org.Dm.eg.db",
                 "Arabidopsis (org.At.tair.db)" = "org.At.tair.db",
                 "Zebrafish (org.Dr.eg.db)" = "org.Dr.eg.db",
                 "Bovine (org.Bt.eg.db)" = "org.Bt.eg.db",
                 "Worm (org.Ce.eg.db)" = "org.Ce.eg.db",
                 "Chicken (org.Gg.eg.db)" = "org.Gg.eg.db",
                 "Canine (org.Cf.eg.db)" = "org.Cf.eg.db",
                 "Pig (org.Ss.eg.db)" = "org.Ss.eg.db",
                 "Rhesus (org.Mmu.eg.db)" = "org.Mmu.eg.db",
                 "E coli strain K12 (org.EcK12.eg.db)" = "org.EcK12.eg.db",
                 "Xenopus (org.Xl.eg.db)" = "org.Xl.eg.db",
                 "Chimp (org.Pt.eg.db)" = "org.Pt.eg.db",
                 "Anopheles (org.Ag.eg.db)" = "org.Ag.eg.db",
                 "Malaria (org.Pf.plasmo.db)" = "org.Pf.plasmo.db",
                 "E coli strain Sakai (org.EcSakai.eg.db)" = 
                    "org.EcSakai.eg.db")
    
    return(orgList)
}

#' User getData utils function
#' 
#' @description To retrieve whole date used for analysis: omic datasets 
#' provided by user, train and test data sets created to build classification
#' models and a concatenated matrices of all selected features used for 
#' network inferences. 
#' 
#' @param dataName Kind of data, by default all generated data
#' \itemize{
#'  \item{"wholeData"}{User's data list}
#'  \item{"dataTrain and dataTest"}{User's data splitted by 
#'  splitDatatoTrainTest() function}
#'  \item{"concatMatrixFeats_diabloNetwork"}{Concatenated matrix of all 
#'  selected features by diablo sPLS-DA method}
#' }
#' @param obj R6 object to wrap all data from different analysis.
#' 
#' @noRd
#' 
#' @return List of data used for analysis. 
getData <- function(obj,
                    dataName = c("wholeData", 
                                 "dataTrain", 
                                 "dataTest", 
                                 "concatMatrixFeats_diabloNetwork"
))
{
    if (exists("obj"))
    {
        return(obj$data[dataName])
      
    }else
    {
        msg <- "Result object 'obj' does not exist"
        stop(msg)
    }
}  

#' User getClassif utils function
#' 
#' @description Provides all results relating to classification models: 
#' models, biosignatures (selected features), performances and feature 
#' details table.  
#' 
#' @param result Kind of result, by default all results.
#' \itemize{
#'  \item{"model"}{Classification models you can use to forecast classes.}
#'  \item{"biosignature"}{Selected features by model, called biosignatures.}
#'  \item{"performance"}{Classification performances, confusion matrix}
#'  \item{"featDetails"}{features tables with computed values
#'  such as class's means or DESEQ2 logF and pvalues if DESEQ2 used}
#'  \item{"design"}{diablo covariance matrix design}
#' }
#' @param obj R6 object to wrap all data from different analysis.
#' 
#' @noRd
#'  
#' @return List of classification results. 
getClassif <- function(obj,
                       method = c("diablo", 
                                  "biosigner"),
                       result = c("model", 
                                  "biosignature",
                                  "performance", 
                                  "featDetails",
                                  "design")
)
{
    if (exists("obj"))
    {
        method <- paste0(method,"Result")
        resList <- lapply(method, function(meth) 
        {
            if (method == "biosigner") 
            {
                result <- result[-5] # no design
            }
            res <- obj$classification[[meth]]
            return(res[[result]])
        })
        return(resList[[1]])
      
    }else
    {
        msg <- "Result object 'obj' does not exist"
        stop(msg)
    }
}  

#' User getEnrichment utils function
#' 
#' @description Provides all results relating to functional enrichment: 
#' features converted to ENTREZ ids, enrichment tables for each omic and 
#' databases selected, multi-omic tables with Stouffer's value and relative
#' enrichment map according user's threshold..
#' 
#' @param sourceFeat enrichment method, source of features used
#' \itemize{
#'  \item{"diablo"}{If features enriched come from diablo selection method}
#'  \item{"deseq"}{If features enriched come from deseq2 analysis and 
#'  threshold}
#'  }
#' @param infoType
#' \itemize{
#'  \item{"pathways"}{Pathways results}
#'  \item{"go"}{Gene Ontology results}
#'  \item{"featureConverted"}{bitr converted feature id to entrezid results}
#'  }
#' @param bioDb
#' \itemize{
#'  \item{"kegg"}
#'  \item{"wikiPathways"}
#'  \item{"reactome"}
#'  \item{"MF"}
#'  \item{"CC"}
#'  \item{"BP"}
#'  }
#' @param obj R6 object to wrap all data from different analysis.
#'   
#' @noRd
#'  
#' @return List of classification results. 
getEnrichment <- function(obj,
                          sourceFeat,
                          infoType,
                          bioDb = NULL)
{
    if (exists("obj"))
    {
        resList <- obj$enrichment[[sourceFeat]][[infoType]]
        if (infoType == "featureConverted")
        {
            resList <- obj$enrichment[[sourceFeat]]$featureConverted
        }else
        {
            if (is.null(bioDb))
            {
                resList <- obj$enrichment[[sourceFeat]][[infoType]]
            }
            else
            {
                resList <- 
                    obj$enrichment[[sourceFeat]][[infoType]][[bioDb]]$result
            }
        }
        
        return(resList)
      
    }else
    {
        msg <- "Result object 'obj' does not exist"
        stop(msg)
    }
}  

#' User getNetwork utils function
#' 
#' @description Provides network graphs  relating to network inference.
#' 
#' @param sourceFeat enrichment method, source of features used
#' \itemize{
#'  \item{"diablo"}{If features come from diablo selection method}
#'  \item{"biosigner"}{If features come from biosigner selection method}
#'  }
#' @param netType
#' \itemize{
#'  \item{"correlation"}{}
#'  \item{"partialCorrelation"}{}
#'  \item{"mutualInformation"}{}
#'  }
#' @param obj R6 object to wrap all data from different analysis.
#'  
#' @noRd
#'  
#' @return List of classification results. 
getNetwork <- function(obj,
                       sourceFeat,
                       netType = NULL)
{
    if (exists("obj"))
    {
        source <- paste0(sourceFeat,"Network")
        if (is.null(netType))
        {
            resList <- obj$networkInference[[source]]$graph
        }else
        {
            resList <- obj$networkInference[[source]]$graph[[netType]]
        }
        
        return(resList)
      
    }else
    {
        msg <- "Result object 'obj' does not exist"
        stop(msg)
    }
}  

#' User getPubmed utils function
#' 
#' @description Provides network graphs  relating to network inference.
#' 
#' @param sourceFeat enrichment method, source of features used
#' \itemize{
#'  \item{"diablo"}{If features come from diablo selection method}
#'  \item{"biosigner"}{If features come from biosigner selection method}
#'  }
#' @param biblioType
#' \itemize{
#'  \item{"request"}{PubMed query giving these articles}
#'  \item{"all"}{All results}
#'  }
#' @param obj R6 object to wrap all data from different analysis.
#'  
#' @noRd
#'  
#' @return List of classification results. 
getPubMed <- function(obj,
                      sourceFeat,
                      biblioType = "request")
{
    if (exists("obj"))
    {
        source <- paste0(sourceFeat,"Network")
        if (biblioType == "request")
        {
            resList <- obj$networkInference[[source]]$biblio$request
        }else if (biblioType == "all")
        {
            resList <- obj$networkInference[[source]]$biblio
        }
    
        return(resList)
  
    }else
    {
        msg <- "Result object 'obj' does not exist"
        stop(msg)
    }
}  