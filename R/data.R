#' Multi-omic data with 2 omics.
#'
#' A dataset containing simulated multi-omic data with 30 samples.
#'
#' @format A list of 2 dataframes and 1 factor vector 30 rows:
#' \describe{
#'   \item{rnaRead}{transcriptomic simulated data}
#'   \item{dnaRead}{genomic simulated data}
#'   \item{Y}{30 samples' classes}
#' }
#' @source MOSim package used to simulate omic data.
"omic2"

#' Diablo results 
#'
#' Diablo results list:
#' \describe{
#'   \item{design}{Covariance matrix design to maximize}
#'   \item{model}{sPLS-DA model}
#'   \item{biosignature}{Selected features of omic data sets}
#' }
#'
#' @source Returned by runSPLSDA() function.
"diabloRes"

#' Biosigner results 
#'
#' Biosigner results list:
#' \describe{
#'   \item{model}{SVM and RF models}
#'   \item{biosignature}{Selected features of omic data sets}
#' }
#'
#' @source Returned by runSVMRFmodels_Biosigner() function.
"biosignerRes"

#' multiSight results
#'
#' enrichRes object obtained with several omic data sets.
#'
#' @source Returned by multiSight runMultiEnrichment() function
#' in enrichTables$pathways$reactome$enrichObj slot
"enrichResList"

#' multiSight results
#'
#' multiOmicEnrichment results object obtained with several omic data sets.
#'
#' @source Returned by multiSight runMultiEnrichment() function
"runMultiEnrichment_result"

#' multiSight results
#'
#' DESeq2 results object obtained with several omic data sets.
#'
#' @source Returned by multiSight runMultiDeseqAnalysis() function
"deseqRes"

