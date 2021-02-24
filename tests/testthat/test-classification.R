##################################################################
##                         DATA to test                         ##
##################################################################
## normal random values matrixes with 5 features correlated to 5 
### omic1
matA <- matrix(rnorm(25),10,10)
colnames(matA) <- LETTERS[seq(1, ncol(matA))]
### omic2
matB <- cbind(matrix(rnorm(50), 10, 10), matA[, seq(1, 5)]*2)
colnames(matB) <- LETTERS[seq(ncol(matA)+1, ncol(matA) + ncol(matB))]

rownames(matA) <- rownames(matB) <- paste0("row_", 
                                           sample(letters)[seq(1, nrow(matA))])

Y_sampleNames <- rownames(matA)
Y_label <- c(rep("A", (nrow(matA)/2)), rep("B", nrow(matA) - (nrow(matB)/2)))
Y <- as.factor(Y_label)
## random multi-blocks data list
matList <- list(matA = matA, matB = matB, Y = Y)
omicdataList <- matList[-3]
## Features List
featuresList <- list(featA = LETTERS[seq(1,10)], featB = LETTERS[seq(18, 25)])

randomMatrix <-
    getDataSelectedFeatures(matList, featuresList)

## Random DESeq table values with only two columns
feats <- featuresList$featA
DESeqSomeColumns <- data.frame(log2FoldChange = rnorm(length(feats), 0, 1.5), 
                               padj =  rbeta(length(feats), 1, 1), 
                               row.names = feats) 
DeseqTables <- list(matA = DESeqSomeColumns, matB = DESeqSomeColumns)

# sPLSDAmodel imported from sysdata.rda 

#################################################################
##                       Functions TESTS                       ##
#################################################################
test_that("splitDataTrainTest",
          {
              freq <- 0.8
              nsample <- sum(vapply(omicdataList, 
                                    nrow, 
                                    FUN.VALUE = numeric(1)))
              dataSplitted <- splitDatatoTrainTest(matList, freq)
              nsample08 <- sum(vapply(dataSplitted$data.train[-3], 
                                      nrow, 
                                      FUN.VALUE = numeric(1)))
              
              expect_equal(nsample * 0.8, nsample08)
          })

test_that("CovarianceDesign_type",
          {
              nOmic <- length(omicdataList)
              designMat <- buildCovarianceDesign(omicdataList)
              expect_type(designMat, "double")
          })

test_that("CovarianceDesign_value",
          {
              nOmic <- length(omicdataList)
              designMat <- buildCovarianceDesign(omicdataList)
              expect_equal(nrow(designMat), nOmic)
              expect_equal(ncol(designMat), nOmic)
              expect_equal(nrow(designMat), ncol(designMat))
          })

test_that("buildFeatTable_class",
          {
              feats <- featuresList$featA
              omicB <- omicdataList$matA
              table <- buildFeatTable(featVec = feats, 
                                      omicBlock =  omicB, 
                                      Y)$x$data
              expect_s3_class(table, "data.frame")
          })

test_that("buildFeatTable_basic",
          {
              feats <- featuresList$featA
              omicB <- omicdataList$matA
              table <- buildFeatTable(featVec = feats, 
                                      omicBlock =  omicB, 
                                      Y)$x$data
              expect_equal(nrow(table), length(feats))
          })

test_that("buildFeatTable_DESeqValues_a",
          {
              feats <- featuresList$featA
              omicB <- omicdataList$matA
              table <- buildFeatTable(featVec = feats, 
                                      omicBlock =  omicB, 
                                      Y,
                                      DeseqTables$matA)$x$data
              expect_equal(nrow(table), length(feats))
              expect_equal(ncol(table), 7)
          })

test_that("buildFeatTable_DESeqValues_b",
          {
              feats <- getSelectedFeatures(sPLSDAmodel)
              total <- sum(vapply(feats, length, FUN.VALUE = numeric(1)))
              omic1Feats <- length(feats$Omic1)
              omic2Feats <- length(feats$Omic2)
              expect_equal(total, 60)
              expect_equal(omic1Feats, 30)
              expect_equal(omic2Feats, 30)
          })
