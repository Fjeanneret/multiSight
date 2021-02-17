##################################################################
##                         DATA to test                         ##
##################################################################
## normal random value matrices with 5 features correlated to 5 
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
test_that("concatenatedMatrix",
{
    matrixWithSelectedFeatures <- 
      getDataSelectedFeatures(matList, featuresList)
    nfeatures <- sum(vapply(featuresList, length, FUN.VALUE = numeric(1)))
    expect_gt(length(matrixWithSelectedFeatures), 0)
    expect_equal(ncol(matrixWithSelectedFeatures), nfeatures)
})

test_that("Correlation network",
{
    corrBuild <- correlationNetworkInference(randomMatrix, 0)
    
    expect_gt(nrow(corrBuild$flattenTable), 0)
    expect_equal(nrow(corrBuild$num.matrix), ncol(corrBuild$num.matrix))
    expect_gt(nrow(corrBuild$graph$x$links), 0)
    expect_gt(nrow(corrBuild$graph$x$nodes), 0)
})

test_that("Partial Correlation network",
{
    pcorBuild <- partialCorrelationNI(randomMatrix, valueThreshold = 0)
    expect_gt(nrow(pcorBuild$flattenTable), 0)
    expect_equal(nrow(pcorBuild$num.matrix), ncol(pcorBuild$num.matrix))
    expect_gt(nrow(pcorBuild$graph$x$links), 0)
    expect_gt(nrow(pcorBuild$graph$x$nodes), 0)
})

test_that("Mutual Information network",
{
    MIBuild <- mutualInformationNI(randomMatrix, 0)

    expect_gt(nrow(MIBuild$flattenTable), 0)
    expect_equal(nrow(MIBuild$num.matrix), ncol(MIBuild$num.matrix))
    expect_gt(nrow(MIBuild$graph$x$links), 0)
    expect_gt(nrow(MIBuild$graph$x$nodes), 0)
})
