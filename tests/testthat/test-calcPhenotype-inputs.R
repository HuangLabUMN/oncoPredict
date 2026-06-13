make_calcPhenotype_inputs <- function(n_genes=8, n_train=6, n_test=3, n_drugs=2) {
  set.seed(1)
  trainingExprData <- matrix(rnorm(n_genes * n_train), nrow=n_genes)
  rownames(trainingExprData) <- paste0("gene", seq_len(n_genes))
  colnames(trainingExprData) <- paste0("train", seq_len(n_train))

  trainingPtype <- matrix(rnorm(n_train * n_drugs), nrow=n_train)
  rownames(trainingPtype) <- colnames(trainingExprData)
  colnames(trainingPtype) <- paste0("drug", seq_len(n_drugs))

  testExprData <- matrix(rnorm(n_genes * n_test), nrow=n_genes)
  rownames(testExprData) <- rownames(trainingExprData)
  colnames(testExprData) <- paste0("test", seq_len(n_test))

  list(
    trainingExprData=trainingExprData,
    trainingPtype=trainingPtype,
    testExprData=testExprData
  )
}

run_calcPhenotype <- function(inputs, minNumSamples=0, pcr=FALSE, ...) {
  calcPhenotype(
    trainingExprData=inputs$trainingExprData,
    trainingPtype=inputs$trainingPtype,
    testExprData=inputs$testExprData,
    batchCorrect="none",
    powerTransformPhenotype=FALSE,
    removeLowVaryingGenes=0,
    minNumSamples=minNumSamples,
    selection=1,
    printOutput=FALSE,
    pcr=pcr,
    removeLowVaringGenesFrom="homogenizeData",
    folder=FALSE,
    ...
  )
}

test_that("calcPhenotype rejects non-matrix inputs", {
  inputs <- make_calcPhenotype_inputs()

  bad_test <- inputs
  bad_test$testExprData <- as.data.frame(bad_test$testExprData)
  expect_error(run_calcPhenotype(bad_test), '"testExprData" must be a matrix', fixed=TRUE)

  bad_training <- inputs
  bad_training$trainingExprData <- as.data.frame(bad_training$trainingExprData)
  expect_error(run_calcPhenotype(bad_training), '"trainingExprData" must be a matrix', fixed=TRUE)

  bad_ptype <- inputs
  bad_ptype$trainingPtype <- as.data.frame(bad_ptype$trainingPtype)
  expect_error(run_calcPhenotype(bad_ptype), '"trainingPtype" must be a matrix', fixed=TRUE)
})

test_that("calcPhenotype accepts matrix subclasses", {
  inputs <- make_calcPhenotype_inputs()
  class(inputs$testExprData) <- c("custom_matrix", class(inputs$testExprData))

  expect_silent(output <- run_calcPhenotype(inputs))
  expect_true(is.matrix(output))
  expect_equal(dim(output), c(3L, 2L))
})

test_that("calcPhenotype preserves one-column trainingPtype as a matrix", {
  inputs <- make_calcPhenotype_inputs(n_drugs=1)

  expect_silent(output <- run_calcPhenotype(inputs))
  expect_true(is.matrix(output))
  expect_equal(dim(output), c(3L, 1L))
})

test_that("calcPhenotype minNumSamples checks samples, not genes", {
  inputs <- make_calcPhenotype_inputs(n_genes=20, n_train=2, n_test=2)

  expect_error(
    run_calcPhenotype(inputs, minNumSamples=3),
    "There are less than 3 samples",
    fixed=TRUE
  )
})

test_that("homogenizeData qn preserves dimnames", {
  inputs <- make_calcPhenotype_inputs()

  output <- homogenizeData(
    testExprMat=inputs$testExprData,
    trainExprMat=inputs$trainingExprData,
    batchCorrect="qn",
    selection=1,
    printOutput=FALSE
  )

  expect_equal(rownames(output$train), rownames(inputs$trainingExprData))
  expect_equal(rownames(output$test), rownames(inputs$testExprData))
  expect_equal(colnames(output$train), colnames(inputs$trainingExprData))
  expect_equal(colnames(output$test), colnames(inputs$testExprData))
})

test_that("calcPhenotype parallel returns the same predictions as serial", {
  skip_on_os("windows")
  inputs <- make_calcPhenotype_inputs(n_drugs=3)

  serial <- run_calcPhenotype(inputs, parallel=FALSE)
  parallel_output <- run_calcPhenotype(inputs, parallel=TRUE, cores=2)

  expect_equal(parallel_output, serial)
})

test_that("calcPhenotype validates cores", {
  inputs <- make_calcPhenotype_inputs()

  expect_error(
    run_calcPhenotype(inputs, parallel=TRUE, cores=0),
    '"cores" must be a positive integer',
    fixed=TRUE
  )
})

test_that("calcPhenotype pcr rsq removes zero-variance genes from train and test", {
  inputs <- make_calcPhenotype_inputs(n_genes=8, n_train=20, n_test=3, n_drugs=1)
  inputs$trainingExprData[1, ] <- 1
  inputs$testExprData[1, ] <- c(1, 2, 3)

  set.seed(1)
  output <- run_calcPhenotype(inputs, pcr=TRUE, rsq=TRUE, percent=50)

  expect_true(is.list(output))
  expect_true("rsq" %in% names(output))
  expect_equal(dim(output$DrugPredictions), c(3L, 1L))
})

run_predictionAccuracybyCV <- function(inputs, trainingPtype=inputs$trainingPtype[, 1, drop=FALSE], ...) {
  predictionAccuracybyCV(
    trainingExprData=inputs$trainingExprData,
    trainingPtype=trainingPtype,
    testExprData=NULL,
    batchCorrect="none",
    powerTransformPhenotype=FALSE,
    removeLowVaryingGenes=0,
    minNumSamples=0,
    selection=1,
    printOutput=FALSE,
    pcr=FALSE,
    removeLowVaringGenesFrom="homogenizeData",
    folder=FALSE,
    ...
  )
}

test_that("predictionAccuracybyCV supports one-drug matrix inputs", {
  inputs <- make_calcPhenotype_inputs(n_drugs=1)

  output <- run_predictionAccuracybyCV(inputs, cvFold=3)

  expect_type(output$cvPtype, "double")
  expect_equal(length(output$cvPtype), ncol(inputs$trainingExprData))
  expect_named(output$cvPtype, colnames(inputs$trainingExprData))
  expect_type(output$realPtype, "double")
  expect_named(output$realPtype, colnames(inputs$trainingExprData))
})

test_that("predictionAccuracybyCV parallel returns the same predictions as serial", {
  skip_on_os("windows")
  inputs <- make_calcPhenotype_inputs(n_drugs=1)
  set.seed(10)
  serial <- run_predictionAccuracybyCV(inputs, cvFold=-1, parallel=FALSE)

  set.seed(10)
  parallel_output <- run_predictionAccuracybyCV(inputs, cvFold=-1, parallel=TRUE, cores=2)

  expect_equal(parallel_output, serial)
})

test_that("predictionAccuracybyCV validates cores", {
  inputs <- make_calcPhenotype_inputs(n_drugs=1)

  expect_error(
    run_predictionAccuracybyCV(inputs, parallel=TRUE, cores=0),
    '"cores" must be a positive integer',
    fixed=TRUE
  )
})
