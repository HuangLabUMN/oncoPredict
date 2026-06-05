#' Remove genes with low variation.
#'
#' This function performs variable selection by removing genes with the lowest variance in the datasets.
#'
#' @param exprMat A matrix of gene expression levels. rownames() are genes, and colnames() are samples.
#' @param removeLowVaryingGenes The proportion of low varying genes to be removed.The default is .2
#' @return A vector of row/genes to keep.
#' @export
doVariableSelection <- function(exprMat, removeLowVaryingGenes=.2)
{
  vars <- apply(exprMat, 1, var)
  return(order(vars, decreasing=TRUE)[seq(1:as.integer(nrow(exprMat)*(1-removeLowVaryingGenes)))])
}
#' Homogenizes two expression matrices
#'
#'This function takes two gene expression matrices (like trainExprMat and testExprMat) and returns homogenized versions of the matrices by employing the homogenization method specified.
#'By default, the Combat method from the sva library is used.
#'In both matrices, genes are row names and samples are column names.
#'It will deal with duplicated gene names, as it subsets and orders the matrices correctly.
#'@param testExprMat A gene expression matrix for samples on which we wish to predict a phenotype.Genes are rows, samples are columns.
#'@param trainExprMat A gene expression matrix for samples for which the phenotype is already known.Genes are rows, samples are columns.
#'@param batchCorrect The type of batch correction to be used. Options are 'eb' for ComBat, 'qn' for quantile normalization, 'standardize' for within-dataset z-score standardization, 'rank', 'rank_then_eb', or 'none'. The 'standardize' option can be useful when using microarray training data to build models on RNA-seq testing data.
#'#The default is 'eb'.
#'@param selection This parameter can be used to specify how duplicates are handled. The default value of -1 means to ask the user.
#'#Other options include '1' to summarize duplicates by their mean, and '2'to discard all duplicated genes.
#'@param printOutput To suppress output, set to false. Default is TRUE.
#'@import sva
#'@importFrom limma normalizeQuantiles
#'@keywords Homogenize gene expression data.
#'@return A list containing two entries $train and $test, which are the homogenized input matrices.
#'@export
homogenizeData<-function (testExprMat, trainExprMat, batchCorrect = "eb", selection = -1, printOutput = TRUE)
{
  #Check the batchCorrect parameter
  if (!(batchCorrect %in% c("eb", "qn", "none",
                            "rank", "rank_then_eb", "standardize")))
    stop("\"batchCorrect\" must be one of \"eb\", \"qn\", \"rank\", \"rank_then_eb\", \"standardize\" or \"none\"")

  #Check if both row and column names have been specified.
  if (is.null(rownames(trainExprMat)) || is.null(rownames(testExprMat))) {
    stop("ERROR: Gene identifiers must be specified as \"rownames()\" on both training and test expression matrices. Both matices must have the same type of gene identifiers.")
  }

  #Check that some of the row names overlap between both datasets (print an error if none overlap)
  if (sum(rownames(trainExprMat) %in% rownames(testExprMat)) == 0) {
    stop("ERROR: The rownames() of the supplied expression matrices do not match. Note that these are case-sensitive.")
  }
  else {
    if (printOutput)
      message(paste("\n", sum(rownames(trainExprMat) %in% rownames(testExprMat)), " gene identifiers overlap between the supplied expression matrices... \n", paste = ""))
  }

  #If there are duplicate gene names, give the option of removing them or summarizing them by their mean.
  if ((sum(duplicated(rownames(trainExprMat))) > 0) || sum(sum(duplicated(rownames(testExprMat))) > 0)) {
    if (selection == -1) {
      message("\nExpression matrix contain duplicated gene identifiers (i.e. duplicate rownames()), how would you like to proceed:")
      message("\n1. Summarize duplicated gene ids by their mean value (acceptable in most cases)")
      message("\n2. Disguard all duplicated genes (recommended if unsure)")
      message("\n3. Abort (if you want to deal with duplicate genes ids manually)\n")
    }
    while (is.na(selection) | selection <= 0 | selection > 3) {
      selection <- readline("Selection: ")
      selection <- ifelse(grepl("[^1-3.]", selection), -1, as.numeric(selection))
    }

    message("\n")

    if (selection == 1) #Summarize duplicates by their mean.
    {
      if ((sum(duplicated(rownames(trainExprMat))) > 0)) {
        trainExprMat <- summarizeGenesByMean(trainExprMat)
      }
      if ((sum(duplicated(rownames(testExprMat))) > 0)) {
        testExprMat <- summarizeGenesByMean(testExprMat)
      }
    }
    else if (selection == 2) #Disguard all duplicated genes.
    {
      if ((sum(duplicated(rownames(trainExprMat))) > 0)) {
        keepGenes <- names(which(table(rownames(trainExprMat)) == 1))
        trainExprMat <- trainExprMat[keepGenes, ]
      }
      if ((sum(duplicated(rownames(testExprMat))) > 0)) {
        keepGenes <- names(which(table(rownames(testExprMat)) == 1))
        testExprMat <- testExprMat[keepGenes, ]
      }
    }
    else {
      stop("Exectution Aborted!")
    }
  }

  #Subset and order gene ids on the expression matrices.
  commonGenesIds <- rownames(trainExprMat)[rownames(trainExprMat) %in%
                                             rownames(testExprMat)]
  trainExprMat <- trainExprMat[commonGenesIds, ]
  testExprMat <- testExprMat[commonGenesIds, ]

  #Subset and order the 2 expression matrices.
  if (batchCorrect == "eb") {
    #Subset to common genes and batch correct using ComBat.
    dataMat <- cbind(trainExprMat, testExprMat)
    mod <- data.frame(`(Intercept)` = rep(1, ncol(dataMat)))
    rownames(mod) <- colnames(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)),
                              rep("test", ncol(testExprMat))))

    # Added
    # Filter out genes with low variances to make sure comBat run correctly
    dataMat <- cbind(trainExprMat, testExprMat)
    gene_vars = apply(dataMat, 1, var)
    genes<-as.vector(gene_vars)

    if (length(which(genes <= 1e-3) != 0)){ #If some genes have low variances (if the variance is not 0), remove them.
      dataMat = dataMat[-(which(genes <= 1e-3)),]
    }

    # End added

    combatout <- ComBat(dataMat, whichbatch, mod = mod)
    return(list(train = combatout[, whichbatch == "train"],
                test = combatout[, whichbatch == "test"], selection = selection))
  }
  else if (batchCorrect == "standardize") #Standardize to mean 0 and variance 1 in each dataset using a non EB based approach.
  {
    for (i in 1:nrow(trainExprMat)) {
      row <- trainExprMat[i, ]
      trainExprMat[i, ] <- ((row - mean(row))/sd(row))
    }
    for (i in 1:nrow(testExprMat)) {
      row <- testExprMat[i, ]
      testExprMat[i, ] <- ((row - mean(row))/sd(row))
    }
    return(list(train = trainExprMat, test = testExprMat,
                selection = selection))
  }
  else if (batchCorrect == "rank") #The random-rank transform approach, that may be better when applying models to RNA-seq data.
  {
    for (i in 1:nrow(trainExprMat)) {
      trainExprMat[i, ] <- rank(trainExprMat[i, ], ties.method = "random")
    }
    for (i in 1:nrow(testExprMat)) {
      testExprMat[i, ] <- rank(testExprMat[i, ], ties.method = "random")
    }
    return(list(train = trainExprMat, test = testExprMat,
                selection = selection))
  }
  else if (batchCorrect == "rank_then_eb") #Rank-transform the RNAseq data, then apply ComBat
  {
    #First, rank transform the RNA-seq data.
    for (i in 1:nrow(testExprMat)) {
      testExprMat[i, ] <- rank(testExprMat[i, ], ties.method = "random")
    }
    #Subset to common genes and batch correct using ComBat.
    dataMat <- cbind(trainExprMat, testExprMat)
    mod <- data.frame(`(Intercept)` = rep(1, ncol(dataMat)))
    rownames(mod) <- colnames(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)),
                              rep("test", ncol(testExprMat))))
    combatout <- ComBat(dataMat, whichbatch, mod = mod)
    return(list(train = combatout[, whichbatch == "train"],
                test = combatout[, whichbatch == "test"], selection = selection))
  }
  else if (batchCorrect == "qn")
  {
    dataMat <- cbind(trainExprMat, testExprMat)
    dataMatNorm <- normalizeQuantiles(dataMat)
    dimnames(dataMatNorm) <- dimnames(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)),
                              rep("test", ncol(testExprMat))))
    return(list(train = dataMatNorm[, whichbatch == "train"],
                test = dataMatNorm[, whichbatch == "test"],
                selection = selection))
  }
  else {
    return(list(train = trainExprMat, test = testExprMat,
                selection = selection))
  }
}
#'Average over duplicate gene values
#'
#'This function takes a gene expression matrix and if duplicate genes are measured, summarizes them by their means.
#'@param exprMat A gene expression matrix with genes as rownames() and samples as colnames().
#'@return A gene expression matrix that does not contain duplicate genes.
#'@keywords Summarize duplicate genes by their mean.
#'@export
summarizeGenesByMean <- function(exprMat)
{
  geneIds <- rownames(exprMat)
  t <- table(geneIds) #How many times is each gene name duplicated.
  allNumDups <- unique(t)
  allNumDups <- allNumDups[-which(allNumDups == 1)]

  #Create a *new* gene expression matrix with everything in the correct order....
  #Starrt by just adding stuff that isn't duplicated
  exprMatUnique <- exprMat[which(geneIds %in% names(t[t == 1])), ]
  gnamesUnique <- geneIds[which(geneIds %in% names(t[t == 1]))]

  #Add all the duplicated genes to the bottom of "exprMatUniqueHuman", summarizing as you go
  for(numDups in allNumDups)
  {
    geneList <- names(which(t == numDups))

    for(i in 1:length(geneList))
    {
      exprMatUnique <- rbind(exprMatUnique, colMeans(exprMat[which(geneIds == geneList[i]), ]))
      gnamesUnique <- c(gnamesUnique, geneList[i])
      # print(i)
    }
  }

  if(is.numeric(exprMatUnique))
  {
    exprMatUnique <- matrix(exprMatUnique, ncol=1)
  }

  rownames(exprMatUnique) <- gnamesUnique
  return(exprMatUnique)
}
#'Generate predicted drug sensitivity scores
#'
#'This function predicts a phenotype (drug sensitivity score) when provided with microarray or bulk RNAseq gene expression data of different platforms.
#'The imputations are performed using ridge regression, training on a gene expression matrix where phenotype is already known.
#'This function integrates training and testing datasets via a user-defined procedure, and power transforming the known phenotype.
#'@param trainingExprData The training data. A matrix of expression levels. rownames() are genes, colnames() are samples (cell line names or cosmic ides, etc.). rownames() must be specified and must contain the same type of gene ids as "testExprData"
#'@param trainingPtype The known phenotype for "trainingExprData". This data must be a matrix with training samples as rows and drugs or phenotypes as columns. This matrix can contain NA values, that is ok (they are removed in the calcPhenotype() function).
#'@param testExprData The test data where the phenotype will be estimated. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
#'@param batchCorrect How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantile normalization, "standardize" for within-dataset z-score standardization, "rank", "rank_then_eb", or "none" for no homogenization.
#'@param powerTransformPhenotype Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
#'@param removeLowVaryingGenes What proportion of low varying genes should be removed? 20 percent be default
#'@param minNumSamples How many training and test samples are required. Print an error if below this threshold
#'@param selection How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
#'@param printOutput Set to FALSE to supress output.
#'@param pcr Indicates whether or not you'd like to use pcr for feature (gene) reduction. Options are 'TRUE' and 'FALSE'. If you indicate 'report_pc=TRUE' you need to also indicate 'pcr=TRUE'
#'@param removeLowVaringGenesFrom Determine method to remove low varying genes. Options are 'homogenizeData' and 'rawData'.
#'@param report_pc Indicates whether you want to output the training principal components. Options are 'TRUE' and 'FALSE'.
#'@param cc Indicate if you want correlation coefficients for biomarker discovery.
#'@param percent Indicate percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is 80 for 80%
#'These are the correlations between a given gene of interest across all samples vs. a given drug response across samples.
#'These correlations can be ranked to obtain a ranked correlation to determine highly correlated drug-gene associations.
#'@param rsq Indicate whether or not you want to output the R^2 values for the data you train on from true and predicted values.
#'These values represent the percentage in which the optimal model accounts for the variance in the training data.
#'Options are 'TRUE' and 'FALSE'.
#'@param folder If TRUE, write calcPhenotype outputs to calcPhenotype_Output in the current working directory. The default is FALSE.
#'@param parallel If TRUE, fit drug models in parallel after the shared homogenization and gene-filtering steps are complete. The default is FALSE.
#'@param cores The number of cores to use when parallel is TRUE. Parallel execution uses forked processes via parallel::mclapply, which is not available for multicore execution on Windows PCs; on Windows, calcPhenotype will warn and run serially.
#'@return A matrix of predicted drug response values. If rsq, cc, or report_pc is TRUE, returns a list containing the predictions and requested optional outputs. If folder is TRUE, the same object is returned invisibly after files are written.
#'@import sva
#'@import ridge
#'@import car
#'@import utils
#'@import stats
#'@import parallel
#'@importFrom pls pcr explvar
#'@keywords predict drug sensitivity and phenotype
#'@export
calcPhenotype<-function (trainingExprData,
                         trainingPtype,
                         testExprData,
                         batchCorrect,
                         powerTransformPhenotype=TRUE,
                         removeLowVaryingGenes=0.2,
                         minNumSamples,
                         selection=1,
                         printOutput,
                         pcr=FALSE,
                         removeLowVaringGenesFrom,
                         report_pc=FALSE,
                         cc=FALSE,
                         percent=80,
                         rsq=FALSE,
                         folder = FALSE,
                         parallel = FALSE,
                         cores = 1)
{

  #Initiate empty lists for each data type you'd like to collect.
  #_______________________________________________________________
  DrugPredictions<-list() #Collects drug predictions.
  rsqs<-list() #Collects R^2 values.
  cors<-list() #Collects correlation coefficient for each gene across all samples vs. each drug across all samples.
  pvalues<-list() #Collects p-values for the correlation coefficients.
  pcs_list<-list() #Collects principal components when report_pc is TRUE.

  #vs=c()

  drugs<-colnames(trainingPtype) #Store all the possible drugs in a vector.

  #Check the supplied data and parameters.
  #_______________________________________________________________
  if (!is.matrix(testExprData))
    stop("\nERROR: \"testExprData\" must be a matrix.", call. = FALSE)
  if (!is.matrix(trainingExprData))
    stop("\nERROR: \"trainingExprData\" must be a matrix.", call. = FALSE)
  if (!is.matrix(trainingPtype))
    stop("\nERROR: \"trainingPtype\" must be a matrix.", call. = FALSE)

  if (report_pc)
    if (pcr == FALSE)
      stop("\nERROR: pcr must be TRUE if report_pc is TRUE")

  if (pcr)
    if (cc)
      stop("\nERROR: pcr must be FALSE if cc is TRUE")

  #Make sure training samples are equivalent in both matrices.
  if (!any(colnames(trainingExprData) %in% rownames(trainingPtype)))
    stop("\nERROR: No Cell Lines Found in Common: Sample names must be consistent in training matrices")

  #Subset and order the training Expr and trainingPtype to the cell lines in common (and order them)

  commonCellLines<-colnames(trainingExprData)[colnames(trainingExprData) %in% rownames(trainingPtype)]

  trainingExprData <- trainingExprData[, commonCellLines, drop=FALSE]
  trainingPtype <- trainingPtype[commonCellLines, , drop=FALSE]

  #Check if an adequate number of training and test samples have been supplied.
  #_______________________________________________________________
  if ((nrow(trainingExprData) < minNumSamples) || (nrow(testExprData) < minNumSamples)) {
    stop(paste("\nThere are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }

  #Get the homogenized data.
  #_______________________________________________________________
  homData <- homogenizeData(testExprMat=testExprData, trainExprMat=trainingExprData, batchCorrect, selection, printOutput)

  #Remove low varying genes.
  #_______________________________________________________________
  #Do variable selection if specified. By default, we remove 20% of least varying genes from the homogenized dataset.
  #We can also remove the intersection of the lowest 20% from both training and test sets (for the gene ids remaining in the homogenized data).
  #Otherwise, keep all genes.

  #Check batchCorrect parameter.
  if (!(removeLowVaringGenesFrom %in% c("homogenizeData", "rawData"))) {
    stop("\nremoveLowVaringGenesFrom\" must be one of \"homogenizeData\", \"rawData\"")
  }

  keepRows <- seq(1:nrow(homData$train)) #By default we will keep all the genes.
  if (removeLowVaryingGenes > 0 && removeLowVaryingGenes < 1) { #If the proportion of variability to keep is between 0 and 1.
    if (removeLowVaringGenesFrom == "homogenizeData") { #If you're filtering based on homogenized data.
      keepRows <- doVariableSelection(cbind(homData$test, homData$train), removeLowVaryingGenes = removeLowVaryingGenes)

      numberGenesRemoved <- nrow(homData$test) - length(keepRows)
      if (printOutput) message(paste("\n", numberGenesRemoved, "low variabilty genes filtered."));
    }
    else if (removeLowVaringGenesFrom == "rawData") { #If we are filtering based on the raw data i.e. the intersection of the things filtered from both datasets.
      evaluabeGenes <- rownames(homData$test)
      keepRowsTrain <- doVariableSelection(trainingExprData[evaluabeGenes,], removeLowVaryingGenes = removeLowVaryingGenes)
      keepRowsTest <- doVariableSelection(testExprData[evaluabeGenes,], removeLowVaryingGenes = removeLowVaryingGenes)
      keepRows <- intersect(keepRowsTrain, keepRowsTest)
      numberGenesRemoved <- nrow(homData$test) - length(keepRows)
      if (printOutput)
        message(paste("\n", numberGenesRemoved, "low variabilty genes filtered."));
    }
  }

  cores <- as.integer(cores)
  if (is.na(cores) || cores < 1) {
    stop("\nERROR: \"cores\" must be a positive integer.", call. = FALSE)
  }
  useParallel <- isTRUE(parallel) && cores > 1
  if (useParallel && .Platform$OS.type == "windows") {
    warning("parallel=TRUE uses forked processes via parallel::mclapply, which is not available for multicore execution on Windows. Running calcPhenotype serially.", call. = FALSE)
    useParallel <- FALSE
  }

  fitDrug <- function(a) {
    drug <- drugs[a]
    trainingPtype2<-trainingPtype[, a, drop=FALSE] #Obtain the response data for the compound of interest.
    NonNAindex <- which(!is.na(trainingPtype2)) #Get the indices of the non NAs.
    samps<-rownames(trainingPtype)[NonNAindex] #Obtain cell lines you have expression and response data for.

    if (length(samps) <= 1){
      return(list(drug=drug, skipped=TRUE,
                  message=paste("\n", drug, "is skipped due to insufficient cell lines to fit the model.")))
    }

    trainingPtype4<-as.numeric(trainingPtype2[NonNAindex])
    offset = 0
    transForm <- 1
    rsq_value <- NA
    pcs <- NULL
    cors_vec <- NULL
    cors_vec2 <- NULL

    if (powerTransformPhenotype){
      if (min(trainingPtype4) < 0){
        offset <- -min(trainingPtype4) + 1
        trainingPtype4 <- trainingPtype4 + offset
      }
      transForm <- powerTransform(trainingPtype4)[[6]]
      trainingPtype4 <- trainingPtype4^transForm
    }

    if (pcr){
      train_x<-(t(homData$train)[samps, keepRows, drop=FALSE])
      train_y<-trainingPtype4
      test_x<-(t(homData$test)[, keepRows, drop=FALSE])

      x<-as.vector(colSums(train_x))
      bad<-which(x == 0)
      if(length(bad) != 0){
        train_x<-train_x[, -bad, drop=FALSE]
        test_x<-data.frame(test_x[, -bad, drop=FALSE])
      }

      variance<-c()
      for(i in 1:ncol(train_x)){
        variance[i]<-var(as.vector(train_x[,i]))
      }
      bi<-which(variance %in% 0)
      if(length(bi) != 0){
        train_x<-train_x[, -bi, drop=FALSE]
        test_x<-data.frame(test_x[, -bi, drop=FALSE])
      }

      trainFrame<-try(data.frame(Resp=train_y, train_x), silent = TRUE)
      if (dim(trainFrame)[1] == 1){
        return(list(drug=drug, skipped=TRUE,
                    message=paste("\n", drug, "is skipped due to insufficient cell lines to fit the model.")))
      }

      pcr_model<-pcr(Resp~., data=trainFrame, validation='CV')
      v=cumsum(explvar(pcr_model))
      ncomp=min(which(v > percent))

      if(printOutput && !useParallel) message("\nCalculating predicted phenotype using pcr...")
      preds<-predict(pcr_model, newdata=test_x, ncomp=ncomp)

      if (rsq){
        if (dim(train_x)[1] < 4){
          if(printOutput && !useParallel) message(paste("\n", drug, 'is skipped for R^2 analysis'))
        }else{
          data<-(cbind(train_x, train_y))
          dt<-sort(sample(nrow(data), nrow(data)*.7))
          train_x<-data[dt,, drop=FALSE]
          ncol<-dim(train_x)[2]
          train_y<-train_x[,ncol]
          train_x<-train_x[,-ncol, drop=FALSE]

          test_x<-data[-dt,, drop=FALSE]
          ncol<-dim(test_x)[2]
          test_y<-test_x[,ncol]
          test_x<-test_x[,-ncol, drop=FALSE]

          x<-as.vector(colSums(train_x))
          bad<-which(x == 0)
          if (length(bad) != 0){
            train_x<-train_x[, -bad, drop=FALSE]
            test_x<-data.frame(test_x[, -bad, drop=FALSE])
          }

          variance<-c()
          for(i in 1:ncol(train_x)){
            variance[i]<-var(as.vector(train_x[,i]))
          }
          bi<-which(variance %in% 0)
          if (length(bi) != 0){
            train_x<-train_x[, -bi, drop=FALSE]
          }

          data<-data.frame(Resp=train_y, train_x)
          pcr_model<-pcr(Resp~., data=data, validation='CV')
          v=cumsum(explvar(pcr_model))
          ncomp=min(which(v > percent))
          pcr_pred<-predict(pcr_model, test_x, ncomp=ncomp)

          if (printOutput && !useParallel) message("\nCalculating R^2...")
          sst<-sum((test_y - mean(test_y))^2)
          sse<-sum((pcr_pred - test_y)^2)
          rsq_value<-1 - sse/sst
        }
      }

      if (report_pc){
        if (printOutput && !useParallel) message("\nObtaining principal components...")
        pcs<-coef(pcr_model, comps = ncomp)
      }

    } else {

      if(printOutput && !useParallel) message("\nFitting Ridge Regression model...");
      expression<-(t(homData$train)[samps, keepRows, drop=FALSE])
      test_x<-(t(homData$test)[, keepRows, drop=FALSE])

      x<-as.vector(colSums(expression))
      bad<-which(x == 0)
      if (length(bad) != 0){
        expression<-expression[, -bad, drop=FALSE]
        test_x<-data.frame(test_x[, -bad, drop=FALSE])
      }

      variance<-c()
      for (i in 1:ncol(expression)){
        variance[i]<-var(as.vector(expression[,i]))
      }
      bi<-which(variance %in% 0)
      if (length(bi) != 0){
        expression<-expression[, -bi, drop=FALSE]
        test_x<-data.frame(test_x[, -bi, drop=FALSE])
      }

      trainFrame<-try(data.frame(Resp=trainingPtype4, expression), silent = TRUE)
      if (dim(trainFrame)[1] == 1){
        return(list(drug=drug, skipped=TRUE,
                    message=paste("\n", drug, "is skipped due to insufficient cell lines to fit the model.")))
      }

      if(printOutput && !useParallel) message("\nCalculating predicted phenotype...")
      rrModel<-linearRidge(Resp ~., data=trainFrame)
      preds<-predict(rrModel, newdata=data.frame(test_x))

      if(rsq){
        if (dim(expression)[1] < 4){
          if(printOutput && !useParallel) message(paste("\n", drug, 'is skipped for R^2 analysis'))
        } else {
          expression<-(cbind(expression, trainingPtype4))
          dt<-sort(sample(nrow(expression), nrow(expression)*.7))

          train_x<-expression[dt,, drop=FALSE]
          ncol<-dim(train_x)[2]
          train_y<-train_x[,ncol]
          train_x<-train_x[,-ncol, drop=FALSE]

          test_x<-expression[-dt,, drop=FALSE]
          ncol<-dim(test_x)[2]
          test_y<-test_x[,ncol]
          test_x<-test_x[,-ncol, drop=FALSE]

          x<-as.vector(colSums(train_x))
          bad<-which(x == 0)
          if (length(bad) != 0){
            train_x<-train_x[, -bad, drop=FALSE]
            test_x<-data.frame(test_x[, -bad, drop=FALSE])
          }

          variance<-c()
          for (i in 1:ncol(train_x)){
            variance[i]<-var(as.vector(train_x[,i]))
          }
          bi<-which(variance %in% 0)
          if (length(bi) != 0){
            train_x<-train_x[, -bi, drop=FALSE]
            test_x<-data.frame(test_x[, -bi, drop=FALSE])
          }

          trainFrame<-data.frame(Resp=train_y, train_x)
          rrModel<-linearRidge(Resp ~., data=trainFrame)
          testFrame<-data.frame(test_x)
          pred<-predict(rrModel, newdata=testFrame)

          if(printOutput && !useParallel) message("\nCalculating R^2...")
          sst<-sum((test_y - mean(test_y))^2)
          sse<-sum((pred - test_y)^2)
          rsq_value<-1 - sse/sst
        }
      }
    }

    if(powerTransformPhenotype) {
      preds <- preds^(1/transForm)
      preds <- preds - offset
    }

    if(cc){
      if(pcr){
        stop('ERROR: pcr must equal FALSE in order to compute correlations')
      }

      if(printOutput && !useParallel) message("\nCalculating correlation coefficients...")
      cors_vec<-c()
      cors_vec2<-c()
      matrix<-homData$test[keepRows,, drop=FALSE]
      for(d in 1:nrow(matrix)){
        cors_vec[d]<-cor.test(as.vector(matrix[d,]), as.vector(preds))$estimate
        cors_vec2[d]<-cor.test(as.vector(matrix[d,]), as.vector(preds))$p.value
      }
    }

    if(printOutput && !useParallel) message(paste("\nDone making prediction for drug", a, "of", ncol(trainingPtype)))

    list(drug=drug, skipped=FALSE, preds=preds, rsq=rsq_value, cors=cors_vec, pvalues=cors_vec2, pcs=pcs)
  }

  if (useParallel && printOutput) {
    message(paste("\nFitting", length(drugs), "drug models using", cores, "cores..."))
  }
  drugResults <- if (useParallel) {
    parallel::mclapply(seq_along(drugs), fitDrug, mc.cores=cores)
  } else {
    lapply(seq_along(drugs), fitDrug)
  }

  skipped <- vapply(drugResults, function(result) isTRUE(result$skipped), logical(1))
  if (any(skipped) && printOutput) {
    for (result in drugResults[skipped]) {
      message(result$message)
    }
  }
  drugResults <- drugResults[!skipped]
  if (!length(drugResults)) {
    stop("\nERROR: No drugs had sufficient cell lines to fit a model.", call. = FALSE)
  }

  drugs <- vapply(drugResults, function(result) result$drug, character(1))
  DrugPredictions<-lapply(drugResults, function(result) result$preds)
  names(DrugPredictions)<-drugs
  DrugPredictions_mat<-do.call(cbind, DrugPredictions)
  colnames(DrugPredictions_mat)<-drugs
  rownames(DrugPredictions_mat)<-colnames(testExprData)

  if(rsq){
    rsqs<-lapply(drugResults, function(result) result$rsq)
  }
  if(cc){
    cors<-lapply(drugResults, function(result) result$cors)
    pvalues<-lapply(drugResults, function(result) result$pvalues)
  }
  if(report_pc){
    pcs_list<-lapply(drugResults, function(result) result$pcs)
    names(pcs_list)<-drugs
  }

  output <- DrugPredictions_mat
  if(rsq){
    names(rsqs)<-drugs
    rsqs_mat<-do.call(cbind, rsqs)
  }
  if(cc){
    names(cors)<-drugs
    cor_mat<-do.call(cbind, cors)
    rownames(cor_mat)<-rownames(homData$test[keepRows,, drop=FALSE])
    colnames(cor_mat)<-drugs

    names(pvalues)<-drugs
    p_mat<-do.call(cbind, pvalues)
    rownames(p_mat)<-rownames(homData$test[keepRows,, drop=FALSE])
    colnames(p_mat)<-drugs
  }
  if(rsq || cc || report_pc){
    output <- list(DrugPredictions=DrugPredictions_mat)
    if(rsq){
      output$rsq <- rsqs_mat
    }
    if(cc){
      output$cors <- cor_mat
      output$pvalues <- p_mat
    }
    if(report_pc){
      output$pcs <- pcs_list
    }
  }

  if(folder){
  dir.create("./calcPhenotype_Output", showWarnings=FALSE)
  write.csv(DrugPredictions_mat, file="./calcPhenotype_Output/DrugPredictions.csv", row.names = TRUE, col.names = TRUE)

  #If rsq=TRUE, save R^2 data.
    if(rsq){
      dir.create("./calcPhenotype_Output", showWarnings=FALSE)
      write.table(rsqs_mat, file="./calcPhenotype_Output/R^2.txt")
    }

  #If CC=TRUE, save correlation coefficient data.
    if(cc){
      dir.create("./calcPhenotype_Output", showWarnings=FALSE)
      write.table(cor_mat, file="./calcPhenotype_Output/cors.txt")

      dir.create("./calcPhenotype_Output", showWarnings=FALSE)
      write.table(p_mat, file="./calcPhenotype_Output/pvalues.txt")
    }
    if(report_pc){
      dir.create("./calcPhenotype_Output", showWarnings=FALSE)
      for(drug in names(pcs_list)){
        pcs <- pcs_list[[drug]]
        path<-paste('./calcPhenotype_Output/', drug, '.RData', sep="")
        save(pcs, file=path)
      }
    }

    return(invisible(output))
  }

  return(output)

  #print(vs)
}

#'Calculate Cross-Validation Scores using OncoPredict
#'
#'This function predicts a phenotype (drug sensitivity score) when provided with microarray or bulk RNAseq gene expression data of different platforms.
#'The imputations are performed using ridge regression, training on a gene expression matrix where phenotype is already known.
#'This function integrates training and testing datasets via a user-defined procedure, and power transforming the known phenotype.
#'@param trainingExprData The training data. A matrix of expression levels. rownames() are genes, colnames() are samples (cell line names or cosmic ides, etc.). rownames() must be specified and must contain the same type of gene ids as "testExprData"
#'@param trainingPtype The known phenotype for "trainingExprData". This can be a one-drug vector with one value per training sample or a matrix with training samples as rows and drugs or phenotypes as columns. This matrix can contain NA values, that is ok (they are removed in the calcPhenotype() function).
#'@param testExprData The test data where the phenotype will be estimated. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
#'@param batchCorrect How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantile normalization, "standardize" for within-dataset z-score standardization, "rank", "rank_then_eb", or "none" for no homogenization.
#'@param powerTransformPhenotype Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
#'@param removeLowVaryingGenes What proportion of low varying genes should be removed? 20 percent be default
#'@param minNumSamples How many training and test samples are required. Print an error if below this threshold
#'@param selection How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
#'@param printOutput Set to FALSE to supress output.
#'@param pcr Indicates whether or not you'd like to use pcr for feature (gene) reduction. Options are 'TRUE' and 'FALSE'. If you indicate 'report_pc=TRUE' you need to also indicate 'pcr=TRUE'
#'@param removeLowVaringGenesFrom Determine method to remove low varying genes. Options are 'homogenizeData' and 'rawData'.
#'@param report_pc Indicates whether you want to output the training principal components. Options are 'TRUE' and 'FALSE'.
#'@param cc Indicate if you want correlation coefficients for biomarker discovery.
#'@param percent Indicate percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is 80 for 80%
#'These are the correlations between a given gene of interest across all samples vs. a given drug response across samples.
#'These correlations can be ranked to obtain a ranked correlation to determine highly correlated drug-gene associations.
#'@param rsq Indicate whether or not you want to output the R^2 values for the data you train on from true and predicted values.
#'These values represent the percentage in which the optimal model accounts for the variance in the training data.
#'Options are 'TRUE' and 'FALSE'.
#'@param folder Retained for compatibility with calcPhenotype arguments. Cross-validation results are returned as an object.
#'@param cvFold Indicate the number of k-folds wanted in the CV calculation. -1 indicates a leave-one-out cross validation
#'@param parallel If TRUE, run cross-validation folds in parallel. The default is FALSE.
#'@param cores The number of cores to use when parallel is TRUE. Parallel execution uses forked processes via parallel::mclapply, which is not available for multicore execution on Windows PCs; on Windows, predictionAccuracybyCV will warn and run serially.

#'@return A list containing cross-validated predicted phenotype values and real phenotype values.
#'@import sva
#'@import ridge
#'@import car
#'@import utils
#'@import stats
#'@import parallel
#'@importFrom pls pcr explvar
#'@keywords predict drug sensitivity and phenotype
#'@export
predictionAccuracybyCV <- function (trainingExprData,
                                    trainingPtype,
                                    testExprData,
                                    batchCorrect,
                                    powerTransformPhenotype=TRUE,
                                    removeLowVaryingGenes=0.2,
                                    minNumSamples,
                                    selection=1,
                                    printOutput,
                                    pcr=FALSE,
                                    removeLowVaringGenesFrom,
                                    report_pc=FALSE,
                                    cc=FALSE,
                                    percent=80,
                                    rsq=FALSE,
                                    folder = FALSE,
                                    cvFold = -1,
                                    parallel = FALSE,
                                    cores = 1)


{
  if (!is.matrix(trainingExprData)) {
    stop("\nERROR: \"trainingExprData\" must be a matrix.", call. = FALSE)
  }

  if ((ncol(trainingExprData) < minNumSamples)) {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for intrinsic difference in your training and test sets and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }

  nTrain <- ncol(trainingExprData)

  if (is.vector(trainingPtype)) {
    if (length(trainingPtype) != nTrain) {
      stop("\nERROR: \"trainingPtype\" must have one value for each training sample.", call. = FALSE)
    }
    trainingPtypeNames <- names(trainingPtype)
    trainingPtype <- matrix(trainingPtype, ncol = 1)
    rownames(trainingPtype) <- if (!is.null(trainingPtypeNames)) trainingPtypeNames else colnames(trainingExprData)
    colnames(trainingPtype) <- "Drug"
  } else if (is.matrix(trainingPtype)) {
    if (nrow(trainingPtype) != nTrain && ncol(trainingPtype) == nTrain) {
      trainingPtype <- t(trainingPtype)
    }
    if (nrow(trainingPtype) != nTrain) {
      stop("\nERROR: \"trainingPtype\" must have one row for each training sample.", call. = FALSE)
    }
    if (is.null(rownames(trainingPtype))) {
      rownames(trainingPtype) <- colnames(trainingExprData)
    }
    if (is.null(colnames(trainingPtype))) {
      colnames(trainingPtype) <- paste0("Drug", seq_len(ncol(trainingPtype)))
    }
  } else {
    stop("\nERROR: \"trainingPtype\" must be a matrix or vector.", call. = FALSE)
  }

  if (!any(colnames(trainingExprData) %in% rownames(trainingPtype))) {
    stop("\nERROR: No Cell Lines Found in Common: Sample names must be consistent in training matrices", call. = FALSE)
  }

  commonCellLines <- colnames(trainingExprData)[colnames(trainingExprData) %in% rownames(trainingPtype)]
  trainingExprData <- trainingExprData[, commonCellLines, drop=FALSE]
  trainingPtype <- trainingPtype[commonCellLines, , drop=FALSE]
  nTrain <- ncol(trainingExprData)

  cores <- as.integer(cores)
  if (is.na(cores) || cores < 1) {
    stop("\nERROR: \"cores\" must be a positive integer.", call. = FALSE)
  }
  useParallel <- isTRUE(parallel) && cores > 1
  if (useParallel && .Platform$OS.type == "windows") {
    warning("parallel=TRUE uses forked processes via parallel::mclapply, which is not available for multicore execution on Windows. Running predictionAccuracybyCV serially.", call. = FALSE)
    useParallel <- FALSE
  }

  if (is.null(testExprData)) {
    homData <- list()
    homData$selection <- selection
    homData$train <- trainingExprData
  }
  else if (!is.null(testExprData)) {
    homData <- homogenizeData(testExprData, trainingExprData,
                              batchCorrect = batchCorrect, selection = selection,
                              printOutput = printOutput)
  }
  predPtype <- matrix(NA, nrow=nTrain, ncol=ncol(trainingPtype))
  rownames(predPtype) <- colnames(trainingExprData)
  colnames(predPtype) <- colnames(trainingPtype)

  runCvTask <- function(testIndex) {
    testCvSet <- homData$train[, testIndex, drop=FALSE]
    trainCvSet <- homData$train[, setdiff(seq_len(nTrain), testIndex), drop=FALSE]
    trainPtypeCv <- trainingPtype[setdiff(seq_len(nTrain), testIndex), , drop=FALSE]

    predictions <- calcPhenotype(trainCvSet,
                                 trainPtypeCv, testCvSet, batchCorrect = "none",
                                 minNumSamples = 0, selection = homData$selection,
                                 removeLowVaryingGenes = removeLowVaryingGenes,
                                 powerTransformPhenotype = powerTransformPhenotype,
                                 printOutput = FALSE, pcr= pcr, percent= percent,
                                 removeLowVaringGenesFrom = removeLowVaringGenesFrom,
                                 report_pc=FALSE, cc=FALSE, rsq=FALSE, folder = FALSE,
                                 parallel = FALSE, cores = 1)

    list(testIndex=testIndex, predictions=predictions)
  }

  if (cvFold == -1) {
    cvTasks <- as.list(seq_len(nTrain))
    cvResults <- if (useParallel) {
      parallel::mclapply(cvTasks, runCvTask, mc.cores=cores)
    } else {
      lapply(cvTasks, runCvTask)
    }

    for (i in seq_along(cvResults)) {
      result <- cvResults[[i]]
      predPtype[result$testIndex, ] <- result$predictions
      if (printOutput && !useParallel && i%%max(1, as.integer(nTrain/5)) == 0) {
        cat(paste(i, "of", nTrain, "iterations complete. \n"))
      }
    }
    if (printOutput && useParallel) {
      cat(paste(nTrain, "of", nTrain, "iterations complete. \n"))
    }
  }
  else if (cvFold > 1) {
    randTestSamplesIndex <- sample(1:nTrain)
    sampleGroup <- rep(cvFold, nTrain)
    groupSize <- as.integer(nTrain/cvFold)
    for (j in 1:(cvFold - 1)) {
      sampleGroup[(((j - 1) * groupSize) + 1):(j * groupSize)] <- rep(j,
                                                                      groupSize)
    }
    cvTasks <- split(randTestSamplesIndex, sampleGroup)
    cvResults <- if (useParallel) {
      parallel::mclapply(cvTasks, runCvTask, mc.cores=cores)
    } else {
      lapply(cvTasks, runCvTask)
    }

    for (j in seq_along(cvResults)) {
      result <- cvResults[[j]]
      predPtype[result$testIndex, ] <- result$predictions
      if (printOutput && !useParallel) {
        cat(paste("\n", j, " of ", cvFold, " iterations complete.",
                  sep = ""))
      }
    }
    if (printOutput && useParallel) {
      cat(paste("\n", cvFold, " of ", cvFold, " iterations complete.",
                sep = ""))
    }
  }
  else {
    stop("Unrecognised value of \"cvFold\"")
  }

  if (ncol(predPtype) == 1) {
    predPtype <- predPtype[, 1]
    trainingPtype <- trainingPtype[, 1]
  }

  finalData <- list(cvPtype = predPtype, realPtype = trainingPtype)
  return(finalData)
}
