#'This function performs an iterative matrix completion algorithm to predict drug response for pre-clinical data when there are missing ('NA') values.
#'@param senMat A matrix of drug sensitivity data with missing ('NA') values. rownames() are samples (e.g. cell lines), and colnames() are drugs.
#'@param nPerms The number of iterations that the EM-algorithm (expectation maximization approach)  run. The default is 50, as previous findings recommend 50 iterations (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1050-9)
#'@return A matrix of drug sensitivity scores without missing values. rownames() are samples, and colnames are drugs.
#'@keywords Drug response prediction.
#'@import glmnet
#'@import stats
#'@import car
#'@import utils
#'@export
completeMatrix <- function(senMat, nPerms=50)
{
  message("\nNumber of iterations:")
  #To initialize the algorithm, all missing values are first imputed to the median.
  numCellLinesNotScreened <- apply(senMat, 2, function(r)return(sum(is.na(r))))
  numDrugsNotScreened <- apply(senMat, 1, function(r)return(sum(is.na(r))))
  indexNotScreened <- apply(senMat, 2, function(r)return(which(is.na(r))))
  drugMed <- apply(senMat, 2, function(col)return(median(na.omit(as.numeric(col)))))
  hundIc50sImpute <- data.matrix(senMat)
  for(i in 1:ncol(senMat))
  {
    datImpute <- as.numeric(senMat[,i])
    datImpute[is.na(datImpute)] <- drugMed[i]
    hundIc50sImpute[,i] <- datImpute
  }

  #Sort the matrix by least number of missing values to most.
  #The rows (e.g. cell lines) of the senMat matrix must be ordered, ascending by the number of missing values (so the drugs toward the bottom have more missing values
  #then the drugs at the top).
  hundIc50sImputeSort <- hundIc50sImpute[, order(numCellLinesNotScreened[colnames(hundIc50sImpute)])]

  #Calculate/estimate the PCs of this matrix X.
  matrix<-matrix(as.numeric(unlist(hundIc50sImputeSort)),nrow=nrow(hundIc50sImputeSort))
  colnames(matrix)<-colnames(hundIc50sImputeSort)
  rownames(matrix)<-rownames(hundIc50sImputeSort)
  hundIc50sImputeSort<-matrix
  pcs <- prcomp(hundIc50sImputeSort)$x

  #Fit a lasso model for all other drugs, for all samples other than those NOT screened by the drug we are predicting for. By default, we will iterate 50 times.
  #For each missing value of X, fit a lasso regression model using the PCs of all other rows of X as predictors of all other values (both measured and imputed values)
  #for that cell line.
  imputeSortList <- list()
  imputeSortList[[1]] <- hundIc50sImputeSort
  medianDistance <- numeric()
  for(j in 1:nPerms)
    #Repeat this procedure iteratively until the total change in A (the sum of the total difference between each of the elements in X and X') converges.
    #This takes approximately 50 iterations, which is why 50 is the default.
  {
    pb <- txtProgressBar(min = 0, max = ncol(hundIc50sImputeSort), style = 3) #Create progress bar.
    for(i in 1:ncol(hundIc50sImputeSort))
    {
      setTxtProgressBar(pb, i) #Update progress bar.
      pcs <- prcomp(hundIc50sImputeSort)$x # calcualte the PCs of the current matrix #hundIc50sImputeSort.
      indexNotScreenedThisDrug <- indexNotScreened[[colnames(hundIc50sImputeSort)[i]]] #Index of "imputed" cell lines for this drug.

      if(length(indexNotScreenedThisDrug) > 0)
      {
        #Make the training (design matrix) for the non-imputed cell lines for this drug.
        pcMat <- pcs[-indexNotScreenedThisDrug, ]
        dmRaw <-  model.matrix(~pcMat)

        #Calculate lambda using CV.
        #The tuning parameter lambda for the lasso regression is selected using cross-validation (cv.glment)
        FitRaw <- cv.glmnet(dmRaw, hundIc50sImputeSort[-indexNotScreenedThisDrug, i], nfolds=20)
        getLambdasRaw <- FitRaw$lambda.min #The optimal lambda value.

        #Fit the model and extract the co-efficients.
        theModRaw <- glmnet(dmRaw, hundIc50sImputeSort[-indexNotScreenedThisDrug, i], lambda=getLambdasRaw)
        coef(theModRaw)[,1][coef(theModRaw)[,1] != 0]
        betas <- coef(theModRaw)[,1][coef(theModRaw)[,1] != 0][-1]
        intercept <- coef(theModRaw)[,1][coef(theModRaw)[,1] != 0][1]
        names(betas) <- substring(names(betas), 6)

        #Use the model to update the estimates for the imputed values for this drug.
        #Apply the model to yield an updated estimate for each missing value.
        #Repeat this procedure for all missing values of X to yield an updated matrix X'
        if(length(indexNotScreenedThisDrug) == 1)
        {
          predictions <- intercept + apply((t(pcs[indexNotScreenedThisDrug, names(betas)]) * betas), 1, sum); # P redict new values
        }
        else
        {
          predictions <- intercept + apply(t(t(pcs[indexNotScreenedThisDrug, names(betas)]) * betas), 1, sum);
        }

        hundIc50sImputeSort[indexNotScreenedThisDrug, i] <- predictions  # Here is the updated matrix
      }

    }
    close(pb)
    imputeSortList[[j + 1]] <- hundIc50sImputeSort
    medianDistance[j] <- mean(as.numeric(imputeSortList[[j]]) - as.numeric(imputeSortList[[j + 1]])) #Estimate A (the sum of the total difference between each of the elements of X and X').
    message(paste("\nIteration: ", j, "\n"), "")
    #plot(medianDistance, main=j, ylab="Median Distance from previous imputed matrix")
  }
  message("\nDone\n")

  write.table(hundIc50sImputeSort[,colnames(senMat)], file='./complete_matrix_output.txt')

  #return(hundIc50sImputeSort[,colnames(senMat)])
}
#'This function determines drug-gene associations for pre-clinical data.
#'@param drugMat A matrix of drug sensitivity data. rownames() are pre-clinical samples, and colnames() are drug names.
#'@param drugRelatedness A matrix in which column 1 contains a list of compounds, and column 2 contains a list of their corresponding target pathways. Given the subjective nature of
#'drug classification, please ensure these pathways are as specific as possible for accurate results.
#'@param markerMat A matrix containing the data for which you are looking for an association with drug sensitivity (e.g. a matrix of somatic mutation data). rownames() are marker names (e.g. gene names), and colnames() are samples.
#'@param threshold Determine the correlation coefficient. Drugs with a correlation coefficient greater than or equal to this number with the drug under scrutiny will be removed from the negative control group.
#'The default is 0.7
#'@param minMuts The minimum number of non-zero entries required so that a p-value can be calculated (e.g. how many somatic mutations must be present). The default is 5.
#'@param additionalCovariateMatrix A matrix containing covariates to be fit in the drug biomarker association models. This could be, for example, tissue of origin or cancer type. Rows are sample names. The default is NULL.
#'@param expression A matrix of expression data. rownames() are genes, and colnames() are the same pre-clinical samples as those in the drugMat (also in the same order).
#'The default is NULL. If expression data is provided, a gene signature will be obtained.
#'@import stats
#'@import utils
#'@import ridge
#'@export
glds <- function(drugMat, drugRelatedness, markerMat, minMuts=5, additionalCovariateMatrix=NULL, expression=NULL, threshold=0.7){

  results_gldsPs <- list()
  results_gldsBetas <- list()
  results_naivePs <- list()
  results_naiveBetas <- list()
  numDrugs <- ncol(drugMat)
  pairCor <- cor(drugMat, method="spearman")
  comNames <- colnames(markerMat)[colnames(markerMat) %in% rownames(drugMat)] # cell lines for which we have both mutation and drug data....

  if(!is.null(additionalCovariateMatrix)) # if additional co variate matrix was provided, then also subset to those samples
  {
    comNames <- comNames[comNames %in% rownames(additionalCovariateMatrix)]
  }

  pb <- txtProgressBar(min = 0, max = ncol(drugMat), style = 3) # create progress bar

  #drugMat=drugMat[,1:20]

  for(i in 1:ncol(drugMat)){
    #index<-(match(colnames(drugMat)[6], drugRelatedness[,1]))

    index<-try(match(colnames(drugMat)[i], drugRelatedness[,1]), silent=TRUE)
    if (is.na(index)){
      drugMat = drugMat[,-i] #Remove that column/drug from the matrix.
      message(paste('\n', colnames(drugMat)[i], 'is skipped because it is not included in your drug relatedness info'))
    }else{
      #Calculate 10 PCs on non-related sets of drugs....
      categoryThisDrug<-(drugRelatedness[,2])[index] #The MOA of the drug under scrutiny.
      thisDrug<-(drugRelatedness[,1])[index] #The drug under scrutiny.

      #Identify all the other drugs in this dataset that are not part of this same drug category...
      negControlDrugs<-c() #Drugs that don't have the same MOA as the drug under scrutiny.
      categoryThisDrug_vec<-scan(text=categoryThisDrug, what="", quiet=TRUE) #The MOA string for the drug of interest will be broken up into elements (each element is a word).
      for (k in 1:length(drugRelatedness[,2])){ #For each MOA...
        s1<-(drugRelatedness[,2])[k] #String one...the MOA for a drug.
        vec<-scan(text = s1, what =  "", quiet=TRUE) #The string is broken up into a vector by spaces such that each element is a word in the string.
        if (FALSE %in% (categoryThisDrug_vec %in% vec)){
          negControlDrugs[k]<-(drugRelatedness[,1])[k]
        }
      }

      #indices<-which(drugRelatedness[,2] != categoryThisDrug) #Indices of categories that aren't the category of the current drug.
      #negControlDrugs<-unique((drugRelatedness[,1])[indices])

      mags<-sort(abs(pairCor[, colnames(drugMat)[i]]), decreasing=FALSE) >= threshold
      pairwiseCorNear<-names(which(mags == "TRUE"))

      negControlDrugs <- setdiff(negControlDrugs, pairwiseCorNear) # remove very highly correlated drugs from "negative controls"

      indices<-match(negControlDrugs, colnames(drugMat))
      indices<-indices[!is.na(indices)]
      controlPCsAll <- prcomp(drugMat[, indices])$x

      controlPCsAllCom <- controlPCsAll[comNames, ]

      #This section determines the gene signature.
      if(!is.null(expression)){ #If expression matrix was provided, then gene signatures can also be obtained...
        spearCorList<-list()
        for(p in 1:10){
          exprComSamples<-rownames(controlPCsAll)[rownames(controlPCsAll) %in% colnames(expression)] #Samples for which you have mutation, expression, and sensitivity data for.
          for(g in 1:nrow(expression)){
            spearCorList[[p]][g]<-cor.test(controlPCsAll[exprComSamples,p],
                                           expression[g,exprComSamples],
                                           method='spearman')$p.value
          }
          names(spearCorList[[p]])<-rownames(expression)
        }
        mdrExprGenesList[[j]] <- unique(as.character(unique(sapply(spearCorList, function(vec)return(names(sort(vec)[1:50])))))) #The full list of genes for this drug.
        names(mdrExprGenesList)<-colnames(drugMat)
        controlGeneFrequency <- table(do.call(c, mdrExprGenesList)) #How often is each gene selected as a negative control
        alwaysControlGene <- names(controlGeneFrequency[controlGeneFrequency == ncol(drugMat)])
      }

      # Calculate the P-values and beta values for each marker for this drug, controlling for GLDS and not controlling for GLDS
      results_gldsPs[[i]] <- numeric()
      results_gldsPs[[i]] <- rep(NA, nrow(markerMat))
      results_gldsBetas[[i]] <- numeric()
      results_gldsBetas[[i]] <- rep(NA, nrow(markerMat))
      results_naivePs[[i]] <- numeric()
      results_naivePs[[i]] <- rep(NA, nrow(markerMat))
      results_naiveBetas[[i]] <- numeric()
      results_naiveBetas[[i]] <- rep(NA, nrow(markerMat))
      for(j in 1:nrow(markerMat)){
        if(sum(markerMat[j, comNames]) > minMuts)
        {
          if(is.null(additionalCovariateMatrix)) # if no additional covariate have been provided....
          {
            theCoefs <- coef(summary(lm(drugMat[comNames, i]~markerMat[j, comNames])))
            results_naivePs[[i]][j] <- theCoefs[2, 4]
            results_naiveBetas[[i]][j] <- theCoefs[2, 1]

            theCoefs <- coef(summary(lm(drugMat[comNames, i]~markerMat[j, comNames]+controlPCsAllCom[,1:10])))
            results_gldsPs[[i]][j] <- theCoefs[2, 4]
            results_gldsBetas[[i]][j] <- theCoefs[2, 1]
          }
          else # if there are other covariates, include them in the models.
          {
            theCoefs <- coef(summary(lm(drugMat[comNames, i]~markerMat[j, comNames]+additionalCovariateMatrix[comNames,])))
            results_naivePs[[i]][j] <- theCoefs[2, 4]
            results_naiveBetas[[i]][j] <- theCoefs[2, 1]

            theCoefs <- coef(summary(lm(drugMat[comNames, i]~markerMat[j, comNames]+controlPCsAllCom[,1:10]+additionalCovariateMatrix[comNames,])))
            results_gldsPs[[i]][j] <- theCoefs[2, 4]
            results_gldsBetas[[i]][j] <- theCoefs[2, 1]
          }
        }
      }
      # cat(paste(i, " ", sep=""))
      names(results_gldsPs[[i]]) <- rownames(markerMat)
      names(results_naivePs[[i]]) <- rownames(markerMat)
      names(results_gldsBetas[[i]]) <- rownames(markerMat)
      names(results_naiveBetas[[i]]) <- rownames(markerMat)

      setTxtProgressBar(pb, i) # update progress bar

    }
  }

  close(pb)
  names(results_gldsPs) <- colnames(drugMat)
  names(results_naivePs) <- colnames(drugMat)
  names(results_gldsBetas) <- colnames(drugMat)
  names(results_naiveBetas) <- colnames(drugMat)

  outList <- list(pGlds=results_gldsPs, betaGlds=results_gldsBetas, pNaive=results_naivePs, betaNaive=results_naiveBetas)

  #return(outList)
  write.csv(results_gldsPs, file="./gldsPs.csv", row.names = TRUE, col.names = TRUE)
  write.csv(results_naivePs, file="./naivePs.csv", row.names = TRUE, col.names = TRUE)
  write.csv(results_gldsBetas, file="./gldsBetas.csv", row.names = TRUE, col.names = TRUE)
  write.csv(results_naiveBetas, file="./naiveBetas.csv", row.names = TRUE, col.names = TRUE)

  if(!is.null(additionalCovariateMatrix)){
    write.csv(alwaysControlGene, file="./gene_signature.txt")
  }

}
