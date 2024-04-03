#'This function maps cnv data to genes. The output of this function is a .RData file called map.RData; this file contains theCnvQuantVecList_mat (rows are genes, and columns are samples) and tumorSamps (indicates which samples are primary tumor samples, 01A).
#'@param Cnvs The cnv data. A table with the following colnames: Sample (named using the TCGA patient barcode), Chromosome, Start, End, Num_Probes, and Segment_Mean.
#'@keywords Map CNV data to genes
#'@return A .RData file called, map.RData, which stores two objects: theCnvQuantVecList_mat (rows are genes, columns are samples), tumorSamps (indicates which samples are primary tumor/01A). This output will serve as the input for test().
#'@import org.Hs.eg.db
#'@import TCGAbiolinks
#'@import GenomicFeatures
#'@import TxDb.Hsapiens.UCSC.hg19.knownGene
#'@import utils
#'@import stats
#'@importFrom BiocGenerics toTable
#'@importFrom GenomicRanges GRanges
#'@importFrom IRanges IRanges subsetByOverlaps findOverlaps countOverlaps
#'@importFrom S4Vectors Rle
#'@export
map_cnv<-function(Cnvs)
{
  #Check colnames() of the cnv data.
  #This is important because depending on what method you use to obtain the cnv data, colnames() might slightly differ (aka Segment.Mean and not Segment_Mean).
  if (('Segment_Mean' %in% colnames(Cnvs)) == FALSE)
    stop("\nERROR: Check colnames() of cnv data. colnames() must include Sample, Chromosome, Start, End, and Segment_Mean")

  #Load the gene ranges for HG19 using.
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  geneRanges <- genes(txdb) #seqnames (chr), ranges, gene id.
  e2s = toTable(org.Hs.egSYMBOL) #gene id, gene symbol.
  syms <- e2s[, "symbol"] #Gene symbols for HG19.
  names(syms) <- e2s[, "gene_id"] #Gene symbols with gene id.
  theGeneSymsOrd <- syms[as.character(geneRanges$gene_id)]

  #path<-paste(theRootDir, '/mut.txt', sep="")
  #Cnvs<-read.table(path, header=TRUE, row.names=1) #The cnv data downloaded with columns sample, chromosome, start, end, num probes, segment mean.
  CnvsList <- split(Cnvs, Cnvs[, "Sample"]) #Put each patient data into list such that each element in the list is for a patient.

  cnvEdGenesList <- list()
  ampGenesList <- list()
  delGenesList <- list()
  numGenesQuantifid <- numeric()
  theCnvQuantVecList <- list()

  for(i in 1:length(CnvsList)){ #Iterate for each patient in the list CnvsList

    chrs <- paste("chr", CnvsList[[i]]$Chromosome, sep="") #chrs contains all the chromosome numbers for these samples in patient i.
    starts <- CnvsList[[i]]$Start #starts contains all the starts.
    ends <- CnvsList[[i]]$End #ends contains all the ends.
    grCnvs <- GRanges(seqnames=Rle(chrs),ranges=IRanges(starts, ends), segMeans=CnvsList[[i]]$Segment_Mean) #chr start-end

    #Amp or del
    grCnvs_ampDel <- grCnvs[grCnvs$segMeans > 1 | grCnvs$segMeans < -1]
    cnvedGenes <- subsetByOverlaps(geneRanges, grCnvs_ampDel, type="within")
    cnvEdGenesList[[i]] <- cnvedGenes$gene_sym

    #Amps
    grCnvs_amp <- grCnvs[grCnvs$segMeans > 1]
    ampedGenes <- subsetByOverlaps(geneRanges, grCnvs_amp, type="within")
    ampGenesList[[i]] <- ampedGenes$gene_sym

    #Dels
    grCnvs_Del <- grCnvs[grCnvs$segMeans < -1]
    deledGenes <- subsetByOverlaps(geneRanges, grCnvs_Del, type="within")
    delGenesList[[i]] <- deledGenes$gene_sym

    #Continuous gene level
    #Use count overlaps to find genes that unambiguously overlap a single peak. Give it an NA it it doesn't overlap a single peak. Assign it the value of the peak if it unambiguously overlaps a peak. PC.
    numOverlaps <- countOverlaps(geneRanges, grCnvs)
    numGenesQuantifid[i] <- sum(numOverlaps == 1)
    inCnv <- which(numOverlaps == 1) # take only gene unambiguously overlaping a peak, this is usually most genes.

    theCnvQuantVec <- rep(NA, length(geneRanges))
    olaps <- findOverlaps(geneRanges, grCnvs, type="within")

    #theCnvQuantVec[olaps@queryHits] <- grCnvs$segMeans[olaps@subjectHits]

    #No slot of name 'subjectHits' for this object of class 'SortedByQueryHits'
    #This problem is caused by a small change in the GenomicAlignments package.
    theCnvQuantVec[olaps@from] <- grCnvs$segMeans[olaps@to]

    theCnvQuantVecList[[i]] <- theCnvQuantVec
    #names(theCnvQuantVecList[[i]]) <- geneRanges$gene_sym
    names(theCnvQuantVecList[[i]]) <- theGeneSymsOrd
  }

  names(theCnvQuantVecList) <- names(CnvsList)
  theCnvQuantVecList_mat <- do.call(rbind, theCnvQuantVecList)
  siteVec <- sapply(strsplit(names(CnvsList), "-"), function(l)return(l[4]))
  tumorSamps <- which(siteVec == "01A")
  save(theCnvQuantVecList_mat, tumorSamps, file="./map.RData") # Save these RData files for use by other scripts.
}
#'This function will test every drug against every CNV or somatic mutation for your cancer type.
#'@param drug_prediction The drug prediction data. Must be a data frame. rownames are samples, colnames are drugs. Make sure sample names are of the same form as the sample names in your cnv or mutation data. e.g. if the rownames() are TCGA barcodes of the form TCGA-##-####-###, make sure your cnv/mutation data also uses samples in the form TCGA-##-####-###
#'@param data The cnv or mutation data. Must be a data frame. If you wish to use cnv data, use the output from map_cnv(), transpose it so that colnames() are samples. Or use data of similar form. If you wish to use mutation data, use the method for downloading mutation data outlined in the vignette, and make sure the TCGA barcodes use '-' instead of '.'; if you use another dataset (and don't download data from TCGA), make sure your data file includes the following columns: 'Variant_Classification', 'Hugo_Symbol', 'Tumor_Sample_Barcode'.
#'@param n The minimum number of samples you want CNVs or mutations to be amplified in. The default is 10 (arbitrarily chosen).
#'@param cnv TRUE or FALSE. Indicate whether or not you would like to test cnv data. If TRUE, you will test cnv data. If FALSE, you will test mutation data.
#'@keywords Test CNV or mutation data to genes.
#'@import org.Hs.eg.db
#'@import TxDb.Hsapiens.UCSC.hg19.knownGene
#'@import GenomicFeatures
#'@import utils
#'@import stats
#'@import ridge
#'@import parallel
#'@return Raw p-value and beta-values for cnv and somatic mutations.
#'@export
idwas<-function(drug_prediction, data, n=10, cnv){
  #Check parameters.
  #_____________________________________________________________________________
  if (!is.data.frame(drug_prediction))
    stop("\nERROR: \"drug_prediction\" must be a data frame")
  if (!is.data.frame(data))
    stop("\nERROR: \"data\" must be a data frame")

  if (cnv == FALSE){
    drug_prediction<-t(drug_prediction)
    #This will make it such that rows are drugs and columns are samples.
  }

  #If TCGA is in my colnames() (as it would if you got cnv data from map_cnv() OR
  #if TCGA is in the column of your mutation data mutation$Tumor_Sample_Barcode, then you have TCGA samples
  #and you want to make sure you only use 01A samples).
  cols<-colnames(data)

  #_____________________________________________________________________________

  if( (sum(grepl('TCGA', cols) == TRUE) > 0) | (TRUE %in% (grepl('TCGA', data$Tumor_Sample_Barcode))) ){ #If TCGA...
    if(cnv){ #If TCGA CNV...
      #If cnv=TRUE, then you have cnv data, so proceed with this testing...
      #Find the indices of the 01A/primary tumor samples and the non-duplicates.
      indices <- which(sapply(strsplit(cols, "-"), function(a)a[4]) == "01A")
      #Create a matrix using only the primary tumor patients as columns.
      matrix <- data[, indices] #Make a matrix of rows/genes and columns/TCGA O1A samples.

      #Collect all the patient names/rows of the drug prediction matrix that are 01A primary tumors.
      drugs_01A <- rownames(drug_prediction)[which(sapply(strsplit(rownames(drug_prediction), "-"), function(a)a[4]) == "01A")]
      #Collect the patient IDs of all the patients in drugs_01A (the patients IDS for the 01A samples you have drug prediction data for). These are the 4 digits after the bacode TCGA-## so TCGA-##-????
      drugs_01Aids <- sapply(strsplit(drugs_01A, "-"), function(a)a[3])
      #Index the drug prediction matrix so that its rows/patients are the primary tumor patients.
      matrix2 <- drug_prediction[drugs_01A,]
      #Rename the rows of the drug prediction matrix to the patient id (the 4 digits).

      #Make sure you only use unique ids (sometimes there are duplicates). If there are duplicates, remove.
      indices<-match(unique(drugs_01Aids), drugs_01Aids)
      matrix2<-matrix2[indices,] #Only keep the rows that have unique ids.
      rownames(matrix2) <- drugs_01Aids[indices]

      #Collect all the patient IDs of the patients in columns of matrix
      ids <- sapply(strsplit(colnames(matrix), "-"), function(a)a[3])

      #Collect the patients that are common to both prediction and cnv/mut data (using their 4 digit ids).
      overlapping_ids<-intersect(drugs_01Aids, ids)

      #Index the drug prediction matrix so that its rows only contains patients it has in common with the cnv/mut data.
      indices<-match(overlapping_ids, drugs_01Aids)
      drug_prediction2<-drug_prediction[indices,]

      indices<-match(overlapping_ids, ids)
      matrix3<-as.matrix(matrix[,indices])

      MatCommonPats_amps <- apply(matrix3, 2, function(theCol)return(as.numeric(theCol > 1))) #Apply this function to each column.
      rownames(MatCommonPats_amps) <- rownames(matrix3)

      MatCommonPats_dels <- apply(matrix3, 2, function(theCol)return(as.numeric(theCol < -1)))

      rownames(MatCommonPats_dels) <- rownames(matrix3)

      theFun <- function(j){
        pVals <- numeric()
        betaVal <- numeric()
        if(sum(na.omit(MatCommonPats_amps[j, ])) > n){ #Make sure the gene is amplifed at least 50 times
          for(i in 1:ncol(drug_prediction2)){
            theMod <- coef(summary(lm(drug_prediction2[,i]~MatCommonPats_amps[j, ]))) #Linear regression betwween each amplified CNV (j/row) and each drug (i/column).
            pVals[i] <- theMod[2,4] #P-values for each model.
            betaVal[i] <- theMod[2,1] #Effect size for each model.
          }
          names(pVals) <- colnames(drug_prediction2)
          names(betaVal) <- colnames(drug_prediction2)
        }
        return(list(pVals, betaVal))
      }

      #mclapply is a parallelized vector of lapply.
      #It's structure is mclapply(X, FUN...) where X is a vector and FUN is the function applied to each element of X.
      #Here, the function is applied to each value of X which is j (each CNV).
      allCors <- mclapply(1:nrow(matrix3), theFun) #You may want to change this number based on the number of available cores.
      #allCors'is a list where each element/drug consists of 2 elements (the p value for that drug's association with the CNV and the beta value).

      #Name each element after each CNV.
      names(allCors) <- rownames(matrix3)

      hasAmps <- apply(MatCommonPats_amps, 1, function(theRow)return(sum(na.omit(theRow)) > n)) #Restrict analysis to CNAs that occur in 50 or more samples.

      allCors_hasAmps <- allCors[hasAmps]
      if(length(allCors_hasAmps) == 0){
        stop((paste("\nERROR: A gene was not mutated in at least", n, "patients. Recommend decreasing the n parameter.", sep=" ")))
      }else{
        pVals <- sapply(allCors_hasAmps, function(item)return(item[[1]]))
        betas <- sapply(allCors_hasAmps, function(item)return(item[[2]]))
        write.csv(pVals, file='./CnvTestOutput_pVals.csv')
        write.csv(betas, file='./CnvTestOutput_betas.csv')
      }
    }else{ #If TCGA MUT...
      #Obtain a list for the patients you have data for and initiate empty lists to fill.
      #_______________________________________
      tcgaIds<-data$Tumor_Sample_Barcode
      unique<-unique(tcgaIds)
      mutsListAll <- list() #A list of all the mutation occurring in each sample.
      proteinChangingMutations <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "In_Frame_Del",
                                    "Frame_Shift_Ins", "In_Frame_Ins", "Nonstop_Mutation", "De_novo_Start_OutOfFrame",
                                    "De_novo_Start_InFrame", "Missense", "Read-through", "Indel")
      genesWithProteinChangeList <- list()

      #Fill those lists.
      #_______________________________________
      for(i in 1:length(unique)){
        #print(paste(i, "of", length(unique)), sep='')
        indices<-unique[i] == tcgaIds

        #Now for each sample, pull out a list of genes with somatic mutations
        variantType <- data[indices,"Variant_Classification"]
        theGenes <- data[indices,"Hugo_Symbol"]
        names(variantType) <- theGenes
        mutsListAll[[i]] <- variantType #A list of all the mutation occurring in each sample.
        genesWithProteinChangeList[[i]] <- unique(names(variantType[variantType %in% proteinChangingMutations]))
      }

      allMutatedGenes <- unique(names(unlist(mutsListAll))) #All the unique genes mutated across all patients.
      mutationTypes <- table(unlist(mutsListAll)) #Frequency of mutation type across all patients.
      mutsListAll_unlist <- unlist(mutsListAll, recursive=F) #Mutation type and gene association.
      genesWithProteinChangeList_unlist <- unlist(genesWithProteinChangeList, recursive=F) #Genes with protein changing mutations across all patients.

      #From mutsListAll we can then create a matrix indicating if the gene has a coding mutation.
      #_______________________________________
      mutMat <- numeric((length(unique)*length(allMutatedGenes)))
      dim(mutMat) <- c(length(allMutatedGenes), length(unique))
      rownames(mutMat) <- allMutatedGenes
      colnames(mutMat) <- unique

      #Now populate this matrix with the relevant information about what kind of mutation each gene has in each sample.
      #_______________________________________
      for(i in 1:length(unique)){
        #print(paste(i, "of", length(unique)), sep='')
        mutMat[genesWithProteinChangeList[[i]], i] <- rep(1, length(genesWithProteinChangeList[[i]]))
      }

      tumorTypeId <- sapply(strsplit(colnames(mutMat), "-", fixed=TRUE), function(l)return(l[4]))

      #Lets remove everything but the "Primary Solid Tumors (i.e. "01")".
      #_______________________________________
      mutMat_only01 <- mutMat[, tumorTypeId == "01A"] #Genes that were mutated in 01A samples.
      theIds <- colnames(mutMat_only01) #Patient ids that were 01A.
      mutIds <- sapply(strsplit(theIds, "-", fixed=T), function(l)return(l[3])) #The TCGA #### part of the barcode.
      colnames(mutMat_only01) <- mutIds

      #Extract the 01a samples from the drug prediction data, i.e. tumor samples.
      #_______________________________________
      all01ASamples <- colnames(drug_prediction)[which(sapply(strsplit(colnames(drug_prediction), "-", fixed=T), function(a)a[4]) == "01A")]
      preds01a <- drug_prediction[, all01ASamples] #Predictions for 01A samples.
      sampIds01a <- sapply(strsplit(all01ASamples, "-", fixed=T), function(l)return(l[3])) #The TCGA #### digit number.
      colnames(preds01a) <- sampIds01a
      inPredAndMutData <- sampIds01a[sampIds01a %in% mutIds] #Samples for which we have both predicted drug response and mutation calls

      if (length(inPredAndMutData) == 0){
        stop((paste("\nERROR: Samples in drug_prediction and data don't overlap. Make sure their sample identifiers are of similar form")))
      }

      #Run the associations between all genes and drugs, for drugs with at least 50 mutations.
      #_______________________________________
      preds01a_filt_ord <- as.matrix(preds01a[, inPredAndMutData]) #The preds for the 01A samples we have both prediction and mutation data for.
      mutMat_nodups_ordFilt <- mutMat_only01[, inPredAndMutData]
      commonMuts <- apply(mutMat_nodups_ordFilt, 1, sum)
      if (length(which(commonMuts >= n)) == 0){
        stop((paste("\nERROR: A gene was not mutated in at least", n, "patients. Recommend decreasing the n parameter.", sep=" ")))
      }
      commonlyMutated <- mutMat_nodups_ordFilt[which(commonMuts >= n), ] #Rows are genes that are mutated in n+ patients. Cols are patients.

      #If there are multiple genes, commonlyMutated will have dimensions.
      #Otherwise, it will be a vector.
      #_______________________________________
      dim<-try(dim(as.vector(unlist(commonlyMutated))), silent=TRUE)
      if(length(dim)){ #If we have a vector (one gene)...
        #Get p values and beta values.
        pValList <- list()
        betaValList <- list()

        suppressWarnings(for(i in 1:nrow(preds01a_filt_ord)){ #For each drug...
          pValList[[i]] <- numeric()
          betaValList[[i]] <- numeric()
          thecoefs <- coef(summary(lm(preds01a_filt_ord[i,]~commonlyMutated)))
          pValList[[i]]<- thecoefs[2,4]
          betaValList[[i]] <- thecoefs[2,1]
        })

        #Get the adjusted p-value for each gene-drug combination, pull out the significant associations
        #and create a supplementary table that lists these for "predictable" drugs?.
        sigPs <- list()
        pAdjListCantype <- list()
        for(i in 1:length(pValList)){
          #names(pValList[[i]]) <- rownames(commonlyMutated)
          #names(betaValList[[i]]) <- rownames(commonlyMutated)
          padj <- p.adjust(pValList[[i]], method="BH")
          sigPs[[i]] <- padj[padj < 0.05]
          pAdjListCantype[[i]] <- padj
        }
      }else{ #Otherwise, we have rows/multiple genes...
        #If there are gene entries with an unknown HUGO ID, remove it.
        #_______________________________________
        if("Unknown" %in% rownames(commonlyMutated)){
          indices<-'Unknown' %in% rownames(commonlyMutated)
          commonlyMutated<-commonlyMutated[-indices,]
        }

        #Get p values and beta values.
        #_______________________________________
        pValList <- list()
        betaValList <- list()

        suppressWarnings(for(i in 1:nrow(preds01a_filt_ord)){ #For each drug...
          pValList[[i]] <- numeric()
          betaValList[[i]] <- numeric()
          for(j in 1:nrow(commonlyMutated)) #For each row/gene that is mutated in n+ patients...
          {
            thecoefs <- coef(summary(lm(preds01a_filt_ord[i,]~commonlyMutated[j,])))
            pValList[[i]][[j]] <- thecoefs[2,4]
            betaValList[[i]][[j]] <- thecoefs[2,1]
          }
        })
        #Get the adjusted p-value for each gene-drug combination, pull out the significant associations
        #and create a supplementary table that lists these for "predictable" drugs?.
        #_______________________________________
        sigPs <- list()
        pAdjListCantype <- list()
        for(i in 1:length(pValList)){
          names(pValList[[i]]) <- rownames(commonlyMutated)
          names(betaValList[[i]]) <- rownames(commonlyMutated)
          padj <- p.adjust(pValList[[i]], method="BH")
          sigPs[[i]] <- padj[padj < 0.05]
          pAdjListCantype[[i]] <- padj
        }
      }

      names(sigPs) <- rownames(preds01a_filt_ord)
      names(pValList) <- rownames(preds01a_filt_ord)
      names(betaValList) <- rownames(preds01a_filt_ord)
      names(pAdjListCantype) <- rownames(preds01a_filt_ord)

      pVal<-unlist(pValList)
      betaVal<-unlist(betaValList)

      #If you only have one gene of interest...make sure the row names are appropriate (drug:gene)
      if(length(dim)){ #If we have a vector (one gene)...
        final_data<-cbind(pVal, betaVal)
        rows<-rownames(final_data)
        gene<-names(which(commonMuts >= n))
        final_rows<-paste(rows, ':', gene, sep='')
        rownames(final_data)<-final_rows
        write.csv(final_data, file='./MutationTestOutput_pVal_and_betaVal.csv')
      }else{
        write.csv(cbind(pVal, betaVal), file='./MutationTestOutput_pVal_and_betaVal.csv')
      }
    }#The end of the first else statement.

    #This code is for when you don't have TCGA barcoded samples (it's similar to above)
  }else{ #If non TCGA...
    if(cnv){ #If non TCGA CNV...
      overlapping_samples<-intersect(rownames(drug_prediction), colnames(data))

      #Index the drug prediction matrix so that it contains the overlapping samples.
      indices<-match(overlapping_samples, rownames(drug_prediction))
      drug_prediction2<-drug_prediction[indices,]
      #Index the cnv.mutation matrix so that it contains the overlapping samples.
      indices<-match(overlapping_samples, colnames(data))
      matrix3<-as.matrix(data[,indices])

      MatCommonPats_amps <- apply(matrix3, 2, function(theCol)return(as.numeric(theCol > 1))) #Apply this function to each column.
      rownames(MatCommonPats_amps) <- rownames(matrix3)

      MatCommonPats_dels <- apply(matrix3, 2, function(theCol)return(as.numeric(theCol < -1)))

      rownames(MatCommonPats_dels) <- rownames(matrix3)

      theFun <- function(j){
        pVals <- numeric()
        betaVal <- numeric()

        if(sum(na.omit(MatCommonPats_amps[j, ])) > n){ # Make sure the gene is amplifed at least 50 times
          for(i in 1:ncol(drug_prediction2)){
            theMod <- coef(summary(lm(drug_prediction2[,i]~MatCommonPats_amps[j, ]))) #Linear regression betwween each amplified CNV (j/row) and each drug (i/column).
            pVals[i] <- theMod[2,4] #P-values for each model.
            betaVal[i] <- theMod[2,1] #Effect size for each model.
          }
          names(pVals) <- colnames(drug_prediction2)
          names(betaVal) <- colnames(drug_prediction2)
        }
        return(list(pVals, betaVal))
      }

      #mclapply is a parallelized vector of lapply.
      #It's structure is mclapply(X, FUN...) where X is a vector and FUN is the function applied to each element of X.
      #Here, the function is applied to each value of X which is j (each CNV).
      allCors <- mclapply(1:nrow(matrix3), theFun) #You may want to change this number based on the number of available cores.
      #allCors'is a list where each element/drug consists of 2 elements (the p value for that drug's association with the CNV and the beta value).

      #Name each element after each CNV.
      names(allCors) <- rownames(matrix3)

      hasAmps <- apply(MatCommonPats_amps, 1, function(theRow)return(sum(na.omit(theRow)) > n)) #Restrict analysis to CNAs that occur in 50 or more samples.

      allCors_hasAmps <- allCors[hasAmps]
      if(length(allCors_hasAmps) == 0){
        stop((paste("\nERROR: A gene was not mutated in at least", n, "patients. Recommend decreasing the n parameter.", sep=" ")))
      }else{
        pVals <- sapply(allCors_hasAmps, function(item)return(item[[1]]))
        betas <- sapply(allCors_hasAmps, function(item)return(item[[2]]))
        write.csv(pVals, file='./CnvTestOutput_pVals.csv')
        write.csv(betas, file='./CnvTestOutput_betas.csv')
      }

    }else{ #If non TCGA Mut...

      drug_prediction<-t(drug_prediction)

      #Obtain a list for the patients you have data for and initiate empty lists to fill.
      #_______________________________________
      sampIds<-data$Tumor_Sample_Barcode
      unique<-unique(sampIds)
      mutsListAll <- list() #A list of all the mutation occurring in each sample.
      proteinChangingMutations <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "In_Frame_Del",
                                    "Frame_Shift_Ins", "In_Frame_Ins", "Nonstop_Mutation", "De_novo_Start_OutOfFrame",
                                    "De_novo_Start_InFrame", "Missense", "Read-through", "Indel")
      genesWithProteinChangeList <- list()

      #Fill those lists.
      #_______________________________________
      for(i in 1:length(unique)){ #For each unique sample in your mut dataset...
        #print(paste(i, "of", length(unique)), sep='')
        indices<-unique[i] == sampIds #Find the index/row of that unique sample in all the samples (in the mut data)

        #Now for each sample, pull out a list of genes with somatic mutations
        variantType <- data[indices,"Variant_Classification"] #Get its variant classification...
        theGenes <- data[indices,"Hugo_Symbol"] #Get its gene
        names(variantType) <- theGenes
        mutsListAll[[i]] <- variantType #A list of all the mutation occurring in each sample.
        genesWithProteinChangeList[[i]] <- unique(names(variantType[variantType %in% proteinChangingMutations]))
      }

      allMutatedGenes <- unique(names(unlist(mutsListAll))) #All the unique genes mutated across all patients.
      mutationTypes <- table(unlist(mutsListAll)) #Frequency of mutation type across all patients.
      mutsListAll_unlist <- unlist(mutsListAll, recursive=F) #Mutation type and gene association.
      genesWithProteinChangeList_unlist <- unlist(genesWithProteinChangeList, recursive=F) #Genes with protein changing mutations across all patients.

      #From mutsListAll we can then create a matrix indicating if the gene has a coding mutation.
      #_______________________________________
      mutMat <- numeric((length(unique)*length(allMutatedGenes)))
      dim(mutMat) <- c(length(allMutatedGenes), length(unique)) #genes/muts/rows by samples.
      rownames(mutMat) <- allMutatedGenes
      colnames(mutMat) <- unique

      #Now populate this matrix with the relevant information about what kind of mutation each gene has in each sample.
      #_______________________________________
      for(i in 1:length(unique)){
        #print(paste(i, "of", length(unique)), sep='')
        mutMat[genesWithProteinChangeList[[i]], i] <- rep(1, length(genesWithProteinChangeList[[i]]))
      }

      #Identify samples that occur in drug and mutation data.
      #_______________________________________
      drug_samps<-rownames(drug_prediction)
      inPredAndMutData <- drug_samps[drug_samps %in% unique] #Samples for which we have both predicted drug response and mutation calls

      if (length(inPredAndMutData) == 0){
        stop((paste("\nERROR: Samples in drug_prediction and data don't overlap. Make sure their sample identifiers are of similar form")))
      }

      #Run the associations between all genes and drugs, for drugs with at least 50 mutations.
      #_______________________________________
      drug_prediction_filt_ord <- as.matrix(drug_prediction[inPredAndMutData,]) #The preds for the 01A samples we have both prediction and mutation data for.
      mutMat_nodups_ordFilt <- mutMat[, inPredAndMutData]
      commonMuts <- apply(mutMat_nodups_ordFilt, 1, sum)
      if (length(which(commonMuts >= n)) == 0){
        stop((paste("\nERROR: A gene was not mutated in at least", n, "patients. Recommend decreasing the n parameter.", sep=" ")))
      }
      #dim(mutMat_nodups_ordFilt) #genes by samples.
      commonlyMutated <- mutMat_nodups_ordFilt[which(commonMuts >= n), ]
      #NOTE: commonlyMutated will become a vector if there is only one mut that is >=n

      #dim(commonlyMutated) #3 208...genes/muts by samples.
      #rownames(commonlyMutated)
      #length(commonlyMutated)

      #If there are multiple genes, commonlyMutated will have dimensions.
      #Otherwise, it will be a vector, and you can calculate the length.
      #_______________________________________
      #length<-try(length(as.vector(unlist(commonlyMutated))), silent=TRUE) #If commonlyMutated is not a vector...(if there are multiple )
      nrow=try(nrow(commonlyMutated))
      if(is.null(nrow)){ #If we have a vector (one gene)...
        #Get p values and beta values.
        pValList <- list()
        betaValList <- list()

        suppressWarnings(for(i in 1:ncol(drug_prediction_filt_ord)){ #For each drug...
          pValList[[i]] <- numeric()
          betaValList[[i]] <- numeric()
          thecoefs <- coef(summary(lm(drug_prediction_filt_ord[,i]~as.vector(unlist(commonlyMutated)))))
          pValList[[i]]<- thecoefs[2,4]
          betaValList[[i]] <- thecoefs[2,1]
        })

        #Get the adjusted p-value for each gene-drug combination, pull out the significant associations
        #and create a supplementary table that lists these for "predictable" drugs?.
        sigPs <- list()
        pAdjListCantype <- list()
        for(i in 1:length(pValList)){
          names(pValList[[i]]) <- rownames(commonlyMutated)
          names(betaValList[[i]]) <- rownames(commonlyMutated)
          padj <- p.adjust(pValList[[i]], method="BH")
          sigPs[[i]] <- padj[padj < 0.05]
          pAdjListCantype[[i]] <- padj
        }

      }else{ #Otherwise, we have rows/multiple genes...

        #If there are gene entries with an unknown HUGO ID, remove it.
        #_______________________________________
        if("Unknown" %in% rownames(commonlyMutated)){
          indices<-'Unknown' %in% rownames(commonlyMutated)
          commonlyMutated<-commonlyMutated[-indices,]
        }

        #Get p values and beta values.
        #_______________________________________
        pValList <- list()
        betaValList <- list()
        #dim(drug_prediction_filt_ord) #samples by drugs.

        suppressWarnings(for(i in 1:ncol(drug_prediction_filt_ord)){ #For each drug...drug_prediction_filt_ord is samples/rows and drugs/columns. 208 11
          pValList[[i]] <- numeric()
          betaValList[[i]] <- numeric()
          for(j in 1:nrow(commonlyMutated)) #For each mut gene...commonlyMutated is mut genes/rows and samples/columns.3 208
          {
            thecoefs <- coef(summary(lm(drug_prediction_filt_ord[,i]~commonlyMutated[j,])))
            pValList[[i]][[j]] <- thecoefs[2,4]
            betaValList[[i]][[j]] <- thecoefs[2,1]
          }
        })

        #Get the adjusted p-value for each gene-drug combination, pull out the significant associations
        #and create a supplementary table that lists these for "predictable" drugs?.
        #_______________________________________
        sigPs <- list()
        pAdjListCantype <- list()
        for(i in 1:length(pValList))
        {
          #names(pValList[[i]]) <- rownames(commonlyMutated)
          #names(betaValList[[i]]) <- rownames(commonlyMutated)
          padj <- p.adjust(pValList[[i]], method="BH")
          sigPs[[i]] <- padj[padj < 0.05]
          pAdjListCantype[[i]] <- padj
        }
      }

      names(sigPs) <- colnames(drug_prediction_filt_ord)
      names(pValList) <- colnames(drug_prediction_filt_ord)
      names(betaValList) <- colnames(drug_prediction_filt_ord)
      names(pAdjListCantype) <- colnames(drug_prediction_filt_ord)

      #If you only have one gene of interest...make sure the row names are appropriate (drug:gene)
      nrow=try(nrow(commonlyMutated))
      if(is.null(nrow)){ #If we have a vector (one gene)...

        pVal<-unlist(pValList)
        betaVal<-unlist(betaValList)

        final_data<-cbind(pVal, betaVal)
        rows<-rownames(final_data)
        gene<-names(which(commonMuts >= n))
        final_rows<-paste(rows, ':', gene, sep='')
        rownames(final_data)<-final_rows
        write.csv(final_data, file='./MutationTestOutput_pVal_and_betaVal.csv')

      }else{

        drugs=names(pValList)
        genes=names(which(commonMuts >= n))
        genes=paste(genes, ":", sep="")
        mut_drug<-do.call(paste0,expand.grid(genes,drugs))
        pvalues=as.vector(unlist(pValList)) #All p values as a vector.
        #names(pvalues)=v1

        drugs=names(betaValList)
        genes=names(which(commonMuts >= n))
        genes=paste(genes, ":", sep="")
        #mut_drug<-do.call(paste0,expand.grid(genes,drugs))
        betavalues=as.vector(unlist(betaValList)) #All p values as a vector.

        df=data.frame(pvalues, betavalues)
        rownames(df)<-mut_drug

        write.csv(df, file='./MutationTestOutput_pVal_and_betaVal.csv')
      }
    }
  }
}
