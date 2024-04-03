# oncoPredict 
(Predict Response from Expression Data and Identify Cell line/Clinical Targets and Trends)

Additional details about this package can be found in our publication [oncoPredict: an R package for predicting in vivo or cancer patient drug response and biomarkers from cell line screening data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8574972/)

An R package for drug response prediction and drug-gene association prediction. The prepared GDSC and CTRP matrices for the calcPhenotype() are located in the oncoPredict OSF.
 *  For drug response prediction, use **calcPhenotype**. 
 *  For pre-clinical biomarker discovery, use **GLDS**. 
 * For clinical biomarker discovery, use **IDWAS** (for CNV or somatic mutation association with drug response) or indicate **cc=TRUE** (for gene expression association with drug response) in calcPhenotype(). 
 * The link to updated CCLE gene expression data is found at [depmap](https://depmap.org/portal/download/). We provide GDSC1/GDSC2 pre-processed expression and response data, as well as CTRP response data and depmap's CCLE expression data (18Q2) [here](https://osf.io/c6tfx/).
 
## R <h2>
 * This directory contains all the R functions included in this package. 

## vignettes <h2> 
  *  This directory contains vignettes which display detailed examples of the functionalities available in this package.
  *  **IDWAS** This directory contains examples of IDWAS code application for clinical drug-gene association prediction. 
      + **cnv.Rmd** Example as to how to download CNV (copy number variation) data from the GDC database, then apply map_cnv() and idwas().
      + **mut.Rmd** Example as to how to download stomatic mutation data from the GDC database, then apply idwas(). 

  * **GLDS** This directory contains examples of GLDS code application for pre-clinical drug-gene association prediction. 
      + **glds_GDSC.Rmd** Example of GLDS application to GDSC data.  

  * **calcPhenotype.Rmd** Example of calcPhenotype() application.

## man <h2>
 * This directory contains .Rd (R documentation) files for each function. These files were automatically generated upon creation of the package. 

## NAMESPACE <h2>
 * This file lists the functions to be imported and exported from this package. 

## DESCRIPTION <h2>
 * This file contains the description documentation and metadata for this package.
 * Dependencies and packages recommended for oncoPredict are listed here. 
  
## Figure 1. 
Flowchart displaying the 3 primary functionalities available through oncoPredict (calcPhenotype, GLDS, IDWAS) as well as the files generated from each function and parameters. Functions and files generated are bold.

![Figure 1_ Overview of pRRophetic_plus (a flow chart or similar diagram to highlight the package’s abilities)   (2)](https://user-images.githubusercontent.com/62571435/114970102-5d471580-9e3f-11eb-8734-a5e40a3d7f41.jpg)



