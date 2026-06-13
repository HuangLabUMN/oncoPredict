test_that("idwas TCGA CNV branch uses filtered 01A drug predictions", {
  drug_prediction <- data.frame(
    drugA=c(-10, 10, 0),
    row.names=c(
      "TCGA-AA-ZZZZ-11A",
      "TCGA-AA-0002-01A",
      "TCGA-AA-0001-01A"
    )
  )

  cnv_data <- data.frame(
    `TCGA-AA-0001-01A`=c(0, 0),
    `TCGA-AA-0002-01A`=c(2, 0),
    row.names=c("geneA", "geneB"),
    check.names=FALSE
  )

  output <- idwas(drug_prediction, cnv_data, n=0, cnv=TRUE, folder=FALSE)

  expect_gt(output$betas["geneA.drugA"], 0)
})

test_that("idwas TCGA mutation branch handles one commonly mutated gene", {
  samples <- paste0("TCGA-AA-000", 1:4, "-01A")
  drug_prediction <- data.frame(
    drugA=c(0, 1, 10, 11),
    row.names=samples
  )

  mutation_data <- data.frame(
    Tumor_Sample_Barcode=c(samples, samples[3], samples[1:2]),
    Variant_Classification=c(
      "Silent",
      "Silent",
      "Missense_Mutation",
      "Missense_Mutation",
      "Missense_Mutation",
      "Missense_Mutation",
      "Missense_Mutation"
    ),
    Hugo_Symbol=c("noise1", "noise2", "geneA", "geneA", "geneB",
                  "Unknown", "Unknown")
  )

  output <- idwas(drug_prediction, mutation_data, n=2, cnv=FALSE, folder=FALSE)

  expect_true("drugA.geneA" %in% rownames(output))
  expect_gt(output["drugA.geneA", "betaVal"], 0)
  expect_false(any(grepl("Unknown", rownames(output), fixed=TRUE)))
})

test_that("idwas non-TCGA mutation branch removes Unknown rows by name", {
  samples <- paste0("sample", seq_len(4))
  drug_prediction <- data.frame(
    drugA=c(0, 1, 10, 11),
    row.names=samples
  )

  mutation_data <- data.frame(
    Tumor_Sample_Barcode=c(samples[3:4], samples[1:2]),
    Variant_Classification=rep("Missense_Mutation", 4),
    Hugo_Symbol=c("geneA", "geneA", "Unknown", "Unknown")
  )

  output <- idwas(drug_prediction, mutation_data, n=2, cnv=FALSE, folder=FALSE)

  expect_true("geneA:drugA" %in% rownames(output))
  expect_false(any(grepl("Unknown", rownames(output), fixed=TRUE)))
})
