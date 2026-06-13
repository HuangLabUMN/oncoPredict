test_that("glds skips unrelated drugs without mutating loop indices", {
  set.seed(2)
  sample_names <- paste0("cell", seq_len(14))
  drug_names <- paste0("drug", seq_len(12))

  drugMat <- matrix(rnorm(length(sample_names) * length(drug_names)),
                    nrow=length(sample_names))
  rownames(drugMat) <- sample_names
  colnames(drugMat) <- drug_names

  drugRelatedness <- cbind(
    drug=drug_names[-2],
    pathway=paste0("path", seq_along(drug_names[-2]))
  )

  markerMat <- matrix(c(rep(c(0, 1), 7), rep(c(1, 0), 7)), nrow=2,
                      byrow=TRUE)
  rownames(markerMat) <- c("marker1", "marker2")
  colnames(markerMat) <- sample_names

  capture.output(
    output <- suppressMessages(
      glds(drugMat, drugRelatedness, markerMat, minMuts=1, threshold=1)
    )
  )

  expect_false("drug2" %in% names(output$pGlds))
  expect_equal(names(output$pGlds), drug_names[-2])
  expect_equal(names(output$pNaive), drug_names[-2])
  expect_named(output$pGlds$drug3, rownames(markerMat))
})

test_that("glds initializes expression gene signature list", {
  set.seed(3)
  sample_names <- paste0("cell", seq_len(14))
  drug_names <- paste0("drug", seq_len(12))

  drugMat <- matrix(rnorm(length(sample_names) * length(drug_names)),
                    nrow=length(sample_names))
  rownames(drugMat) <- sample_names
  colnames(drugMat) <- drug_names

  drugRelatedness <- cbind(
    drug=drug_names,
    pathway=paste0("path", seq_along(drug_names))
  )

  markerMat <- matrix(c(rep(c(0, 1), 7), rep(c(1, 0), 7)), nrow=2,
                      byrow=TRUE)
  rownames(markerMat) <- c("marker1", "marker2")
  colnames(markerMat) <- sample_names

  expression <- matrix(rnorm(60 * length(sample_names)), nrow=60)
  rownames(expression) <- paste0("gene", seq_len(nrow(expression)))
  colnames(expression) <- sample_names

  expect_no_error(
    capture.output(
      glds(drugMat, drugRelatedness, markerMat, minMuts=1,
           expression=expression, threshold=1)
    )
  )
})

test_that("glds uses available control PCs when fewer than 10 exist", {
  set.seed(4)
  sample_names <- paste0("cell", seq_len(14))
  drug_names <- paste0("drug", seq_len(6))

  drugMat <- matrix(rnorm(length(sample_names) * length(drug_names)),
                    nrow=length(sample_names))
  rownames(drugMat) <- sample_names
  colnames(drugMat) <- drug_names

  drugRelatedness <- cbind(
    drug=drug_names,
    pathway=paste0("path", seq_along(drug_names))
  )

  markerMat <- matrix(c(rep(c(0, 1), 7), rep(c(1, 0), 7)), nrow=2,
                      byrow=TRUE)
  rownames(markerMat) <- c("marker1", "marker2")
  colnames(markerMat) <- sample_names

  expect_no_error(
    capture.output(
      glds(drugMat, drugRelatedness, markerMat, minMuts=1, threshold=1)
    )
  )
})

test_that("glds gives a clear error when drug and marker samples do not overlap", {
  drugMat <- matrix(rnorm(14 * 6), nrow=14)
  rownames(drugMat) <- paste0("drug_sample", seq_len(14))
  colnames(drugMat) <- paste0("drug", seq_len(6))

  drugRelatedness <- cbind(
    drug=colnames(drugMat),
    pathway=paste0("path", seq_len(ncol(drugMat)))
  )

  markerMat <- matrix(c(rep(c(0, 1), 7), rep(c(1, 0), 7)), nrow=2,
                      byrow=TRUE)
  rownames(markerMat) <- c("marker1", "marker2")
  colnames(markerMat) <- paste0("marker_sample", seq_len(14))

  expect_error(
    glds(drugMat, drugRelatedness, markerMat, minMuts=1, threshold=1),
    "No overlapping samples"
  )
})
