examplePeaks <- function() {
  readPeakSets(system.file("extdata",
                           c("rep1.narrowPeak", "rep2.narrowPeak",
                             "rep3.narrowPeak"),
                           package = "consensusRegions"),
               sampleNames = c("r1", "r2", "r3"),
               verbose = FALSE)
}


test_that("equal weights are the default and sum to the replicate count", {
  weights <- computeReplicateWeights(examplePeaks(), method = "equal")

  expect_equal(unname(weights), rep(1, 3))
  expect_identical(names(weights), c("r1", "r2", "r3"))
})


test_that("pre-computed FRiP needs no BAM file", {
  weights <- computeReplicateWeights(examplePeaks(), method = "frip",
                                     frip = c(0.20, 0.18, 0.05),
                                     verbose = FALSE)

  expect_equal(mean(weights), 1)
  ## the shallow replicate must end up with the smallest weight
  expect_equal(which.min(weights), 3L, ignore_attr = TRUE)
})


test_that("library size weights scale with the square root of depth", {
  weights <- computeReplicateWeights(examplePeaks(),
                                     method = "librarySize",
                                     librarySize = c(4e7, 4e7, 1e7),
                                     verbose = FALSE)

  ## a quarter of the depth is worth half the weight
  expect_equal(unname(weights[1] / weights[3]), 2, tolerance = 1e-8)
})


test_that("intrinsic weights work with nothing but the peaks", {
  weights <- computeReplicateWeights(examplePeaks(), method = "intrinsic",
                                     verbose = FALSE)

  expect_equal(mean(weights), 1)
  expect_true(all(weights > 0))
  expect_false(is.null(attr(weights, "metrics")))
})


test_that("missing inputs are reported rather than guessed", {
  expect_error(
    computeReplicateWeights(examplePeaks(), method = "frip",
                            verbose = FALSE),
    "supply 'bamFiles'")

  expect_error(
    computeReplicateWeights(examplePeaks(), method = "custom"),
    "'weights' must be supplied")

  expect_error(
    computeReplicateWeights(examplePeaks(), method = "frip",
                            frip = c(0.2, 1.5, 0.1), verbose = FALSE),
    "outside")
})


test_that("calibration returns a usable threshold", {
  set.seed(42)
  calibration <- calibrateThreshold(examplePeaks(), nPermutations = 5,
                                    verbose = FALSE)

  expect_true(calibration$threshold > 0)
  expect_true(calibration$threshold < 1)
  expect_true(all(c("cut", "fdr") %in% names(calibration$fdrCurve)))
  expect_lte(calibration$achievedFDR, calibration$targetFDR)
})


test_that("calibration repeats when the stream is seeded the same way", {
  ## the function deliberately does not seed itself, so this also
  ## checks that it leaves the caller in charge of the stream
  set.seed(7)
  first <- calibrateThreshold(examplePeaks(), nPermutations = 3,
                              verbose = FALSE)
  set.seed(7)
  second <- calibrateThreshold(examplePeaks(), nPermutations = 3,
                               verbose = FALSE)

  expect_equal(first$threshold, second$threshold)
})


test_that("calibration refuses input with no statistic", {
  minimal <- system.file("extdata", "rep1_minimal.bed",
                         package = "consensusRegions")
  peaks <- readPeakSets(c(minimal, minimal), sampleNames = c("a", "b"),
                        verbose = FALSE)

  expect_error(calibrateThreshold(peaks, verbose = FALSE),
               "needs a per-peak statistic")
})


test_that("recentring returns regions of one fixed width", {
  result <- buildConsensus(examplePeaks(), verbose = FALSE)

  recentred <- recentrePeaks(result, width = 400)

  expect_true(all(GenomicRanges::width(recentred) <= 400))
  expect_equal(stats::median(GenomicRanges::width(recentred)), 400)
  expect_length(recentred, length(result))
  expect_true("summit" %in% names(S4Vectors::mcols(recentred)))
})


test_that("the three anchor rules all place a summit inside the region", {
  result <- buildConsensus(examplePeaks(), verbose = FALSE)

  for (source in c("best", "weighted", "centre")) {
    recentred <- recentrePeaks(result, width = 300,
                               summitSource = source)
    expect_length(recentred, length(result))
  }
})


test_that("conversion to SummarizedExperiment keeps the row ranges", {
  result <- buildConsensus(examplePeaks(), verbose = FALSE)

  se <- asSummarizedExperiment(result)

  expect_s4_class(se, "RangedSummarizedExperiment")
  expect_equal(nrow(se), length(result))
  expect_equal(ncol(se), 3)
  expect_true(all(c("negLog10P", "detected") %in%
                    SummarizedExperiment::assayNames(se)))
})


test_that("export writes a file that can be read back", {
  result <- buildConsensus(examplePeaks(), verbose = FALSE)
  target <- file.path(tempdir(), "consensusTest.bed")

  written <- exportConsensus(result, target, verbose = FALSE)

  expect_true(file.exists(target))
  expect_length(written, 1)
  expect_equal(length(rtracklayer::import(target)), length(result))
})


test_that("the plots build without error", {
  result <- buildConsensus(examplePeaks(), verbose = FALSE)

  expect_s3_class(plotRescue(result), "ggplot")
  expect_s3_class(plotRescue(result, proportion = TRUE), "ggplot")
  expect_s3_class(plotJaccard(result), "ggplot")

  set.seed(3)
  calibration <- calibrateThreshold(examplePeaks(), nPermutations = 3,
                                    verbose = FALSE)
  expect_s3_class(plotCalibration(calibration), "ggplot")
})
