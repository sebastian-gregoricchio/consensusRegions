## Regression tests for the inferential fixes. Each one starts from a
## case that the earlier implementation got wrong.

examplePeaks <- function() {
    readPeakSets(system.file("extdata",
                             c("rep1.narrowPeak", "rep2.narrowPeak",
                               "rep3.narrowPeak"),
                             package = "consensusRegions"),
                 sampleNames = c("r1", "r2", "r3"),
                 verbose = FALSE)
}


## Two replicates hold the same peaks and the third overlaps nothing, so
## no peak ever has more than one supporter.
lopsidedPeaks <- function() {
    shared <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(start = seq(1000, 100000, by = 5000),
                                 width = 500))
    apart <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(start = seq(500000, 599000, by = 5000),
                                 width = 500))
    readPeakSets(list(A = shared, B = shared, C = apart), verbose = FALSE)
}


test_that("a heavy weight cannot stand in for a missing replicate", {
    peaks <- lopsidedPeaks()
    weights <- computeReplicateWeights(peaks, method = "custom",
                                       weights = c(0.3, 2.4, 0.3))

    result <- suppressWarnings(
        buildConsensus(peaks, weights = weights, minSupport = 2,
                       verbose = FALSE))

    confirmed <- BiocGenerics::unlist(peakSets(result), use.names = FALSE)
    confirmed <- confirmed[
        S4Vectors::mcols(confirmed)$status == "confirmed"]

    ## replicate B carries a weight above minSupport on its own; that
    ## must not be enough to confirm a peak only B supports
    expect_true(all(S4Vectors::mcols(confirmed)$nSupport >= 2))
    expect_equal(length(confirmed), 0L)
})


test_that("minSupportWeight applies on top of the replicate count", {
    peaks <- examplePeaks()
    weights <- computeReplicateWeights(peaks, method = "custom",
                                       weights = c(1.5, 1.0, 0.5))

    without <- buildConsensus(peaks, weights = weights, verbose = FALSE)
    with <- suppressWarnings(
        buildConsensus(peaks, weights = weights, minSupportWeight = 2.4,
                       verbose = FALSE))

    ## an extra requirement can only remove regions, never add them
    expect_lte(length(with), length(without))
})


test_that("unequal weights do not unravel the recursive confirmation", {
    peaks <- examplePeaks()
    weights <- computeReplicateWeights(peaks, method = "custom",
                                       weights = c(0.3, 2.4, 0.3))

    recursive <- buildConsensus(peaks, weights = weights, recursive = TRUE,
                                verbose = FALSE)
    singlePass <- buildConsensus(peaks, weights = weights,
                                 recursive = FALSE, verbose = FALSE)

    ## eligibility now turns on reproducibility rather than on a peak's
    ## own score, so recursion settles instead of cascading to nothing
    expect_gt(length(recursive), 0)
    expect_equal(length(recursive), length(singlePass))
})


test_that("a region counts each replicate once however fragmented", {
    peaks <- examplePeaks()
    result <- buildConsensus(peaks, verbose = FALSE)
    regions <- consensusRanges(result)
    regionData <- S4Vectors::mcols(regions)

    ## the chained loci in the example data put five peaks per replicate
    ## inside one region, which is the case that used to inflate the score
    fragmented <- which(regionData$nPeaks > regionData$nReplicates)
    expect_gt(length(fragmented), 0)

    ## a fragmented region must not outscore what its replicates support
    perReplicate <- regionData$nPeaks / regionData$nReplicates
    correlation <- suppressWarnings(
        stats::cor(perReplicate, regionData$combinedNegLog10P,
                   use = "complete.obs", method = "spearman"))
    expect_lt(abs(correlation), 0.5)
})


test_that("the BH family defaults to every peak that was tested", {
    peaks <- examplePeaks()

    tested <- buildConsensus(peaks, verbose = FALSE)
    mspcStyle <- buildConsensus(peaks, adjustmentFamily = "confirmed",
                                verbose = FALSE)

    expect_identical(analysisParameters(tested)$adjustmentFamily, "tested")

    ## correcting across the whole family is the more conservative of the
    ## two, so it can never return more regions
    expect_lte(length(tested), length(mspcStyle))
})


test_that("calibration can choose a cut below the median statistic", {
    peaks <- examplePeaks()
    calibration <- calibrateThreshold(peaks, nPermutations = 5, seed = 4,
                                      verbose = FALSE)

    medianObserved <- stats::median(calibration$observed, na.rm = TRUE)
    expect_lt(min(calibration$fdrCurve$cut), medianObserved)
})


test_that("a finite permutation never reports an FDR of exactly zero", {
    peaks <- examplePeaks()
    calibration <- calibrateThreshold(peaks, nPermutations = 5, seed = 5,
                                      verbose = FALSE)

    expect_true(all(calibration$fdrCurve$expectedNull > 0))
    expect_true(all(calibration$fdrCurve$fdr > 0))

    ## the floor is one expected null peak in nPermutations + 1
    expect_gte(min(calibration$fdrCurve$expectedNull), 1 / 6)
})


test_that("calibration uses the same settings as the analysis", {
    peaks <- examplePeaks()

    ## recursive is TRUE by default in buildConsensus, so a calibration
    ## hard-coded to FALSE would be describing a different procedure
    expect_silent(
        calibration <- calibrateThreshold(peaks, nPermutations = 3,
                                          recursive = TRUE, seed = 6,
                                          verbose = FALSE))
    expect_true(calibration$threshold > 0)

    withFraction <- calibrateThreshold(peaks, nPermutations = 3,
                                       minOverlapFraction = 0.5, seed = 6,
                                       verbose = FALSE)
    expect_true(withFraction$threshold > 0)
})


test_that("minReplicates counts the peak's own replicate as MSPC does", {
    peaks <- examplePeaks()

    ## two of three replicates is one supporting replicate
    byTotal <- buildConsensus(peaks, minReplicates = 2, verbose = FALSE)
    bySupport <- buildConsensus(peaks, minSupport = 1, verbose = FALSE)
    expect_equal(length(byTotal), length(bySupport))

    ## all three is two supporting replicates
    allThree <- buildConsensus(peaks, minReplicates = 3, verbose = FALSE)
    twoOthers <- buildConsensus(peaks, minSupport = 2, verbose = FALSE)
    expect_equal(length(allThree), length(twoOthers))
    expect_lte(length(allThree), length(byTotal))
})


test_that("proportions and percentages resolve the same way", {
    ## rounded up, so two thirds of three replicates asks for two
    expect_equal(consensusRegions:::.parseReplicateCount(0.66, 3), 2L)
    expect_equal(consensusRegions:::.parseReplicateCount("66%", 3), 2L)
    expect_equal(consensusRegions:::.parseReplicateCount("75%", 8), 6L)
    expect_equal(consensusRegions:::.parseReplicateCount(0.75, 8), 6L)

    ## anything at or above one is taken literally
    expect_equal(consensusRegions:::.parseReplicateCount(5, 8), 5L)
    expect_equal(consensusRegions:::.parseReplicateCount("5", 8), 5L)

    ## a percentage stays a share even at the top of the range
    expect_equal(consensusRegions:::.parseReplicateCount("100%", 3), 3L)
    expect_equal(consensusRegions:::.parseReplicateCount("50%", 2), 1L)

    expect_error(consensusRegions:::.parseReplicateCount("150%", 3),
                 "not a percentage")
})


test_that("a percentage runs end to end", {
    peaks <- examplePeaks()

    ## 60% of three replicates rounds up to two, so one supporter
    byPercent <- buildConsensus(peaks, minReplicates = "60%",
                                verbose = FALSE)
    expect_equal(length(byPercent),
                 length(buildConsensus(peaks, minSupport = 1,
                                       verbose = FALSE)))
})


test_that("contradictory or impossible requirements are refused", {
    peaks <- examplePeaks()

    expect_error(
        buildConsensus(peaks, minSupport = 1, minReplicates = 2,
                       verbose = FALSE),
        "not both")

    expect_error(
        buildConsensus(peaks, minReplicates = 1, verbose = FALSE),
        "not a consensus")

    expect_error(
        buildConsensus(peaks, minReplicates = 5, verbose = FALSE),
        "5 replicates in total but only 3")

    expect_error(
        buildConsensus(peaks, minReplicates = "many", verbose = FALSE),
        "could not read")
})
