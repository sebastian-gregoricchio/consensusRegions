examplePeaks <- function(...) {
    readPeakSets(system.file("extdata",
                             c("rep1.narrowPeak", "rep2.narrowPeak",
                               "rep3.narrowPeak"),
                             package = "consensusRegions"),
                 sampleNames = c("r1", "r2", "r3"),
                 verbose = FALSE, ...)
}


test_that("a plain run returns a populated object", {
    result <- buildConsensus(examplePeaks(), verbose = FALSE)

    expect_s4_class(result, "ConsensusRegions")
    expect_gt(length(result), 0)
    expect_identical(nrow(consensusStats(result)), 3L)
    expect_true(all(c("nReplicates", "nPeaks", "combinedNegLog10P") %in%
                        names(S4Vectors::mcols(consensusRanges(result)))))
})


test_that("weak peaks with agreement are rescued", {
    result <- buildConsensus(examplePeaks(), verbose = FALSE)
    statistics <- consensusStats(result)

    ## the example data was built with a set of positions that only ever
    ## reach the weak tier, so some rescue must happen
    expect_gt(sum(statistics$nRescued), 0)
    expect_true(all(statistics$nRescued <= statistics$nWeak))
})


test_that("asking for more support yields fewer regions", {
    permissive <- buildConsensus(examplePeaks(), minSupport = 1,
                                 verbose = FALSE)
    strict <- buildConsensus(examplePeaks(), replicateType = "technical",
                             verbose = FALSE)

    expect_lte(length(strict), length(permissive))
})


test_that("a stricter combined threshold yields fewer regions", {
    loose <- buildConsensus(examplePeaks(), combinedThreshold = 1e-6,
                            verbose = FALSE)
    tight <- buildConsensus(examplePeaks(), combinedThreshold = 1e-20,
                            verbose = FALSE)

    expect_lt(length(tight), length(loose))
})


test_that("down-weighting a replicate is carried through the analysis", {
    peaks <- examplePeaks()

    weighted <- buildConsensus(
        peaks,
        weights = computeReplicateWeights(peaks, method = "custom",
                                          weights = c(1, 1, 0.2)),
        combinationMethod = "stouffer",
        verbose = FALSE)

    storedWeights <- replicateWeights(weighted)

    expect_equal(mean(storedWeights), 1)
    expect_lt(storedWeights[["r3"]], storedWeights[["r1"]])
    expect_gt(length(weighted), 0)
})


test_that("every combination scheme runs end to end", {
    for (method in c("fisher", "stouffer", "lancaster")) {
        result <- buildConsensus(examplePeaks(),
                                 combinationMethod = method,
                                 verbose = FALSE)
        expect_gt(length(result), 0)
    }

    ## The rank product sits on its own scale and needs its own cut. On a
    ## peak set this small its ranks are also too coarse to survive a
    ## correction across the whole family, so this checks it against the
    ## MSPC-style family instead.
    result <- buildConsensus(examplePeaks(),
                             combinationMethod = "rankProduct",
                             combinedThreshold = 0.05,
                             adjustmentFamily = "confirmed",
                             verbose = FALSE)
    expect_gt(length(result), 0)
})


test_that("an unreachable rank-product threshold is refused up front", {
    expect_error(
        buildConsensus(examplePeaks(), combinationMethod = "rankProduct",
                       combinedThreshold = 1e-8, verbose = FALSE),
        "cannot reach a combined p-value")
})


test_that("BED3 input takes the presence route", {
    minimal <- system.file("extdata", "rep1_minimal.bed",
                           package = "consensusRegions")
    peaks <- readPeakSets(c(minimal, minimal, minimal),
                          sampleNames = c("a", "b", "c"),
                          verbose = FALSE)

    result <- buildConsensus(peaks, verbose = FALSE)

    expect_true(analysisParameters(result)$presenceOnly)
    expect_true(is.na(analysisParameters(result)$combinationMethod))
    ## the three sets are identical, so everything is supported
    expect_gt(length(result), 0)
})


test_that("the iterative merge keeps regions no wider than reduce", {
    peaks <- examplePeaks()

    reduced <- buildConsensus(peaks, mergeMethod = "reduce",
                              verbose = FALSE)
    seeded <- buildConsensus(peaks, mergeMethod = "iterative",
                             verbose = FALSE)

    expect_lte(max(GenomicRanges::width(consensusRanges(seeded))),
               max(GenomicRanges::width(consensusRanges(reduced))))
})


test_that("thresholds the wrong way round are refused", {
    expect_error(
        buildConsensus(examplePeaks(), stringencyThreshold = 1e-4,
                       weakThreshold = 1e-8, verbose = FALSE),
        "more permissive")
})


test_that("asking for more support than exists is refused", {
    expect_error(
        buildConsensus(examplePeaks(), minSupport = 5, verbose = FALSE),
        "supporting replicates")
})


test_that("excluded regions are removed from the consensus", {
    result <- buildConsensus(examplePeaks(), verbose = FALSE)
    blacklist <- consensusRanges(result)[seq_len(10)]

    filtered <- buildConsensus(examplePeaks(), excludeRegions = blacklist,
                               verbose = FALSE)

    expect_lt(length(filtered), length(result))
    expect_equal(
        length(IRanges::subsetByOverlaps(consensusRanges(filtered),
                                         blacklist)),
        0)
})


test_that("a single-pass run can differ from the recursive one", {
    recursive <- buildConsensus(examplePeaks(), recursive = TRUE,
                                verbose = FALSE)
    singlePass <- buildConsensus(examplePeaks(), recursive = FALSE,
                                 verbose = FALSE)

    ## recursion can only remove support, never add it
    expect_lte(length(recursive), length(singlePass))
})
