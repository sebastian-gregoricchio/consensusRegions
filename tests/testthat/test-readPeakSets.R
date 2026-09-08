exampleFiles <- function(which = c("rep1.narrowPeak", "rep2.narrowPeak",
                                   "rep3.narrowPeak")) {
    system.file("extdata", which, package = "consensusRegions")
}


test_that("narrowPeak input is recognised as carrying p-values", {
    peaks <- readPeakSets(exampleFiles(),
                          sampleNames = c("r1", "r2", "r3"),
                          verbose = FALSE)

    expect_s4_class(peaks, "GRangesList")
    expect_length(peaks, 3)
    expect_identical(S4Vectors::metadata(peaks)$scoreType, "log10pvalue")
    expect_true("negLog10P" %in% names(S4Vectors::mcols(peaks[[1]])))
    expect_true(all(S4Vectors::mcols(peaks[[1]])$negLog10P > 0))
})


test_that("BED3 input falls back to no statistic at all", {
    minimal <- system.file("extdata", "rep1_minimal.bed",
                           package = "consensusRegions")

    peaks <- readPeakSets(c(minimal, minimal),
                          sampleNames = c("a", "b"),
                          verbose = FALSE)

    expect_identical(S4Vectors::metadata(peaks)$scoreType, "none")
    expect_true(all(is.na(S4Vectors::mcols(peaks[[1]])$negLog10P)))
})


test_that("a score column is rank-transformed onto the p-value scale", {
    peaks <- readPeakSets(exampleFiles(),
                          sampleNames = c("r1", "r2", "r3"),
                          scoreType = "score",
                          verbose = FALSE)

    transformed <- S4Vectors::mcols(peaks[[1]])$negLog10P

    expect_true(all(transformed > 0))
    ## the best peak of the replicate must also be the most significant
    expect_equal(which.max(transformed),
                 which.max(S4Vectors::mcols(peaks[[1]])$score))
})


test_that("a single replicate is refused", {
    expect_error(readPeakSets(exampleFiles()[1], verbose = FALSE),
                 "at least two replicates")
})


test_that("mismatched or duplicated names are caught", {
    expect_error(
        readPeakSets(exampleFiles(), sampleNames = c("a", "b"),
                     verbose = FALSE),
        "'sampleNames' has length")

    expect_error(
        readPeakSets(exampleFiles(), sampleNames = c("a", "a", "b"),
                     verbose = FALSE),
        "must be unique")
})


test_that("a missing file is reported by name", {
    expect_error(
        readPeakSets(c(exampleFiles()[1], "no_such_file.narrowPeak"),
                     verbose = FALSE),
        "file not found")
})


test_that("a list of GRanges is accepted directly", {
    peaks <- readPeakSets(exampleFiles(),
                          sampleNames = c("r1", "r2", "r3"),
                          verbose = FALSE)

    fromList <- readPeakSets(as.list(peaks), verbose = FALSE)

    expect_length(fromList, 3)
    expect_identical(names(fromList), c("r1", "r2", "r3"))
})
