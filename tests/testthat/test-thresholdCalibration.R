examplePeaks <- function() {
    readPeakSets(
        system.file("extdata",
            c(
                "rep1.narrowPeak", "rep2.narrowPeak",
                "rep3.narrowPeak"
            ),
            package = "consensusRegions"
        ),
        sampleNames = c("r1", "r2", "r3"),
        verbose = FALSE
    )
}


exampleCalibration <- function() {
    calibrateThreshold(examplePeaks(),
        nPermutations = 5, seed = 11,
        verbose = FALSE
    )
}


test_that("calibrateThreshold returns a ThresholdCalibration", {
    calibration <- exampleCalibration()

    expect_s4_class(calibration, "ThresholdCalibration")
    expect_true(methods::validObject(calibration))
})


test_that("the accessors return what the slots hold", {
    calibration <- exampleCalibration()

    expect_length(calibrationThreshold(calibration), 1L)
    expect_gt(calibrationThreshold(calibration), 0)
    expect_lte(calibrationThreshold(calibration), 1)

    expect_equal(
        calibrationThresholdNegLog10(calibration),
        -log10(calibrationThreshold(calibration))
    )
    expect_lte(
        calibrationAchievedFDR(calibration),
        calibrationTargetFDR(calibration)
    )

    expect_s3_class(calibrationCurve(calibration), "data.frame")
    expect_true(all(c("cut", "fdr") %in% names(calibrationCurve(calibration))))

    expect_equal(calibrationSeed(calibration), 11)
    expect_true(length(nullStatistics(calibration)) > 0)
    expect_true(length(observedStatistics(calibration)) > 0)
})


test_that("the dollar shorthand still reaches the slots", {
    calibration <- exampleCalibration()

    expect_identical(
        calibration$threshold,
        calibrationThreshold(calibration)
    )
    expect_identical(calibration$null, nullStatistics(calibration))
    expect_identical(calibration$seedSource, "seed")
    expect_error(calibration$notASlot, "no 'notASlot'")
})


test_that("an unseeded calibration keeps a NULL seed", {
    calibration <- calibrateThreshold(examplePeaks(),
        nPermutations = 5,
        verbose = FALSE
    )

    expect_null(calibrationSeed(calibration))
    expect_identical(calibration$seedSource, "unseeded")
})


test_that("validity rejects an impossible threshold", {
    expect_error(
        methods::new("ThresholdCalibration",
            threshold = 2,
            thresholdNegLog10 = -log10(2),
            achievedFDR = 0.01,
            targetFDR = 0.05,
            fdrCurve = data.frame(cut = 1, fdr = 0.01),
            observed = 1, null = 1,
            nPermutations = 5,
            seed = NULL, seedSource = "unseeded"
        ),
        "p-value in"
    )
})


test_that("validity rejects a curve without the columns it needs", {
    expect_error(
        methods::new("ThresholdCalibration",
            threshold = 0.01,
            thresholdNegLog10 = 2,
            achievedFDR = 0.01,
            targetFDR = 0.05,
            fdrCurve = data.frame(x = 1),
            observed = 1, null = 1,
            nPermutations = 5,
            seed = NULL, seedSource = "unseeded"
        ),
        "fdrCurve"
    )
})


test_that("buildConsensus refuses a calibration that is a bare list", {
    expect_error(
        buildConsensus(examplePeaks(),
            calibration = list(threshold = 0.01),
            verbose = FALSE
        ),
        "output of calibrateThreshold"
    )
})


test_that("a consensus built without calibration carries none", {
    result <- buildConsensus(examplePeaks(), verbose = FALSE)

    expect_null(calibrationResult(result))
})


test_that("a calibrated consensus keeps the calibration object", {
    calibration <- exampleCalibration()
    result <- buildConsensus(examplePeaks(),
        calibration = calibration,
        verbose = FALSE
    )

    expect_s4_class(calibrationResult(result), "ThresholdCalibration")
    expect_equal(calibrationSeed(calibrationResult(result)), 11)
})


test_that("the replicate summary is a DataFrame", {
    result <- buildConsensus(examplePeaks(), verbose = FALSE)

    expect_s4_class(consensusStats(result), "DataFrame")
    expect_identical(nrow(consensusStats(result)), 3L)
})
