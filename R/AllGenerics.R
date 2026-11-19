# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' Accessors for ConsensusRegions objects
#'
#' @description
#' Retrieve the components of a [ConsensusRegions-class] object.
#'
#' @param object A [ConsensusRegions-class] object.
#'
#' @return
#' `consensusRanges()` returns a `GRanges` with the consensus regions.
#' `peakSets()` returns the annotated per-replicate peaks as a
#' `GRangesList`. `replicateWeights()` returns the named numeric vector of
#' weights. `analysisParameters()` returns the list of settings used.
#' `consensusStats()` returns a per-replicate summary `data.frame`.
#' `calibrationResult()` returns the [ThresholdCalibration-class] object
#' the consensus was built with, or `NULL` when it was built without one.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setGeneric
#'
#' @examples
#' peakFiles <- system.file("extdata",
#'     c(
#'         "rep1.narrowPeak", "rep2.narrowPeak",
#'         "rep3.narrowPeak"
#'     ),
#'     package = "consensusRegions"
#' )
#' peaks <- readPeakSets(peakFiles, sampleNames = c("r1", "r2", "r3"))
#' result <- buildConsensus(peaks, verbose = FALSE)
#'
#' consensusRanges(result)
#' consensusStats(result)
#'
#' @name ConsensusRegions-accessors
NULL


#' @rdname ConsensusRegions-accessors
#' @export
methods::setGeneric(
    "consensusRanges",
    function(object) standardGeneric("consensusRanges")
)

#' @rdname ConsensusRegions-accessors
#' @export
methods::setGeneric(
    "peakSets",
    function(object) standardGeneric("peakSets")
)

#' @rdname ConsensusRegions-accessors
#' @export
methods::setGeneric(
    "replicateWeights",
    function(object) standardGeneric("replicateWeights")
)

#' @rdname ConsensusRegions-accessors
#' @export
methods::setGeneric(
    "analysisParameters",
    function(object) standardGeneric("analysisParameters")
)

#' @rdname ConsensusRegions-accessors
#' @export
methods::setGeneric(
    "consensusStats",
    function(object) standardGeneric("consensusStats")
)

#' @rdname ConsensusRegions-accessors
#' @export
methods::setGeneric(
    "calibrationResult",
    function(object) standardGeneric("calibrationResult")
)


#' Accessors for ThresholdCalibration objects
#'
#' @description
#' Retrieve the components of a [ThresholdCalibration-class] object, the
#' result of [calibrateThreshold()].
#'
#' @param object A [ThresholdCalibration-class] object.
#' @param x A [ThresholdCalibration-class] object.
#' @param name Slot name, for the `$` shorthand.
#'
#' @return
#' `calibrationThreshold()` returns the combined p-value threshold, and
#' `calibrationThresholdNegLog10()` the same value on the -log10 scale.
#' `calibrationCurve()` returns the FDR curve as a `data.frame`.
#' `calibrationTargetFDR()` and `calibrationAchievedFDR()` return the FDR
#' that was asked for and the one that was reached. `calibrationSeed()`
#' returns the seed of the permutations, or `NULL`. `nullStatistics()` and
#' `observedStatistics()` return the statistics the curve was built from.
#'
#' `$` reaches the slots by name, which the result of this function
#' carried as a plain list up to version 0.99.2.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setGeneric
#'
#' @examples
#' peakFiles <- system.file("extdata",
#'     c(
#'         "rep1.narrowPeak", "rep2.narrowPeak",
#'         "rep3.narrowPeak"
#'     ),
#'     package = "consensusRegions"
#' )
#' peaks <- readPeakSets(peakFiles,
#'     sampleNames = c("r1", "r2", "r3"),
#'     verbose = FALSE
#' )
#' calibration <- calibrateThreshold(peaks,
#'     nPermutations = 5, seed = 42,
#'     verbose = FALSE
#' )
#'
#' calibrationThreshold(calibration)
#' calibrationSeed(calibration)
#'
#' @name ThresholdCalibration-accessors
NULL


#' @rdname ThresholdCalibration-accessors
#' @export
methods::setGeneric(
    "calibrationThreshold",
    function(object) standardGeneric("calibrationThreshold")
)

#' @rdname ThresholdCalibration-accessors
#' @export
methods::setGeneric(
    "calibrationThresholdNegLog10",
    function(object) {
        standardGeneric("calibrationThresholdNegLog10")
    }
)

#' @rdname ThresholdCalibration-accessors
#' @export
methods::setGeneric(
    "calibrationCurve",
    function(object) standardGeneric("calibrationCurve")
)

#' @rdname ThresholdCalibration-accessors
#' @export
methods::setGeneric(
    "calibrationTargetFDR",
    function(object) standardGeneric("calibrationTargetFDR")
)

#' @rdname ThresholdCalibration-accessors
#' @export
methods::setGeneric(
    "calibrationAchievedFDR",
    function(object) {
        standardGeneric("calibrationAchievedFDR")
    }
)

#' @rdname ThresholdCalibration-accessors
#' @export
methods::setGeneric(
    "calibrationSeed",
    function(object) standardGeneric("calibrationSeed")
)

#' @rdname ThresholdCalibration-accessors
#' @export
methods::setGeneric(
    "nullStatistics",
    function(object) standardGeneric("nullStatistics")
)

#' @rdname ThresholdCalibration-accessors
#' @export
methods::setGeneric(
    "observedStatistics",
    function(object) standardGeneric("observedStatistics")
)
