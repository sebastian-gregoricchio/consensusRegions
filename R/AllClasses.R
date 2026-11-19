# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.


#' ThresholdCalibration class
#'
#' @description
#' The result of [calibrateThreshold()]. It carries the threshold itself,
#' the FDR curve it was read off, the observed and null statistics behind
#' that curve, and the seed the permutations were drawn with, so that a
#' calibration can be inspected and repeated rather than merely trusted.
#'
#' Until version 0.99.2 this was a plain list, and `$` still reaches the
#' slots so that code written against it keeps working.
#'
#' @slot threshold The combined p-value threshold.
#' @slot thresholdNegLog10 The same threshold on the -log10 scale.
#' @slot achievedFDR The empirical FDR reached at that threshold.
#' @slot targetFDR The FDR that was asked for.
#' @slot fdrCurve A `data.frame` with one row per candidate cut, holding
#'   `cut`, `nObserved`, `expectedNull` and `fdr`.
#' @slot observed The observed combined statistics.
#' @slot null The combined statistics of the shuffled peak sets.
#' @slot nPermutations How many shuffles the null came from.
#' @slot seed The seed the permutations were drawn with, or `NULL`.
#' @slot seedSource `"seed"` when the seed came from the argument,
#'   `"BPPARAM"` when it came from the backend, `"unseeded"` when there
#'   was none.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setClass setValidity setClassUnion
#'
#' @name ThresholdCalibration-class
#' @rdname ThresholdCalibration-class
#' @exportClass ThresholdCalibration
methods::setClass(
    "ThresholdCalibration",
    slots = c(
        threshold = "numeric",
        thresholdNegLog10 = "numeric",
        achievedFDR = "numeric",
        targetFDR = "numeric",
        fdrCurve = "data.frame",
        observed = "numeric",
        null = "numeric",
        nPermutations = "numeric",
        seed = "ANY",
        seedSource = "character"
    )
)


## buildConsensus() reads the threshold out of this object and thresholds
## every peak against it, so a missing or impossible value here would
## surface much later as an empty consensus rather than as an error
methods::setValidity("ThresholdCalibration", function(object) {
    problems <- character(0)

    counts <- c(
        threshold = length(object@threshold),
        thresholdNegLog10 = length(object@thresholdNegLog10),
        achievedFDR = length(object@achievedFDR),
        targetFDR = length(object@targetFDR)
    )
    if (any(counts != 1)) {
        problems <- c(
            problems,
            paste(
                "these must hold exactly one value:",
                paste(names(counts)[counts != 1],
                    collapse = ", "
                )
            )
        )
    }

    if (counts[["threshold"]] == 1 && (is.na(object@threshold) ||
        object@threshold <= 0 ||
        object@threshold > 1)) {
        problems <- c(problems, "'threshold' must be a p-value in (0, 1]")
    }
    if (counts[["targetFDR"]] == 1 && (is.na(object@targetFDR) ||
        object@targetFDR <= 0 ||
        object@targetFDR >= 1)) {
        problems <- c(
            problems,
            "'targetFDR' must lie strictly between 0 and 1"
        )
    }
    if (!all(c("cut", "fdr") %in% names(object@fdrCurve))) {
        problems <- c(
            problems,
            "'fdrCurve' must have 'cut' and 'fdr' columns"
        )
    }

    ## NULL is the unseeded case and has to stay representable
    if (!is.null(object@seed) && (!is.numeric(object@seed) ||
        length(object@seed) != 1)) {
        problems <- c(problems, "'seed' must be a single number, or NULL")
    }
    if (length(object@seedSource) != 1 ||
        !object@seedSource %in% c("seed", "BPPARAM", "unseeded")) {
        problems <- c(
            problems,
            "'seedSource' must be 'seed', 'BPPARAM' or 'unseeded'"
        )
    }

    if (length(problems) == 0) TRUE else problems
})


## a consensus built without calibration has none, and NULL says that
## more plainly than an empty object would
methods::setClassUnion(
    "ThresholdCalibrationOrNULL",
    c("ThresholdCalibration", "NULL")
)


#' ConsensusRegions class
#'
#' @description
#' Container for the result of [buildConsensus()]. It keeps the peaks as
#' they were submitted, the per-peak verdict of the combination step, the
#' consensus regions themselves, and the settings used to get there, so
#' that an analysis can be traced back from its output.
#'
#' @slot peaks A `GRangesList`, one element per replicate, holding the
#'   peaks after classification. Each element carries the metadata columns
#'   `negLog10P`, `class`, `nSupport`, `supportWeight`,
#'   `combinedNegLog10P`, `combinedNegLog10Padj` and `status`.
#' @slot consensus A `GRanges` of consensus regions.
#' @slot weights A named numeric vector of replicate weights.
#' @slot parameters A list with the call arguments.
#' @slot stats A `DataFrame` summarising each replicate.
#' @slot calibration A [ThresholdCalibration-class] object when the
#'   consensus was built with one, `NULL` otherwise.
#'
#' @author Sebastian Gregoricchio
#'
#' @importClassesFrom GenomicRanges GRanges GRangesList
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom methods setClass validObject new is prototype
#'
#' @name ConsensusRegions-class
#' @rdname ConsensusRegions-class
#' @exportClass ConsensusRegions
methods::setClass(
    "ConsensusRegions",
    slots = c(
        peaks = "GRangesList",
        consensus = "GRanges",
        weights = "numeric",
        parameters = "list",
        stats = "DataFrame",
        calibration = "ThresholdCalibrationOrNULL"
    ),
    prototype = methods::prototype(calibration = NULL)
)


## the checks below catch objects assembled by hand rather than by
## buildConsensus, where slots can easily fall out of sync
methods::setValidity("ConsensusRegions", function(object) {
    problems <- character(0)

    ## every replicate must be named, since weights and stats are matched
    ## to the peak list by name and not by position
    replicateNames <- names(object@peaks)
    if (is.null(replicateNames) || any(replicateNames == "")) {
        problems <- c(problems, "all elements of 'peaks' must be named")
    }

    ## weights are looked up by replicate name during the combination
    if (length(object@weights) > 0 && !is.null(replicateNames)) {
        if (!setequal(names(object@weights), replicateNames)) {
            problems <- c(
                problems,
                "names of 'weights' must match names of 'peaks'"
            )
        }
        if (any(object@weights <= 0, na.rm = TRUE)) {
            problems <- c(problems, "'weights' must be strictly positive")
        }
    }

    if (length(problems) == 0) TRUE else problems
})
