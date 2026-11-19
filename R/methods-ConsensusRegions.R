# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @rdname ConsensusRegions-accessors
#' @importFrom methods setMethod
#' @export
methods::setMethod("consensusRanges", "ConsensusRegions", function(object) {
    object@consensus
})


#' @rdname ConsensusRegions-accessors
#' @export
methods::setMethod("peakSets", "ConsensusRegions", function(object) {
    object@peaks
})


#' @rdname ConsensusRegions-accessors
#' @export
methods::setMethod("replicateWeights", "ConsensusRegions", function(object) {
    object@weights
})


#' @rdname ConsensusRegions-accessors
#' @export
methods::setMethod(
    "analysisParameters", "ConsensusRegions",
    function(object) {
        object@parameters
    }
)


#' @rdname ConsensusRegions-accessors
#' @export
methods::setMethod("consensusStats", "ConsensusRegions", function(object) {
    object@stats
})


#' @rdname ConsensusRegions-accessors
#' @export
methods::setMethod(
    "calibrationResult", "ConsensusRegions",
    function(object) {
        object@calibration
    }
)


#' Replace the consensus regions of an object
#'
#' @param object A [ConsensusRegions-class] object.
#' @param value A `GRanges` of consensus regions.
#'
#' @return The object with its consensus replaced.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods validObject
#'
#' @keywords internal
#' @noRd
.setConsensus <- function(object, value) {
    object@consensus <- value
    methods::validObject(object)
    object
}


#' Record one setting in an object
#'
#' @param object A [ConsensusRegions-class] object.
#' @param name Name of the setting.
#' @param value Value to store under that name.
#'
#' @return The object with the setting added.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal
#' @noRd
.setParameter <- function(object, name, value) {
    object@parameters[[name]] <- value
    object
}


#' Number of consensus regions
#'
#' @param x A [ConsensusRegions-class] object.
#'
#' @return Integer, the number of consensus regions.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setMethod
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
#' length(buildConsensus(peaks, verbose = FALSE))
#'
#' @rdname ConsensusRegions-length
#' @export
methods::setMethod("length", "ConsensusRegions", function(x) {
    length(x@consensus)
})


#' Display a ConsensusRegions object
#'
#' @param object A [ConsensusRegions-class] object.
#'
#' @return Invisibly `NULL`, called for the printed output.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setMethod show
#' @importFrom GenomicRanges width
#' @importFrom stats median
#'
#' @rdname ConsensusRegions-show
#' @export
methods::setMethod("show", "ConsensusRegions", function(object) {
    settings <- object@parameters

    cat("ConsensusRegions\n")
    cat(
        "  replicates      :", length(object@peaks), "(",
        paste(names(object@peaks), collapse = ", "), ")\n"
    )
    cat("  score type      :", settings$scoreType, "\n")

    if (isTRUE(settings$presenceOnly)) {
        cat("  combination     : weighted presence rule\n")
    } else {
        cat("  combination     :", settings$combinationMethod, "\n")
        cat("  combined cut    :", settings$combinedThreshold, "\n")
    }

    cat(
        "  weights         :",
        paste(names(object@weights), round(object@weights, 3),
            sep = "=", collapse = ", "
        ), "\n"
    )
    cat("  consensus       :", length(object@consensus), "regions\n")

    ## the width distribution is the first thing to check after a run,
    ## since runaway chaining shows up here before anywhere else
    if (length(object@consensus) > 0) {
        regionWidth <- GenomicRanges::width(object@consensus)
        cat("  width (median)  :", stats::median(regionWidth), "bp\n")
        cat("  width (max)     :", max(regionWidth), "bp\n")
    }

    calibration <- calibrationResult(object)
    if (!is.null(calibration)) {
        ## calibrations made before the seed was recorded carry none
        calibrationSeed <- calibration$seed
        cat(
            "  calibrated at   : FDR", calibration$targetFDR,
            if (is.null(calibrationSeed)) {
                "(no seed)"
            } else {
                paste0("(seed ", calibrationSeed, ")")
            }, "\n"
        )
    }

    invisible(NULL)
})


## ---------------------------------------------------------------------
## ThresholdCalibration
## ---------------------------------------------------------------------

#' @rdname ThresholdCalibration-accessors
#' @export
methods::setMethod(
    "calibrationThreshold", "ThresholdCalibration",
    function(object) {
        object@threshold
    }
)


#' @rdname ThresholdCalibration-accessors
#' @export
methods::setMethod(
    "calibrationThresholdNegLog10", "ThresholdCalibration",
    function(object) {
        object@thresholdNegLog10
    }
)


#' @rdname ThresholdCalibration-accessors
#' @export
methods::setMethod(
    "calibrationCurve", "ThresholdCalibration",
    function(object) {
        object@fdrCurve
    }
)


#' @rdname ThresholdCalibration-accessors
#' @export
methods::setMethod(
    "calibrationTargetFDR", "ThresholdCalibration",
    function(object) {
        object@targetFDR
    }
)


#' @rdname ThresholdCalibration-accessors
#' @export
methods::setMethod(
    "calibrationAchievedFDR", "ThresholdCalibration",
    function(object) {
        object@achievedFDR
    }
)


#' @rdname ThresholdCalibration-accessors
#' @export
methods::setMethod(
    "calibrationSeed", "ThresholdCalibration",
    function(object) {
        object@seed
    }
)


#' @rdname ThresholdCalibration-accessors
#' @export
methods::setMethod(
    "nullStatistics", "ThresholdCalibration",
    function(object) {
        object@null
    }
)


#' @rdname ThresholdCalibration-accessors
#' @export
methods::setMethod(
    "observedStatistics", "ThresholdCalibration",
    function(object) {
        object@observed
    }
)


## the calibration was a plain list up to 0.99.2, and `$` keeps the code
## written against it working
#' @rdname ThresholdCalibration-accessors
#' @importFrom methods slot slotNames
#' @export
methods::setMethod("$", "ThresholdCalibration", function(x, name) {
    if (!name %in% methods::slotNames(x)) {
        stop("a ThresholdCalibration has no '", name, "'")
    }
    methods::slot(x, name)
})


#' Summary of a threshold calibration
#'
#' @param object A [ThresholdCalibration-class] object.
#'
#' @return `NULL`, invisibly. Called for the summary it prints.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setMethod show
#'
#' @export
methods::setMethod("show", "ThresholdCalibration", function(object) {
    cat("ThresholdCalibration\n")
    cat(
        "  threshold       :", format(object@threshold, digits = 3),
        "( -log10", format(object@thresholdNegLog10, digits = 4), ")\n"
    )
    cat(
        "  FDR             :", object@achievedFDR, "achieved,",
        object@targetFDR, "asked for\n"
    )
    cat("  permutations    :", object@nPermutations, "\n")
    cat(
        "  statistics      :", length(object@observed), "observed,",
        length(object@null), "null\n"
    )

    ## an unseeded calibration cannot be repeated, which is worth saying
    ## where the seed would otherwise be printed
    cat(
        "  seed            :",
        if (is.null(object@seed)) {
            "none (not repeatable)"
        } else {
            paste0(object@seed, " (from ", object@seedSource, ")")
        }, "\n"
    )

    invisible(NULL)
})
