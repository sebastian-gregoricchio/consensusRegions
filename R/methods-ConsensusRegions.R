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
methods::setMethod("analysisParameters", "ConsensusRegions",
                   function(object) {
    object@parameters
})


#' @rdname ConsensusRegions-accessors
#' @export
methods::setMethod("consensusStats", "ConsensusRegions", function(object) {
    object@stats
})


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
#'                          c("rep1.narrowPeak", "rep2.narrowPeak",
#'                            "rep3.narrowPeak"),
#'                          package = "consensusRegions")
#' peaks <- readPeakSets(peakFiles, sampleNames = c("r1", "r2", "r3"),
#'                       verbose = FALSE)
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
    cat("  replicates      :", length(object@peaks), "(",
        paste(names(object@peaks), collapse = ", "), ")\n")
    cat("  score type      :", settings$scoreType, "\n")

    if (isTRUE(settings$presenceOnly)) {
        cat("  combination     : weighted presence rule\n")
    } else {
        cat("  combination     :", settings$combinationMethod, "\n")
        cat("  combined cut    :", settings$combinedThreshold, "\n")
    }

    cat("  weights         :",
        paste(names(object@weights), round(object@weights, 3),
              sep = "=", collapse = ", "), "\n")
    cat("  consensus       :", length(object@consensus), "regions\n")

    ## the width distribution is the first thing to check after a run,
    ## since runaway chaining shows up here before anywhere else
    if (length(object@consensus) > 0) {
        regionWidth <- GenomicRanges::width(object@consensus)
        cat("  width (median)  :", stats::median(regionWidth), "bp\n")
        cat("  width (max)     :", max(regionWidth), "bp\n")
    }

    if (length(object@calibration) > 0) {
        cat("  calibrated at   : FDR", object@calibration$targetFDR, "\n")
    }

    invisible(NULL)
})
