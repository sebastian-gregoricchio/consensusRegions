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
#'   `negLog10P`, `class`, `nSupport`, `combinedNegLog10P`,
#'   `combinedPadj` and `status`.
#' @slot consensus A `GRanges` of consensus regions.
#' @slot weights A named numeric vector of replicate weights.
#' @slot parameters A list with the call arguments.
#' @slot stats A `data.frame` summarising each replicate.
#' @slot calibration A list holding the output of [calibrateThreshold()]
#'   when one was supplied, empty otherwise.
#'
#' @author Sebastian Gregoricchio
#'
#' @importClassesFrom GenomicRanges GRanges GRangesList
#' @importFrom methods setClass validObject new is slot
#'
#' @name ConsensusRegions-class
#' @rdname ConsensusRegions-class
#' @exportClass ConsensusRegions
methods::setClass(
    "ConsensusRegions",
    slots = c(peaks = "GRangesList",
              consensus = "GRanges",
              weights = "numeric",
              parameters = "list",
              stats = "data.frame",
              calibration = "list")
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
            problems <- c(problems,
                          "names of 'weights' must match names of 'peaks'")
        }
        if (any(object@weights <= 0, na.rm = TRUE)) {
            problems <- c(problems, "'weights' must be strictly positive")
        }
    }

    if (length(problems) == 0) TRUE else problems
})
