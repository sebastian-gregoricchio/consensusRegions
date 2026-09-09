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
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setGeneric
#'
#' @examples
#' peakFiles <- system.file("extdata",
#'                          c("rep1.narrowPeak", "rep2.narrowPeak",
#'                            "rep3.narrowPeak"),
#'                          package = "consensusRegions")
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
methods::setGeneric("consensusRanges",
                    function(object) standardGeneric("consensusRanges"))

#' @rdname ConsensusRegions-accessors
#' @export
methods::setGeneric("peakSets",
                    function(object) standardGeneric("peakSets"))

#' @rdname ConsensusRegions-accessors
#' @export
methods::setGeneric("replicateWeights",
                    function(object) standardGeneric("replicateWeights"))

#' @rdname ConsensusRegions-accessors
#' @export
methods::setGeneric("analysisParameters",
                    function(object) standardGeneric("analysisParameters"))

#' @rdname ConsensusRegions-accessors
#' @export
methods::setGeneric("consensusStats",
                    function(object) standardGeneric("consensusStats"))
