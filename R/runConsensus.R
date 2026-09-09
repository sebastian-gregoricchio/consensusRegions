# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' Run the whole analysis in one call
#'
#' @description
#' Chains the steps that a consensus analysis almost always performs in
#' the same order: read the peaks, work out the replicate weights,
#' optionally calibrate the threshold, build the consensus, optionally
#' recentre it on summits, and optionally write it out.
#'
#' The individual functions remain the way to go when a step needs
#' inspecting before the next one runs, and calibration in particular is
#' worth looking at rather than trusting blind. This is for the case where
#' the settings are already decided and the analysis is being repeated
#' across experiments.
#'
#' @param peaks Passed to [readPeakSets()]: file paths, a `GRangesList`,
#'   or a list of `GRanges`.
#' @param sampleNames Character vector naming the replicates.
#'   Default: \code{NULL}.
#' @param weightMethod Passed to [computeReplicateWeights()]. `"equal"`
#'   leaves every replicate counting the same, which is the default and
#'   the reproducible choice. Default: \code{"equal"}.
#' @param bamFiles,frip,librarySize Passed to
#'   [computeReplicateWeights()] when `weightMethod` needs them.
#'   Default: \code{NULL}.
#' @param calibrate Estimate the combined threshold with
#'   [calibrateThreshold()] instead of using `combinedThreshold`.
#'   Default: \code{FALSE}.
#' @param nPermutations Number of permutations, passed to
#'   [calibrateThreshold()]. Default: \code{50L}.
#' @param targetFDR Empirical false discovery rate to aim for, passed to
#'   [calibrateThreshold()]. Default: \code{0.05}.
#' @param recentre Reduce the consensus to fixed-width regions around
#'   summits with [recentrePeaks()]. Sensible for transcription factors
#'   and ATAC, wrong for broad domains. Default: \code{FALSE}.
#' @param width Width of the recentred regions. Default: \code{400L}.
#' @param outputFile Optional path. When given the consensus is written
#'   there with [exportConsensus()]. Default: \code{NULL}.
#' @param seqlevelsStyle Chromosome naming style, passed to
#'   [readPeakSets()]. Default: \code{"UCSC"}.
#' @param excludeRegions Optional `GRanges` of blacklisted positions.
#'   Default: \code{NULL}.
#' @param BPPARAM Either the number of cores to use, or a
#'   `BiocParallelParam` object. Only the calibration uses it.
#'   Default: \code{1}.
#' @param verbose Report progress through the steps. Default: \code{TRUE}.
#' @param ... Further arguments passed to [buildConsensus()], for example
#'   `combinationMethod`, `minReplicates` or `mergeMethod`.
#'
#' @return A [ConsensusRegions-class] object. When `recentre = TRUE` the
#'   consensus it carries is the fixed-width version.
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso [readPeakSets()], [computeReplicateWeights()],
#'   [calibrateThreshold()], [buildConsensus()], [recentrePeaks()]
#'
#' @importFrom methods slot slot<-
#'
#' @examples
#' peakFiles <- system.file("extdata",
#'                          c("rep1.narrowPeak", "rep2.narrowPeak",
#'                            "rep3.narrowPeak"),
#'                          package = "consensusRegions")
#'
#' result <- runConsensus(peakFiles,
#'                        sampleNames = c("r1", "r2", "r3"),
#'                        minReplicates = 2,
#'                        verbose = FALSE)
#' result
#'
#' ## the same with measured weights and a calibrated threshold
#' \donttest{
#' set.seed(42)
#' calibrated <- runConsensus(peakFiles,
#'                            sampleNames = c("r1", "r2", "r3"),
#'                            weightMethod = "frip",
#'                            frip = c(0.21, 0.19, 0.06),
#'                            calibrate = TRUE,
#'                            nPermutations = 10,
#'                            verbose = FALSE)
#' }
#'
#' @export
runConsensus <- function(peaks,
                         sampleNames = NULL,
                         weightMethod = c("equal", "frip", "librarySize",
                                          "intrinsic"),
                         bamFiles = NULL,
                         frip = NULL,
                         librarySize = NULL,
                         calibrate = FALSE,
                         nPermutations = 50L,
                         targetFDR = 0.05,
                         recentre = FALSE,
                         width = 400L,
                         outputFile = NULL,
                         seqlevelsStyle = "UCSC",
                         excludeRegions = NULL,
                         BPPARAM = 1,
                         verbose = TRUE,
                         ...) {
    weightMethod <- match.arg(weightMethod)

    ## ---- read ---------------------------------------------------------
    .messageIf(verbose, "1/", if (isTRUE(calibrate)) "4" else "3",
               " Reading peaks")
    peakList <- readPeakSets(peaks,
                             sampleNames = sampleNames,
                             seqlevelsStyle = seqlevelsStyle,
                             verbose = verbose)

    ## ---- weights ------------------------------------------------------
    .messageIf(verbose, "2/", if (isTRUE(calibrate)) "4" else "3",
               " Deriving replicate weights (", weightMethod, ")")
    weights <- computeReplicateWeights(peakList,
                                       method = weightMethod,
                                       bamFiles = bamFiles,
                                       frip = frip,
                                       librarySize = librarySize,
                                       verbose = verbose)

    ## ---- calibration --------------------------------------------------
    calibration <- NULL
    if (isTRUE(calibrate)) {
        .messageIf(verbose, "3/4 Calibrating the combined threshold")
        calibration <- calibrateThreshold(peakList,
                                          weights = weights,
                                          nPermutations = nPermutations,
                                          targetFDR = targetFDR,
                                          excludeRegions = excludeRegions,
                                          BPPARAM = BPPARAM,
                                          verbose = verbose)
    }

    ## ---- consensus ----------------------------------------------------
    .messageIf(verbose, if (isTRUE(calibrate)) "4/4" else "3/3",
               " Building the consensus")
    result <- buildConsensus(peakList,
                             weights = weights,
                             calibration = calibration,
                             excludeRegions = excludeRegions,
                             verbose = verbose,
                             ...)

    ## ---- optional post-processing -------------------------------------
    if (isTRUE(recentre)) {
        .messageIf(verbose, "Recentring on summits at ", width, " bp")
        methods::slot(result, "consensus") <- recentrePeaks(result,
                                                            width = width)
        result@parameters$recentred <- width
    }

    if (!is.null(outputFile)) {
        exportConsensus(result, file = outputFile, verbose = verbose)
    }

    result
}
