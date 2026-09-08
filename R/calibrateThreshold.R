#' Calibrate the combined significance threshold
#'
#' @description
#' The threshold on the combined p-value is normally set by hand, and
#' `1e-8` has been carried around for long enough that it now looks like a
#' property of the method rather than a guess. It is a guess, and the
#' right value depends on how many replicates there are, how they were
#' called, and how densely the peaks sit.
#'
#' This function estimates it instead. Peak positions are shuffled within
#' each chromosome while their widths, counts and statistics are kept, so
#' any overlap between replicates in the shuffled sets happens by chance.
#' Running the combination on those gives a null distribution, and the
#' threshold is the point where the expected number of null peaks reaching
#' it falls to `targetFDR` of the observed number.
#'
#' @param peakList A `GRangesList` from [readPeakSets()].
#' @param weights Named numeric vector of replicate weights, or `NULL`.
#' @param combinationMethod Passed to [combineEvidence()].
#' @param nPermutations Number of shuffled replicate sets to generate.
#' @param targetFDR Empirical false discovery rate to aim for.
#' @param stringencyThreshold Stringent p-value cut, as in
#'   [buildConsensus()].
#' @param weakThreshold Background p-value cut.
#' @param minSupport Minimum supporting replicates.
#' @param minOverlap Minimum overlap in base pairs.
#' @param multipleIntersections Resolution rule for several overlapping
#'   peaks from one replicate.
#' @param chromosomeLengths Named integer vector. Taken from the input, or
#'   inferred from the furthest peak, when left `NULL`.
#' @param excludeRegions Optional `GRanges` the shuffled peaks must avoid.
#' @param seed Optional integer for reproducibility.
#' @param verbose Report progress.
#'
#' @return A list with the recommended `threshold` on the p-value scale,
#'   the `fdrCurve` it was read off, and the observed and null combined
#'   statistics.
#'
#' @details
#' Shuffling within a chromosome preserves peak density but not the
#' relationship between peaks and genes, so the null is a little
#' optimistic wherever peaks cluster around promoters. Passing
#' `excludeRegions` with a blacklist, and restricting the analysis to
#' mappable chromosomes beforehand, both help.
#'
#' Fifty permutations are usually enough to place the threshold within a
#' factor of two, which is as much precision as the choice deserves. Push
#' it higher only if the curve looks ragged near `targetFDR`.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges GRangesList mcols mcols<-
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom S4Vectors metadata
#' @importFrom dplyr tibble filter arrange slice_head pull mutate
#' @importFrom rlang .data
#' @importFrom stats quantile
#'
#' @examples
#' peakFiles <- system.file("extdata",
#'                          c("rep1.narrowPeak", "rep2.narrowPeak",
#'                            "rep3.narrowPeak"),
#'                          package = "consensusRegions")
#' peaks <- readPeakSets(peakFiles, sampleNames = c("r1", "r2", "r3"),
#'                       verbose = FALSE)
#'
#' calibration <- calibrateThreshold(peaks, nPermutations = 5, seed = 1,
#'                                   verbose = FALSE)
#' calibration$threshold
#'
#' result <- buildConsensus(peaks,
#'                          combinedThreshold = calibration$threshold,
#'                          verbose = FALSE)
#'
#' @export
calibrateThreshold <- function(peakList,
                               weights = NULL,
                               combinationMethod = c("stouffer", "fisher",
                                                     "lancaster",
                                                     "rankProduct"),
                               nPermutations = 50L,
                               targetFDR = 0.05,
                               stringencyThreshold = 1e-8,
                               weakThreshold = 1e-4,
                               minSupport = 1L,
                               minOverlap = 1L,
                               multipleIntersections = c("lowest", "highest"),
                               chromosomeLengths = NULL,
                               excludeRegions = NULL,
                               seed = NULL,
                               verbose = TRUE) {
    combinationMethod <- match.arg(combinationMethod)
    multipleIntersections <- match.arg(multipleIntersections)

    scoreType <- S4Vectors::metadata(peakList)$scoreType
    if (identical(scoreType, "none")) {
        stop("calibration needs a per-peak statistic; with BED3 input ",
             "there is no combined p-value to calibrate")
    }
    if (targetFDR <= 0 || targetFDR >= 1) {
        stop("'targetFDR' must lie strictly between 0 and 1")
    }
    if (nPermutations < 1) {
        stop("'nPermutations' must be at least 1")
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }

    replicateNames <- names(peakList)
    if (is.null(weights)) {
        weights <- rep(1, length(peakList))
        names(weights) <- replicateNames
    }
    weights <- .normaliseWeights(weights[replicateNames])

    ## chromosome bounds for the shuffling
    if (is.null(chromosomeLengths)) {
        chromosomeLengths <- .inferSeqlengths(
            BiocGenerics::unlist(peakList, use.names = FALSE))
    }

    settings <- list(weights = weights,
                     combinationMethod = combinationMethod,
                     stringencyThreshold = stringencyThreshold,
                     weakThreshold = weakThreshold,
                     minSupport = minSupport,
                     minOverlap = minOverlap,
                     multipleIntersections = multipleIntersections)

    ## observed statistics, computed once
    .messageIf(verbose, "Computing observed statistics")
    observed <- .combinedStatisticsOnly(peakList, settings)

    ## null statistics, one round per permutation
    nullStatistics <- vector("list", nPermutations)
    for (permutation in seq_len(nPermutations)) {
        .messageIf(verbose, "Permutation ", permutation, " of ",
                   nPermutations)
        shuffled <- .shufflePeakList(peakList,
                                     chromosomeLengths = chromosomeLengths,
                                     excludeRegions = excludeRegions)
        nullStatistics[[permutation]] <-
            .combinedStatisticsOnly(shuffled, settings)
    }
    nullStatistics <- unlist(nullStatistics, use.names = FALSE)

    ## sweep candidate cuts across the observed range
    curve <- .empiricalFdrCurve(observed = observed,
                                nullValues = nullStatistics,
                                nPermutations = nPermutations)

    acceptable <- dplyr::filter(curve, .data$fdr <= targetFDR)
    if (nrow(acceptable) == 0) {
        warning("no cut reached an FDR of ", targetFDR,
                "; returning the most stringent cut examined")
        chosen <- dplyr::slice_head(dplyr::arrange(curve, .data$fdr), n = 1)
    } else {
        ## the lowest cut still meeting the target keeps the most peaks
        chosen <- dplyr::slice_head(
            dplyr::arrange(acceptable, .data$cut), n = 1)
    }

    thresholdP <- 10^(-dplyr::pull(chosen, .data$cut))

    .messageIf(verbose, "Threshold at FDR ", targetFDR, ": ", thresholdP)

    list(threshold = thresholdP,
         thresholdNegLog10 = dplyr::pull(chosen, .data$cut),
         achievedFDR = dplyr::pull(chosen, .data$fdr),
         targetFDR = targetFDR,
         fdrCurve = curve,
         observed = observed,
         null = nullStatistics,
         nPermutations = nPermutations)
}


#' Combined statistics without the rest of the analysis
#'
#' @param peakList A `GRangesList`.
#' @param settings List of analysis settings.
#'
#' @return Numeric vector of combined significance for the peaks that had
#'   enough support.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom BiocGenerics unlist
#' @importFrom S4Vectors mcols mcols<-
#'
#' @keywords internal
#' @noRd
.combinedStatisticsOnly <- function(peakList, settings) {
    flatPeaks <- BiocGenerics::unlist(peakList, use.names = FALSE)
    S4Vectors::mcols(flatPeaks)$replicate <-
        rep(names(peakList), lengths(peakList))
    S4Vectors::mcols(flatPeaks)$weight <-
        unname(settings$weights[S4Vectors::mcols(flatPeaks)$replicate])

    flatPeaks <- .classifyPeaks(
        flatPeaks,
        stringencyThreshold = settings$stringencyThreshold,
        weakThreshold = settings$weakThreshold,
        presenceOnly = FALSE)

    retained <- flatPeaks[S4Vectors::mcols(flatPeaks)$class != "background"]
    if (length(retained) == 0) {
        return(numeric(0))
    }

    S4Vectors::mcols(retained)$rho <- .relativeRanks(
        negLog10P = S4Vectors::mcols(retained)$negLog10P,
        replicate = S4Vectors::mcols(retained)$replicate,
        presenceOnly = FALSE)

    overlapTable <- .buildOverlapTable(retained,
                                       minOverlap = settings$minOverlap,
                                       minOverlapFraction = NULL)

    verdict <- .runConfirmation(
        retained = retained,
        overlapTable = overlapTable,
        weights = settings$weights,
        combinationMethod = settings$combinationMethod,
        ## nothing is filtered here, the whole distribution is wanted
        combinedThreshold = 1,
        requiredSupport = settings$minSupport,
        multipleIntersections = settings$multipleIntersections,
        recursive = FALSE,
        maxIterations = 1L,
        presenceOnly = FALSE,
        verbose = FALSE)

    ## unsupported peaks never reach the combination in a real run either
    verdict$combinedNegLog10P[verdict$nSupport >= settings$minSupport]
}


#' Shuffle every replicate independently
#'
#' @param peakList A `GRangesList`.
#' @param chromosomeLengths Named integer vector.
#' @param excludeRegions Optional `GRanges` to avoid.
#'
#' @return A `GRangesList` with the same metadata but shuffled positions.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges GRangesList
#' @importFrom S4Vectors mcols mcols<- metadata metadata<-
#'
#' @keywords internal
#' @noRd
.shufflePeakList <- function(peakList, chromosomeLengths,
                             excludeRegions = NULL) {
    shuffled <- lapply(peakList, function(gr) {
        moved <- .shuffleRanges(gr,
                                chromosomeLengths = chromosomeLengths,
                                excludeRegions = excludeRegions)
        ## the statistics travel with the peaks, only the positions move
        S4Vectors::mcols(moved) <- S4Vectors::mcols(gr)
        moved
    })

    result <- GenomicRanges::GRangesList(shuffled)
    S4Vectors::metadata(result) <- S4Vectors::metadata(peakList)
    result
}


#' Empirical FDR across candidate cuts
#'
#' @param observed Numeric vector of observed statistics.
#' @param nullValues Numeric vector pooled over permutations.
#' @param nPermutations Number of permutations that produced `nullValues`.
#'
#' @return A `tibble` with `cut`, `nObserved`, `expectedNull` and `fdr`.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr tibble filter arrange
#' @importFrom stats quantile
#' @importFrom rlang .data
#'
#' @keywords internal
#' @noRd
.empiricalFdrCurve <- function(observed, nullValues, nPermutations) {
    if (length(observed) == 0) {
        stop("no observed statistics were produced, nothing to calibrate")
    }

    ## a grid over the observed range is enough; finer steps would only
    ## add noise given how few permutations are usually affordable
    candidateCuts <- unique(stats::quantile(
        observed, probs = seq(0.5, 0.999, length.out = 100),
        na.rm = TRUE, names = FALSE))

    nObserved <- vapply(candidateCuts,
                        function(cut) sum(observed >= cut, na.rm = TRUE),
                        numeric(1))
    expectedNull <- vapply(
        candidateCuts,
        function(cut) sum(nullValues >= cut, na.rm = TRUE) / nPermutations,
        numeric(1))

    dplyr::arrange(
        dplyr::tibble(cut = candidateCuts,
                      nObserved = nObserved,
                      expectedNull = expectedNull,
                      fdr = expectedNull / pmax(nObserved, 1)),
        .data$cut)
}
