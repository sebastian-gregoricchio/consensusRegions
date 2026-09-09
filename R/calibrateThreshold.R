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
#'   Default: \code{NULL}.
#' @param combinationMethod Passed to [combineEvidence()].
#'   Default: \code{"stouffer"}.
#' @param nPermutations Number of shuffled replicate sets to generate.
#'   Default: \code{50L}.
#' @param targetFDR Empirical false discovery rate to aim for.
#'   Default: \code{0.05}.
#' @param stringencyThreshold Stringent p-value cut, as in
#'   [buildConsensus()]. Default: \code{1e-08}.
#' @param weakThreshold Background p-value cut. Default: \code{1e-04}.
#' @param minSupport Minimum supporting replicates. Give this or
#'   `minReplicates`, not both. Default: \code{NULL}.
#' @param minReplicates Total replicates that must hold the peak, counting
#'   its own. Accepts a count, a proportion or a percentage string, as in
#'   [buildConsensus()]. Must match the value used there.
#'   Default: \code{NULL}.
#' @param minOverlap Minimum overlap in base pairs. Default: \code{1L}.
#' @param minOverlapFraction Optional fractional overlap requirement.
#'   Must match the value used in [buildConsensus()]. Default: \code{NULL}.
#' @param recursive Whether the confirmation is re-run after pruning
#'   unsupported peaks. Must match the value used in [buildConsensus()].
#'   Default: \code{TRUE}.
#' @param maxIterations Cap on the recursive rounds. Default: \code{10L}.
#' @param multipleIntersections Resolution rule for several overlapping
#'   peaks from one replicate. Default: \code{"lowest"}.
#' @param chromosomeLengths Named integer vector. Taken from the input, or
#'   inferred from the furthest peak, when left `NULL`. Default: \code{NULL}.
#' @param excludeRegions Optional `GRanges` the shuffled peaks must avoid.
#'   Default: \code{NULL}.
#' @param BPPARAM Either the number of cores to use, or a
#'   `BiocParallelParam` object for finer control. The permutations are
#'   independent of one another and are where nearly all the time goes, so
#'   raising this is worth it on a large peak set: a single round takes
#'   around half a minute on 100,000 peaks per replicate. Default: \code{1}.
#' @param verbose Report progress. Default: \code{TRUE}.
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
#' The peak positions are drawn at random, so call [base::set.seed()]
#' beforehand if you need the same threshold back. The function does not
#' set the seed itself, since doing so would silently reset the random
#' number stream the rest of your session is drawing from.
#'
#' This is also why the default runs on a single core. Workers draw from
#' their own random streams, so `set.seed()` no longer governs the result
#' once `BPPARAM` is above one, and the seed has to travel to the workers
#' instead: `BiocParallel::MulticoreParam(workers = 8, RNGseed = 42)`.
#' Reach for that when a parallel run has to be reproducible.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom BiocParallel bplapply bpnworkers SerialParam
#'   MulticoreParam SnowParam
#' @importFrom GenomicRanges GRangesList mcols mcols<-
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom S4Vectors metadata
#' @importFrom dplyr tibble filter arrange slice_head pull mutate desc
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
#' set.seed(42)
#' calibration <- calibrateThreshold(peaks, nPermutations = 5,
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
                               minSupport = NULL,
                               minReplicates = NULL,
                               minOverlap = 1L,
                               minOverlapFraction = NULL,
                               recursive = TRUE,
                               maxIterations = 10L,
                               multipleIntersections = c("lowest", "highest"),
                               chromosomeLengths = NULL,
                               excludeRegions = NULL,
                               BPPARAM = 1,
                               verbose = TRUE) {
    combinationMethod <- match.arg(combinationMethod)
    multipleIntersections <- match.arg(multipleIntersections)
    BPPARAM <- .resolveBPPARAM(BPPARAM)

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
    ## the null has to be filtered by the same support rule as the data
    minSupport <- .resolveRequiredSupport(minSupport = minSupport,
                                          minReplicates = minReplicates,
                                          nReplicates = length(peakList))

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

    ## Every operation that shapes the observed statistic has to shape
    ## the null the same way, so the settings are carried across whole
    ## rather than partly re-specified inside the permutation loop.
    settings <- list(weights = weights,
                     combinationMethod = combinationMethod,
                     stringencyThreshold = stringencyThreshold,
                     weakThreshold = weakThreshold,
                     minSupport = minSupport,
                     minOverlap = minOverlap,
                     minOverlapFraction = minOverlapFraction,
                     recursive = recursive,
                     maxIterations = maxIterations,
                     multipleIntersections = multipleIntersections)

    ## observed statistics, computed once
    .messageIf(verbose, "Computing observed statistics")
    observed <- .combinedStatisticsOnly(peakList, settings)

    ## Null statistics, one round per permutation. Each round is
    ## self-contained, so they are handed to BiocParallel rather than run
    ## in sequence: on a genome-scale peak set a single round takes tens
    ## of seconds and fifty of them is most of an afternoon.
    .messageIf(verbose, "Running ", nPermutations, " permutations on ",
               BiocParallel::bpnworkers(BPPARAM), " worker(s)")

    nullStatistics <- BiocParallel::bplapply(
        seq_len(nPermutations),
        function(permutation) {
            shuffled <- .shufflePeakList(
                peakList,
                chromosomeLengths = chromosomeLengths,
                excludeRegions = excludeRegions)
            .combinedStatisticsOnly(shuffled, settings)
        },
        BPPARAM = BPPARAM)

    nullStatistics <- unlist(nullStatistics, use.names = FALSE)

    ## sweep candidate cuts across the observed range
    curve <- .empiricalFdrCurve(observed = observed,
                                nullValues = nullStatistics,
                                nPermutations = nPermutations)

    acceptable <- dplyr::filter(curve, .data$fdr <= targetFDR)
    if (nrow(acceptable) == 0) {
        warning("no cut reached an FDR of ", targetFDR,
                "; returning the cut with the lowest estimated FDR")
        chosen <- dplyr::slice_head(
            dplyr::arrange(curve, .data$fdr, dplyr::desc(.data$cut)), n = 1)
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

    overlapTable <- .buildOverlapTable(
        retained,
        minOverlap = settings$minOverlap,
        minOverlapFraction = settings$minOverlapFraction)

    verdict <- .runConfirmation(
        retained = retained,
        overlapTable = overlapTable,
        weights = settings$weights,
        combinationMethod = settings$combinationMethod,
        ## nothing is filtered here, the whole distribution is wanted
        combinedThreshold = 1,
        requiredSupport = settings$minSupport,
        multipleIntersections = settings$multipleIntersections,
        recursive = settings$recursive,
        maxIterations = settings$maxIterations,
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
#' @importFrom dplyr tibble arrange
#' @importFrom stats quantile
#' @importFrom rlang .data
#'
#' @keywords internal
#' @noRd
.empiricalFdrCurve <- function(observed, nullValues, nPermutations) {
    if (length(observed) == 0) {
        stop("no observed statistics were produced, nothing to calibrate")
    }

    ## The grid has to span the whole observed range. Starting it at the
    ## median, as an earlier version did, put a floor under the answer:
    ## no threshold below the median observed statistic could ever be
    ## chosen, however clean the null turned out to be.
    candidateCuts <- unique(stats::quantile(
        observed, probs = seq(0, 0.999, length.out = 200),
        na.rm = TRUE, names = FALSE))

    nObserved <- vapply(candidateCuts,
                        function(cut) sum(observed >= cut, na.rm = TRUE),
                        numeric(1))

    ## Seeing no null peak above a cut in a finite number of permutations
    ## is not evidence that none exists, and dividing straight through
    ## would report an FDR of exactly zero. The added pseudocount puts a
    ## floor of 1 / (nPermutations + 1) on the expected null count, which
    ## is the usual correction for permutation p-values.
    expectedNull <- vapply(
        candidateCuts,
        function(cut) {
            (1 + sum(nullValues >= cut, na.rm = TRUE)) / (nPermutations + 1)
        },
        numeric(1))

    dplyr::arrange(
        dplyr::tibble(cut = candidateCuts,
                      nObserved = nObserved,
                      expectedNull = expectedNull,
                      fdr = expectedNull / pmax(nObserved, 1)),
        .data$cut)
}
