#' Rescue and discard rates per replicate
#'
#' @description
#' Shows what happened to the peaks of each replicate. The bar to watch is
#' the rescued one: peaks that were only weak on their own and were kept
#' because the other replicates agreed. That is the whole point of the
#' method, but a rescued fraction far above what the replicates otherwise
#' share is a sign that `weakThreshold` was set too permissively and the
#' combination is confirming noise.
#'
#' @param object A [ConsensusRegions-class] object.
#' @param proportion Scale the bars to one instead of showing counts.
#'
#' @return A `ggplot` object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom ggplot2 ggplot aes geom_col labs theme_bw scale_fill_manual
#'   position_stack theme element_text
#' @importFrom dplyr tibble mutate select
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @examples
#' peakFiles <- system.file("extdata",
#'                          c("rep1.narrowPeak", "rep2.narrowPeak",
#'                            "rep3.narrowPeak"),
#'                          package = "consensusRegions")
#' peaks <- readPeakSets(peakFiles, sampleNames = c("r1", "r2", "r3"),
#'                       verbose = FALSE)
#' result <- buildConsensus(peaks, verbose = FALSE)
#'
#' plotRescue(result)
#'
#' @export
plotRescue <- function(object, proportion = FALSE) {
    if (!methods::is(object, "ConsensusRegions")) {
        stop("'object' must be a ConsensusRegions object")
    }

    ## split the confirmed peaks into those that stood on their own and
    ## those that needed the other replicates
    summary <- dplyr::mutate(
        tibble::as_tibble(object@stats),
        confirmedStringent = .data$nConfirmed - .data$nRescued)

    plotData <- tidyr::pivot_longer(
        dplyr::select(summary, "replicate", "confirmedStringent",
                      "nRescued", "nFalsePositive", "nDiscarded"),
        cols = -"replicate",
        names_to = "outcome",
        values_to = "count")

    plotData <- dplyr::mutate(
        plotData,
        outcome = factor(.data$outcome,
                         levels = c("confirmedStringent", "nRescued",
                                    "nFalsePositive", "nDiscarded"),
                         labels = c("confirmed", "rescued",
                                    "failed correction", "discarded")))

    barPosition <- if (isTRUE(proportion)) "fill" else "stack"

    ggplot2::ggplot(plotData,
                    ggplot2::aes(x = .data$replicate, y = .data$count,
                                 fill = .data$outcome)) +
        ggplot2::geom_col(position = barPosition, width = 0.7) +
        ggplot2::scale_fill_manual(
            values = c(confirmed = "#2C6E91",
                       rescued = "#78B4C8",
                       `failed correction` = "#E0B25F",
                       discarded = "#BDBDBD")) +
        ggplot2::labs(x = NULL,
                      y = if (isTRUE(proportion)) "fraction" else "peaks",
                      fill = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                           hjust = 1))
}


#' Pairwise agreement between replicates
#'
#' @description
#' Jaccard index over the peaks entering the analysis. A replicate that
#' agrees poorly with all the others is a candidate for down-weighting;
#' two blocks of mutually agreeing replicates usually mean a batch effect
#' rather than a quality problem.
#'
#' @param object A [ConsensusRegions-class] object.
#' @param showValues Print the index inside each tile.
#'
#' @return A `ggplot` object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient
#'   labs theme_minimal theme element_text element_blank
#' @importFrom dplyr tibble mutate
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @examples
#' peakFiles <- system.file("extdata",
#'                          c("rep1.narrowPeak", "rep2.narrowPeak",
#'                            "rep3.narrowPeak"),
#'                          package = "consensusRegions")
#' peaks <- readPeakSets(peakFiles, sampleNames = c("r1", "r2", "r3"),
#'                       verbose = FALSE)
#' plotJaccard(buildConsensus(peaks, verbose = FALSE))
#'
#' @export
plotJaccard <- function(object, showValues = TRUE) {
    if (!methods::is(object, "ConsensusRegions")) {
        stop("'object' must be a ConsensusRegions object")
    }

    jaccard <- .pairwiseJaccard(object@peaks)

    plotData <- tidyr::pivot_longer(
        dplyr::mutate(as.data.frame(jaccard),
                      replicateA = rownames(jaccard)),
        cols = -"replicateA",
        names_to = "replicateB",
        values_to = "jaccard")

    plot <- ggplot2::ggplot(
        plotData,
        ggplot2::aes(x = .data$replicateA, y = .data$replicateB,
                     fill = .data$jaccard)) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::scale_fill_gradient(low = "#F5F5F5", high = "#2C6E91",
                                     limits = c(0, 1)) +
        ggplot2::labs(x = NULL, y = NULL, fill = "Jaccard") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    if (isTRUE(showValues)) {
        plot <- plot +
            ggplot2::geom_text(
                ggplot2::aes(label = sprintf("%.2f", .data$jaccard)),
                size = 3, colour = "grey20")
    }

    plot
}


#' Which replicates contribute to each consensus region
#'
#' @description
#' An UpSet plot of the replicate combinations found across the consensus
#' regions. Regions supported by every replicate should dominate; a large
#' partial set points at a replicate that is pulling in a different
#' direction.
#'
#' @param object A [ConsensusRegions-class] object.
#' @param minSize Drop intersections smaller than this.
#'
#' @return A `ggplot` object built by `ComplexUpset`.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges findOverlaps mcols
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom methods is
#'
#' @examples
#' peakFiles <- system.file("extdata",
#'                          c("rep1.narrowPeak", "rep2.narrowPeak",
#'                            "rep3.narrowPeak"),
#'                          package = "consensusRegions")
#' peaks <- readPeakSets(peakFiles, sampleNames = c("r1", "r2", "r3"),
#'                       verbose = FALSE)
#' result <- buildConsensus(peaks, verbose = FALSE)
#'
#' if (requireNamespace("ComplexUpset", quietly = TRUE)) {
#'     plotUpset(result)
#' }
#'
#' @export
plotUpset <- function(object, minSize = 0) {
    if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
        stop("plotUpset() needs the ComplexUpset package, install it with ",
             "install.packages('ComplexUpset')")
    }
    if (!methods::is(object, "ConsensusRegions")) {
        stop("'object' must be a ConsensusRegions object")
    }

    ## a membership matrix is all ComplexUpset needs
    membership <- .membershipMatrix(object)
    if (nrow(membership) == 0) {
        stop("the consensus is empty, there is nothing to plot")
    }

    ComplexUpset::upset(
        as.data.frame(membership),
        intersect = colnames(membership),
        min_size = minSize,
        name = "replicates")
}


#' The calibration curve
#'
#' @description
#' Plots the empirical false discovery rate against the cut on the
#' combined statistic, with the chosen threshold marked. A curve that
#' never descends to the target says the replicates do not agree beyond
#' chance often enough to support a consensus at that stringency.
#'
#' @param calibration Output of [calibrateThreshold()].
#'
#' @return A `ggplot` object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_vline labs
#'   theme_bw scale_y_continuous
#' @importFrom rlang .data
#'
#' @examples
#' peakFiles <- system.file("extdata",
#'                          c("rep1.narrowPeak", "rep2.narrowPeak",
#'                            "rep3.narrowPeak"),
#'                          package = "consensusRegions")
#' peaks <- readPeakSets(peakFiles, sampleNames = c("r1", "r2", "r3"),
#'                       verbose = FALSE)
#' calibration <- calibrateThreshold(peaks, nPermutations = 5, seed = 1,
#'                                   verbose = FALSE)
#'
#' plotCalibration(calibration)
#'
#' @export
plotCalibration <- function(calibration) {
    if (!is.list(calibration) || is.null(calibration$fdrCurve)) {
        stop("'calibration' must be the output of calibrateThreshold()")
    }

    ggplot2::ggplot(calibration$fdrCurve,
                    ggplot2::aes(x = .data$cut, y = .data$fdr)) +
        ggplot2::geom_line(colour = "#2C6E91", linewidth = 0.7) +
        ggplot2::geom_hline(yintercept = calibration$targetFDR,
                            linetype = "dashed", colour = "grey40") +
        ggplot2::geom_vline(xintercept = calibration$thresholdNegLog10,
                            linetype = "dotted", colour = "#B4553C") +
        ggplot2::labs(x = "combined significance cut (-log10 p)",
                      y = "empirical FDR") +
        ggplot2::theme_bw()
}


#' Replicate membership of the consensus regions
#'
#' @param object A [ConsensusRegions-class] object.
#'
#' @return A logical matrix, regions by replicates.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits mcols
#'
#' @keywords internal
#' @noRd
.membershipMatrix <- function(object) {
    consensus <- object@consensus
    replicateNames <- names(object@peaks)

    membership <- matrix(FALSE, nrow = length(consensus),
                         ncol = length(replicateNames),
                         dimnames = list(NULL, replicateNames))
    if (length(consensus) == 0) {
        return(membership)
    }

    for (thisReplicate in replicateNames) {
        peaks <- object@peaks[[thisReplicate]]
        peaks <- peaks[S4Vectors::mcols(peaks)$status == "confirmed"]
        if (length(peaks) == 0) {
            next
        }
        hits <- GenomicRanges::findOverlaps(peaks, consensus,
                                            ignore.strand = TRUE)
        membership[unique(S4Vectors::subjectHits(hits)), thisReplicate] <-
            TRUE
    }

    membership
}
