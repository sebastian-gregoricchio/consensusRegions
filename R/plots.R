#' A theme in the style of ggpubr::theme_pubr
#'
#' @description
#' Reassembles the look of `ggpubr::theme_pubr()` without taking on the
#' dependency: no grid, black axis text, a solid axis line, and the legend
#' on the right rather than on top. Passing `border = TRUE` swaps the axis
#' lines for a closed panel border, which suits a heatmap.
#'
#' @param baseSize Base font size in points.
#' @param legendPosition Where to put the legend, or `"none"`.
#' @param border Draw a panel border instead of axis lines.
#' @param borderWidth Line width of the panel border.
#'
#' @return A `ggplot2` theme object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom ggplot2 theme_bw theme element_text element_line
#'   element_blank element_rect
#'
#' @keywords internal
#' @noRd
.themePubrLike <- function(baseSize = 12,
                           legendPosition = "right",
                           border = FALSE,
                           borderWidth = 1) {
  ## a heatmap wants a closed panel, everything else wants open axes
  panelBorder <- if (isTRUE(border)) {
    ggplot2::element_rect(colour = "black", fill = NA,
                          linewidth = borderWidth)
  } else {
    ggplot2::element_blank()
  }
  axisLine <- if (isTRUE(border)) {
    ggplot2::element_blank()
  } else {
    ggplot2::element_line(colour = "black", linewidth = 0.5)
  }

  ggplot2::theme_bw(base_size = baseSize) +
    ggplot2::theme(
      panel.border = panelBorder,
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = axisLine,
      axis.text = ggplot2::element_text(colour = "black",
                                        size = baseSize - 1),
      axis.title = ggplot2::element_text(colour = "black",
                                         size = baseSize),
      axis.ticks = ggplot2::element_line(colour = "black",
                                         linewidth = 0.5),
      legend.position = legendPosition,
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(colour = "black",
                                          size = baseSize - 1),
      legend.title = ggplot2::element_text(colour = "black",
                                           size = baseSize),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(colour = "black",
                                         size = baseSize),
      plot.title = ggplot2::element_text(colour = "black",
                                         size = baseSize + 2,
                                         face = "bold", hjust = 0.5))
}


#' Pick a label colour that stays readable on its tile
#'
#' @description
#' Interpolates the fill that will actually be drawn for each value, then
#' decides from its perceived brightness: white lettering on a dark tile,
#' black on a light one. Deciding from the value alone would break as soon
#' as the palette changed.
#'
#' @param value Numeric vector of the values being labelled.
#' @param lowColour Colour at the bottom of the fill scale.
#' @param highColour Colour at the top of the fill scale.
#' @param limits Numeric of length two, the fill scale limits.
#'
#' @return Character vector of colours, one per value.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom grDevices colorRamp
#'
#' @keywords internal
#' @noRd
.contrastingLabelColour <- function(value,
                                    lowColour,
                                    highColour,
                                    limits = c(0, 1)) {
  scaled <- (value - limits[1]) / diff(limits)
  scaled[is.na(scaled)] <- 0
  scaled <- pmin(pmax(scaled, 0), 1)

  channels <- grDevices::colorRamp(c(lowColour, highColour))(scaled)

  ## the usual weighted sum for perceived brightness: the eye is far
  ## more sensitive to green than it is to blue
  luminance <- (channels[, 1] * 0.299 +
                  channels[, 2] * 0.587 +
                  channels[, 3] * 0.114) / 255

  ifelse(luminance > 0.5, "black", "white")
}


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
#' @param baseSize Base font size in points.
#'
#' @return A `ggplot` object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom ggplot2 ggplot aes geom_col labs scale_fill_manual theme
#'   element_blank element_text
#' @importFrom dplyr mutate select
#' @importFrom tibble as_tibble
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
plotRescue <- function(object, proportion = FALSE, baseSize = 12) {
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
    .themePubrLike(baseSize = baseSize, legendPosition = "right") +
    ## the replicate names sit under the bars already, so the ticks
    ## only add clutter
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank())
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
#' @param showValues Print the index inside each tile. The lettering
#'   switches between black and white so that it stays legible whatever
#'   the tile is filled with.
#' @param lowColour Colour for a Jaccard index of zero.
#' @param highColour Colour for a Jaccard index of one.
#' @param baseSize Base font size in points.
#'
#' @return A `ggplot` object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient
#'   scale_colour_identity scale_x_discrete scale_y_discrete expansion
#'   labs theme
#' @importFrom dplyr mutate
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
plotJaccard <- function(object,
                        showValues = TRUE,
                        lowColour = "#F5F5F5",
                        highColour = "#2C6E91",
                        baseSize = 12) {
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

  ## the lettering follows the fill rather than the value, so it keeps
  ## working if the palette is changed
  plotData <- dplyr::mutate(
    plotData,
    labelColour = .contrastingLabelColour(.data$jaccard,
                                          lowColour = lowColour,
                                          highColour = highColour,
                                          limits = c(0, 1)))

  plot <- ggplot2::ggplot(
    plotData,
    ggplot2::aes(x = .data$replicateA, y = .data$replicateB,
                 fill = .data$jaccard)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_gradient(low = lowColour, high = highColour,
                                 limits = c(0, 1)) +
    ## half a category of padding puts the panel edge exactly on the
    ## outer tile boundary: no gap, and no tile cut in half
    ggplot2::scale_x_discrete(
      expand = ggplot2::expansion(mult = 0, add = 0.5)) +
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(mult = 0, add = 0.5)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Jaccard") +
    .themePubrLike(baseSize = baseSize, legendPosition = "right",
                   border = TRUE, borderWidth = 1) +
    ggplot2::theme(aspect.ratio = 1,
                   axis.ticks = element_blank())

  if (isTRUE(showValues)) {
    plot <- plot +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", .data$jaccard),
                     colour = .data$labelColour),
        size = baseSize / 4) +
      ggplot2::scale_colour_identity()
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
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits mcols
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
#' @param baseSize Base font size in points.
#'
#' @return A `ggplot` object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_vline labs
#'   theme
#' @importFrom ggtext element_markdown
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
plotCalibration <- function(calibration, baseSize = 12) {
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
    ## the base of the logarithm belongs in subscript, which plotmath
    ## cannot mix with running text the way markdown can
    ggplot2::labs(
      x = "combined significance cut (-log<sub>10</sub> p)",
      y = "empirical FDR") +
    .themePubrLike(baseSize = baseSize, legendPosition = "none") +
    ggplot2::theme(
      axis.title.x = ggtext::element_markdown(colour = "black",
                                              size = baseSize))
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
#' @importFrom S4Vectors subjectHits mcols
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
