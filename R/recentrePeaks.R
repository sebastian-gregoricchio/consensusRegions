# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' Recover summits and return fixed-width regions
#'
#' @description
#' Consensus regions inherit their edges from whichever peaks happened to
#' be merged, so their widths vary a great deal. That variation biases
#' every count-based and enrichment-based analysis downstream, since a
#' wider region collects more reads and more motif matches for reasons
#' that have nothing to do with the biology.
#'
#' Recentring fixes it. Each region is reduced to a single position and
#' re-expanded to a common width, which is what the ATAC-seq peak atlases
#' do and what motif analysis wants.
#'
#' @param object A [ConsensusRegions-class] object.
#' @param width Width of the returned regions in base pairs.
#'   Default: \code{400L}.
#' @param summitSource How to place the anchor. `"best"` takes the summit
#'   of the most significant member peak, `"weighted"` averages the member
#'   summits weighted by their significance, and `"centre"` uses the
#'   midpoint of the merged region. Default: \code{"best"}.
#' @param chromosomeLengths Named integer vector used to trim regions that
#'   would run off a chromosome end. Taken from the object when `NULL`.
#'   Default: \code{NULL}.
#'
#' @return A `GRanges` of fixed-width regions carrying the summit position
#'   and the metadata of the consensus it came from.
#'
#' @details
#' `"best"` needs a summit offset, which narrowPeak files carry in their
#' tenth column and other formats do not. When it is missing the midpoint
#' of the member peak is used instead, and a warning says so.
#'
#' Widths between 200 and 500 bp suit transcription factors and ATAC
#' peaks. Broad histone domains should not be recentred at all, since
#' collapsing a domain to a point throws away the thing being measured.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges GRanges start end width findOverlaps
#'   seqnames
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom GenomeInfoDb seqlengths seqinfo
#' @importFrom BiocGenerics unlist
#' @importFrom dplyr tibble group_by summarise mutate
#' @importFrom rlang .data
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
#' recentrePeaks(result, width = 400)
#'
#' @export
recentrePeaks <- function(object,
                          width = 400L,
                          summitSource = c("best", "weighted", "centre"),
                          chromosomeLengths = NULL) {
    summitSource <- match.arg(summitSource)

    if (!methods::is(object, "ConsensusRegions")) {
        stop("'object' must be a ConsensusRegions object")
    }
    if (width < 1) {
        stop("'width' must be a positive number of base pairs")
    }

    consensus <- object@consensus
    if (length(consensus) == 0) {
        stop("the consensus is empty, there is nothing to recentre")
    }

    ## the midpoint route needs no member peaks at all
    if (summitSource == "centre") {
        anchor <- GenomicRanges::start(consensus) +
            round(GenomicRanges::width(consensus) / 2)
        return(.expandAroundAnchor(consensus, anchor, width,
                                   chromosomeLengths))
    }

    confirmedPeaks <- .confirmedPeaks(object)
    if (length(confirmedPeaks) == 0) {
        stop("no confirmed peaks are stored, cannot locate summits")
    }

    peakMetadata <- S4Vectors::mcols(confirmedPeaks)

    ## narrowPeak keeps the summit offset in a column named peak; without
    ## it the best that can be done is the middle of the peak
    if ("peak" %in% names(peakMetadata) &&
        any(peakMetadata$peak >= 0, na.rm = TRUE)) {
        summitPosition <- GenomicRanges::start(confirmedPeaks) +
            pmax(peakMetadata$peak, 0)
    } else {
        warning("no summit offsets found in the peak files, using peak ",
                "midpoints instead")
        summitPosition <- GenomicRanges::start(confirmedPeaks) +
            round(GenomicRanges::width(confirmedPeaks) / 2)
    }

    ## map every confirmed peak onto the consensus region holding it
    hits <- GenomicRanges::findOverlaps(confirmedPeaks, consensus,
                                        ignore.strand = TRUE)
    members <- dplyr::tibble(
        region = S4Vectors::subjectHits(hits),
        summit = summitPosition[S4Vectors::queryHits(hits)],
        negLog10P = peakMetadata$negLog10P[S4Vectors::queryHits(hits)])

    ## a missing statistic must not silently drop out of a weighted mean
    members <- dplyr::mutate(
        members,
        negLog10P = ifelse(is.na(.data$negLog10P), 1, .data$negLog10P))

    anchorTable <- if (summitSource == "best") {
        dplyr::summarise(
            dplyr::group_by(members, .data$region),
            anchor = .data$summit[which.max(.data$negLog10P)],
            .groups = "drop")
    } else {
        dplyr::summarise(
            dplyr::group_by(members, .data$region),
            anchor = round(sum(.data$summit * .data$negLog10P) /
                               sum(.data$negLog10P)),
            .groups = "drop")
    }

    ## regions with no member peak fall back to their own midpoint
    anchor <- GenomicRanges::start(consensus) +
        round(GenomicRanges::width(consensus) / 2)
    anchor[anchorTable$region] <- anchorTable$anchor

    .expandAroundAnchor(consensus, anchor, width, chromosomeLengths)
}


#' Grow a fixed-width window around each anchor
#'
#' @param consensus `GRanges` of consensus regions.
#' @param anchor Integer vector of anchor positions.
#' @param width Requested width.
#' @param chromosomeLengths Optional named integer vector.
#'
#' @return A `GRanges` of fixed-width regions.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges GRanges seqnames start end
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlengths seqinfo
#'
#' @keywords internal
#' @noRd
.expandAroundAnchor <- function(consensus, anchor, width,
                                chromosomeLengths = NULL) {
    halfWidth <- floor(width / 2)

    newStart <- pmax(anchor - halfWidth, 1)
    newEnd <- newStart + width - 1L

    ## a window near a chromosome end would otherwise stick out past it,
    ## so trim before the ranges are built rather than patching after
    if (is.null(chromosomeLengths)) {
        chromosomeLengths <- GenomeInfoDb::seqlengths(consensus)
    }
    if (!all(is.na(chromosomeLengths))) {
        bound <- chromosomeLengths[
            as.character(GenomicRanges::seqnames(consensus))]
        trimmable <- !is.na(bound) & newEnd > bound
        newEnd[trimmable] <- bound[trimmable]
    }

    recentred <- GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(consensus),
        ranges = IRanges::IRanges(start = newStart, end = newEnd),
        seqinfo = GenomeInfoDb::seqinfo(consensus))

    S4Vectors::mcols(recentred) <- S4Vectors::mcols(consensus)
    S4Vectors::mcols(recentred)$summit <- anchor

    recentred
}


#' Confirmed peaks pulled out of the object
#'
#' @param object A [ConsensusRegions-class] object.
#'
#' @return A flattened `GRanges` of confirmed peaks.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom BiocGenerics unlist
#' @importFrom S4Vectors mcols
#'
#' @keywords internal
#' @noRd
.confirmedPeaks <- function(object) {
    allPeaks <- BiocGenerics::unlist(object@peaks, use.names = FALSE)
    allPeaks[S4Vectors::mcols(allPeaks)$status == "confirmed"]
}
