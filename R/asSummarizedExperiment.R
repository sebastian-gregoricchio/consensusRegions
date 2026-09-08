#' Hand the consensus over to a differential analysis
#'
#' @description
#' Builds a `RangedSummarizedExperiment` whose rows are the consensus
#' regions and whose columns are the replicates. The assay holds the
#' per-replicate significance at each region, zero where the replicate had
#' no peak, together with a presence matrix.
#'
#' The point of this is the row ranges. Counting reads over a fixed,
#' shared set of regions is what csaw, DiffBind and edgeR all want, and
#' this puts the consensus in the shape they expect without asking anyone
#' to write out a BED file and read it back.
#'
#' @param object A [ConsensusRegions-class] object.
#' @param assayName Name given to the significance assay.
#'
#' @return A `RangedSummarizedExperiment`.
#'
#' @details
#' The assay is not a count matrix and must not be used as one. It records
#' how convincing each replicate found each region, which is useful for
#' clustering and for spotting a replicate that disagrees with the rest,
#' but a differential test needs actual reads counted over
#' `rowRanges(se)`.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges findOverlaps mcols
#' @importFrom S4Vectors queryHits subjectHits mcols DataFrame
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
#' se <- asSummarizedExperiment(result)
#' se
#'
#' @export
asSummarizedExperiment <- function(object, assayName = "negLog10P") {
    if (!methods::is(object, "ConsensusRegions")) {
        stop("'object' must be a ConsensusRegions object")
    }

    consensus <- object@consensus
    if (length(consensus) == 0) {
        stop("the consensus is empty, there is nothing to convert")
    }

    replicateNames <- names(object@peaks)
    significance <- matrix(0, nrow = length(consensus),
                           ncol = length(replicateNames),
                           dimnames = list(NULL, replicateNames))
    presence <- matrix(FALSE, nrow = length(consensus),
                       ncol = length(replicateNames),
                       dimnames = list(NULL, replicateNames))

    ## fill one column at a time, keeping the strongest peak when a
    ## replicate has several inside the same consensus region
    for (thisReplicate in replicateNames) {
        peaks <- object@peaks[[thisReplicate]]
        peaks <- peaks[S4Vectors::mcols(peaks)$status == "confirmed"]
        if (length(peaks) == 0) {
            next
        }

        hits <- GenomicRanges::findOverlaps(peaks, consensus,
                                            ignore.strand = TRUE)
        rowIndex <- S4Vectors::subjectHits(hits)
        value <- S4Vectors::mcols(peaks)$negLog10P[
            S4Vectors::queryHits(hits)]
        value[is.na(value)] <- 1

        significance[, thisReplicate] <- .maxByIndex(
            rowIndex, value, length(consensus))
        presence[unique(rowIndex), thisReplicate] <- TRUE
    }

    columnData <- S4Vectors::DataFrame(
        replicate = replicateNames,
        weight = unname(object@weights[replicateNames]),
        row.names = replicateNames)

    assayList <- list(significance, presence)
    names(assayList) <- c(assayName, "detected")

    SummarizedExperiment::SummarizedExperiment(
        assays = assayList,
        rowRanges = consensus,
        colData = columnData,
        metadata = object@parameters)
}


#' Largest value per row index
#'
#' @param index Integer vector of row positions.
#' @param value Numeric vector of the same length.
#' @param n Number of rows in the result.
#'
#' @return Numeric vector of length `n`.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal
#' @noRd
.maxByIndex <- function(index, value, n) {
    result <- numeric(n)
    ## sorting ascending means the last write per index is the largest
    ordering <- order(value, decreasing = FALSE)
    result[index[ordering]] <- value[ordering]
    result
}
