#' Write the results to disk
#'
#' @description
#' Writes the consensus regions, and optionally the annotated per-replicate
#' peaks, as BED files. The combined significance is carried in the score
#' column, capped at the 1000 that the BED specification allows, so the
#' output loads into a genome browser without complaint.
#'
#' @param object A [ConsensusRegions-class] object.
#' @param file Path for the consensus file.
#' @param format `"bed"` or `"narrowPeak"`.
#' @param perReplicate Also write one file per replicate holding every
#'   peak with its verdict. The names are derived from `file`.
#' @param verbose Report what was written.
#'
#' @return Invisibly, a character vector of the paths written.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom rtracklayer export
#' @importFrom GenomicRanges width
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom S4Vectors mcols
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
#' target <- file.path(tempdir(), "consensus.bed")
#' exportConsensus(result, target, verbose = FALSE)
#'
#' @export
exportConsensus <- function(object,
                            file,
                            format = c("bed", "narrowPeak"),
                            perReplicate = FALSE,
                            verbose = TRUE) {
    format <- match.arg(format)

    if (!methods::is(object, "ConsensusRegions")) {
        stop("'object' must be a ConsensusRegions object")
    }
    if (length(object@consensus) == 0) {
        stop("the consensus is empty, there is nothing to write")
    }

    written <- character(0)

    ## the browser expects a name and a bounded score in those columns
    consensus <- object@consensus
    S4Vectors::mcols(consensus)$name <-
        paste0("consensus_", seq_along(consensus))
    S4Vectors::mcols(consensus)$score <-
        .scoreForExport(S4Vectors::mcols(consensus)$combinedNegLog10P)

    ## narrowPeak is a fixed ten-column format, so fill what it expects
    ## rather than letting the export fail on a missing column
    if (format == "narrowPeak") {
        combined <- S4Vectors::mcols(consensus)$combinedNegLog10P
        if (is.null(combined)) {
            combined <- rep(NA_real_, length(consensus))
        }
        S4Vectors::mcols(consensus)$signalValue <-
            ifelse(is.na(combined), 0, combined)
        S4Vectors::mcols(consensus)$pValue <-
            ifelse(is.na(combined), -1, combined)
        S4Vectors::mcols(consensus)$qValue <- -1
        S4Vectors::mcols(consensus)$peak <-
            round(GenomicRanges::width(consensus) / 2)
    }

    rtracklayer::export(consensus, con = file, format = format)
    written <- c(written, file)
    .messageIf(verbose, "Consensus regions written to ", file)

    if (isTRUE(perReplicate)) {
        stem <- tools::file_path_sans_ext(file)
        extension <- tools::file_ext(file)

        for (thisReplicate in names(object@peaks)) {
            replicateFile <- paste0(stem, "_", thisReplicate, ".",
                                    extension)
            peaks <- object@peaks[[thisReplicate]]
            S4Vectors::mcols(peaks)$name <-
                S4Vectors::mcols(peaks)$status
            S4Vectors::mcols(peaks)$score <-
                .scoreForExport(S4Vectors::mcols(peaks)$combinedNegLog10P)
            rtracklayer::export(peaks, con = replicateFile, format = "bed")
            written <- c(written, replicateFile)
        }
        .messageIf(verbose, "Per-replicate files written: ",
                   length(written) - 1)
    }

    invisible(written)
}


#' Squeeze a -log10 p-value into the BED score range
#'
#' @param negLog10P Numeric vector.
#'
#' @return Integer vector between 0 and 1000.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal
#' @noRd
.scoreForExport <- function(negLog10P) {
    if (is.null(negLog10P) || all(is.na(negLog10P))) {
        return(rep(0L, max(length(negLog10P), 1L)))
    }
    ## scale by ten so a p-value of 1e-100 lands at the ceiling rather
    ## than everything above 1e-1000 collapsing onto it
    capped <- pmin(round(negLog10P * 10), 1000)
    capped[is.na(capped)] <- 0
    as.integer(pmax(capped, 0))
}
