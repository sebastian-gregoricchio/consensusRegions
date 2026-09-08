#' Derive per-replicate weights
#'
#' @description
#' Fisher's method assumes every replicate is worth the same. In practice
#' one of them is usually shallower, or has a worse signal-to-noise ratio,
#' and the usual responses are to drop it or to let it drag the combined
#' statistic down. Weighting is the middle option: the replicate still
#' contributes, just less.
#'
#' Weights can come from three places, in decreasing order of how much
#' extra material they require. FRiP and library size need the BAM files
#' or the numbers already computed from them. The intrinsic route needs
#' nothing beyond the peaks themselves and scores each replicate by how
#' well it agrees with the others, which is the only option left when all
#' that survives of an old analysis is a folder of BED files.
#'
#' @param peakList A `GRangesList`, typically from [readPeakSets()].
#' @param method One of `"equal"`, `"frip"`, `"librarySize"`,
#'   `"intrinsic"` or `"custom"`.
#' @param bamFiles Character vector of indexed BAM files, in the same
#'   order as `peakList`. Needed by `"frip"` and `"librarySize"` unless
#'   the values are supplied directly.
#' @param frip Numeric vector of pre-computed fractions of reads in
#'   peaks. Skips the BAM pass.
#' @param librarySize Numeric vector of pre-computed mapped read counts.
#' @param weights Numeric vector used as is when `method` is `"custom"`.
#' @param minMapq Minimum mapping quality when counting from BAM.
#' @param verbose Report progress.
#'
#' @return A named numeric vector with mean one, carrying the raw metrics
#'   as the `metrics` attribute.
#'
#' @details
#' For FRiP the weight is proportional to the fraction itself. For library
#' size it is proportional to the square root of the depth, which is how
#' the information in a z-score scales, so doubling the depth is worth
#' rather less than twice as much.
#'
#' The intrinsic weight is the mean Jaccard index between one replicate
#' and each of the others. A replicate that shares little with the rest is
#' either the odd one out biologically or the noisy one, and in both cases
#' letting it drive a consensus is unwise. This measure is circular by
#' construction, so prefer FRiP whenever the alignments are available.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom Rsamtools idxstatsBam countBam ScanBamParam BamFile
#' @importFrom GenomicRanges reduce
#' @importFrom dplyr tibble filter pull
#' @importFrom rlang .data
#'
#' @examples
#' peakFiles <- system.file("extdata",
#'                          c("rep1.narrowPeak", "rep2.narrowPeak",
#'                            "rep3.narrowPeak"),
#'                          package = "consensusRegions")
#' peaks <- readPeakSets(peakFiles, sampleNames = c("r1", "r2", "r3"),
#'                       verbose = FALSE)
#'
#' ## no BAM files at hand
#' computeReplicateWeights(peaks, method = "intrinsic")
#'
#' ## FRiP measured elsewhere
#' computeReplicateWeights(peaks, method = "frip",
#'                         frip = c(0.21, 0.19, 0.06))
#'
#' @export
computeReplicateWeights <- function(peakList,
                                    method = c("equal", "frip",
                                               "librarySize", "intrinsic",
                                               "custom"),
                                    bamFiles = NULL,
                                    frip = NULL,
                                    librarySize = NULL,
                                    weights = NULL,
                                    minMapq = 0L,
                                    verbose = TRUE) {
    method <- match.arg(method)
    replicateNames <- names(peakList)

    if (is.null(replicateNames)) {
        stop("'peakList' must be named")
    }

    ## equal weights keep the analysis reproducible for anyone who does
    ## not opt into weighting
    if (method == "equal") {
        result <- rep(1, length(peakList))
        names(result) <- replicateNames
        return(result)
    }

    if (method == "custom") {
        if (is.null(weights)) {
            stop("'weights' must be supplied when method is 'custom'")
        }
        if (length(weights) != length(peakList)) {
            stop("'weights' has length ", length(weights), " but ",
                 length(peakList), " replicates were supplied")
        }
        result <- .normaliseWeights(weights)
        names(result) <- replicateNames
        return(result)
    }

    if (method == "intrinsic") {
        .messageIf(verbose, "Scoring replicates by mutual agreement")
        jaccard <- .pairwiseJaccard(peakList)
        diag(jaccard) <- NA_real_

        ## the mean overlap with the other replicates stands in for
        ## quality when nothing better is available
        agreement <- rowMeans(jaccard, na.rm = TRUE)
        if (any(agreement <= 0)) {
            stop("at least one replicate shares no peaks with the others, ",
                 "intrinsic weights cannot be computed")
        }
        result <- .normaliseWeights(agreement)
        names(result) <- replicateNames
        attr(result, "metrics") <- dplyr::tibble(replicate = replicateNames,
                                                 jaccard = agreement)
        return(result)
    }

    ## from here on the metrics either come pre-computed or are read off
    ## the alignments
    measured <- .collectBamMetrics(peakList = peakList,
                                   required = method,
                                   bamFiles = bamFiles,
                                   frip = frip,
                                   librarySize = librarySize,
                                   minMapq = minMapq,
                                   verbose = verbose)

    if (method == "frip") {
        usable <- dplyr::filter(measured, !is.na(.data$frip))
        if (nrow(usable) != nrow(measured)) {
            stop("FRiP is missing for: ",
                 paste(dplyr::pull(dplyr::filter(measured, is.na(.data$frip)),
                                   .data$replicate), collapse = ", "))
        }
        result <- .normaliseWeights(measured$frip)
    } else {
        usable <- dplyr::filter(measured, !is.na(.data$librarySize))
        if (nrow(usable) != nrow(measured)) {
            stop("library size is missing for at least one replicate")
        }
        ## information in a z-score grows with the square root of depth
        result <- .normaliseWeights(sqrt(measured$librarySize))
    }

    names(result) <- replicateNames
    attr(result, "metrics") <- measured
    result
}


#' Gather FRiP and depth, reading the BAM files only when needed
#'
#' @param peakList A `GRangesList`.
#' @param required Which metric the caller needs, `"frip"` or
#'   `"librarySize"`.
#' @param bamFiles Character vector of BAM paths or `NULL`.
#' @param frip Pre-computed FRiP values or `NULL`.
#' @param librarySize Pre-computed depths or `NULL`.
#' @param minMapq Minimum mapping quality.
#' @param verbose Report progress.
#'
#' @return A `tibble` with one row per replicate.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom Rsamtools idxstatsBam countBam ScanBamParam BamFile
#' @importFrom GenomicRanges reduce
#' @importFrom dplyr tibble
#'
#' @keywords internal
#' @noRd
.collectBamMetrics <- function(peakList,
                               required = c("frip", "librarySize"),
                               bamFiles,
                               frip,
                               librarySize,
                               minMapq = 0L,
                               verbose = TRUE) {
    required <- match.arg(required)
    replicateNames <- names(peakList)
    nReplicates <- length(peakList)

    ## supplied numbers win over anything that would have to be measured
    if (!is.null(frip) && length(frip) != nReplicates) {
        stop("'frip' has length ", length(frip), " but ", nReplicates,
             " replicates were supplied")
    }
    if (!is.null(librarySize) && length(librarySize) != nReplicates) {
        stop("'librarySize' has length ", length(librarySize), " but ",
             nReplicates, " replicates were supplied")
    }

    fripValues <- if (is.null(frip)) rep(NA_real_, nReplicates) else frip
    depthValues <- if (is.null(librarySize)) {
        rep(NA_real_, nReplicates)
    } else {
        librarySize
    }

    ## a FRiP-weighted run has no use for the library size, so asking for
    ## BAM files just because that column is empty would be wrong
    missingFrip <- required == "frip" && any(is.na(fripValues))
    missingDepth <- required == "librarySize" && any(is.na(depthValues))
    needsBam <- missingFrip || missingDepth

    if (needsBam) {
        if (is.null(bamFiles)) {
            stop("supply 'bamFiles', or pass 'frip' and 'librarySize' ",
                 "directly if the alignments are no longer available")
        }
        if (length(bamFiles) != nReplicates) {
            stop("'bamFiles' has length ", length(bamFiles), " but ",
                 nReplicates, " replicates were supplied")
        }
        missingFiles <- bamFiles[!file.exists(bamFiles)]
        if (length(missingFiles) > 0) {
            stop("BAM file not found: ",
                 paste(missingFiles, collapse = ", "))
        }

        for (i in seq_len(nReplicates)) {
            .messageIf(verbose, "Counting reads for ", replicateNames[i])

            ## the index carries the mapped totals, so the whole file
            ## never has to be traversed for the denominator; FRiP needs
            ## it too, as the denominator of the fraction
            if (is.na(depthValues[i]) &&
                (missingDepth || is.na(fripValues[i]))) {
                indexStats <- Rsamtools::idxstatsBam(bamFiles[i])
                depthValues[i] <- sum(indexStats$mapped)
            }

            if (missingFrip && is.na(fripValues[i])) {
                ## flattening the peaks first stops a read that spans two
                ## overlapping calls from being counted twice
                targets <- GenomicRanges::reduce(peakList[[i]],
                                                 ignore.strand = TRUE)
                countParam <- Rsamtools::ScanBamParam(which = targets,
                                                      mapqFilter = minMapq)
                counted <- Rsamtools::countBam(
                    Rsamtools::BamFile(bamFiles[i]), param = countParam)
                fripValues[i] <- sum(counted$records) / depthValues[i]
            }
        }
    }

    outOfRange <- which(!is.na(fripValues) &
                            (fripValues <= 0 | fripValues > 1))
    if (length(outOfRange) > 0) {
        stop("FRiP outside (0, 1] for: ",
             paste(replicateNames[outOfRange], collapse = ", "))
    }

    dplyr::tibble(replicate = replicateNames,
                  frip = fripValues,
                  librarySize = depthValues)
}
