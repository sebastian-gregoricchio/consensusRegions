# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

## Internal helpers. None of these are exported.


#' Emit a message only when the caller asked for it
#'
#' @param verbose Single logical.
#' @param ... Passed on to [base::message()].
#'
#' @return Invisibly `NULL`, called for the side effect.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal
#' @noRd
.messageIf <- function(verbose, ...) {
    if (isTRUE(verbose)) {
        message(...)
    }
    invisible(NULL)
}


#' Turn arbitrary peak scores into -log10 p-values
#'
#' @description
#' Ranks the peaks of one replicate from best to worst and maps rank `r`
#' out of `n` onto the uniform quantile `r / (n + 1)`. Any monotonic score
#' therefore becomes a p-value on the scale the combination functions
#' expect. The absolute values mean little; what carries over is the order,
#' which is all a rank-based combination uses anyway.
#'
#' @param score Numeric vector of peak scores, larger meaning better.
#'
#' @return Numeric vector of -log10 transformed empirical p-values.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal
#' @noRd
.negLog10FromScore <- function(score) {
    if (length(score) == 0) {
        return(numeric(0))
    }
    if (all(is.na(score))) {
        stop("all peak scores are missing, cannot rank-transform them")
    }

    ## ties get the average rank so that equally scored peaks keep
    ## identical evidence rather than being split arbitrarily
    peakRank <- rank(-score, na.last = "keep", ties.method = "average")
    -log10(peakRank / (length(score) + 1))
}


#' Rescale weights to a mean of one
#'
#' @param weights Numeric vector, optionally named.
#'
#' @return Numeric vector of the same length with mean one.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal
#' @noRd
.normaliseWeights <- function(weights) {
    if (any(!is.finite(weights)) || any(weights <= 0)) {
        stop("weights must be finite and strictly positive")
    }
    ## centring on one keeps the weighted statistics on the same scale as
    ## their unweighted counterparts, so thresholds stay comparable
    weights / mean(weights)
}


#' Pairwise Jaccard index between peak sets
#'
#' @param peakList A `GRangesList`.
#'
#' @return A square numeric matrix of Jaccard indices.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges intersect union reduce width
#' @importFrom BiocGenerics unlist
#'
#' @keywords internal
#' @noRd
.pairwiseJaccard <- function(peakList) {
    n <- length(peakList)
    replicateNames <- names(peakList)
    jaccard <- matrix(NA_real_, nrow = n, ncol = n,
                      dimnames = list(replicateNames, replicateNames))

    ## collapsing each set first stops overlapping peaks inside a single
    ## replicate from inflating the covered base count
    flattened <- lapply(peakList, function(x) {
        GenomicRanges::reduce(x, ignore.strand = TRUE)
    })

    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            shared <- GenomicRanges::intersect(flattened[[i]], flattened[[j]],
                                               ignore.strand = TRUE)
            combined <- GenomicRanges::union(flattened[[i]], flattened[[j]],
                                             ignore.strand = TRUE)
            totalWidth <- sum(as.numeric(GenomicRanges::width(combined)))
            jaccard[i, j] <- if (totalWidth == 0) {
                0
            } else {
                sum(as.numeric(GenomicRanges::width(shared))) / totalWidth
            }
        }
    }

    jaccard
}


#' Guess chromosome lengths when the input carries none
#'
#' @param gr A `GRanges`.
#'
#' @return A named integer vector of sequence lengths.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomeInfoDb seqlengths seqnames
#' @importFrom GenomicRanges end
#' @importFrom S4Vectors split
#'
#' @keywords internal
#' @noRd
.inferSeqlengths <- function(gr) {
    known <- GenomeInfoDb::seqlengths(gr)
    if (!any(is.na(known))) {
        return(known)
    }

    ## without a genome to lean on the furthest peak end is the only
    ## available bound; it underestimates the chromosome, which makes the
    ## permutation slightly conservative
    lastEnd <- vapply(
        S4Vectors::split(GenomicRanges::end(gr),
                         as.character(GenomeInfoDb::seqnames(gr))),
        function(x) if (length(x) == 0) NA_real_ else max(x),
        numeric(1)
    )
    known[names(lastEnd)] <- as.integer(lastEnd)
    known
}


#' Shuffle intervals within their own chromosome
#'
#' @description
#' Keeps the widths and the per-chromosome counts of the original set and
#' redraws the start positions uniformly. Regions listed in
#' `excludeRegions` are avoided by redrawing the offending intervals a
#' bounded number of times.
#'
#' @param gr A `GRanges` to shuffle.
#' @param chromosomeLengths Named integer vector of chromosome lengths.
#' @param excludeRegions Optional `GRanges` of positions to avoid.
#' @param maxAttempts Number of redraw rounds allowed.
#'
#' @return A shuffled `GRanges`.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges GRanges width start seqnames countOverlaps
#' @importFrom IRanges IRanges
#' @importFrom stats runif
#'
#' @keywords internal
#' @noRd
.shuffleRanges <- function(gr,
                           chromosomeLengths,
                           excludeRegions = NULL,
                           maxAttempts = 20L) {
    peakWidth <- GenomicRanges::width(gr)
    chromosome <- as.character(GenomicRanges::seqnames(gr))
    upperBound <- chromosomeLengths[chromosome]

    ## fall back to the observed span when a chromosome has no length
    upperBound[is.na(upperBound)] <- max(peakWidth) + 1L

    drawStarts <- function(index) {
        span <- pmax(upperBound[index] - peakWidth[index], 1L)
        as.integer(stats::runif(length(index), min = 1, max = span))
    }

    newStart <- drawStarts(seq_along(gr))

    buildRanges <- function(starts) {
        GenomicRanges::GRanges(
            seqnames = chromosome,
            ranges = IRanges::IRanges(start = starts, width = peakWidth))
    }
    shuffled <- buildRanges(newStart)

    ## redraw only the intervals that landed on an excluded region, and
    ## give up quietly rather than looping forever on a crowded genome
    if (!is.null(excludeRegions)) {
        for (attempt in seq_len(maxAttempts)) {
            bad <- which(GenomicRanges::countOverlaps(
                shuffled, excludeRegions, ignore.strand = TRUE) > 0)
            if (length(bad) == 0) {
                break
            }
            newStart[bad] <- drawStarts(bad)
            shuffled <- buildRanges(newStart)
        }
    }

    shuffled
}


#' Convert a GRanges into a tibble of the columns used downstream
#'
#' @param gr A `GRanges`.
#'
#' @return A `tibble`.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom tibble as_tibble
#' @importFrom S4Vectors mcols
#'
#' @keywords internal
#' @noRd
.mcolsAsTibble <- function(gr) {
    tibble::as_tibble(as.data.frame(S4Vectors::mcols(gr)))
}


#' Accept a core count or a BiocParallel object
#'
#' @description
#' Naming a number of cores is how most people think about this, so a bare
#' number is accepted and turned into the right backend for the platform.
#' A `BiocParallelParam` object is passed through untouched, which is what
#' you need for anything beyond the core count, a fixed random seed above
#' all.
#'
#' @param BPPARAM A number of cores or a `BiocParallelParam` object.
#'
#' @return A `BiocParallelParam` object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom BiocParallel SerialParam MulticoreParam SnowParam
#' @importFrom methods is
#'
#' @keywords internal
#' @noRd
.resolveBPPARAM <- function(BPPARAM) {
    if (methods::is(BPPARAM, "BiocParallelParam")) {
        return(BPPARAM)
    }

    if (!is.numeric(BPPARAM) || length(BPPARAM) != 1 ||
        !is.finite(BPPARAM) || BPPARAM < 1) {
        stop("'BPPARAM' must be a number of cores, or a ",
             "BiocParallelParam object")
    }

    nCores <- as.integer(BPPARAM)
    if (nCores == 1L) {
        return(BiocParallel::SerialParam())
    }

    ## forking is unavailable on Windows, where a socket cluster is the
    ## equivalent arrangement
    if (.Platform$OS.type == "windows") {
        BiocParallel::SnowParam(workers = nCores)
    } else {
        BiocParallel::MulticoreParam(workers = nCores)
    }
}
