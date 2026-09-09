#' Build consensus regions from replicated peak calls
#'
#' @description
#' Runs the whole analysis. Each peak is classified as stringent, weak or
#' background by two p-value thresholds; the peaks overlapping it in the
#' other replicates are collected; their evidence is combined; and the
#' result is tested against a threshold on the combined significance and
#' then corrected within each replicate. Peaks that survive are merged
#' into consensus regions.
#'
#' Feed this permissive input. Calling peaks at `q < 0.05` and running the
#' consensus on the survivors leaves nothing to rescue and turns the whole
#' exercise into an intersection with extra steps. Something around
#' `p < 1e-3` gives the combination room to work.
#'
#' @param peakList A `GRangesList` from [readPeakSets()].
#' @param weights Named numeric vector of replicate weights. See
#'   [computeReplicateWeights()]. Default: \code{NULL}, meaning that every
#'   replicate counts the same.
#' @param combinationMethod Passed to [combineEvidence()].
#'   Default: \code{"stouffer"}.
#' @param replicateType `"biological"` asks for `minSupport` supporting
#'   replicates; `"technical"` asks for all of them.
#'   Default: \code{"biological"}.
#' @param stringencyThreshold P-value below which a peak is stringent.
#'   Default: \code{1e-08}.
#' @param weakThreshold P-value above which a peak is background and takes
#'   no further part. Default: \code{1e-04}.
#' @param combinedThreshold Threshold on the combined p-value.
#'   [calibrateThreshold()] can estimate it instead. The scale differs
#'   between combination schemes, so a value carried over from one will
#'   not mean the same under another. The rank product in particular is
#'   bounded by the number of peaks per replicate and usually needs a far
#'   more permissive threshold. Default: \code{NULL}, meaning that
#'   \code{stringencyThreshold} is used.
#' @param alpha Level for the within-replicate Benjamini-Hochberg step.
#'   Default: \code{0.05}.
#' @param minSupport Number of *other* replicates that must hold an
#'   overlapping peak. This is a count of replicates and weights never
#'   substitute for it. Give this or `minReplicates`, not both.
#'   Default: \code{NULL}, meaning that one supporting replicate is
#'   enough.
#' @param minReplicates How many replicates in total must hold the peak,
#'   counting the one it came from. This is the way the requirement is
#'   usually spoken about, and it matches MSPC's `-c`. Accepts a count
#'   (`5`), a proportion (`0.75`) or a percentage (`"75%"`), the last two
#'   rounded up. `minReplicates = 2` and `minSupport = 1` are the same
#'   requirement. Default: \code{NULL}.
#' @param minSupportWeight Optional additional requirement on the summed
#'   weight of the supporting replicates. Applied on top of `minSupport`,
#'   never instead of it. Only meaningful when weights are not equal.
#'   Default: \code{NULL}.
#' @param adjustmentFamily Which peaks form the family for the
#'   Benjamini-Hochberg step. `"tested"` corrects across every peak that
#'   entered the analysis; `"confirmed"` corrects only across the peaks
#'   that cleared `combinedThreshold`, which is what MSPC does.
#'   Default: \code{"tested"}.
#' @param minOverlap Minimum overlap in base pairs. Default: \code{1L}.
#' @param minOverlapFraction Optional minimum overlap as a fraction of the
#'   shorter of the two peaks. Applied on top of `minOverlap`.
#'   Default: \code{NULL}.
#' @param multipleIntersections How to pick among several overlapping
#'   peaks from the same replicate: `"lowest"` takes the smallest p-value,
#'   `"highest"` the largest. Default: \code{"lowest"}.
#' @param recursive Re-run the confirmation after dropping discarded peaks
#'   from the pool of eligible supporters. Default: \code{TRUE}.
#' @param maxIterations Cap on the recursive rounds. Default: \code{10L}.
#' @param mergeMethod `"reduce"` merges anything that touches;
#'   `"iterative"` seeds on the most significant peak and removes what it
#'   overlaps, then repeats. Default: \code{"reduce"}.
#' @param maxConsensusWidth Optional cap in base pairs. Merged regions
#'   wider than this are rebuilt with the iterative rule.
#'   Default: \code{NULL}.
#' @param calibration Optional output of [calibrateThreshold()]. When
#'   supplied its threshold is used, overriding `combinedThreshold`, and
#'   the whole calibration is kept with the result. Default: \code{NULL}.
#' @param excludeRegions Optional `GRanges` of blacklisted positions,
#'   removed from the consensus at the end. Default: \code{NULL}.
#' @param verbose Report progress. Default: \code{TRUE}.
#'
#' @return A [ConsensusRegions-class] object.
#'
#' @details
#' Two behaviours are worth knowing about.
#'
#' Transitive merging can run away. Peak A overlaps B, B overlaps C, and
#' A never touches C, yet `"reduce"` puts all three in one region. On
#' permissive input this occasionally produces regions tens of kilobases
#' long. `maxConsensusWidth` catches them and rebuilds only the offenders.
#'
#' Recursive confirmation makes a peak's fate depend on peaks that were
#' themselves discarded, which is arguably the more defensible reading of
#' the method but is not what the original implementation does. Set
#' `recursive = FALSE` to reproduce a single-pass analysis.
#'
#' When the peak sets carry no usable statistic the combination is skipped
#' and a weighted presence rule is applied instead: a peak is kept when
#' the weights of the replicates supporting it sum to at least
#' `minSupport`. With equal weights this is the familiar k-of-n rule.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges findOverlaps reduce pintersect width mcols
#'   GRanges seqnames countOverlaps
#' @importFrom S4Vectors queryHits subjectHits mcols mcols<- metadata
#' @importFrom BiocGenerics unlist
#' @importFrom IRanges subsetByOverlaps
#' @importFrom dplyr tibble filter group_by summarise mutate select
#'   slice_max slice_min left_join bind_rows arrange n_distinct pull
#'   distinct ungroup case_when
#' @importFrom rlang .data
#' @importFrom methods new
#'
#' @examples
#' peakFiles <- system.file("extdata",
#'                          c("rep1.narrowPeak", "rep2.narrowPeak",
#'                            "rep3.narrowPeak"),
#'                          package = "consensusRegions")
#' peaks <- readPeakSets(peakFiles, sampleNames = c("r1", "r2", "r3"),
#'                       verbose = FALSE)
#'
#' result <- buildConsensus(peaks, verbose = FALSE)
#' consensusRanges(result)
#'
#' ## down-weighting a replicate known to be poor
#' weighted <- buildConsensus(
#'     peaks,
#'     weights = computeReplicateWeights(peaks, method = "frip",
#'                                       frip = c(0.2, 0.18, 0.05)),
#'     combinationMethod = "stouffer",
#'     verbose = FALSE)
#' consensusStats(weighted)
#'
#' @export
buildConsensus <- function(peakList,
                           weights = NULL,
                           combinationMethod = c("stouffer", "fisher",
                                                 "lancaster", "rankProduct"),
                           replicateType = c("biological", "technical"),
                           stringencyThreshold = 1e-8,
                           weakThreshold = 1e-4,
                           combinedThreshold = NULL,
                           alpha = 0.05,
                           minSupport = NULL,
                           minReplicates = NULL,
                           minSupportWeight = NULL,
                           adjustmentFamily = c("tested", "confirmed"),
                           minOverlap = 1L,
                           minOverlapFraction = NULL,
                           multipleIntersections = c("lowest", "highest"),
                           recursive = TRUE,
                           maxIterations = 10L,
                           mergeMethod = c("reduce", "iterative"),
                           maxConsensusWidth = NULL,
                           calibration = NULL,
                           excludeRegions = NULL,
                           verbose = TRUE) {
    combinationMethod <- match.arg(combinationMethod)
    replicateType <- match.arg(replicateType)
    adjustmentFamily <- match.arg(adjustmentFamily)
    multipleIntersections <- match.arg(multipleIntersections)
    mergeMethod <- match.arg(mergeMethod)

    ## ---- settings ----------------------------------------------------
    if (!methods::is(peakList, "GRangesList")) {
        stop("'peakList' must be a GRangesList, see readPeakSets()")
    }
    if (length(peakList) < 2) {
        stop("consensus needs at least two replicates")
    }
    if (weakThreshold < stringencyThreshold) {
        stop("'weakThreshold' must be larger than 'stringencyThreshold'; ",
             "they are p-values, so the weak cut is the more permissive one")
    }
    ## a calibrated threshold was derived from these peak sets, so it
    ## takes precedence over anything carried over by hand
    if (!is.null(calibration)) {
        if (!is.list(calibration) || is.null(calibration$threshold)) {
            stop("'calibration' must be the output of calibrateThreshold()")
        }
        combinedThreshold <- calibration$threshold
        .messageIf(verbose, "Using the calibrated threshold ",
                   signif(combinedThreshold, 3))
    }
    if (is.null(combinedThreshold)) {
        combinedThreshold <- stringencyThreshold
    }

    nReplicates <- length(peakList)
    replicateNames <- names(peakList)

    requiredSupport <- if (replicateType == "technical") {
        ## technical replicates are expected to agree everywhere
        nReplicates - 1L
    } else {
        .resolveRequiredSupport(minSupport = minSupport,
                                minReplicates = minReplicates,
                                nReplicates = nReplicates)
    }
    if (requiredSupport > nReplicates - 1L) {
        stop("the requested support asks for ", requiredSupport + 1L,
             " replicates in total but only ", nReplicates,
             " were supplied")
    }

    ## ---- weights -----------------------------------------------------
    if (is.null(weights)) {
        weights <- rep(1, nReplicates)
        names(weights) <- replicateNames
    }
    if (!setequal(names(weights), replicateNames)) {
        stop("names of 'weights' do not match the replicate names")
    }
    weights <- .normaliseWeights(weights[replicateNames])

    ## ---- score type --------------------------------------------------
    scoreType <- S4Vectors::metadata(peakList)$scoreType
    if (is.null(scoreType)) {
        scoreType <- .detectScoreType(as.list(peakList))
    }
    presenceOnly <- identical(scoreType, "none")

    if (presenceOnly) {
        .messageIf(verbose,
                   "No usable statistic found, falling back to the ",
                   "weighted presence rule")
    }

    ## ---- flatten -----------------------------------------------------
    flatPeaks <- BiocGenerics::unlist(peakList, use.names = FALSE)
    if (length(flatPeaks) == 0) {
        stop("the peak sets are empty")
    }
    S4Vectors::mcols(flatPeaks)$replicate <-
        rep(replicateNames, lengths(peakList))
    S4Vectors::mcols(flatPeaks)$weight <-
        unname(weights[S4Vectors::mcols(flatPeaks)$replicate])

    ## ---- classification ----------------------------------------------
    flatPeaks <- .classifyPeaks(flatPeaks,
                                stringencyThreshold = stringencyThreshold,
                                weakThreshold = weakThreshold,
                                presenceOnly = presenceOnly)

    retained <- flatPeaks[S4Vectors::mcols(flatPeaks)$class != "background"]
    if (length(retained) == 0) {
        stop("every peak fell below 'weakThreshold'; the input is likely ",
             "already filtered, call peaks more permissively")
    }
    .messageIf(verbose, "Peaks entering the analysis: ", length(retained),
               " of ", length(flatPeaks))

    ## within-replicate relative ranks, needed by the rank product
    S4Vectors::mcols(retained)$rho <- .relativeRanks(
        negLog10P = S4Vectors::mcols(retained)$negLog10P,
        replicate = S4Vectors::mcols(retained)$replicate,
        presenceOnly = presenceOnly)

    ## the rank product is bounded by the number of peaks per replicate,
    ## so the usual 1e-8 is often unreachable and would silently return
    ## nothing at all
    if (!presenceOnly && combinationMethod == "rankProduct") {
        reachable <- .rankProductCeiling(retained, requiredSupport)
        if (-log10(combinedThreshold) > reachable) {
            stop("the rank product cannot reach a combined p-value of ",
                 combinedThreshold, " with these peak sets: the smallest ",
                 "attainable is ", signif(10^(-reachable), 3),
                 ", because the statistic is bounded by the number of ",
                 "peaks per replicate. Lower 'combinedThreshold', or let ",
                 "calibrateThreshold() choose one")
        }
    }

    ## ---- overlap graph -----------------------------------------------
    .messageIf(verbose, "Collecting cross-replicate overlaps")
    overlapTable <- .buildOverlapTable(
        retained,
        minOverlap = minOverlap,
        minOverlapFraction = minOverlapFraction)

    ## ---- confirmation ------------------------------------------------
    verdict <- .runConfirmation(
        retained = retained,
        overlapTable = overlapTable,
        weights = weights,
        combinationMethod = combinationMethod,
        combinedThreshold = combinedThreshold,
        requiredSupport = requiredSupport,
        minSupportWeight = minSupportWeight,
        multipleIntersections = multipleIntersections,
        recursive = recursive,
        maxIterations = maxIterations,
        presenceOnly = presenceOnly,
        verbose = verbose)

    S4Vectors::mcols(retained)$nSupport <- verdict$nSupport
    S4Vectors::mcols(retained)$supportWeight <- verdict$supportWeight
    S4Vectors::mcols(retained)$combinedNegLog10P <- verdict$combinedNegLog10P

    ## ---- within-replicate correction ---------------------------------
    adjusted <- .adjustWithinReplicates(
        combinedNegLog10P = verdict$combinedNegLog10P,
        passed = verdict$passed,
        replicate = S4Vectors::mcols(retained)$replicate,
        presenceOnly = presenceOnly,
        adjustmentFamily = adjustmentFamily)

    S4Vectors::mcols(retained)$combinedNegLog10Padj <- adjusted
    keptPeak <- verdict$passed &
        (presenceOnly | adjusted >= -log10(alpha))

    S4Vectors::mcols(retained)$status <- dplyr::case_when(
        keptPeak ~ "confirmed",
        verdict$passed ~ "falsePositive",
        .default = "discarded")

    .messageIf(verbose, "Peaks confirmed: ", sum(keptPeak))

    ## ---- consensus ---------------------------------------------------
    consensus <- .assembleConsensus(
        retained[keptPeak],
        combinationMethod = combinationMethod,
        weights = weights,
        mergeMethod = mergeMethod,
        maxConsensusWidth = maxConsensusWidth,
        multipleIntersections = multipleIntersections,
        presenceOnly = presenceOnly)

    ## blacklisting belongs here and not earlier: removing regions before
    ## the combination would distort the rank distributions it relies on
    if (!is.null(excludeRegions)) {
        before <- length(consensus)
        consensus <- IRanges::subsetByOverlaps(consensus, excludeRegions,
                                               invert = TRUE,
                                               ignore.strand = TRUE)
        .messageIf(verbose, "Excluded regions removed: ",
                   before - length(consensus))
    }

    .messageIf(verbose, "Consensus regions: ", length(consensus))

    ## ---- assemble the object -----------------------------------------
    annotatedPeaks <- S4Vectors::split(retained,
                                       S4Vectors::mcols(retained)$replicate)
    annotatedPeaks <- annotatedPeaks[replicateNames]

    methods::new(
        "ConsensusRegions",
        peaks = annotatedPeaks,
        consensus = consensus,
        weights = weights,
        parameters = list(
            scoreType = scoreType,
            combinationMethod = if (presenceOnly) NA_character_ else
                combinationMethod,
            replicateType = replicateType,
            stringencyThreshold = stringencyThreshold,
            weakThreshold = weakThreshold,
            combinedThreshold = combinedThreshold,
            alpha = alpha,
            minSupport = requiredSupport,
            minReplicates = requiredSupport + 1L,
            minSupportWeight = minSupportWeight,
            adjustmentFamily = adjustmentFamily,
            minOverlap = minOverlap,
            minOverlapFraction = minOverlapFraction,
            multipleIntersections = multipleIntersections,
            recursive = recursive,
            mergeMethod = mergeMethod,
            maxConsensusWidth = maxConsensusWidth,
            presenceOnly = presenceOnly),
        stats = .replicateSummary(retained, replicateNames),
        calibration = if (is.null(calibration)) list() else calibration
    )
}


#' Work out how many supporting replicates are required
#'
#' @description
#' The requirement can be given either way round. `minSupport` counts the
#' other replicates, which is what the confirmation step works with;
#' `minReplicates` counts all of them including the peak's own, which is
#' how people describe it out loud and how MSPC's `-c` is defined. This
#' resolves whichever was supplied to the former.
#'
#' @param minSupport Supporting replicates, or `NULL`.
#' @param minReplicates Total replicates, a count, a proportion or a
#'   percentage string, or `NULL`.
#' @param nReplicates Number of replicates in the experiment.
#'
#' @return Integer, the number of supporting replicates required.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal
#' @noRd
.resolveRequiredSupport <- function(minSupport, minReplicates,
                                    nReplicates) {
    if (!is.null(minSupport) && !is.null(minReplicates)) {
        stop("give either 'minSupport' or 'minReplicates', not both: ",
             "they state the same requirement, one counting the other ",
             "replicates and one counting all of them")
    }

    if (is.null(minReplicates)) {
        return(as.integer(if (is.null(minSupport)) 1L else minSupport))
    }

    total <- .parseReplicateCount(minReplicates, nReplicates)
    if (total < 2L) {
        stop("'minReplicates' resolved to ", total,
             "; a peak supported by no other replicate is not a consensus")
    }

    ## the peak's own replicate is one of the total, the rest must support
    total - 1L
}


#' Read a replicate requirement given as a count, a proportion or a
#' percentage
#'
#' @param value The user supplied requirement.
#' @param nReplicates Number of replicates in the experiment.
#'
#' @return Integer count of replicates.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal
#' @noRd
.parseReplicateCount <- function(value, nReplicates) {
    if (is.character(value)) {
        trimmed <- trimws(value)

        ## check the shape first rather than converting and inspecting
        ## the damage, which would need the warning to be muffled
        if (!grepl("^[0-9]*[.]?[0-9]+[%]?$", trimmed)) {
            stop("could not read '", value,
                 "' as a count, a proportion or a percentage")
        }
        number <- as.numeric(sub("%$", "", trimmed))

        ## a percentage is always a share of the experiment, including
        ## when it is the whole of it; falling through to the branch
        ## below would read "100%" as a single replicate
        if (grepl("%$", trimmed)) {
            if (number <= 0 || number > 100) {
                stop("'", value, "' is not a percentage between 0 and 100")
            }
            return(as.integer(ceiling(number / 100 * nReplicates)))
        }
    } else {
        number <- as.numeric(value)
    }

    if (!is.finite(number) || number <= 0) {
        stop("the replicate requirement must be a positive number")
    }

    ## anything below one is read as a share of the experiment, rounded
    ## up so that 0.75 of four replicates asks for three and not two
    if (number < 1) {
        return(as.integer(ceiling(number * nReplicates)))
    }
    as.integer(round(number))
}


#' Split peaks into stringent, weak and background
#'
#' @param gr Flattened `GRanges` of peaks.
#' @param stringencyThreshold Stringent p-value cut.
#' @param weakThreshold Background p-value cut.
#' @param presenceOnly Skip the classification entirely.
#'
#' @return The `GRanges` with a `class` metadata column.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom dplyr case_when
#'
#' @keywords internal
#' @noRd
.classifyPeaks <- function(gr, stringencyThreshold, weakThreshold,
                           presenceOnly) {
    if (presenceOnly) {
        ## nothing to threshold on, every peak stays in
        S4Vectors::mcols(gr)$class <- rep("stringent", length(gr))
        return(gr)
    }

    negLog10P <- S4Vectors::mcols(gr)$negLog10P
    if (any(is.na(negLog10P))) {
        stop("missing values in 'negLog10P'; some peaks carry no statistic")
    }

    S4Vectors::mcols(gr)$class <- dplyr::case_when(
        negLog10P >= -log10(stringencyThreshold) ~ "stringent",
        negLog10P >= -log10(weakThreshold) ~ "weak",
        .default = "background")
    gr
}


#' Relative within-replicate ranks
#'
#' @param negLog10P Numeric vector.
#' @param replicate Character vector of replicate labels.
#' @param presenceOnly Return `NA` when there is nothing to rank.
#'
#' @return Numeric vector inside `(0, 1]`.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr tibble group_by mutate ungroup pull n
#' @importFrom rlang .data
#'
#' @keywords internal
#' @noRd
.relativeRanks <- function(negLog10P, replicate, presenceOnly) {
    if (presenceOnly) {
        return(rep(NA_real_, length(replicate)))
    }

    ranked <- dplyr::tibble(negLog10P = negLog10P, replicate = replicate)
    ranked <- dplyr::group_by(ranked, .data$replicate)
    ranked <- dplyr::mutate(
        ranked,
        rho = rank(-.data$negLog10P, ties.method = "average") /
            (dplyr::n() + 1))
    dplyr::pull(dplyr::ungroup(ranked), .data$rho)
}


#' Smallest combined p-value the rank product could ever produce
#'
#' @description
#' The rank product works on relative ranks, so the best a peak can do is
#' come first in every replicate. That bound depends on how many peaks
#' each replicate holds, which is why a threshold borrowed from a
#' p-value-based scheme often cannot be met.
#'
#' @param retained Flattened `GRanges` of retained peaks.
#' @param requiredSupport Minimum supporting replicates.
#'
#' @return The attainable maximum on the -log10 scale.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom S4Vectors mcols split
#' @importFrom stats pgamma
#'
#' @keywords internal
#' @noRd
.rankProductCeiling <- function(retained, requiredSupport) {
    metadataColumns <- S4Vectors::mcols(retained)
    bestPerReplicate <- vapply(
        split(metadataColumns$rho, metadataColumns$replicate),
        min, numeric(1))

    ## a peak plus the replicates that have to support it
    nMembers <- min(requiredSupport + 1L, length(bestPerReplicate))
    bestMembers <- sort(bestPerReplicate)[seq_len(nMembers)]

    -stats::pgamma(sum(-log(bestMembers)), shape = nMembers, rate = 1,
                   lower.tail = FALSE, log.p = TRUE) / log(10)
}


#' Cross-replicate overlaps as a table
#'
#' @param retained Flattened `GRanges` of retained peaks.
#' @param minOverlap Minimum overlap in base pairs.
#' @param minOverlapFraction Optional fractional requirement.
#'
#' @return A `tibble` with `queryHits`, `subjectHits` and the subject
#'   replicate, weight and statistics.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges findOverlaps pintersect width
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom dplyr tibble filter mutate
#' @importFrom rlang .data
#'
#' @keywords internal
#' @noRd
.buildOverlapTable <- function(retained, minOverlap, minOverlapFraction) {
    hits <- GenomicRanges::findOverlaps(retained, retained,
                                        minoverlap = max(minOverlap, 1L),
                                        ignore.strand = TRUE)

    metadataColumns <- S4Vectors::mcols(retained)
    queryIndex <- S4Vectors::queryHits(hits)
    subjectIndex <- S4Vectors::subjectHits(hits)

    overlaps <- dplyr::tibble(
        queryHits = queryIndex,
        subjectHits = subjectIndex,
        queryReplicate = metadataColumns$replicate[queryIndex],
        subjectReplicate = metadataColumns$replicate[subjectIndex],
        subjectNegLog10P = metadataColumns$negLog10P[subjectIndex],
        subjectRho = metadataColumns$rho[subjectIndex],
        subjectWeight = metadataColumns$weight[subjectIndex])

    ## a peak cannot support itself, and neither can its neighbours from
    ## the same replicate
    overlaps <- dplyr::filter(
        overlaps, .data$queryReplicate != .data$subjectReplicate)

    if (!is.null(minOverlapFraction)) {
        if (minOverlapFraction <= 0 || minOverlapFraction > 1) {
            stop("'minOverlapFraction' must lie inside (0, 1]")
        }
        ## measured against the shorter peak, so a narrow summit inside a
        ## broad domain still counts as supported
        sharedWidth <- GenomicRanges::width(
            GenomicRanges::pintersect(retained[overlaps$queryHits],
                                      retained[overlaps$subjectHits],
                                      ignore.strand = TRUE))
        shorterWidth <- pmin(
            GenomicRanges::width(retained[overlaps$queryHits]),
            GenomicRanges::width(retained[overlaps$subjectHits]))
        overlaps <- dplyr::filter(
            dplyr::mutate(overlaps,
                          overlapFraction = sharedWidth / shorterWidth),
            .data$overlapFraction >= minOverlapFraction)
    }

    overlaps
}


#' Confirmation loop
#'
#' @param retained Flattened `GRanges`.
#' @param overlapTable Output of `.buildOverlapTable`.
#' @param weights Named numeric vector.
#' @param combinationMethod Combination scheme.
#' @param combinedThreshold Threshold on the combined p-value.
#' @param requiredSupport Minimum supporting replicates.
#' @param minSupportWeight Optional extra requirement on summed support
#'   weight.
#' @param multipleIntersections Resolution rule.
#' @param recursive Re-evaluate after discarding.
#' @param maxIterations Cap on rounds.
#' @param presenceOnly Skip the combination.
#' @param verbose Report progress.
#'
#' @return A list with `nSupport`, `supportWeight`, `combinedNegLog10P`
#'   and `passed`.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr filter group_by slice_max slice_min ungroup summarise
#'   tibble bind_rows left_join select mutate n_distinct
#' @importFrom S4Vectors mcols
#' @importFrom rlang .data
#'
#' @keywords internal
#' @noRd
.runConfirmation <- function(retained,
                             overlapTable,
                             weights,
                             combinationMethod,
                             combinedThreshold,
                             requiredSupport,
                             minSupportWeight = NULL,
                             multipleIntersections,
                             recursive,
                             maxIterations,
                             presenceOnly,
                             verbose) {
    nPeaks <- length(retained)
    metadataColumns <- S4Vectors::mcols(retained)
    combinedCut <- -log10(combinedThreshold)

    ## every peak starts out eligible to support the others
    eligible <- rep(TRUE, nPeaks)
    previousPassed <- NULL
    result <- NULL

    for (iteration in seq_len(max(maxIterations, 1L))) {
        currentOverlaps <- dplyr::filter(overlapTable,
                                         eligible[.data$subjectHits])

        ## keep one supporter per replicate, the most or the least
        ## significant depending on the rule asked for
        selected <- if (presenceOnly) {
            ## nothing to rank on, so any one overlapping peak per
            ## replicate carries the same information
            dplyr::distinct(currentOverlaps, .data$queryHits,
                            .data$subjectReplicate, .keep_all = TRUE)
        } else {
            grouped <- dplyr::group_by(currentOverlaps, .data$queryHits,
                                       .data$subjectReplicate)
            picked <- if (multipleIntersections == "lowest") {
                dplyr::slice_max(grouped, order_by = .data$subjectNegLog10P,
                                 n = 1, with_ties = FALSE)
            } else {
                dplyr::slice_min(grouped, order_by = .data$subjectNegLog10P,
                                 n = 1, with_ties = FALSE)
            }
            dplyr::ungroup(picked)
        }

        support <- dplyr::summarise(
            dplyr::group_by(selected, .data$queryHits),
            nSupport = dplyr::n_distinct(.data$subjectReplicate),
            supportWeight = sum(.data$subjectWeight),
            .groups = "drop")

        nSupport <- integer(nPeaks)
        supportWeight <- numeric(nPeaks)
        nSupport[support$queryHits] <- support$nSupport
        supportWeight[support$queryHits] <- support$supportWeight

        ## presence mode never combines anything, it only counts weight
        combined <- if (presenceOnly) {
            rep(NA_real_, nPeaks)
        } else {
            .combineOverPeaks(retained = retained,
                              selected = selected,
                              combinationMethod = combinationMethod,
                              nPeaks = nPeaks)
        }

        ## the replicate count is the reproducibility criterion and a
        ## large weight never stands in for a missing replicate
        enoughReplicates <- nSupport >= requiredSupport
        enoughWeight <- if (is.null(minSupportWeight)) {
            rep(TRUE, nPeaks)
        } else {
            supportWeight >= minSupportWeight
        }

        passed <- if (presenceOnly) {
            enoughReplicates & enoughWeight
        } else {
            enoughReplicates & enoughWeight & combined >= combinedCut
        }

        result <- list(nSupport = nSupport,
                       supportWeight = supportWeight,
                       combinedNegLog10P = combined,
                       passed = passed)

        if (!recursive) {
            break
        }
        if (!is.null(previousPassed) && identical(previousPassed, passed)) {
            .messageIf(verbose, "Confirmation stable after ", iteration,
                       " rounds")
            break
        }
        if (iteration == maxIterations) {
            warning("confirmation did not stabilise within ", maxIterations,
                    " rounds; consider recursive = FALSE")
            break
        }

        previousPassed <- passed
        ## Eligibility to support others turns on reproducibility, not on
        ## the peak's own significance. Overlap is symmetric, so this
        ## shrinks monotonically and settles; keying it on `passed`
        ## instead lets unequal weights confirm a peak while refusing its
        ## own supporter, and the whole analysis then unravels to nothing.
        eligible <- enoughReplicates
    }

    result
}


#' Combine each peak with the supporters selected for it
#'
#' @param retained Flattened `GRanges`.
#' @param selected Table of retained supporters.
#' @param combinationMethod Combination scheme.
#' @param nPeaks Number of peaks.
#'
#' @return Numeric vector of combined significance, -log10 scale.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr tibble bind_rows group_by summarise mutate select
#' @importFrom S4Vectors mcols
#' @importFrom rlang .data
#'
#' @keywords internal
#' @noRd
.combineOverPeaks <- function(retained, selected, combinationMethod,
                              nPeaks) {
    metadataColumns <- S4Vectors::mcols(retained)

    ## each peak contributes its own evidence alongside its supporters
    members <- dplyr::bind_rows(
        dplyr::tibble(queryHits = seq_len(nPeaks),
                      negLog10P = metadataColumns$negLog10P,
                      rho = metadataColumns$rho,
                      weight = metadataColumns$weight),
        dplyr::tibble(queryHits = selected$queryHits,
                      negLog10P = selected$subjectNegLog10P,
                      rho = selected$subjectRho,
                      weight = selected$subjectWeight))

    ## all four schemes are sums over members, so the whole table can be
    ## reduced in one grouped pass rather than peak by peak
    terms <- .transformMembers(negLog10P = members$negLog10P,
                               weight = members$weight,
                               method = combinationMethod,
                               rho = members$rho)
    members <- dplyr::mutate(members,
                             term = terms$term,
                             scale = terms$scale)

    aggregated <- dplyr::summarise(
        dplyr::group_by(members, .data$queryHits),
        sumTerm = sum(.data$term),
        sumScale = sum(.data$scale),
        nMembers = dplyr::n(),
        .groups = "drop")

    combined <- numeric(nPeaks)
    combined[aggregated$queryHits] <- .closeCombination(
        sumTerm = aggregated$sumTerm,
        sumScale = aggregated$sumScale,
        nMembers = aggregated$nMembers,
        method = combinationMethod)

    combined
}


#' Benjamini-Hochberg correction inside each replicate
#'
#' @param combinedNegLog10P Numeric vector of combined significance.
#' @param passed Logical vector marking the peaks under test.
#' @param replicate Character vector of replicate labels.
#' @param presenceOnly Skip the correction.
#' @param adjustmentFamily `"tested"` or `"confirmed"`.
#'
#' @return Numeric vector of adjusted significance, -log10 scale.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom S4Vectors split
#'
#' @keywords internal
#' @noRd
.adjustWithinReplicates <- function(combinedNegLog10P, passed, replicate,
                                    presenceOnly,
                                    adjustmentFamily = "tested") {
    adjusted <- rep(NA_real_, length(combinedNegLog10P))
    if (presenceOnly) {
        return(adjusted)
    }

    ## Correcting only across the peaks that already cleared the combined
    ## threshold means the family was chosen for being significant, and
    ## the usual Benjamini-Hochberg guarantee does not survive that. The
    ## default therefore corrects across every peak that was tested. The
    ## alternative reproduces what MSPC does and is kept for comparison.
    for (thisReplicate in unique(replicate)) {
        index <- if (adjustmentFamily == "confirmed") {
            which(replicate == thisReplicate & passed)
        } else {
            which(replicate == thisReplicate & is.finite(combinedNegLog10P))
        }
        if (length(index) == 0) {
            next
        }
        adjusted[index] <- .bhAdjustLog(combinedNegLog10P[index])
    }
    adjusted
}


#' Benjamini-Hochberg on the -log10 scale
#'
#' @description
#' Peak callers report p-values well below the smallest representable
#' double, so the usual correction has to be done without ever forming
#' the p-value itself.
#'
#' @param negLog10P Numeric vector of -log10 p-values.
#'
#' @return Numeric vector of adjusted values on the same scale.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal
#' @noRd
.bhAdjustLog <- function(negLog10P) {
    n <- length(negLog10P)
    if (n == 0) {
        return(numeric(0))
    }

    ## mirrors stats::p.adjust(method = "BH"): walk from the least to the
    ## most significant peak, keeping a running best
    ordering <- order(negLog10P, decreasing = FALSE)
    rankIndex <- n:1L
    stepped <- cummax(negLog10P[ordering] + log10(rankIndex) - log10(n))

    adjusted <- numeric(n)
    adjusted[ordering] <- pmax(stepped, 0)
    adjusted
}


#' Merge confirmed peaks into consensus regions
#'
#' @param confirmed `GRanges` of confirmed peaks.
#' @param combinationMethod Combination scheme.
#' @param weights Named numeric vector.
#' @param mergeMethod `"reduce"` or `"iterative"`.
#' @param maxConsensusWidth Optional width cap.
#' @param multipleIntersections Rule for choosing among several peaks
#'   contributed by one replicate.
#' @param presenceOnly Skip the re-combination.
#'
#' @return A `GRanges` of consensus regions.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges reduce findOverlaps width GRanges mcols
#'   mcols<-
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom dplyr tibble group_by summarise n_distinct slice_max
#'   slice_min ungroup
#' @importFrom rlang .data
#'
#' @keywords internal
#' @noRd
.assembleConsensus <- function(confirmed,
                               combinationMethod,
                               weights,
                               mergeMethod,
                               maxConsensusWidth,
                               multipleIntersections = "lowest",
                               presenceOnly) {
    if (length(confirmed) == 0) {
        warning("no peak survived the analysis, the consensus is empty")
        return(GenomicRanges::GRanges())
    }

    merged <- if (mergeMethod == "reduce") {
        GenomicRanges::reduce(confirmed, ignore.strand = TRUE)
    } else {
        .iterativeMerge(confirmed)
    }

    ## transitive merging occasionally chains a whole neighbourhood into
    ## one region; rebuild only those with the seeded rule
    if (!is.null(maxConsensusWidth) && mergeMethod == "reduce") {
        tooWide <- GenomicRanges::width(merged) > maxConsensusWidth
        if (any(tooWide)) {
            warning(sum(tooWide), " merged regions exceeded ",
                    maxConsensusWidth, " bp and were rebuilt iteratively")
            rebuilt <- .iterativeMerge(
                IRanges::subsetByOverlaps(confirmed, merged[tooWide],
                                          ignore.strand = TRUE))
            merged <- c(merged[!tooWide], rebuilt)
            merged <- sort(merged)
        }
    }

    ## describe each region by the peaks that fall inside it
    hits <- GenomicRanges::findOverlaps(confirmed, merged,
                                        ignore.strand = TRUE)
    peakMetadata <- S4Vectors::mcols(confirmed)

    membership <- dplyr::tibble(
        region = S4Vectors::subjectHits(hits),
        replicate = peakMetadata$replicate[S4Vectors::queryHits(hits)],
        negLog10P = peakMetadata$negLog10P[S4Vectors::queryHits(hits)],
        rho = peakMetadata$rho[S4Vectors::queryHits(hits)],
        weight = peakMetadata$weight[S4Vectors::queryHits(hits)])

    summarised <- dplyr::summarise(
        dplyr::group_by(membership, .data$region),
        nReplicates = dplyr::n_distinct(.data$replicate),
        nPeaks = dplyr::n(),
        replicates = paste(sort(unique(.data$replicate)), collapse = ","),
        .groups = "drop")

    S4Vectors::mcols(merged)$nReplicates <- 0L
    S4Vectors::mcols(merged)$nPeaks <- 0L
    S4Vectors::mcols(merged)$replicates <- NA_character_
    S4Vectors::mcols(merged)$nReplicates[summarised$region] <-
        summarised$nReplicates
    S4Vectors::mcols(merged)$nPeaks[summarised$region] <- summarised$nPeaks
    S4Vectors::mcols(merged)$replicates[summarised$region] <-
        summarised$replicates

    ## Recombine the member peaks so the region carries a significance of
    ## its own. A replicate contributes once, exactly as it did when the
    ## peaks were confirmed: without this a replicate that fragmented a
    ## domain into five calls would count five times, and the region
    ## would look more significant the worse the peak calling was.
    if (!presenceOnly) {
        grouped <- dplyr::group_by(membership, .data$region,
                                   .data$replicate)
        representative <- if (multipleIntersections == "lowest") {
            dplyr::slice_max(grouped, order_by = .data$negLog10P,
                             n = 1, with_ties = FALSE)
        } else {
            dplyr::slice_min(grouped, order_by = .data$negLog10P,
                             n = 1, with_ties = FALSE)
        }
        representative <- dplyr::ungroup(representative)

        terms <- .transformMembers(negLog10P = representative$negLog10P,
                                   weight = representative$weight,
                                   method = combinationMethod,
                                   rho = representative$rho)
        representative$term <- terms$term
        representative$scale <- terms$scale

        regionStat <- dplyr::summarise(
            dplyr::group_by(representative, .data$region),
            sumTerm = sum(.data$term),
            sumScale = sum(.data$scale),
            nMembers = dplyr::n(),
            .groups = "drop")

        S4Vectors::mcols(merged)$combinedNegLog10P <- NA_real_
        S4Vectors::mcols(merged)$combinedNegLog10P[regionStat$region] <-
            .closeCombination(sumTerm = regionStat$sumTerm,
                              sumScale = regionStat$sumScale,
                              nMembers = regionStat$nMembers,
                              method = combinationMethod)
    }

    merged
}


#' Seeded merging that avoids transitive chaining
#'
#' @description
#' Takes the most significant peak, claims everything it overlaps as one
#' region, removes them, and repeats. Regions therefore stay close to the
#' size of a real peak instead of growing along a chain of partial
#' overlaps.
#'
#' @param gr `GRanges` of confirmed peaks.
#'
#' @return A `GRanges` of merged regions.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges findOverlaps reduce
#' @importFrom S4Vectors queryHits subjectHits mcols split
#'
#' @keywords internal
#' @noRd
.iterativeMerge <- function(gr) {
    if (length(gr) == 0) {
        return(GenomicRanges::GRanges())
    }

    score <- S4Vectors::mcols(gr)$combinedNegLog10P
    if (is.null(score) || all(is.na(score))) {
        score <- S4Vectors::mcols(gr)$negLog10P
    }
    if (is.null(score) || all(is.na(score))) {
        ## nothing to rank on, so width decides which peak seeds a region
        score <- GenomicRanges::width(gr)
    }

    hits <- GenomicRanges::findOverlaps(gr, gr, ignore.strand = TRUE)
    neighbours <- S4Vectors::split(S4Vectors::subjectHits(hits),
                                   S4Vectors::queryHits(hits))

    available <- rep(TRUE, length(gr))
    clusters <- vector("list", length(gr))
    nClusters <- 0L

    ## walk the peaks best first; each one absorbs its overlappers
    for (index in order(score, decreasing = TRUE)) {
        if (!available[index]) {
            next
        }
        members <- as.integer(neighbours[[index]])
        members <- members[available[members]]
        available[members] <- FALSE
        nClusters <- nClusters + 1L
        clusters[[nClusters]] <- members
    }

    clusters <- clusters[seq_len(nClusters)]
    merged <- lapply(clusters, function(members) {
        GenomicRanges::reduce(gr[members], ignore.strand = TRUE)
    })

    sort(unlist(GenomicRanges::GRangesList(merged)))
}


#' Per-replicate summary table
#'
#' @param retained Flattened `GRanges` after confirmation.
#' @param replicateNames Character vector of replicate names.
#'
#' @return A `data.frame`.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr tibble group_by summarise mutate arrange
#' @importFrom S4Vectors mcols
#' @importFrom rlang .data
#'
#' @keywords internal
#' @noRd
.replicateSummary <- function(retained, replicateNames) {
    peakTable <- .mcolsAsTibble(retained)

    summarised <- dplyr::summarise(
        dplyr::group_by(peakTable, .data$replicate),
        nTested = dplyr::n(),
        nStringent = sum(.data$class == "stringent"),
        nWeak = sum(.data$class == "weak"),
        nConfirmed = sum(.data$status == "confirmed"),
        nRescued = sum(.data$status == "confirmed" & .data$class == "weak"),
        nFalsePositive = sum(.data$status == "falsePositive"),
        nDiscarded = sum(.data$status == "discarded"),
        .groups = "drop")

    ## the rescue rate is the number to look at first: a value far above
    ## what the replicates share otherwise usually means the weak
    ## threshold was set too permissively
    summarised <- dplyr::mutate(
        summarised,
        rescueRate = .data$nRescued / pmax(.data$nConfirmed, 1))

    summarised <- summarised[match(replicateNames, summarised$replicate), ]
    as.data.frame(summarised)
}
