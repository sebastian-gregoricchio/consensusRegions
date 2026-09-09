# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' Read peak files into a GRangesList
#'
#' @description
#' Imports peak calls and works out which statistic can actually be used.
#' A narrowPeak carries -log10 p-values in its eighth column and gets the
#' full treatment. A BED with a score column is rank-transformed. A bare
#' BED3 has nothing to combine, so the analysis later falls back to a
#' presence rule. The detection is reported so the choice is never silent.
#'
#' @param peaks Character vector of file paths, or a `GRanges`,
#'   `GRangesList`, or list of `GRanges`.
#' @param sampleNames Character vector naming the replicates. Taken from
#'   the file names or the list names when left `NULL`. Default: \code{NULL}.
#' @param scoreType One of `"auto"`, `"log10pvalue"`, `"pvalue"`,
#'   `"score"` or `"none"`. `"auto"` inspects the input.
#'   Default: \code{"auto"}.
#' @param scoreColumn Name of the metadata column holding the statistic.
#'   Only needed when the input is not a standard peak format.
#'   Default: \code{NULL}.
#' @param keepStandardChromosomes Drop scaffolds and patches.
#'   Default: \code{TRUE}.
#' @param seqlevelsStyle String indicating the chromosome naming style to
#'   apply to every replicate, one among `"UCSC"` (chr1), `"Ensembl"` (1)
#'   or `"NCBI"`. When set to `NULL` the names are left as they are and
#'   the loading is interrupted if the replicates use different styles.
#'   Default: `"UCSC"`.
#' @param genome Genome identifier passed to the importer, for instance
#'   `"hg38"`. Fills in the chromosome lengths used by
#'   [calibrateThreshold()].
#' @param verbose Report what was detected. Default: \code{TRUE}.
#'
#' @return A `GRangesList`, one element per replicate, with a `negLog10P`
#'   metadata column when a statistic was available. The resolved score
#'   type is stored in the object metadata.
#'
#' @details
#' A narrowPeak whose p-value column is filled with `-1`, which is what
#' several callers write when they have nothing to report, is treated as
#' score-only rather than as a set of p-values equal to `10`.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRangesList GRanges mcols
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom S4Vectors metadata mcols mcols<- metadata<-
#' @importFrom methods is
#' @importFrom dplyr case_when
#'
#' @examples
#' peakFiles <- system.file("extdata",
#'                          c("rep1.narrowPeak", "rep2.narrowPeak",
#'                            "rep3.narrowPeak"),
#'                          package = "consensusRegions")
#' peaks <- readPeakSets(peakFiles, sampleNames = c("r1", "r2", "r3"))
#' peaks
#'
#' @export
readPeakSets <- function(peaks,
                         sampleNames = NULL,
                         scoreType = c("auto", "log10pvalue", "pvalue",
                                       "score", "none"),
                         scoreColumn = NULL,
                         keepStandardChromosomes = TRUE,
                         seqlevelsStyle = "UCSC",
                         genome = NA,
                         verbose = TRUE) {
    scoreType <- match.arg(scoreType)

    if (!is.null(seqlevelsStyle)) {
        seqlevelsStyle <- unname(
            c(ucsc = "UCSC", ensembl = "Ensembl",
              ncbi = "NCBI")[tolower(seqlevelsStyle[1])])
        if (is.na(seqlevelsStyle)) {
            stop("The 'seqlevelsStyle' parameter must be one among ",
                 "'UCSC', 'Ensembl' or 'NCBI'.", call. = FALSE)
        }
    }

    ## bring every accepted input shape down to a plain list of GRanges
    peakList <- .coercePeakInput(peaks, genome = genome)

    if (length(peakList) < 2) {
        stop("at least two replicates are needed, ", length(peakList),
             " were supplied")
    }

    ## names drive the weight lookup and every table produced later
    if (is.null(sampleNames)) {
        sampleNames <- names(peakList)
    }
    if (is.null(sampleNames) || any(is.na(sampleNames))) {
        stop("replicate names could not be determined, supply 'sampleNames'")
    }
    if (length(sampleNames) != length(peakList)) {
        stop("'sampleNames' has length ", length(sampleNames),
             " but ", length(peakList), " peak sets were supplied")
    }
    if (anyDuplicated(sampleNames) > 0) {
        stop("'sampleNames' must be unique")
    }
    names(peakList) <- sampleNames

    ## Peak files rarely agree on whether a chromosome is called chr1 or
    ## 1, and two replicates naming the same chromosome differently share
    ## no overlap at all: every peak would be discarded for lack of
    ## support and the analysis would look like a biological result.
    peakList <- .applySeqlevelsStyle(peakList,
                                     seqlevelsStyle = seqlevelsStyle,
                                     verbose = verbose)

    if (isTRUE(keepStandardChromosomes)) {
        peakList <- lapply(peakList, function(x) {
            GenomeInfoDb::keepStandardChromosomes(x, pruning.mode = "coarse")
        })
    }

    ## the score type is settled once for the whole experiment, since
    ## combining a p-value from one replicate with a rank from another
    ## would not mean anything
    resolvedType <- if (scoreType == "auto") {
        .detectScoreType(peakList, scoreColumn = scoreColumn)
    } else {
        scoreType
    }

    .messageIf(verbose, "Score type in use: ", resolvedType)

    ## attach the statistic every replicate will be judged on
    peakList <- lapply(peakList, .attachNegLog10P,
                       scoreType = resolvedType,
                       scoreColumn = scoreColumn)

    emptySets <- names(peakList)[vapply(peakList, length, integer(1)) == 0]
    if (length(emptySets) > 0) {
        warning("no peaks left in: ", paste(emptySets, collapse = ", "))
    }

    result <- GenomicRanges::GRangesList(peakList)
    S4Vectors::metadata(result) <- list(scoreType = resolvedType)
    result
}


#' Put every replicate on the same chromosome naming style
#'
#' @description
#' `GenomeInfoDb` knows the mapping for the assemblies it ships, and falls
#' over on anything it does not recognise. The prefix rule below is the
#' fallback for those cases, which is the same arrangement used by
#' `RegionSetDE::loadRegions()`.
#'
#' @param peakList Named list of `GRanges`.
#' @param seqlevelsStyle `"UCSC"`, `"Ensembl"`, `"NCBI"` or `NULL`.
#' @param verbose Report what was changed.
#'
#' @return The list with harmonised sequence levels.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqlevelsStyle
#'   seqlevelsStyle<-
#'
#' @keywords internal
#' @noRd
.applySeqlevelsStyle <- function(peakList, seqlevelsStyle, verbose = TRUE) {
    currentStyle <- function(gr) {
        style <- tryCatch(GenomeInfoDb::seqlevelsStyle(gr),
                          error = function(e) NA_character_)
        style[1]
    }

    ## nothing to impose, but replicates that disagree cannot be compared
    if (is.null(seqlevelsStyle)) {
        styles <- vapply(peakList, currentStyle, character(1))
        prefixed <- vapply(peakList, function(gr) {
            any(grepl("^chr", GenomeInfoDb::seqlevels(gr)))
        }, logical(1))

        if (length(unique(prefixed)) > 1) {
            stop("the replicates use different chromosome naming styles (",
                 paste(names(peakList)[prefixed], collapse = ", "),
                 " carry a 'chr' prefix and the others do not); set ",
                 "'seqlevelsStyle' to harmonise them", call. = FALSE)
        }
        if (length(unique(stats::na.omit(styles))) > 1) {
            stop("the replicates use different chromosome naming styles: ",
                 paste(unique(stats::na.omit(styles)), collapse = ", "),
                 "; set 'seqlevelsStyle' to harmonise them", call. = FALSE)
        }
        return(peakList)
    }

    before <- vapply(peakList, function(gr) {
        paste(utils::head(GenomeInfoDb::seqlevels(gr), 1), collapse = "")
    }, character(1))

    harmonised <- lapply(peakList, function(gr) {
        tryCatch({
            GenomeInfoDb::seqlevelsStyle(gr) <- seqlevelsStyle
            gr
        }, error = function(e) .harmonizeSeqlevels(gr, seqlevelsStyle))
    })

    after <- vapply(harmonised, function(gr) {
        paste(utils::head(GenomeInfoDb::seqlevels(gr), 1), collapse = "")
    }, character(1))

    changed <- names(peakList)[before != after]
    if (length(changed) > 0) {
        .messageIf(verbose, "Chromosome names set to ", seqlevelsStyle,
                   " style in: ", paste(changed, collapse = ", "))
    }

    harmonised
}


#' Add or strip the chr prefix when the assembly is not recognised
#'
#' @param gr A `GRanges`.
#' @param style `"UCSC"`, `"Ensembl"` or `"NCBI"`.
#'
#' @return The `GRanges` with renamed sequence levels.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#'
#' @keywords internal
#' @noRd
.harmonizeSeqlevels <- function(gr, style) {
    newLevels <- GenomeInfoDb::seqlevels(gr)

    if (style == "UCSC") {
        newLevels <- ifelse(grepl("^chr", newLevels), newLevels,
                            paste0("chr", newLevels))
        newLevels <- gsub("^chrMT$", "chrM", newLevels)
    } else {
        newLevels <- gsub("^chr", "", newLevels)
        newLevels <- gsub("^M$", "MT", newLevels)
    }

    if (any(duplicated(newLevels))) {
        stop("The chromosome names could not be converted to the ",
             "requested style without collisions.", call. = FALSE)
    }

    GenomeInfoDb::seqlevels(gr) <- newLevels
    gr
}


#' Accept files, a GRanges, a GRangesList or a list of GRanges
#'
#' @param peaks The user supplied input.
#' @param genome Genome identifier handed to the importer.
#'
#' @return A named list of `GRanges`.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom rtracklayer import
#' @importFrom methods is
#' @importFrom tools file_path_sans_ext
#'
#' @keywords internal
#' @noRd
.coercePeakInput <- function(peaks, genome = NA) {
    if (is.character(peaks)) {
        missingFiles <- peaks[!file.exists(peaks)]
        if (length(missingFiles) > 0) {
            stop("file not found: ", paste(missingFiles, collapse = ", "))
        }
        imported <- lapply(peaks, .importOnePeakFile, genome = genome)
        names(imported) <- tools::file_path_sans_ext(basename(peaks))
        return(imported)
    }

    if (methods::is(peaks, "GRangesList")) {
        return(as.list(peaks))
    }

    if (methods::is(peaks, "GRanges")) {
        stop("a single GRanges was supplied; consensus needs several ",
             "replicates, pass a list or a GRangesList")
    }

    if (is.list(peaks) &&
        all(vapply(peaks, methods::is, logical(1), "GRanges"))) {
        return(peaks)
    }

    stop("'peaks' must be file paths, a GRangesList, or a list of GRanges")
}


#' Import one peak file, choosing the format from its extension
#'
#' @param path Path to the file.
#' @param genome Genome identifier handed to the importer.
#'
#' @return A `GRanges`.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom rtracklayer import
#' @importFrom tools file_ext
#'
#' @keywords internal
#' @noRd
.importOnePeakFile <- function(path, genome = NA) {
    extension <- tolower(tools::file_ext(path))

    peakFormat <- dplyr::case_when(
        extension %in% c("narrowpeak", "np") ~ "narrowPeak",
        extension %in% c("broadpeak", "bp") ~ "broadPeak",
        extension == "gappedpeak" ~ "bed",
        .default = "bed"
    )

    ## a mislabelled extension is common enough that a hard failure here
    ## is unhelpful; retry as plain BED before giving up
    imported <- tryCatch(
        rtracklayer::import(path, format = peakFormat, genome = genome),
        error = function(e) {
            tryCatch(
                rtracklayer::import(path, format = "bed", genome = genome),
                error = function(e2) {
                    stop("could not read ", path, ": ", conditionMessage(e2))
                }
            )
        }
    )

    imported
}


#' Work out which statistic the peak sets can support
#'
#' @param peakList List of `GRanges`.
#' @param scoreColumn Optional metadata column name.
#'
#' @return One of `"log10pvalue"`, `"score"` or `"none"`.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom S4Vectors mcols
#'
#' @keywords internal
#' @noRd
.detectScoreType <- function(peakList, scoreColumn = NULL) {
    columnNames <- lapply(peakList, function(x) names(S4Vectors::mcols(x)))
    shared <- Reduce(intersect, columnNames)

    if (!is.null(scoreColumn)) {
        if (!scoreColumn %in% shared) {
            stop("column '", scoreColumn,
                 "' is not present in every peak set")
        }
        return("score")
    }

    ## narrowPeak and broadPeak keep -log10 p in a column named pValue;
    ## callers that have none write -1 there, which must not be mistaken
    ## for a p-value of 0.1
    if ("pValue" %in% shared) {
        usable <- vapply(peakList, function(x) {
            pv <- S4Vectors::mcols(x)$pValue
            length(pv) > 0 && any(pv > 0, na.rm = TRUE)
        }, logical(1))
        if (all(usable)) {
            return("log10pvalue")
        }
    }

    ## a score column is enough to rank the peaks, provided it varies
    if ("score" %in% shared) {
        informative <- vapply(peakList, function(x) {
            sc <- S4Vectors::mcols(x)$score
            length(sc) > 0 &&
                length(unique(sc[!is.na(sc)])) > 1
        }, logical(1))
        if (all(informative)) {
            return("score")
        }
    }

    "none"
}


#' Attach the negLog10P column used by every downstream step
#'
#' @param gr A `GRanges`.
#' @param scoreType Resolved score type.
#' @param scoreColumn Optional metadata column name.
#'
#' @return The `GRanges` with a `negLog10P` metadata column.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom S4Vectors mcols mcols<-
#'
#' @keywords internal
#' @noRd
.attachNegLog10P <- function(gr, scoreType, scoreColumn = NULL) {
    if (length(gr) == 0) {
        S4Vectors::mcols(gr)$negLog10P <- numeric(0)
        return(gr)
    }

    metadataColumns <- S4Vectors::mcols(gr)

    negLog10P <- if (scoreType == "log10pvalue") {
        ## already on the right scale, nothing to do
        as.numeric(metadataColumns$pValue)
    } else if (scoreType == "pvalue") {
        column <- if (is.null(scoreColumn)) "pValue" else scoreColumn
        rawP <- as.numeric(metadataColumns[[column]])
        if (any(rawP < 0 | rawP > 1, na.rm = TRUE)) {
            stop("values outside [0, 1] found where p-values were expected")
        }
        -log10(pmax(rawP, .Machine$double.xmin))
    } else if (scoreType == "score") {
        column <- if (is.null(scoreColumn)) "score" else scoreColumn
        .negLog10FromScore(as.numeric(metadataColumns[[column]]))
    } else {
        ## nothing usable; the presence rule takes over later
        rep(NA_real_, length(gr))
    }

    S4Vectors::mcols(gr)$negLog10P <- negLog10P
    gr
}
