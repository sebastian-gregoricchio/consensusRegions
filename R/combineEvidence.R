#' Combine the evidence carried by overlapping peaks
#'
#' @description
#' Pools a set of p-values into a single one. Four schemes are offered and
#' they behave differently enough that the choice matters.
#'
#' Fisher's method is what MSPC uses. It treats every replicate as
#' equally informative and is the right default when they really are
#' comparable. Its weakness shows up as the number of replicates grows,
#' because the combined statistic keeps accumulating and a region detected
#' faintly in many samples starts to look highly significant.
#'
#' Stouffer's Z accepts weights, so a shallow replicate can be told to
#' count for less instead of being dropped or trusted blindly. This is the
#' default here.
#'
#' Lancaster's method also takes weights, folding them into the degrees of
#' freedom of a chi-squared sum. It is closer to Fisher in spirit and
#' tends to react more strongly to a single very small p-value.
#'
#' The rank product ignores the p-values and works on within-replicate
#' ranks, which makes it the natural option when the input carries only a
#' score. It is also the least sensitive to a replicate whose p-values are
#' poorly calibrated.
#'
#' @param negLog10P Numeric vector of -log10 transformed p-values.
#' @param weights Numeric vector of replicate weights, recycled to the
#'   length of `negLog10P`. Ignored by `"fisher"` and `"rankProduct"`.
#'   Defaults to equal weights.
#' @param method One of `"stouffer"`, `"fisher"`, `"lancaster"` or
#'   `"rankProduct"`.
#' @param rho Numeric vector of relative within-replicate ranks in
#'   `(0, 1]`, required by `"rankProduct"` and ignored otherwise.
#'
#' @return A single numeric value, the combined significance on the
#'   -log10 scale.
#'
#' @details
#' Everything is computed in log space. Peak callers routinely report
#' p-values far below the smallest representable double, and exponentiating
#' them first would collapse the whole set to zero.
#'
#' @author Sebastian Gregoricchio
#'
#' @references
#' Lancaster H.O. (1961). The combination of probabilities: an application
#' of orthonormal functions. *Australian Journal of Statistics* 3:20-33.
#'
#' Breitling R., Armengaud P., Amtmann A., Herzyk P. (2004). Rank products:
#' a simple, yet powerful, new method to detect differentially regulated
#' genes. *FEBS Letters* 573(1-3):83-92.
#'
#' @examples
#' ## three replicates, none convincing alone
#' combineEvidence(c(3, 3.2, 2.8), method = "fisher")
#'
#' ## the same peaks, but the third replicate is known to be poor
#' combineEvidence(c(3, 3.2, 2.8), weights = c(1.3, 1.3, 0.4),
#'                 method = "stouffer")
#'
#' @export
combineEvidence <- function(negLog10P,
                            weights = NULL,
                            method = c("stouffer", "fisher", "lancaster",
                                       "rankProduct"),
                            rho = NULL) {
    method <- match.arg(method)

    if (length(negLog10P) == 0) {
        stop("'negLog10P' is empty, there is nothing to combine")
    }
    if (is.null(weights)) {
        weights <- rep(1, length(negLog10P))
    }
    weights <- rep_len(weights, length(negLog10P))

    terms <- .transformMembers(negLog10P = negLog10P,
                               weight = weights,
                               method = method,
                               rho = rho)

    .closeCombination(sumTerm = sum(terms$term),
                      sumScale = sum(terms$scale),
                      nMembers = length(negLog10P),
                      method = method)
}


#' Per-member contribution to a combined statistic
#'
#' @description
#' All four schemes reduce to a sum over members plus a count, which is
#' what lets the whole peak table be combined with a single grouped
#' aggregation instead of one call per peak.
#'
#' @param negLog10P Numeric vector of -log10 p-values.
#' @param weight Numeric vector of weights.
#' @param method Combination scheme.
#' @param rho Numeric vector of relative ranks, used by `"rankProduct"`.
#'
#' @return A `tibble` with columns `term` and `scale`.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom stats qnorm qchisq
#' @importFrom tibble tibble
#'
#' @keywords internal
#' @noRd
.transformMembers <- function(negLog10P, weight, method, rho = NULL) {
    ## a p-value of exactly one sends the normal quantile to -Inf, so nudge
    ## it just below; the shift is far smaller than any real signal
    negLog10P <- pmax(negLog10P, 1e-10)
    logP <- -negLog10P * log(10)

    if (method == "fisher") {
        ## Fisher only needs the sum of -log10 p; the constants are applied
        ## once when the sum is closed
        return(tibble::tibble(term = negLog10P,
                              scale = rep(1, length(negLog10P))))
    }

    if (method == "stouffer") {
        zScore <- stats::qnorm(logP, lower.tail = FALSE, log.p = TRUE)
        return(tibble::tibble(term = weight * zScore,
                              scale = weight^2))
    }

    if (method == "lancaster") {
        ## the weight becomes the degrees of freedom of each contribution
        chiTerm <- stats::qchisq(logP, df = weight,
                                 lower.tail = FALSE, log.p = TRUE)
        return(tibble::tibble(term = chiTerm, scale = weight))
    }

    ## rank product
    if (is.null(rho)) {
        stop("'rho' is required when method is 'rankProduct'")
    }
    if (any(rho <= 0 | rho > 1, na.rm = TRUE)) {
        stop("'rho' must hold relative ranks inside (0, 1]")
    }
    tibble::tibble(term = -log(rho), scale = rep(1, length(rho)))
}


#' Close a combined statistic into a significance value
#'
#' @param sumTerm Sum of the per-member terms.
#' @param sumScale Sum of the per-member scales.
#' @param nMembers Number of members contributing.
#' @param method Combination scheme.
#'
#' @return Numeric vector of combined significance, -log10 scale.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom stats pchisq pnorm pgamma
#'
#' @keywords internal
#' @noRd
.closeCombination <- function(sumTerm, sumScale, nMembers, method) {
    logTen <- log(10)

    combined <- switch(
        method,
        ## X2 = -2 * sum(ln p) = 2 * ln(10) * sum(-log10 p)
        fisher = stats::pchisq(2 * logTen * sumTerm,
                               df = 2 * nMembers,
                               lower.tail = FALSE, log.p = TRUE),
        stouffer = stats::pnorm(sumTerm / sqrt(sumScale),
                                lower.tail = FALSE, log.p = TRUE),
        lancaster = stats::pchisq(sumTerm, df = sumScale,
                                  lower.tail = FALSE, log.p = TRUE),
        rankProduct = stats::pgamma(sumTerm, shape = nMembers, rate = 1,
                                    lower.tail = FALSE, log.p = TRUE)
    )

    ## report on the same -log10 scale as the input
    -combined / logTen
}
