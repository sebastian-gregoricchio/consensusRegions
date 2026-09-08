#' consensusRegions: consensus peaks across replicated experiments
#'
#' @description
#' Replicated ChIP-seq and ATAC-seq experiments rarely agree peak for peak.
#' A region can be convincing in two replicates and just below threshold in
#' the third, and discarding it wastes real signal. `consensusRegions`
#' looks at every peak together with the peaks that overlap it in the other
#' replicates and combines their p-values, so that repeated weak evidence
#' can compensate for the lack of a single strong call.
#'
#' The approach follows the one introduced by MSPC (Jalili and colleagues,
#' 2015) and adds three things: replicate weights, so a shallow or noisy
#' sample counts for less; a permutation null, so the threshold on the
#' combined significance does not have to be guessed; and tolerance for
#' incomplete input, since not everyone keeps the BAM files or the
#' p-value column.
#'
#' @section Entry points:
#' \describe{
#'   \item{[readPeakSets()]}{reads narrowPeak, broadPeak, BED6 or BED3
#'     files and decides which statistics are usable.}
#'   \item{[computeReplicateWeights()]}{derives weights from BAM files,
#'     from pre-computed FRiP values, or from the peak sets alone.}
#'   \item{[buildConsensus()]}{runs the analysis and returns a
#'     [ConsensusRegions-class] object.}
#'   \item{[calibrateThreshold()]}{estimates the combined threshold from
#'     shuffled peak sets.}
#'   \item{[recentrePeaks()]}{recovers summits and returns fixed-width
#'     regions.}
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @references
#' Jalili V., Matteucci M., Masseroli M., Morelli M.J. (2015). Using
#' combined evidence from replicates to evaluate ChIP-seq peaks.
#' *Bioinformatics* 31(17):2761-2769.
#'
#' Corces M.R. *et al.* (2018). The chromatin accessibility landscape of
#' primary human cancers. *Science* 362(6413):eaav1898.
#'
#' @keywords internal
"_PACKAGE"
