# consensusRegions 0.99.0

First release.

`consensusRegions` combines the p-values of overlapping peaks across
replicated ChIP-seq, ATAC-seq and CUT&RUN experiments. A region that
convinces two samples and falls short in the third survives, where an
intersection would drop it. The method follows MSPC (Jalili *et al.*,
2015) and runs in R, with no .NET runtime to install.

## Analysis

* `runConsensus()` takes the whole workflow in one call. The individual
  functions stay available for the runs you want to inspect partway
  through.
* Four combination schemes. Fisher reproduces MSPC. Stouffer and
  Lancaster accept per-replicate weights, so a shallow sample counts for
  less rather than being dropped or trusted whole. The rank product works
  on ranks and needs no p-values at all.
* `computeReplicateWeights()` reads those weights from BAM files, from
  FRiP values you measured elsewhere, or from the peak sets themselves
  once the alignments are gone.
* `calibrateThreshold()` estimates the cut on the combined statistic from
  shuffled peak sets, instead of leaving it at the customary `1e-8`. Give
  `BPPARAM` a number of cores to spread the permutations.

## Input

* narrowPeak p-values go straight in. `readPeakSets()` rank-transforms a
  score column when that is all a file carries, and falls back to a
  weighted k-of-n presence rule for BED3.
* `readPeakSets()` puts every replicate on one chromosome naming style.
  Two replicates calling the same chromosome `chr1` and `1` share no
  overlap, and the analysis would report that as biology.

## Output

* `recentrePeaks()` recovers summits and returns fixed-width regions,
  which is what motif analysis and count matrices want.
* `asSummarizedExperiment()` hands the consensus to csaw, DiffBind or
  edgeR. `exportConsensus()` writes BED or narrowPeak.
* `plotRescue()`, `plotJaccard()`, `plotUpset()` and `plotCalibration()`
  cover the quality control.

## Multiple testing

`buildConsensus()` corrects within each replicate across every peak it
tested, not across the subset that cleared the combined threshold.
Correcting a family picked for being significant is not
Benjamini-Hochberg, whatever the column is called. Set
`adjustmentFamily = "confirmed"` to reproduce MSPC.

Neither setting delivers exact false discovery rate control, and the
package does not claim it. Each replicate contributes its most
significant overlapping peak, and a maximum over a selected set does not
behave like an unselected p-value. Overlapping peaks share supporters
too, so their combined values are not independent. Reach for
`calibrateThreshold()` when you need a defensible error rate rather than
a ranking.
