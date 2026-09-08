# consensusRegions 0.99.0

First release, submitted to Bioconductor.

* `buildConsensus()` combines the evidence of overlapping peaks across
  replicates, following the approach of MSPC: two p-value thresholds, a
  minimum number of supporting replicates, a test on the combined
  statistic, and a Benjamini-Hochberg correction inside each replicate.
* Four combination schemes. Fisher reproduces the original behaviour;
  Stouffer and Lancaster accept per-replicate weights; the rank product
  works on ranks alone and needs no p-values.
* `computeReplicateWeights()` derives weights from BAM files, from FRiP
  values measured elsewhere, or from the peak sets themselves when
  neither is available.
* `calibrateThreshold()` estimates the threshold on the combined
  statistic from shuffled peak sets rather than leaving it at a fixed
  1e-8.
* Input degrades gracefully. narrowPeak p-values are used directly, a
  plain score column is rank-transformed, and BED3 input falls back to a
  weighted k-of-n presence rule.
* `recentrePeaks()` recovers summits and returns fixed-width regions,
  which is what motif analysis and count matrices want.
* `plotRescue()`, `plotJaccard()`, `plotUpset()` and `plotCalibration()`
  for quality control; `asSummarizedExperiment()` to hand the consensus
  to csaw, DiffBind or edgeR.
