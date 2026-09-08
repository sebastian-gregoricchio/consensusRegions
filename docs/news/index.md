# Changelog

## consensusRegions 0.99.0

First release, submitted to Bioconductor.

- [`buildConsensus()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/buildConsensus.md)
  combines the evidence of overlapping peaks across replicates,
  following the approach of MSPC: two p-value thresholds, a minimum
  number of supporting replicates, a test on the combined statistic, and
  a Benjamini-Hochberg correction inside each replicate.
- Four combination schemes. Fisher reproduces the original behaviour;
  Stouffer and Lancaster accept per-replicate weights; the rank product
  works on ranks alone and needs no p-values.
- [`computeReplicateWeights()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/computeReplicateWeights.md)
  derives weights from BAM files, from FRiP values measured elsewhere,
  or from the peak sets themselves when neither is available.
- [`calibrateThreshold()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/calibrateThreshold.md)
  estimates the threshold on the combined statistic from shuffled peak
  sets rather than leaving it at a fixed 1e-8.
- Input degrades gracefully. narrowPeak p-values are used directly, a
  plain score column is rank-transformed, and BED3 input falls back to a
  weighted k-of-n presence rule.
- [`recentrePeaks()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/recentrePeaks.md)
  recovers summits and returns fixed-width regions, which is what motif
  analysis and count matrices want.
- [`plotRescue()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/plotRescue.md),
  [`plotJaccard()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/plotJaccard.md),
  [`plotUpset()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/plotUpset.md)
  and
  [`plotCalibration()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/plotCalibration.md)
  for quality control;
  [`asSummarizedExperiment()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/asSummarizedExperiment.md)
  to hand the consensus to csaw, DiffBind or edgeR.
