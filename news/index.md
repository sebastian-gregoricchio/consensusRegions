# Changelog

## consensusRegions 0.99.3

### New

- `ThresholdCalibration` class.
  [`calibrateThreshold()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/calibrateThreshold.md)
  returned a plain list, which
  [`buildConsensus()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/buildConsensus.md)
  then had to check by looking for a `threshold` element. It now returns
  an S4 object with a validity method, so a malformed calibration is
  refused where it is made rather than surfacing later as an empty
  consensus. Accessors are
  [`calibrationThreshold()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/ThresholdCalibration-accessors.md),
  [`calibrationThresholdNegLog10()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/ThresholdCalibration-accessors.md),
  [`calibrationCurve()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/ThresholdCalibration-accessors.md),
  [`calibrationTargetFDR()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/ThresholdCalibration-accessors.md),
  [`calibrationAchievedFDR()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/ThresholdCalibration-accessors.md),
  [`calibrationSeed()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/ThresholdCalibration-accessors.md),
  [`nullStatistics()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/ThresholdCalibration-accessors.md)
  and
  [`observedStatistics()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/ThresholdCalibration-accessors.md),
  and there is a `show()` method. `$` still reaches the components by
  name, so code written against the list keeps working.
- [`calibrationResult()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/ConsensusRegions-accessors.md)
  returns the calibration of a `ConsensusRegions` object, so that every
  slot now has an accessor. A consensus built without calibration
  carries `NULL` rather than an empty list.

### Changes

- The components of a `ConsensusRegions` object are read through their
  accessors throughout, so `@` appears only in the class definition and
  in the methods themselves.
- The `stats` slot is an
  [`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  rather than a `data.frame`.
- `BPPARAM` defaults to `NULL`, which builds a
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  seeded from `seed`, so nothing runs in parallel unasked. A
  `BiocParallelParam` is used untouched and a number of cores is still
  accepted as a shorthand.
- [`buildConsensus()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/buildConsensus.md)
  hands its argument checking to `.resolveBuildSettings()`, which can be
  read and tested on its own.
- [`plotCalibration()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/plotCalibration.md)
  labels its axis with a plotmath expression, and `ggtext` is no longer
  a dependency.
- `.pairwiseJaccard()` computes the upper triangle only and mirrors it,
  rather than every ordered pair. A set still scores 1 against itself,
  and 0 when it covers no bases at all.
- The package man page is no longer marked internal, so
  [`?consensusRegions`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/consensusRegions-package.md)
  is reachable from the help index, and points at the vignette.
- The vignette has an installation section, says how the approach
  relates to MSPC and to the consensus options already in Bioconductor,
  and evaluates the parallel and blacklist examples instead of showing
  them unevaluated.
- `LICENSE` is no longer shipped in the tarball, since `GPL-3` needs no
  copy of its text.
- The chromosome-name fallback credits `RegionSetDE`, by the same
  author, and links to it.

## consensusRegions 0.99.2

### Reproducibility

- [`calibrateThreshold()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/calibrateThreshold.md)
  and
  [`runConsensus()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/runConsensus.md)
  take `seed`, which seeds the random positions of the shuffled peaks.
  It goes to the `BiocParallel` backend built from `BPPARAM` as its
  `RNGseed`, so the same seed returns the same threshold on one core or
  on many, and the random numbers of the session are left as they were.
  The calibration keeps the `seed` and its `seedSource`, and the object
  prints it.
- The documentation said that
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  brought the same threshold back on a single core. It does not:
  `BiocParallel` draws from random numbers of its own, which
  [`set.seed()`](https://rdrr.io/r/base/Random.html) does not reach, so
  two identical calls could calibrate two different thresholds. The
  examples, the vignette and the tests use `seed` now. A
  `BiocParallelParam` passed as `BPPARAM` is still used untouched, and
  carries its own `RNGseed`; `seed` is ignored with it, with a warning.

## consensusRegions 0.99.1

### Performance

- The confirmation step of
  [`buildConsensus()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/buildConsensus.md)
  is about four times faster on genome-scale peak sets. It picked the
  supporting peak of every replicate with a grouped `slice_max()` and
  summed the evidence with a grouped `summarise()`, which dplyr
  evaluates one group at a time, tens of thousands of times per
  iteration. A single sort followed by `distinct()`, and
  [`split()`](https://rdrr.io/r/base/split.html) with
  [`sum()`](https://rdrr.io/r/base/sum.html), give the same peaks and
  the same sums to the last bit, so the consensus is unchanged. Three
  replicates of about 12,000 peaks each went from 16 s to 4.4 s. The
  calibration gains as much, since every permutation runs the same step.

## consensusRegions 0.99.0

First release.

`consensusRegions` combines the p-values of overlapping peaks across
replicated ChIP-seq, ATAC-seq and CUT&RUN experiments. A region that
convinces two samples and falls short in the third survives, where an
intersection would drop it. The method follows MSPC (Jalili *et al.*,
2015) and runs in R, with no .NET runtime to install.

### Analysis

- [`runConsensus()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/runConsensus.md)
  takes the whole workflow in one call. The individual functions stay
  available for the runs you want to inspect partway through.
- Four combination schemes. Fisher reproduces MSPC. Stouffer and
  Lancaster accept per-replicate weights, so a shallow sample counts for
  less rather than being dropped or trusted whole. The rank product
  works on ranks and needs no p-values at all.
- [`computeReplicateWeights()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/computeReplicateWeights.md)
  reads those weights from BAM files, from FRiP values you measured
  elsewhere, or from the peak sets themselves once the alignments are
  gone.
- [`calibrateThreshold()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/calibrateThreshold.md)
  estimates the cut on the combined statistic from shuffled peak sets,
  instead of leaving it at the customary `1e-8`. Give `BPPARAM` a number
  of cores to spread the permutations.

### Input

- narrowPeak p-values go straight in.
  [`readPeakSets()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/readPeakSets.md)
  rank-transforms a score column when that is all a file carries, and
  falls back to a weighted k-of-n presence rule for BED3.
- [`readPeakSets()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/readPeakSets.md)
  puts every replicate on one chromosome naming style. Two replicates
  calling the same chromosome `chr1` and `1` share no overlap, and the
  analysis would report that as biology.

### Output

- [`recentrePeaks()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/recentrePeaks.md)
  recovers summits and returns fixed-width regions, which is what motif
  analysis and count matrices want.
- [`asSummarizedExperiment()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/asSummarizedExperiment.md)
  hands the consensus to csaw, DiffBind or edgeR.
  [`exportConsensus()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/exportConsensus.md)
  writes BED or narrowPeak.
- [`plotRescue()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/plotRescue.md),
  [`plotJaccard()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/plotJaccard.md),
  [`plotUpset()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/plotUpset.md)
  and
  [`plotCalibration()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/plotCalibration.md)
  cover the quality control.

### Multiple testing

[`buildConsensus()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/buildConsensus.md)
corrects within each replicate across every peak it tested, not across
the subset that cleared the combined threshold. Correcting a family
picked for being significant is not Benjamini-Hochberg, whatever the
column is called. Set `adjustmentFamily = "confirmed"` to reproduce
MSPC.

Neither setting delivers exact false discovery rate control, and the
package does not claim it. Each replicate contributes its most
significant overlapping peak, and a maximum over a selected set does not
behave like an unselected p-value. Overlapping peaks share supporters
too, so their combined values are not independent. Reach for
[`calibrateThreshold()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/calibrateThreshold.md)
when you need a defensible error rate rather than a ranking.
