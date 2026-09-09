# Calibrate the combined significance threshold

The threshold on the combined p-value is normally set by hand, and
\`1e-8\` has been carried around for long enough that it now looks like
a property of the method rather than a guess. It is a guess, and the
right value depends on how many replicates there are, how they were
called, and how densely the peaks sit.

This function estimates it instead. Peak positions are shuffled within
each chromosome while their widths, counts and statistics are kept, so
any overlap between replicates in the shuffled sets happens by chance.
Running the combination on those gives a null distribution, and the
threshold is the point where the expected number of null peaks reaching
it falls to \`targetFDR\` of the observed number.

## Usage

``` r
calibrateThreshold(
  peakList,
  weights = NULL,
  combinationMethod = c("stouffer", "fisher", "lancaster", "rankProduct"),
  nPermutations = 50L,
  targetFDR = 0.05,
  stringencyThreshold = 1e-08,
  weakThreshold = 1e-04,
  minSupport = 1L,
  minOverlap = 1L,
  minOverlapFraction = NULL,
  recursive = TRUE,
  maxIterations = 10L,
  multipleIntersections = c("lowest", "highest"),
  chromosomeLengths = NULL,
  excludeRegions = NULL,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- peakList:

  A \`GRangesList\` from \[readPeakSets()\].

- weights:

  Named numeric vector of replicate weights, or \`NULL\`.

- combinationMethod:

  Passed to \[combineEvidence()\].

- nPermutations:

  Number of shuffled replicate sets to generate.

- targetFDR:

  Empirical false discovery rate to aim for.

- stringencyThreshold:

  Stringent p-value cut, as in \[buildConsensus()\].

- weakThreshold:

  Background p-value cut.

- minSupport:

  Minimum supporting replicates.

- minOverlap:

  Minimum overlap in base pairs.

- minOverlapFraction:

  Optional fractional overlap requirement. Must match the value used in
  \[buildConsensus()\].

- recursive:

  Whether the confirmation is re-run after pruning unsupported peaks.
  Must match the value used in \[buildConsensus()\].

- maxIterations:

  Cap on the recursive rounds.

- multipleIntersections:

  Resolution rule for several overlapping peaks from one replicate.

- chromosomeLengths:

  Named integer vector. Taken from the input, or inferred from the
  furthest peak, when left \`NULL\`.

- excludeRegions:

  Optional \`GRanges\` the shuffled peaks must avoid.

- seed:

  Optional integer for reproducibility.

- verbose:

  Report progress.

## Value

A list with the recommended \`threshold\` on the p-value scale, the
\`fdrCurve\` it was read off, and the observed and null combined
statistics.

## Details

Shuffling within a chromosome preserves peak density but not the
relationship between peaks and genes, so the null is a little optimistic
wherever peaks cluster around promoters. Passing \`excludeRegions\` with
a blacklist, and restricting the analysis to mappable chromosomes
beforehand, both help.

Fifty permutations are usually enough to place the threshold within a
factor of two, which is as much precision as the choice deserves. Push
it higher only if the curve looks ragged near \`targetFDR\`.

## Author

Sebastian Gregoricchio

## Examples

``` r
peakFiles <- system.file("extdata",
                         c("rep1.narrowPeak", "rep2.narrowPeak",
                           "rep3.narrowPeak"),
                         package = "consensusRegions")
peaks <- readPeakSets(peakFiles, sampleNames = c("r1", "r2", "r3"),
                      verbose = FALSE)

calibration <- calibrateThreshold(peaks, nPermutations = 5, seed = 1,
                                  verbose = FALSE)
calibration$threshold
#> [1] 4.242908e-11

result <- buildConsensus(peaks,
                         combinedThreshold = calibration$threshold,
                         verbose = FALSE)
```
