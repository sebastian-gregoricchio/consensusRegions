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
  minSupport = NULL,
  minReplicates = NULL,
  minOverlap = 1L,
  minOverlapFraction = NULL,
  recursive = TRUE,
  maxIterations = 10L,
  multipleIntersections = c("lowest", "highest"),
  chromosomeLengths = NULL,
  excludeRegions = NULL,
  BPPARAM = 1,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- peakList:

  A \`GRangesList\` from \[readPeakSets()\].

- weights:

  Named numeric vector of replicate weights, or \`NULL\`. Default:
  `NULL`.

- combinationMethod:

  Passed to \[combineEvidence()\]. Default: `"stouffer"`.

- nPermutations:

  Number of shuffled replicate sets to generate. Default: `50L`.

- targetFDR:

  Empirical false discovery rate to aim for. Default: `0.05`.

- stringencyThreshold:

  Stringent p-value cut, as in \[buildConsensus()\]. Default: `1e-08`.

- weakThreshold:

  Background p-value cut. Default: `1e-04`.

- minSupport:

  Minimum supporting replicates. Give this or \`minReplicates\`, not
  both. Default: `NULL`.

- minReplicates:

  Total replicates that must hold the peak, counting its own. Accepts a
  count, a proportion or a percentage string, as in
  \[buildConsensus()\]. Must match the value used there. Default:
  `NULL`.

- minOverlap:

  Minimum overlap in base pairs. Default: `1L`.

- minOverlapFraction:

  Optional fractional overlap requirement. Must match the value used in
  \[buildConsensus()\]. Default: `NULL`.

- recursive:

  Whether the confirmation is re-run after pruning unsupported peaks.
  Must match the value used in \[buildConsensus()\]. Default: `TRUE`.

- maxIterations:

  Cap on the recursive rounds. Default: `10L`.

- multipleIntersections:

  Resolution rule for several overlapping peaks from one replicate.
  Default: `"lowest"`.

- chromosomeLengths:

  Named integer vector. Taken from the input, or inferred from the
  furthest peak, when left \`NULL\`. Default: `NULL`.

- excludeRegions:

  Optional \`GRanges\` the shuffled peaks must avoid. Default: `NULL`.

- BPPARAM:

  Either the number of cores to use, or a \`BiocParallelParam\` object
  for finer control. The permutations are independent of one another and
  are where nearly all the time goes, so raising this is worth it on a
  large peak set: a single round takes around half a minute on 100,000
  peaks per replicate. Default: `1`.

- seed:

  Number seeding the random positions of the shuffled peaks. It goes to
  the backend built from \`BPPARAM\` as its \`RNGseed\`, so the same
  seed returns the same threshold on any number of cores. Ignored, with
  a warning, when \`BPPARAM\` is a \`BiocParallelParam\` object: give
  the seed to it as \`RNGseed\` instead. Default: `NULL`, no seed.

- verbose:

  Report progress. Default: `TRUE`.

## Value

A list with the recommended \`threshold\` on the p-value scale, the
\`fdrCurve\` it was read off, the observed and null combined statistics,
the \`seed\` the permutations were drawn with and its \`seedSource\`:
\`"seed"\` for the argument, \`"BPPARAM"\` for the \`RNGseed\` of a
backend passed as \`BPPARAM\`, \`"unseeded"\` when there was none, in
which case \`seed\` is \`NULL\`.

## Details

Shuffling within a chromosome preserves peak density but not the
relationship between peaks and genes, so the null is a little optimistic
wherever peaks cluster around promoters. Passing \`excludeRegions\` with
a blacklist, and restricting the analysis to mappable chromosomes
beforehand, both help.

Fifty permutations are usually enough to place the threshold within a
factor of two, which is as much precision as the choice deserves. Push
it higher only if the curve looks ragged near \`targetFDR\`.

The peak positions are drawn at random inside \`BiocParallel\`, which
keeps its random numbers apart from those of the session:
\[base::set.seed()\] before the call does not reach them, on one core or
on many. \`seed\` is what brings the same threshold back, whatever the
number of cores. The function never calls \`set.seed()\`, so the random
numbers the rest of the session draws are left as they were. A
\`BiocParallelParam\` object passed as \`BPPARAM\` is used as it is and
carries its own seed, for instance
\`BiocParallel::MulticoreParam(workers = 8, RNGseed = 42)\`.

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

calibration <- calibrateThreshold(peaks, nPermutations = 5, seed = 42,
                                  verbose = FALSE)
calibration$threshold
#> [1] 4.242908e-11

result <- buildConsensus(peaks,
                         combinedThreshold = calibration$threshold,
                         verbose = FALSE)
```
