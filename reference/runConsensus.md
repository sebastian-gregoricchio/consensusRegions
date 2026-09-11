# Run the whole analysis in one call

Chains the steps that a consensus analysis almost always performs in the
same order: read the peaks, work out the replicate weights, optionally
calibrate the threshold, build the consensus, optionally recentre it on
summits, and optionally write it out.

The individual functions remain the way to go when a step needs
inspecting before the next one runs, and calibration in particular is
worth looking at rather than trusting blind. This is for the case where
the settings are already decided and the analysis is being repeated
across experiments.

## Usage

``` r
runConsensus(
  peaks,
  sampleNames = NULL,
  weightMethod = c("equal", "frip", "librarySize", "intrinsic"),
  bamFiles = NULL,
  frip = NULL,
  librarySize = NULL,
  calibrate = FALSE,
  nPermutations = 50L,
  targetFDR = 0.05,
  recentre = FALSE,
  width = 400L,
  outputFile = NULL,
  seqlevelsStyle = "UCSC",
  excludeRegions = NULL,
  BPPARAM = 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- peaks:

  Passed to \[readPeakSets()\]: file paths, a \`GRangesList\`, or a list
  of \`GRanges\`.

- sampleNames:

  Character vector naming the replicates. Default: `NULL`.

- weightMethod:

  Passed to \[computeReplicateWeights()\]. \`"equal"\` leaves every
  replicate counting the same, which is the default and the reproducible
  choice. Default: `"equal"`.

- bamFiles, frip, librarySize:

  Passed to \[computeReplicateWeights()\] when \`weightMethod\` needs
  them. Default: `NULL`.

- calibrate:

  Estimate the combined threshold with \[calibrateThreshold()\] instead
  of using \`combinedThreshold\`. Default: `FALSE`.

- nPermutations:

  Number of permutations, passed to \[calibrateThreshold()\]. Default:
  `50L`.

- targetFDR:

  Empirical false discovery rate to aim for, passed to
  \[calibrateThreshold()\]. Default: `0.05`.

- recentre:

  Reduce the consensus to fixed-width regions around summits with
  \[recentrePeaks()\]. Sensible for transcription factors and ATAC,
  wrong for broad domains. Default: `FALSE`.

- width:

  Width of the recentred regions. Default: `400L`.

- outputFile:

  Optional path. When given the consensus is written there with
  \[exportConsensus()\]. Default: `NULL`.

- seqlevelsStyle:

  Chromosome naming style, passed to \[readPeakSets()\]. Default:
  `"UCSC"`.

- excludeRegions:

  Optional \`GRanges\` of blacklisted positions. Default: `NULL`.

- BPPARAM:

  Either the number of cores to use, or a \`BiocParallelParam\` object.
  Only the calibration uses it. Default: `1`.

- verbose:

  Report progress through the steps. Default: `TRUE`.

- ...:

  Further arguments passed to \[buildConsensus()\], for example
  \`combinationMethod\`, \`minReplicates\` or \`mergeMethod\`.

## Value

A \[ConsensusRegions-class\] object. When \`recentre = TRUE\` the
consensus it carries is the fixed-width version.

## See also

\[readPeakSets()\], \[computeReplicateWeights()\],
\[calibrateThreshold()\], \[buildConsensus()\], \[recentrePeaks()\]

## Author

Sebastian Gregoricchio

## Examples

``` r
peakFiles <- system.file("extdata",
    c("rep1.narrowPeak", "rep2.narrowPeak",
        "rep3.narrowPeak"),
    package = "consensusRegions")

result <- runConsensus(peakFiles,
    sampleNames = c("r1", "r2", "r3"),
    minReplicates = 2,
    verbose = FALSE)
result
#> ConsensusRegions
#>   replicates      : 3 ( r1, r2, r3 )
#>   score type      : log10pvalue 
#>   combination     : stouffer 
#>   combined cut    : 1e-08 
#>   weights         : r1=1, r2=1, r3=1 
#>   consensus       : 292 regions
#>   width (median)  : 786.5 bp
#>   width (max)     : 2200 bp

## the same with measured weights and a calibrated threshold. Five
## permutations are enough to show the shape of it; fifty is the
## working number, and the positions are drawn at random, so the seed
## is what brings the same threshold back.
set.seed(42)
calibrated <- runConsensus(peakFiles,
    sampleNames = c("r1", "r2", "r3"),
    weightMethod = "frip",
    frip = c(0.21, 0.19, 0.06),
    calibrate = TRUE,
    nPermutations = 5,
    verbose = FALSE)
calibrated
#> ConsensusRegions
#>   replicates      : 3 ( r1, r2, r3 )
#>   score type      : log10pvalue 
#>   combination     : stouffer 
#>   combined cut    : 1.294522e-09 
#>   weights         : r1=1.37, r2=1.239, r3=0.391 
#>   consensus       : 292 regions
#>   width (median)  : 786.5 bp
#>   width (max)     : 2200 bp
#>   calibrated at   : FDR 0.05 
```
