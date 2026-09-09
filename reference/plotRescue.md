# Rescue and discard rates per replicate

Shows what happened to the peaks of each replicate. The bar to watch is
the rescued one: peaks that were only weak on their own and were kept
because the other replicates agreed. That is the whole point of the
method, but a rescued fraction far above what the replicates otherwise
share is a sign that \`weakThreshold\` was set too permissively and the
combination is confirming noise.

## Usage

``` r
plotRescue(object, proportion = FALSE, baseSize = 12)
```

## Arguments

- object:

  A \[ConsensusRegions-class\] object.

- proportion:

  Scale the bars to one instead of showing counts. Default: `FALSE`.

- baseSize:

  Base font size in points. Default: `12`.

## Value

A \`ggplot\` object.

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
result <- buildConsensus(peaks, verbose = FALSE)

plotRescue(result)

```
