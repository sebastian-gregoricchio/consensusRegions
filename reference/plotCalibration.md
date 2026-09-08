# The calibration curve

Plots the empirical false discovery rate against the cut on the combined
statistic, with the chosen threshold marked. A curve that never descends
to the target says the replicates do not agree beyond chance often
enough to support a consensus at that stringency.

## Usage

``` r
plotCalibration(calibration, baseSize = 12)
```

## Arguments

- calibration:

  Output of \[calibrateThreshold()\].

- baseSize:

  Base font size in points.

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
calibration <- calibrateThreshold(peaks, nPermutations = 5, seed = 1,
                                  verbose = FALSE)

plotCalibration(calibration)

```
