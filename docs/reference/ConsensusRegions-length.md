# Number of consensus regions

Number of consensus regions

## Usage

``` r
# S4 method for class 'ConsensusRegions'
length(x)
```

## Arguments

- x:

  A \[ConsensusRegions-class\] object.

## Value

Integer, the number of consensus regions.

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
length(buildConsensus(peaks, verbose = FALSE))
#> [1] 279
```
