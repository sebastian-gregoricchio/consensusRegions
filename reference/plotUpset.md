# Which replicates contribute to each consensus region

An UpSet plot of the replicate combinations found across the consensus
regions. Regions supported by every replicate should dominate; a large
partial set points at a replicate that is pulling in a different
direction.

## Usage

``` r
plotUpset(object, minSize = 0)
```

## Arguments

- object:

  A \[ConsensusRegions-class\] object.

- minSize:

  Drop intersections smaller than this.

## Value

A \`ggplot\` object built by \`ComplexUpset\`.

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

if (requireNamespace("ComplexUpset", quietly = TRUE)) {
    plotUpset(result)
}
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the ComplexUpset package.
#>   Please report the issue at
#>   <https://github.com/krassowski/complex-upset/issues>.

```
