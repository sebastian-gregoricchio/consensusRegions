# Write the results to disk

Writes the consensus regions, and optionally the annotated per-replicate
peaks, as BED files. The combined significance is carried in the score
column, capped at the 1000 that the BED specification allows, so the
output loads into a genome browser without complaint.

## Usage

``` r
exportConsensus(
  object,
  file,
  format = c("bed", "narrowPeak"),
  perReplicate = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A \[ConsensusRegions-class\] object.

- file:

  Path for the consensus file.

- format:

  \`"bed"\` or \`"narrowPeak"\`. Default: `"bed"`.

- perReplicate:

  Also write one file per replicate holding every peak with its verdict.
  The names are derived from \`file\`. Default: `FALSE`.

- verbose:

  Report what was written. Default: `TRUE`.

## Value

Invisibly, a character vector of the paths written.

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

target <- file.path(tempdir(), "consensus.bed")
exportConsensus(result, target, verbose = FALSE)
```
