# Hand the consensus over to a differential analysis

Builds a \`RangedSummarizedExperiment\` whose rows are the consensus
regions and whose columns are the replicates. The assay holds the
per-replicate significance at each region, zero where the replicate had
no peak, together with a presence matrix.

The point of this is the row ranges. Counting reads over a fixed, shared
set of regions is what csaw, DiffBind and edgeR all want, and this puts
the consensus in the shape they expect without asking anyone to write
out a BED file and read it back.

## Usage

``` r
asSummarizedExperiment(object, assayName = "negLog10P")
```

## Arguments

- object:

  A \[ConsensusRegions-class\] object.

- assayName:

  Name given to the significance assay.

## Value

A \`RangedSummarizedExperiment\`.

## Details

The assay is not a count matrix and must not be used as one. It records
how convincing each replicate found each region, which is useful for
clustering and for spotting a replicate that disagrees with the rest,
but a differential test needs actual reads counted over
\`rowRanges(se)\`.

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

se <- asSummarizedExperiment(result)
se
#> class: RangedSummarizedExperiment 
#> dim: 292 3 
#> metadata(18): scoreType combinationMethod ... maxConsensusWidth
#>   presenceOnly
#> assays(2): negLog10P detected
#> rownames: NULL
#> rowData names(4): nReplicates nPeaks replicates combinedNegLog10P
#> colnames(3): r1 r2 r3
#> colData names(2): replicate weight
```
