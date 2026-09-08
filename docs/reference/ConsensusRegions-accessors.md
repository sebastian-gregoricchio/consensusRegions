# Accessors for ConsensusRegions objects

Retrieve the components of a \[ConsensusRegions-class\] object.

## Usage

``` r
consensusRanges(object)

peakSets(object)

replicateWeights(object)

analysisParameters(object)

consensusStats(object)

# S4 method for class 'ConsensusRegions'
consensusRanges(object)

# S4 method for class 'ConsensusRegions'
peakSets(object)

# S4 method for class 'ConsensusRegions'
replicateWeights(object)

# S4 method for class 'ConsensusRegions'
analysisParameters(object)

# S4 method for class 'ConsensusRegions'
consensusStats(object)
```

## Arguments

- object:

  A \[ConsensusRegions-class\] object.

## Value

\`consensusRanges()\` returns a \`GRanges\` with the consensus regions.
\`peakSets()\` returns the annotated per-replicate peaks as a
\`GRangesList\`. \`replicateWeights()\` returns the named numeric vector
of weights. \`analysisParameters()\` returns the list of settings used.
\`consensusStats()\` returns a per-replicate summary \`data.frame\`.

## Author

Sebastian Gregoricchio

## Examples

``` r
peakFiles <- system.file("extdata",
                         c("rep1.narrowPeak", "rep2.narrowPeak",
                           "rep3.narrowPeak"),
                         package = "consensusRegions")
peaks <- readPeakSets(peakFiles, sampleNames = c("r1", "r2", "r3"))
#> Score type in use: log10pvalue
result <- buildConsensus(peaks, verbose = FALSE)

consensusRanges(result)
#> GRanges object with 279 ranges and 4 metadata columns:
#>         seqnames          ranges strand | nReplicates    nPeaks  replicates
#>            <Rle>       <IRanges>  <Rle> |   <integer> <integer> <character>
#>     [1]     chr1     88346-88911      * |           3         3    r1,r2,r3
#>     [2]     chr1   123906-124747      * |           3         3    r1,r2,r3
#>     [3]     chr1   193385-193791      * |           3         3    r1,r2,r3
#>     [4]     chr1   406763-407340      * |           3         3    r1,r2,r3
#>     [5]     chr1   479769-480472      * |           3         3    r1,r2,r3
#>     ...      ...             ...    ... .         ...       ...         ...
#>   [275]     chr2 8866200-8866881      * |           3         3    r1,r2,r3
#>   [276]     chr2 8895424-8896197      * |           3         3    r1,r2,r3
#>   [277]     chr2 8942749-8943594      * |           3         3    r1,r2,r3
#>   [278]     chr2 8971717-8972598      * |           3         3    r1,r2,r3
#>   [279]     chr2 8988237-8989040      * |           3         3    r1,r2,r3
#>         combinedNegLog10P
#>                 <numeric>
#>     [1]           13.0846
#>     [2]           25.1267
#>     [3]           11.8195
#>     [4]           25.8863
#>     [5]           26.6986
#>     ...               ...
#>   [275]           22.7642
#>   [276]           20.6173
#>   [277]           24.0368
#>   [278]           29.3225
#>   [279]           25.5599
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
consensusStats(result)
#>   replicate nTested nStringent nWeak nConfirmed nRescued nFalsePositive
#> 1        r1     370        180   190        280      100              0
#> 2        r2     370        190   180        281       91              0
#> 3        r3     370         31   339        281      250              0
#>   nDiscarded rescueRate
#> 1         90  0.3571429
#> 2         89  0.3238434
#> 3         89  0.8896797
```
