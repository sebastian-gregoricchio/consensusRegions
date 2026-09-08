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
#> GRanges object with 292 ranges and 4 metadata columns:
#>         seqnames          ranges strand | nReplicates    nPeaks  replicates
#>            <Rle>       <IRanges>  <Rle> |   <integer> <integer> <character>
#>     [1]     chr1     57439-58142      * |           3         3    r1,r2,r3
#>     [2]     chr1   102191-103030      * |           3         3    r1,r2,r3
#>     [3]     chr1   145851-146446      * |           3         3    r1,r2,r3
#>     [4]     chr1   235338-235885      * |           3         3    r1,r2,r3
#>     [5]     chr1   240015-240714      * |           3         3    r1,r2,r3
#>     ...      ...             ...    ... .         ...       ...         ...
#>   [288]     chr2 8631832-8632669      * |           3         3    r1,r2,r3
#>   [289]     chr2 8742442-8743024      * |           3         3    r1,r2,r3
#>   [290]     chr2 8745833-8746585      * |           3         3    r1,r2,r3
#>   [291]     chr2 8772694-8773607      * |           3         3    r1,r2,r3
#>   [292]     chr2 8862514-8863473      * |           3         3    r1,r2,r3
#>         combinedNegLog10P
#>                 <numeric>
#>     [1]           25.2442
#>     [2]           22.3973
#>     [3]           26.9232
#>     [4]           13.5687
#>     [5]           21.6011
#>     ...               ...
#>   [288]           21.2867
#>   [289]           12.7704
#>   [290]           30.0325
#>   [291]           31.6482
#>   [292]           25.1038
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
consensusStats(result)
#>   replicate nTested nStringent nWeak nConfirmed nRescued nFalsePositive
#> 1        r1     430        241   189        340       99              0
#> 2        r2     430        240   190        340      100              0
#> 3        r3     430         30   400        340      310              0
#>   nDiscarded rescueRate
#> 1         90  0.2911765
#> 2         90  0.2941176
#> 3         90  0.9117647
```
