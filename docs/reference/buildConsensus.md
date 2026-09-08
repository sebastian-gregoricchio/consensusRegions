# Build consensus regions from replicated peak calls

Runs the whole analysis. Each peak is classified as stringent, weak or
background by two p-value thresholds; the peaks overlapping it in the
other replicates are collected; their evidence is combined; and the
result is tested against a threshold on the combined significance and
then corrected within each replicate. Peaks that survive are merged into
consensus regions.

Feed this permissive input. Calling peaks at \`q \< 0.05\` and running
the consensus on the survivors leaves nothing to rescue and turns the
whole exercise into an intersection with extra steps. Something around
\`p \< 1e-3\` gives the combination room to work.

## Usage

``` r
buildConsensus(
  peakList,
  weights = NULL,
  combinationMethod = c("stouffer", "fisher", "lancaster", "rankProduct"),
  replicateType = c("biological", "technical"),
  stringencyThreshold = 1e-08,
  weakThreshold = 1e-04,
  combinedThreshold = NULL,
  alpha = 0.05,
  minSupport = 1L,
  minOverlap = 1L,
  minOverlapFraction = NULL,
  multipleIntersections = c("lowest", "highest"),
  recursive = TRUE,
  maxIterations = 10L,
  mergeMethod = c("reduce", "iterative"),
  maxConsensusWidth = NULL,
  excludeRegions = NULL,
  verbose = TRUE
)
```

## Arguments

- peakList:

  A \`GRangesList\` from \[readPeakSets()\].

- weights:

  Named numeric vector of replicate weights, or \`NULL\` for equal
  weights. See \[computeReplicateWeights()\].

- combinationMethod:

  Passed to \[combineEvidence()\].

- replicateType:

  \`"biological"\` asks for \`minSupport\` supporting replicates;
  \`"technical"\` asks for all of them.

- stringencyThreshold:

  P-value below which a peak is stringent.

- weakThreshold:

  P-value above which a peak is background and takes no further part.

- combinedThreshold:

  Threshold on the combined p-value. Defaults to
  \`stringencyThreshold\`; \[calibrateThreshold()\] can estimate it
  instead. The scale differs between combination schemes, so a value
  carried over from one will not mean the same under another. The rank
  product in particular is bounded by the number of peaks per replicate
  and usually needs a far more permissive threshold.

- alpha:

  Level for the within-replicate Benjamini-Hochberg step.

- minSupport:

  Number of \*other\* replicates that must hold an overlapping peak.

- minOverlap:

  Minimum overlap in base pairs.

- minOverlapFraction:

  Optional minimum overlap as a fraction of the shorter of the two
  peaks. Applied on top of \`minOverlap\`.

- multipleIntersections:

  How to pick among several overlapping peaks from the same replicate:
  \`"lowest"\` takes the smallest p-value, \`"highest"\` the largest.

- recursive:

  Re-run the confirmation after dropping discarded peaks from the pool
  of eligible supporters.

- maxIterations:

  Cap on the recursive rounds.

- mergeMethod:

  \`"reduce"\` merges anything that touches; \`"iterative"\` seeds on
  the most significant peak and removes what it overlaps, then repeats.

- maxConsensusWidth:

  Optional cap in base pairs. Merged regions wider than this are rebuilt
  with the iterative rule.

- excludeRegions:

  Optional \`GRanges\` of blacklisted positions, removed from the
  consensus at the end.

- verbose:

  Report progress.

## Value

A \[ConsensusRegions-class\] object.

## Details

Two behaviours are worth knowing about.

Transitive merging can run away. Peak A overlaps B, B overlaps C, and A
never touches C, yet \`"reduce"\` puts all three in one region. On
permissive input this occasionally produces regions tens of kilobases
long. \`maxConsensusWidth\` catches them and rebuilds only the
offenders.

Recursive confirmation makes a peak's fate depend on peaks that were
themselves discarded, which is arguably the more defensible reading of
the method but is not what the original implementation does. Set
\`recursive = FALSE\` to reproduce a single-pass analysis.

When the peak sets carry no usable statistic the combination is skipped
and a weighted presence rule is applied instead: a peak is kept when the
weights of the replicates supporting it sum to at least \`minSupport\`.
With equal weights this is the familiar k-of-n rule.

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

## down-weighting a replicate known to be poor
weighted <- buildConsensus(
    peaks,
    weights = computeReplicateWeights(peaks, method = "frip",
                                      frip = c(0.2, 0.18, 0.05)),
    combinationMethod = "stouffer",
    verbose = FALSE)
consensusStats(weighted)
#>   replicate nTested nStringent nWeak nConfirmed nRescued nFalsePositive
#> 1        r1     370        180   190        280      100              0
#> 2        r2     370        190   180        281       91              0
#> 3        r3     370         31   339        281      250              0
#>   nDiscarded rescueRate
#> 1         90  0.3571429
#> 2         89  0.3238434
#> 3         89  0.8896797
```
