# Recover summits and return fixed-width regions

Consensus regions inherit their edges from whichever peaks happened to
be merged, so their widths vary a great deal. That variation biases
every count-based and enrichment-based analysis downstream, since a
wider region collects more reads and more motif matches for reasons that
have nothing to do with the biology.

Recentring fixes it. Each region is reduced to a single position and
re-expanded to a common width, which is what the ATAC-seq peak atlases
do and what motif analysis wants.

## Usage

``` r
recentrePeaks(
  object,
  width = 400L,
  summitSource = c("best", "weighted", "centre"),
  chromosomeLengths = NULL
)
```

## Arguments

- object:

  A \[ConsensusRegions-class\] object.

- width:

  Width of the returned regions in base pairs.

- summitSource:

  How to place the anchor. \`"best"\` takes the summit of the most
  significant member peak, \`"weighted"\` averages the member summits
  weighted by their significance, and \`"centre"\` uses the midpoint of
  the merged region.

- chromosomeLengths:

  Named integer vector used to trim regions that would run off a
  chromosome end. Taken from the object when \`NULL\`.

## Value

A \`GRanges\` of fixed-width regions carrying the summit position and
the metadata of the consensus it came from.

## Details

\`"best"\` needs a summit offset, which narrowPeak files carry in their
tenth column and other formats do not. When it is missing the midpoint
of the member peak is used instead, and a warning says so.

Widths between 200 and 500 bp suit transcription factors and ATAC peaks.
Broad histone domains should not be recentred at all, since collapsing a
domain to a point throws away the thing being measured.

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

recentrePeaks(result, width = 400)
#> GRanges object with 279 ranges and 5 metadata columns:
#>         seqnames          ranges strand | nReplicates    nPeaks  replicates
#>            <Rle>       <IRanges>  <Rle> |   <integer> <integer> <character>
#>     [1]     chr1     88410-88809      * |           3         3    r1,r2,r3
#>     [2]     chr1   124000-124399      * |           3         3    r1,r2,r3
#>     [3]     chr1   193359-193758      * |           3         3    r1,r2,r3
#>     [4]     chr1   406839-407238      * |           3         3    r1,r2,r3
#>     [5]     chr1   479963-480362      * |           3         3    r1,r2,r3
#>     ...      ...             ...    ... .         ...       ...         ...
#>   [275]     chr2 8866334-8866733      * |           3         3    r1,r2,r3
#>   [276]     chr2 8895637-8896036      * |           3         3    r1,r2,r3
#>   [277]     chr2 8942945-8943344      * |           3         3    r1,r2,r3
#>   [278]     chr2 8971831-8972230      * |           3         3    r1,r2,r3
#>   [279]     chr2 8988388-8988787      * |           3         3    r1,r2,r3
#>         combinedNegLog10P    summit
#>                 <numeric> <numeric>
#>     [1]           13.0846     88610
#>     [2]           25.1267    124200
#>     [3]           11.8195    193559
#>     [4]           25.8863    407039
#>     [5]           26.6986    480163
#>     ...               ...       ...
#>   [275]           22.7642   8866534
#>   [276]           20.6173   8895837
#>   [277]           24.0368   8943145
#>   [278]           29.3225   8972031
#>   [279]           25.5599   8988588
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
