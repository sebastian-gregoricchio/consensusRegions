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
#> GRanges object with 292 ranges and 5 metadata columns:
#>         seqnames          ranges strand | nReplicates    nPeaks  replicates
#>            <Rle>       <IRanges>  <Rle> |   <integer> <integer> <character>
#>     [1]     chr1     57617-58016      * |           3         3    r1,r2,r3
#>     [2]     chr1   102389-102788      * |           3         3    r1,r2,r3
#>     [3]     chr1   145949-146348      * |           3         3    r1,r2,r3
#>     [4]     chr1   235462-235861      * |           3         3    r1,r2,r3
#>     [5]     chr1   239978-240377      * |           3         3    r1,r2,r3
#>     ...      ...             ...    ... .         ...       ...         ...
#>   [288]     chr2 8631856-8632255      * |           3         3    r1,r2,r3
#>   [289]     chr2 8742403-8742802      * |           3         3    r1,r2,r3
#>   [290]     chr2 8746053-8746452      * |           3         3    r1,r2,r3
#>   [291]     chr2 8772995-8773394      * |           3         3    r1,r2,r3
#>   [292]     chr2 8862888-8863287      * |           3         3    r1,r2,r3
#>         combinedNegLog10P    summit
#>                 <numeric> <numeric>
#>     [1]           25.2442     57817
#>     [2]           22.3973    102589
#>     [3]           26.9232    146149
#>     [4]           13.5687    235662
#>     [5]           21.6011    240178
#>     ...               ...       ...
#>   [288]           21.2867   8632056
#>   [289]           12.7704   8742603
#>   [290]           30.0325   8746253
#>   [291]           31.6482   8773195
#>   [292]           25.1038   8863088
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
