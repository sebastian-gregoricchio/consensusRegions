# Read peak files into a GRangesList

Imports peak calls and works out which statistic can actually be used. A
narrowPeak carries -log10 p-values in its eighth column and gets the
full treatment. A BED with a score column is rank-transformed. A bare
BED3 has nothing to combine, so the analysis later falls back to a
presence rule. The detection is reported so the choice is never silent.

## Usage

``` r
readPeakSets(
  peaks,
  sampleNames = NULL,
  scoreType = c("auto", "log10pvalue", "pvalue", "score", "none"),
  scoreColumn = NULL,
  keepStandardChromosomes = TRUE,
  genome = NA,
  verbose = TRUE
)
```

## Arguments

- peaks:

  Character vector of file paths, or a \`GRanges\`, \`GRangesList\`, or
  list of \`GRanges\`.

- sampleNames:

  Character vector naming the replicates. Taken from the file names or
  the list names when left \`NULL\`.

- scoreType:

  One of \`"auto"\`, \`"log10pvalue"\`, \`"pvalue"\`, \`"score"\` or
  \`"none"\`. \`"auto"\` inspects the input.

- scoreColumn:

  Name of the metadata column holding the statistic. Only needed when
  the input is not a standard peak format.

- keepStandardChromosomes:

  Drop scaffolds and patches.

- genome:

  Genome identifier passed to the importer, for instance \`"hg38"\`.
  Fills in the chromosome lengths used by \[calibrateThreshold()\].

- verbose:

  Report what was detected.

## Value

A \`GRangesList\`, one element per replicate, with a \`negLog10P\`
metadata column when a statistic was available. The resolved score type
is stored in the object metadata.

## Details

A narrowPeak whose p-value column is filled with \`-1\`, which is what
several callers write when they have nothing to report, is treated as
score-only rather than as a set of p-values equal to \`10\`.

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
peaks
#> GRangesList object of length 3:
#> $r1
#> GRanges object with 370 ranges and 7 metadata columns:
#>         seqnames          ranges strand |        name     score signalValue
#>            <Rle>       <IRanges>  <Rle> | <character> <numeric>   <numeric>
#>     [1]     chr1     88346-88808      * |      peak_1        54     1.50060
#>     [2]     chr1   114261-114654      * |      peak_2        64     1.79430
#>     [3]     chr1   123906-124747      * |      peak_3       109     3.03999
#>     [4]     chr1   136167-136426      * |      peak_4        57     1.60768
#>     [5]     chr1   193385-193784      * |      peak_5        52     1.46941
#>     ...      ...             ...    ... .         ...       ...         ...
#>   [366]     chr2 8866207-8866552      * |    peak_366        96     2.68352
#>   [367]     chr2 8895500-8896058      * |    peak_367        55     1.53920
#>   [368]     chr2 8942805-8943535      * |    peak_368       146     4.05697
#>   [369]     chr2 8971717-8972165      * |    peak_369       130     3.62981
#>   [370]     chr2 8988240-8988910      * |    peak_370       139     3.87212
#>            pValue    qValue      peak negLog10P
#>         <numeric> <numeric> <integer> <numeric>
#>     [1]   4.50179   3.50179       261   4.50179
#>     [2]   5.38290   4.38290       180   5.38290
#>     [3]   9.11996   8.11996       408   9.11996
#>     [4]   4.82304   3.82304       118   4.82304
#>     [5]   4.40824   3.40824       201   4.40824
#>     ...       ...       ...       ...       ...
#>   [366]   8.05055   7.05055       184   8.05055
#>   [367]   4.61759   3.61759       252   4.61759
#>   [368]  12.17092  11.17092       340  12.17092
#>   [369]  10.88943   9.88943       251  10.88943
#>   [370]  11.61636  10.61636       348  11.61636
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
#> 
#> $r2
#> GRanges object with 370 ranges and 7 metadata columns:
#>         seqnames          ranges strand |        name     score signalValue
#>            <Rle>       <IRanges>  <Rle> | <character> <numeric>   <numeric>
#>     [1]     chr1     51278-51539      * |      peak_1        69     1.94252
#>     [2]     chr1     88369-88861      * |      peak_2        79     2.19907
#>     [3]     chr1   124022-124387      * |      peak_3       144     4.00577
#>     [4]     chr1   193461-193708      * |      peak_4        62     1.74570
#>     [5]     chr1   406843-407195      * |      peak_5       174     4.85535
#>     ...      ...             ...    ... .         ...       ...         ...
#>   [366]     chr2 8866207-8866881      * |    peak_366       155     4.32219
#>   [367]     chr2 8895424-8896197      * |    peak_367       133     3.69461
#>   [368]     chr2 8942806-8943594      * |    peak_368        74     2.06403
#>   [369]     chr2 8971768-8972291      * |    peak_369       145     4.02792
#>   [370]     chr2 8988238-8988669      * |    peak_370       135     3.76715
#>            pValue    qValue      peak negLog10P
#>         <numeric> <numeric> <integer> <numeric>
#>     [1]   5.82755   4.82755       138   5.82755
#>     [2]   6.59722   5.59722       241   6.59722
#>     [3]  12.01732  11.01732       178  12.01732
#>     [4]   5.23709   4.23709        98   5.23709
#>     [5]  14.56604  13.56604       196  14.56604
#>     ...       ...       ...       ...       ...
#>   [366]   12.9666   11.9666       327   12.9666
#>   [367]   11.0838   10.0838       413   11.0838
#>   [368]    6.1921    5.1921       405    6.1921
#>   [369]   12.0838   11.0838       263   12.0838
#>   [370]   11.3015   10.3015       186   11.3015
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
#> 
#> $r3
#> GRanges object with 370 ranges and 7 metadata columns:
#>         seqnames          ranges strand |        name     score signalValue
#>            <Rle>       <IRanges>  <Rle> | <character> <numeric>   <numeric>
#>     [1]     chr1     31937-32346      * |      peak_1        48     1.36094
#>     [2]     chr1     88383-88911      * |      peak_2        48     1.35000
#>     [3]     chr1   123972-124259      * |      peak_3        78     2.19156
#>     [4]     chr1   193460-193791      * |      peak_4        48     1.35000
#>     [5]     chr1   213986-214335      * |      peak_5        48     1.34916
#>     ...      ...             ...    ... .         ...       ...         ...
#>   [366]     chr2 8895441-8895769      * |    peak_366        92     2.57831
#>   [367]     chr2 8942749-8943558      * |    peak_367       100     2.78221
#>   [368]     chr2 8956938-8957360      * |    peak_368        55     1.53327
#>   [369]     chr2 8971761-8972598      * |    peak_369       104     2.90552
#>   [370]     chr2 8988237-8989040      * |    peak_370        66     1.85235
#>            pValue    qValue      peak negLog10P
#>         <numeric> <numeric> <integer> <numeric>
#>     [1]   4.08282   3.08282       191   4.08282
#>     [2]   4.05000   3.05000       279   4.05000
#>     [3]   6.57467   5.57467       123   6.57467
#>     [4]   4.05000   3.05000       185   4.05000
#>     [5]   4.04747   3.04747       170   4.04747
#>     ...       ...       ...       ...       ...
#>   [366]   7.73494   6.73494       163   7.73494
#>   [367]   8.34663   7.34663       417   8.34663
#>   [368]   4.59980   3.59980       240   4.59980
#>   [369]   8.71656   7.71656       418   8.71656
#>   [370]   5.55706   4.55706       419   5.55706
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
#> 
```
