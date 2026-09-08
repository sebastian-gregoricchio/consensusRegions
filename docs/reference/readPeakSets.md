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
#> GRanges object with 430 ranges and 7 metadata columns:
#>         seqnames          ranges strand |        name     score signalValue
#>            <Rle>       <IRanges>  <Rle> | <character> <numeric>   <numeric>
#>     [1]     chr1     57515-58142      * |      peak_1        98     2.72737
#>     [2]     chr1   102197-103030      * |      peak_2       124     3.46213
#>     [3]     chr1   145913-146341      * |      peak_3       192     5.35816
#>     [4]     chr1   217889-218140      * |      peak_4        58     1.63109
#>     [5]     chr1   235442-235885      * |      peak_5        85     2.37401
#>     ...      ...             ...    ... .         ...       ...         ...
#>   [426]     chr2 8631942-8632279      * |    peak_426        98     2.74174
#>   [427]     chr2 8742508-8742933      * |    peak_427        61     1.71393
#>   [428]     chr2 8745954-8746585      * |    peak_428       212     5.90034
#>   [429]     chr2 8772842-8773607      * |    peak_429       204     5.66718
#>   [430]     chr2 8862653-8863473      * |    peak_430       170     4.73122
#>            pValue    qValue      peak negLog10P
#>         <numeric> <numeric> <integer> <numeric>
#>     [1]   8.18212   7.18212       295   8.18212
#>     [2]  10.38638   9.38638       392  10.38638
#>     [3]  16.07448  15.07448       236  16.07448
#>     [4]   4.89327   3.89327       112   4.89327
#>     [5]   7.12203   6.12203       220   7.12203
#>     ...       ...       ...       ...       ...
#>   [426]   8.22521   7.22521       188   8.22521
#>   [427]   5.14178   4.14178       233   5.14178
#>   [428]  17.70101  16.70101       299  17.70101
#>   [429]  17.00155  16.00155       353  17.00155
#>   [430]  14.19366  13.19366       435  14.19366
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
#> 
#> $r2
#> GRanges object with 430 ranges and 7 metadata columns:
#>         seqnames          ranges strand |        name     score signalValue
#>            <Rle>       <IRanges>  <Rle> | <character> <numeric>   <numeric>
#>     [1]     chr1     57552-58092      * |      peak_1       141     3.93745
#>     [2]     chr1   102200-102630      * |      peak_2       118     3.28321
#>     [3]     chr1   145909-146277      * |      peak_3       103     2.87270
#>     [4]     chr1   201963-202312      * |      peak_4        67     1.87344
#>     [5]     chr1   235377-235868      * |      peak_5        54     1.51754
#>     ...      ...             ...    ... .         ...       ...         ...
#>   [426]     chr2 8651848-8652239      * |    peak_426        56     1.57910
#>   [427]     chr2 8742442-8742750      * |    peak_427        66     1.83676
#>   [428]     chr2 8745997-8746366      * |    peak_428       108     3.02740
#>   [429]     chr2 8772694-8773478      * |    peak_429       140     3.89196
#>   [430]     chr2 8862514-8863166      * |    peak_430       114     3.18792
#>            pValue    qValue      peak negLog10P
#>         <numeric> <numeric> <integer> <numeric>
#>     [1]  11.81234  10.81234       265  11.81234
#>     [2]   9.84962   8.84962       214   9.84962
#>     [3]   8.61810   7.61810       161   8.61810
#>     [4]   5.62032   4.62032       170   5.62032
#>     [5]   4.55263   3.55263       239   4.55263
#>     ...       ...       ...       ...       ...
#>   [426]   4.73729   3.73729       219   4.73729
#>   [427]   5.51028   4.51028       161   5.51028
#>   [428]   9.08219   8.08219       175   9.08219
#>   [429]  11.67587  10.67587       402  11.67587
#>   [430]   9.56377   8.56377       325   9.56377
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
#> 
#> $r3
#> GRanges object with 430 ranges and 7 metadata columns:
#>         seqnames          ranges strand |        name     score signalValue
#>            <Rle>       <IRanges>  <Rle> | <character> <numeric>   <numeric>
#>     [1]     chr1     57439-57955      * |      peak_1        92     2.55995
#>     [2]     chr1   102191-102497      * |      peak_2        59     1.64347
#>     [3]     chr1   145851-146446      * |      peak_3        69     1.93957
#>     [4]     chr1   235338-235729      * |      peak_4        48     1.35000
#>     [5]     chr1   240015-240714      * |      peak_5        66     1.84047
#>     ...      ...             ...    ... .         ...       ...         ...
#>   [426]     chr2 8742492-8743024      * |    peak_426        48     1.35000
#>   [427]     chr2 8745833-8746478      * |    peak_427        82     2.30257
#>   [428]     chr2 8772831-8773465      * |    peak_428        78     2.17616
#>   [429]     chr2 8862620-8863418      * |    peak_429        57     1.60899
#>   [430]     chr2 8986598-8986983      * |    peak_430        55     1.55471
#>            pValue    qValue      peak negLog10P
#>         <numeric> <numeric> <integer> <numeric>
#>     [1]   7.67984   6.67984       234   7.67984
#>     [2]   4.93041   3.93041       125   4.93041
#>     [3]   5.81871   4.81871       304   5.81871
#>     [4]   4.05000   3.05000       181   4.05000
#>     [5]   5.52140   4.52140       369   5.52140
#>     ...       ...       ...       ...       ...
#>   [426]   4.05000   3.05000       264   4.05000
#>   [427]   6.90770   5.90770       293   6.90770
#>   [428]   6.52848   5.52848       324   6.52848
#>   [429]   4.82697   3.82697       379   4.82697
#>   [430]   4.66414   3.66414       180   4.66414
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
#> 
```
