# Derive per-replicate weights

Fisher's method assumes every replicate is worth the same. In practice
one of them is usually shallower, or has a worse signal-to-noise ratio,
and the usual responses are to drop it or to let it drag the combined
statistic down. Weighting is the middle option: the replicate still
contributes, just less.

Weights can come from three places, in decreasing order of how much
extra material they require. FRiP and library size need the BAM files or
the numbers already computed from them. The intrinsic route needs
nothing beyond the peaks themselves and scores each replicate by how
well it agrees with the others, which is the only option left when all
that survives of an old analysis is a folder of BED files.

## Usage

``` r
computeReplicateWeights(
  peakList,
  method = c("equal", "frip", "librarySize", "intrinsic", "custom"),
  bamFiles = NULL,
  frip = NULL,
  librarySize = NULL,
  weights = NULL,
  minMapq = 0L,
  verbose = TRUE
)
```

## Arguments

- peakList:

  A \`GRangesList\`, typically from \[readPeakSets()\].

- method:

  One of \`"equal"\`, \`"frip"\`, \`"librarySize"\`, \`"intrinsic"\` or
  \`"custom"\`. Default: `"equal"`.

- bamFiles:

  Character vector of indexed BAM files, in the same order as
  \`peakList\`. Needed by \`"frip"\` and \`"librarySize"\` unless the
  values are supplied directly. Default: `NULL`.

- frip:

  Numeric vector of pre-computed fractions of reads in peaks. Skips the
  BAM pass. Default: `NULL`.

- librarySize:

  Numeric vector of pre-computed mapped read counts. Default: `NULL`.

- weights:

  Numeric vector used as is when \`method\` is \`"custom"\`. Default:
  `NULL`.

- minMapq:

  Minimum mapping quality when counting from BAM. Default: `0L`.

- verbose:

  Report progress. Default: `TRUE`.

## Value

A named numeric vector with mean one, carrying the raw metrics as the
\`metrics\` attribute.

## Details

For FRiP the weight is proportional to the fraction itself. For library
size it is proportional to the square root of the depth, which is how
the information in a z-score scales, so doubling the depth is worth
rather less than twice as much. Bear in mind that peak-caller p-values
already respond to sequencing depth, so weighting by depth on top of
them counts it twice to some degree; FRiP is the safer choice when both
are available.

When \`minMapq\` is above zero both the reads in peaks and the library
total are counted under that filter, which means the whole file is read
rather than just the index.

The intrinsic weight is the mean Jaccard index between one replicate and
each of the others. A replicate that shares little with the rest is
either the odd one out biologically or the noisy one, and in both cases
letting it drive a consensus is unwise. This measure is circular by
construction, so prefer FRiP whenever the alignments are available.

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

## no BAM files at hand
computeReplicateWeights(peaks, method = "intrinsic")
#> Scoring replicates by mutual agreement
#>        r1        r2        r3 
#> 0.9925346 1.0069021 1.0005633 
#> attr(,"metrics")
#> # A tibble: 3 × 2
#>   replicate jaccard
#>   <chr>       <dbl>
#> 1 r1          0.525
#> 2 r2          0.533
#> 3 r3          0.529

## FRiP measured elsewhere
computeReplicateWeights(peaks, method = "frip",
                        frip = c(0.21, 0.19, 0.06))
#>        r1        r2        r3 
#> 1.3695652 1.2391304 0.3913043 
#> attr(,"metrics")
#> # A tibble: 3 × 3
#>   replicate  frip librarySize
#>   <chr>     <dbl>       <dbl>
#> 1 r1         0.21          NA
#> 2 r2         0.19          NA
#> 3 r3         0.06          NA
```
