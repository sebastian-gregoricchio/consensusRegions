# Consensus regions across replicated experiments

------------------------------------------------------------------------

## **Introduction**

Replicates disagree. A region can be convincing in two samples and just
below threshold in the third, and the usual response, taking the
intersection of the peak sets, throws it away. Taking the union instead
keeps it along with every piece of noise the worst replicate produced.

`consensusRegions` does neither. Each peak is examined together with the
peaks that overlap it in the other replicates, and their p-values are
combined into one. Repeated weak evidence then adds up to something that
survives, while a strong peak nobody else saw does not. The idea comes
from [MSPC](https://genometric.github.io/MSPC/docs/) ([Jalili *et al.*,
2015](https://doi.org/10.1093/bioinformatics/btv293)); what is added
here is the ability to weight replicates, to calibrate the decision
threshold rather than fix it, and to work with peak files that have lost
their statistics.

``` r

library(consensusRegions)
```

------------------------------------------------------------------------

  

## **Reading peaks**

The example data is three replicates over two chromosomes, with a
deliberately shallow third sample.

``` r
peakFiles <- system.file("extdata",
                         c("rep1.narrowPeak",
                           "rep2.narrowPeak",
                           "rep3.narrowPeak"),
                         package = "consensusRegions")

peaks <- readPeakSets(peakFiles, sampleNames = c("rep1", "rep2", "rep3"))
peaks
> GRangesList object of length 3:
> $rep1
> GRanges object with 370 ranges and 7 metadata columns:
>         seqnames          ranges strand |        name     score signalValue
>            <Rle>       <IRanges>  <Rle> | <character> <numeric>   <numeric>
>     [1]     chr1     88346-88808      * |      peak_1        54     1.50060
>     [2]     chr1   114261-114654      * |      peak_2        64     1.79430
>     [3]     chr1   123906-124747      * |      peak_3       109     3.03999
>     [4]     chr1   136167-136426      * |      peak_4        57     1.60768
>     [5]     chr1   193385-193784      * |      peak_5        52     1.46941
>     ...      ...             ...    ... .         ...       ...         ...
>   [366]     chr2 8866207-8866552      * |    peak_366        96     2.68352
>   [367]     chr2 8895500-8896058      * |    peak_367        55     1.53920
>   [368]     chr2 8942805-8943535      * |    peak_368       146     4.05697
>   [369]     chr2 8971717-8972165      * |    peak_369       130     3.62981
>   [370]     chr2 8988240-8988910      * |    peak_370       139     3.87212
>            pValue    qValue      peak negLog10P
>         <numeric> <numeric> <integer> <numeric>
>     [1]   4.50179   3.50179       261   4.50179
>     [2]   5.38290   4.38290       180   5.38290
>     [3]   9.11996   8.11996       408   9.11996
>     [4]   4.82304   3.82304       118   4.82304
>     [5]   4.40824   3.40824       201   4.40824
>     ...       ...       ...       ...       ...
>   [366]   8.05055   7.05055       184   8.05055
>   [367]   4.61759   3.61759       252   4.61759
>   [368]  12.17092  11.17092       340  12.17092
>   [369]  10.88943   9.88943       251  10.88943
>   [370]  11.61636  10.61636       348  11.61636
>   -------
>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
> 
> $rep2
> GRanges object with 370 ranges and 7 metadata columns:
>         seqnames          ranges strand |        name     score signalValue
>            <Rle>       <IRanges>  <Rle> | <character> <numeric>   <numeric>
>     [1]     chr1     51278-51539      * |      peak_1        69     1.94252
>     [2]     chr1     88369-88861      * |      peak_2        79     2.19907
>     [3]     chr1   124022-124387      * |      peak_3       144     4.00577
>     [4]     chr1   193461-193708      * |      peak_4        62     1.74570
>     [5]     chr1   406843-407195      * |      peak_5       174     4.85535
>     ...      ...             ...    ... .         ...       ...         ...
>   [366]     chr2 8866207-8866881      * |    peak_366       155     4.32219
>   [367]     chr2 8895424-8896197      * |    peak_367       133     3.69461
>   [368]     chr2 8942806-8943594      * |    peak_368        74     2.06403
>   [369]     chr2 8971768-8972291      * |    peak_369       145     4.02792
>   [370]     chr2 8988238-8988669      * |    peak_370       135     3.76715
>            pValue    qValue      peak negLog10P
>         <numeric> <numeric> <integer> <numeric>
>     [1]   5.82755   4.82755       138   5.82755
>     [2]   6.59722   5.59722       241   6.59722
>     [3]  12.01732  11.01732       178  12.01732
>     [4]   5.23709   4.23709        98   5.23709
>     [5]  14.56604  13.56604       196  14.56604
>     ...       ...       ...       ...       ...
>   [366]   12.9666   11.9666       327   12.9666
>   [367]   11.0838   10.0838       413   11.0838
>   [368]    6.1921    5.1921       405    6.1921
>   [369]   12.0838   11.0838       263   12.0838
>   [370]   11.3015   10.3015       186   11.3015
>   -------
>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
> 
> $rep3
> GRanges object with 370 ranges and 7 metadata columns:
>         seqnames          ranges strand |        name     score signalValue
>            <Rle>       <IRanges>  <Rle> | <character> <numeric>   <numeric>
>     [1]     chr1     31937-32346      * |      peak_1        48     1.36094
>     [2]     chr1     88383-88911      * |      peak_2        48     1.35000
>     [3]     chr1   123972-124259      * |      peak_3        78     2.19156
>     [4]     chr1   193460-193791      * |      peak_4        48     1.35000
>     [5]     chr1   213986-214335      * |      peak_5        48     1.34916
>     ...      ...             ...    ... .         ...       ...         ...
>   [366]     chr2 8895441-8895769      * |    peak_366        92     2.57831
>   [367]     chr2 8942749-8943558      * |    peak_367       100     2.78221
>   [368]     chr2 8956938-8957360      * |    peak_368        55     1.53327
>   [369]     chr2 8971761-8972598      * |    peak_369       104     2.90552
>   [370]     chr2 8988237-8989040      * |    peak_370        66     1.85235
>            pValue    qValue      peak negLog10P
>         <numeric> <numeric> <integer> <numeric>
>     [1]   4.08282   3.08282       191   4.08282
>     [2]   4.05000   3.05000       279   4.05000
>     [3]   6.57467   5.57467       123   6.57467
>     [4]   4.05000   3.05000       185   4.05000
>     [5]   4.04747   3.04747       170   4.04747
>     ...       ...       ...       ...       ...
>   [366]   7.73494   6.73494       163   7.73494
>   [367]   8.34663   7.34663       417   8.34663
>   [368]   4.59980   3.59980       240   4.59980
>   [369]   8.71656   7.71656       418   8.71656
>   [370]   5.55706   4.55706       419   5.55706
>   -------
>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

[`readPeakSets()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/readPeakSets.md)
inspects the files and reports which statistic it can use. narrowPeak
keeps -log10 p-values in its eighth column, so those go straight in. A
BED with only a score gets rank-transformed. A bare BED3 has nothing at
all, and the analysis will fall back to counting replicates instead of
combining evidence.

  

### Call peaks permissively

This is the mistake that ruins the whole exercise, so it is worth saying
plainly. If you call peaks at `q < 0.05` and run the consensus on what
survives, there is nothing left to rescue and you have written an
intersection with extra steps. Call at something like `p < 1e-3` and let
[`buildConsensus()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/buildConsensus.md)
do the thresholding. The two-tier design assumes exactly this.

  

------------------------------------------------------------------------

## A first run

``` r
result <- buildConsensus(peaks, verbose = FALSE)
result
> ConsensusRegions
>   replicates      : 3 ( rep1, rep2, rep3 )
>   score type      : log10pvalue 
>   combination     : stouffer 
>   combined cut    : 1e-08 
>   weights         : rep1=1, rep2=1, rep3=1 
>   consensus       : 279 regions
>   width (median)  : 779 bp
>   width (max)     : 1436 bp
```

The per-replicate summary is where to look first.

``` r
consensusStats(result)
>   replicate nTested nStringent nWeak nConfirmed nRescued nFalsePositive
> 1      rep1     370        180   190        280      100              0
> 2      rep2     370        190   180        281       91              0
> 3      rep3     370         31   339        281      250              0
>   nDiscarded rescueRate
> 1         90  0.3571429
> 2         89  0.3238434
> 3         89  0.8896797
```

`nRescued` counts peaks that were only weak on their own and were kept
because the others agreed. That number is the reason to use this method,
and also the thing to be suspicious of. If it dwarfs everything else,
`weakThreshold` is too permissive and the combination is confirming
noise.

``` r

plotRescue(result)
```

![](consensusRegions.vignette_files/figure-html/plotRescue-1.png)

The consensus regions carry the number of contributing replicates and a
significance recomputed from their member peaks.

``` r
head(consensusRanges(result), 3)
> GRanges object with 3 ranges and 4 metadata columns:
>       seqnames        ranges strand | nReplicates    nPeaks     replicates
>          <Rle>     <IRanges>  <Rle> |   <integer> <integer>    <character>
>   [1]     chr1   88346-88911      * |           3         3 rep1,rep2,rep3
>   [2]     chr1 123906-124747      * |           3         3 rep1,rep2,rep3
>   [3]     chr1 193385-193791      * |           3         3 rep1,rep2,rep3
>       combinedNegLog10P
>               <numeric>
>   [1]           13.0846
>   [2]           25.1267
>   [3]           11.8195
>   -------
>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

  

------------------------------------------------------------------------

## **Weighting the replicates**

Fisher’s method treats every replicate as equally informative. Yours
probably are not: one is shallower, or has a worse antibody, or was
prepared on a bad day. The choices are usually to drop it, which wastes
data, or to keep it and let it drag the combined statistic around.
Weighting is the third option.

When the BAM files are available, FRiP and library size can be measured
directly:

``` r

weights <- computeReplicateWeights(
    peaks,
    method = "frip",
    bamFiles = c("rep1.bam", "rep2.bam", "rep3.bam"))
```

More often the alignments are long gone and only the numbers survive in
an old QC report. Pass them in:

``` r
weights <- computeReplicateWeights(peaks, method = "frip",
                                   frip = c(0.21, 0.19, 0.06),
                                   verbose = FALSE)
weights
>      rep1      rep2      rep3 
> 1.3695652 1.2391304 0.3913043 
> attr(,"metrics")
>  [38;5;246m# A tibble: 3 × 3 [39m
>   replicate  frip librarySize
>    [3m [38;5;246m<chr> [39m [23m      [3m [38;5;246m<dbl> [39m [23m        [3m [38;5;246m<dbl> [39m [23m
>  [38;5;250m1 [39m rep1       0.21           [31mNA [39m
>  [38;5;250m2 [39m rep2       0.19           [31mNA [39m
>  [38;5;250m3 [39m rep3       0.06           [31mNA [39m
```

And when there is nothing but the peaks themselves, each replicate can
be scored by how well it agrees with the others:

``` r
computeReplicateWeights(peaks, method = "intrinsic", verbose = FALSE)
>      rep1      rep2      rep3 
> 0.9958619 1.0109825 0.9931556 
> attr(,"metrics")
>  [38;5;246m# A tibble: 3 × 2 [39m
>   replicate jaccard
>    [3m [38;5;246m<chr> [39m [23m        [3m [38;5;246m<dbl> [39m [23m
>  [38;5;250m1 [39m rep1        0.474
>  [38;5;250m2 [39m rep2        0.482
>  [38;5;250m3 [39m rep3        0.473
```

That last measure is circular, since agreement is also what the
consensus is testing. Use it when there is no alternative, not in
preference to FRiP.

Weights only reach the combination through Stouffer’s Z or Lancaster’s
method. Fisher and the rank product ignore them.

``` r
weighted <- buildConsensus(peaks, weights = weights,
                           combinationMethod = "stouffer",
                           verbose = FALSE)
consensusStats(weighted)
>   replicate nTested nStringent nWeak nConfirmed nRescued nFalsePositive
> 1      rep1     370        180   190        280      100              0
> 2      rep2     370        190   180        281       91              0
> 3      rep3     370         31   339        281      250              0
>   nDiscarded rescueRate
> 1         90  0.3571429
> 2         89  0.3238434
> 3         89  0.8896797
```

``` r

plotJaccard(weighted)
```

![](consensusRegions.vignette_files/figure-html/jaccard-1.png)

  

------------------------------------------------------------------------

## **Choosing the combination**

| Method | Weights | Use it when |
|----|----|----|
| `fisher` | no | Reproducing MSPC, or replicates really are comparable |
| `stouffer` | yes | Default. Replicate quality varies |
| `lancaster` | yes | Weights wanted, Fisher-like sensitivity to one strong peak |
| `rankProduct` | no | Only scores available, or p-values you do not trust |

One caveat about Fisher that matters more than it usually gets credit
for: the combined statistic keeps accumulating as replicates are added.
With three samples this is fine. With twenty donors, a region detected
faintly in all of them starts to look overwhelmingly significant, which
is not what anyone means by reproducible. Stouffer and the rank product
degrade more gracefully at that scale.

The threshold does not carry across methods. Fisher, Stouffer and
Lancaster all sit on roughly the same scale, so a cut chosen for one is
a reasonable starting point for the others. The rank product does not.
It works on relative ranks, so the best a peak can possibly do is come
first in every replicate, and that bound depends on how many peaks each
replicate holds. With the few hundred peaks in this example the smallest
attainable combined p-value is around `1e-5`, and asking for `1e-8`
would return nothing at all.
[`buildConsensus()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/buildConsensus.md)
checks for this and refuses rather than handing back an empty result:

``` r
buildConsensus(peaks, combinationMethod = "rankProduct", verbose = FALSE)
>  [1m [33mError [39m in `buildConsensus()`: [22m
>  [33m! [39m the rank product cannot reach a combined p-value of 1e-08 with these peak sets: the smallest attainable is 9.32e-05, because the statistic is bounded by the number of peaks per replicate. Lower 'combinedThreshold', or let calibrateThreshold() choose one
```

Give it a threshold on its own scale, or let the calibration below pick
one:

``` r
ranked <- buildConsensus(peaks, combinationMethod = "rankProduct",
                         combinedThreshold = 0.01, verbose = FALSE)
length(ranked)
> [1] 5
```

  

------------------------------------------------------------------------

## **Calibrating the threshold**

The cut on the combined p-value is conventionally `1e-8`. That number is
a guess, and the right one depends on the number of replicates, how the
peaks were called, and how densely they sit.

[`calibrateThreshold()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/calibrateThreshold.md)
estimates it. Peak positions are shuffled within each chromosome,
keeping widths, counts and statistics, so any overlap between the
shuffled replicates arises by chance. Running the combination on those
gives a null, and the threshold is read off where the expected number of
null peaks falls to a chosen fraction of the observed count.

``` r
calibration <- calibrateThreshold(peaks, nPermutations = 10,
                                  targetFDR = 0.05, seed = 1,
                                  verbose = FALSE)

calibration$threshold
> [1] 6.462449e-24
```

``` r

plotCalibration(calibration)
```

![](consensusRegions.vignette_files/figure-html/plotCalibration-1.png)

``` r
calibrated <- buildConsensus(peaks,
                             combinedThreshold = calibration$threshold,
                             verbose = FALSE)
length(calibrated)
> [1] 141
```

Ten permutations is enough for a demonstration; fifty is a reasonable
working number. Passing a blacklist through `excludeRegions` makes the
null more honest, because shuffled peaks otherwise land in artefact
regions where real peaks cluster too.

  

------------------------------------------------------------------------

## **Input without statistics**

Sometimes all that survives of an old analysis is a folder of BED3
files. There is no p-value to combine, so the package says so and
switches to a weighted presence rule: a region is kept when the weights
of the replicates supporting it sum to at least `minSupport`. With equal
weights that is the familiar k-of-n rule.

``` r
minimalFile <- system.file("extdata", "rep1_minimal.bed",
                           package = "consensusRegions")

minimal <- readPeakSets(rep(minimalFile, 3),
                        sampleNames = c("a", "b", "c"))

minimalResult <- buildConsensus(minimal, verbose = FALSE)
analysisParameters(minimalResult)$presenceOnly
> [1] TRUE
```

This is a fallback, not an equivalent. Without per-peak statistics
nothing can be rescued, and the result is reproducibility filtering
rather than combined evidence.

------------------------------------------------------------------------

  

## **Merging, and relative issues**

Peak A overlaps B, B overlaps C, and A never touches C. Transitive
merging puts all three in one region, and on permissive input this
occasionally produces regions tens of kilobases long that mean nothing.

Two defences. `maxConsensusWidth` catches oversized regions and rebuilds
only those. `mergeMethod = "iterative"` avoids the problem entirely by
seeding on the most significant peak, claiming what it overlaps, and
repeating, which is the approach the ATAC-seq peak atlases use.

``` r
seeded <- buildConsensus(peaks, mergeMethod = "iterative",
                         verbose = FALSE)

summary(GenomicRanges::width(consensusRanges(result)))
>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
>   341.0   655.5   779.0   757.4   872.0  1436.0
summary(GenomicRanges::width(consensusRanges(seeded)))
>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
>   341.0   652.0   778.5   755.0   872.0  1081.0
```

  

------------------------------------------------------------------------

## **Fixed-width regions**

Consensus regions inherit their edges from whichever peaks were merged,
so their widths vary. That variation biases everything downstream: a
wider region collects more reads and more motif matches for reasons that
have nothing to do with biology. For transcription factors and ATAC
data, recentre on the summit and use a common width.

``` r
fixed <- recentrePeaks(result, width = 400, summitSource = "best")
head(fixed, 3)
> GRanges object with 3 ranges and 5 metadata columns:
>       seqnames        ranges strand | nReplicates    nPeaks     replicates
>          <Rle>     <IRanges>  <Rle> |   <integer> <integer>    <character>
>   [1]     chr1   88410-88809      * |           3         3 rep1,rep2,rep3
>   [2]     chr1 124000-124399      * |           3         3 rep1,rep2,rep3
>   [3]     chr1 193359-193758      * |           3         3 rep1,rep2,rep3
>       combinedNegLog10P    summit
>               <numeric> <numeric>
>   [1]           13.0846     88610
>   [2]           25.1267    124200
>   [3]           11.8195    193559
>   -------
>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

Do not do this to broad histone domains. Collapsing a domain to a point
discards the thing being measured.

  

------------------------------------------------------------------------

## **Handing over to a differential analysis**

If the consensus is a step on the way to differential binding rather
than the endpoint, convert it and count reads over the row ranges with
csaw, DiffBind or edgeR.

``` r
se <- asSummarizedExperiment(result)
se
> class: RangedSummarizedExperiment 
> dim: 279 3 
> metadata(15): scoreType combinationMethod ... maxConsensusWidth
>   presenceOnly
> assays(2): negLog10P detected
> rownames: NULL
> rowData names(4): nReplicates nPeaks replicates combinedNegLog10P
> colnames(3): rep1 rep2 rep3
> colData names(2): replicate weight
```

The assay holds per-replicate significance, not counts. It is useful for
spotting a replicate that disagrees with the others; it is not a
substitute for counting reads.

And a word on where to spend your effort. If the consensus is an
intermediate, do not over-engineer it. Aggressive filtering upfront
mostly discards peaks that a proper count-based model would have handled
correctly anyway. The stringency of this package earns its keep when the
peak set is the deliverable.

  

------------------------------------------------------------------------

## Blacklisting

Filter after the consensus, not before. Removing regions upfront
distorts the rank distributions the combination relies on, which is why
`excludeRegions` is applied at the end of
[`buildConsensus()`](https://sebastian-gregoricchio.github.io/consensusRegions/reference/buildConsensus.md).

``` r

result <- buildConsensus(peaks, excludeRegions = blacklistGRanges)
```

  

------------------------------------------------------------------------

## **Checking the result**

Nothing in this package can tell you whether your consensus set is any
good. Validate it against something orthogonal: central motif enrichment
for a transcription factor, TSS enrichment and cCRE overlap for ATAC,
correlation with expression for activating marks. If the rescued peaks
are motif-poor and TSS-depleted, the settings are too permissive
regardless of what the statistics say.

============================

## **Session info**

    > R version 4.6.1 (2026-06-24)
    > Platform: x86_64-pc-linux-gnu
    > Running under: Ubuntu 24.04.4 LTS
    > 
    > Matrix products: default
    > BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    > LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    > 
    > locale:
    >  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    >  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    >  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    >  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    >  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    > [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    > 
    > time zone: Europe/Amsterdam
    > tzcode source: system (glibc)
    > 
    > attached base packages:
    > [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    > [8] base     
    > 
    > other attached packages:
    > [1] consensusRegions_0.99.0 GenomicRanges_1.64.0    Seqinfo_1.2.0          
    > [4] IRanges_2.46.0          S4Vectors_0.50.2        BiocGenerics_0.58.1    
    > [7] generics_0.1.4          BiocStyle_2.40.0       
    > 
    > loaded via a namespace (and not attached):
    >  [1] SummarizedExperiment_1.42.0 gtable_0.3.6               
    >  [3] rjson_0.2.23                xfun_0.60                  
    >  [5] bslib_0.12.0                ggplot2_4.0.3              
    >  [7] htmlwidgets_1.6.4           Biobase_2.72.0             
    >  [9] lattice_0.23-1              vctrs_0.7.3                
    > [11] tools_4.6.1                 bitops_1.1-0               
    > [13] curl_8.0.0                  parallel_4.6.1             
    > [15] tibble_3.3.1                pkgconfig_2.0.3            
    > [17] Matrix_1.7-6                RColorBrewer_1.1-3         
    > [19] cigarillo_1.2.1             S7_0.2.2                   
    > [21] desc_1.4.3                  lifecycle_1.0.5            
    > [23] farver_2.1.2                compiler_4.6.1             
    > [25] Rsamtools_2.28.0            textshaping_1.0.5          
    > [27] Biostrings_2.80.2           codetools_0.2-20           
    > [29] GenomeInfoDb_1.48.0         htmltools_0.5.9            
    > [31] sass_0.4.10                 RCurl_1.98-1.20            
    > [33] yaml_2.3.12                 tidyr_1.3.2                
    > [35] pkgdown_2.2.1               pillar_1.11.1              
    > [37] crayon_1.5.3                jquerylib_0.1.4            
    > [39] BiocParallel_1.46.0         cachem_1.1.0               
    > [41] DelayedArray_0.38.2         abind_1.4-8                
    > [43] tidyselect_1.2.1            digest_0.6.39              
    > [45] purrr_1.2.2                 restfulr_0.0.17            
    > [47] dplyr_1.2.1                 bookdown_0.48              
    > [49] labeling_0.4.3              fastmap_1.2.0              
    > [51] grid_4.6.1                  cli_3.6.6                  
    > [53] SparseArray_1.12.2          magrittr_2.0.5             
    > [55] S4Arrays_1.12.0             utf8_1.2.6                 
    > [57] XML_3.99-0.24               dichromat_2.0-1            
    > [59] withr_3.0.3                 UCSC.utils_1.8.0           
    > [61] scales_1.4.0                rmarkdown_2.32             
    > [63] XVector_0.52.0              httr_1.4.9                 
    > [65] matrixStats_1.5.0           otel_0.2.0                 
    > [67] ragg_1.5.2                  evaluate_1.0.5             
    > [69] knitr_1.51                  BiocIO_1.22.0              
    > [71] rtracklayer_1.72.0          rlang_1.3.0                
    > [73] glue_1.8.1                  BiocManager_1.30.27        
    > [75] rstudioapi_0.19.0           jsonlite_2.0.0             
    > [77] R6_2.6.1                    GenomicAlignments_1.48.0   
    > [79] MatrixGenerics_1.24.0       systemfonts_1.3.2          
    > [81] fs_2.1.0
