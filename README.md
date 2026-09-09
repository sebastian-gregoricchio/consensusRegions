<!-- badges: start -->
![release](https://img.shields.io/github/v/release/sebastian-gregoricchio/consensusRegions?sort=semver)
[![license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://sebastian-gregoricchio.github.io/consensusRegions/LICENSE.md/LICENSE)
[![R-CMD-check-bioc](https://github.com/sebastian-gregoricchio/consensusRegions/actions/workflows/R-CMD-check-bioc.yaml/badge.svg?branch=devel-tools)](https://github.com/sebastian-gregoricchio/consensusRegions/actions/workflows/R-CMD-check-bioc.yaml)
[![pkgdown](https://github.com/sebastian-gregoricchio/consensusRegions/actions/workflows/pkgdown.yaml/badge.svg?branch=devel-tools)](https://github.com/sebastian-gregoricchio/consensusRegions/actions/workflows/pkgdown.yaml)
[![test-coverage](https://github.com/sebastian-gregoricchio/consensusRegions/actions/workflows/test-coverage.yaml/badge.svg?branch=devel-tools)](https://github.com/sebastian-gregoricchio/consensusRegions/actions/workflows/test-coverage.yaml)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/consensusRegions?style=social)](https://github.com/sebastian-gregoricchio/consensusRegions/fork)

<!-- badges: end -->

# consensusRegions

Consensus peak regions from replicated ChIP-seq, ATAC-seq and CUT&RUN experiments.

Replicates disagree. A region can be convincing in two samples and just
below threshold in the third, and the usual response, intersecting the
peak sets, throws it away. `consensusRegions` combines the p-values of
overlapping peaks across replicates instead, so repeated weak evidence
adds up to something that survives while a strong peak nobody else saw
does not.

The approach follows MSPC ([Jalili *et al.*, 2015](https://doi.org/10.1093/bioinformatics/btv293)).
What is added here:

* **Replicate weights.** Stouffer and Lancaster combinations accept
  per-replicate weights from FRiP, library size, or the peak sets
  themselves when the BAM files are gone.
* **A calibrated threshold.** The cut on the combined p-value is
  estimated from shuffled peak sets rather than left at a fixed `1e-8`.
* **Tolerant input.** narrowPeak p-values are used directly, a plain
  score column is rank-transformed, and BED3 falls back to a weighted
  k-of-n presence rule.
* **Native R.** No .NET runtime, and the output is a `GRanges` that goes
  straight into csaw, DiffBind or ChIPseeker.

<br>


----------------------------------------


## Installation
```r
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("consensusRegions")
```

The development version:

```r
if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
remotes::install_github("sebastian-gregoricchio/consensusRegions",
                        build_manual = TRUE,
                        build_vignettes = TRUE)
```

------------------------

<br>

## Quick start
```r
library(consensusRegions)

peaks <- readPeakSets(c("rep1.narrowPeak", "rep2.narrowPeak",
                        "rep3.narrowPeak"),
                      sampleNames = c("rep1", "rep2", "rep3"))

weights <- computeReplicateWeights(peaks, method = "frip",
                                   bamFiles = c("rep1.bam", "rep2.bam",
                                                "rep3.bam"))

calibration <- calibrateThreshold(peaks, weights = weights,
                                  nPermutations = 50, targetFDR = 0.05)

result <- buildConsensus(peaks,
                         weights = weights,
                         combinationMethod = "stouffer",
                         combinedThreshold = calibration$threshold,
                         excludeRegions = blacklist)

consensusRanges(result)
consensusStats(result)
plotRescue(result)
```

### One thing worth knowing before you start
Feed this permissive input. Calling peaks at `q < 0.05` and running the
consensus on the survivors leaves nothing to rescue and turns the whole
exercise into an intersection with extra steps. Call at around
`p < 1e-3` and let `buildConsensus()` do the thresholding.

<br>

-----------------------------------------

## Documentation
`browseVignettes("consensusRegions")` after installation, or the
[web manual](https://sebastian-gregoricchio.github.io/consensusRegions/).

<br>

------------------------

## Citation
If the package is useful, please cite MSPC alongside it, since the core
approach is theirs:

> Jalili V., Matteucci M., Masseroli M., Morelli M.J. (2015). Using
> combined evidence from replicates to evaluate ChIP-seq peaks.
> *Bioinformatics* 31(17):2761-2769.


<br>

## Issues
Bug reports and suggestions in the
[issues tab](https://github.com/sebastian-gregoricchio/consensusRegions/issues).

