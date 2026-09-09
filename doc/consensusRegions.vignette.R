## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = ">",
                      warning = FALSE, message = FALSE, fig.align = "center",
                      dev = "png", dpi = 96, fig.width = 7, fig.height = 4.5)

## ----load---------------------------------------------------------------------
library(consensusRegions)

## ----read---------------------------------------------------------------------
peakFiles <- system.file("extdata",
                         c("rep1.narrowPeak",
                           "rep2.narrowPeak",
                           "rep3.narrowPeak"),
                         package = "consensusRegions")

peaks <- readPeakSets(peakFiles, sampleNames = c("rep1", "rep2", "rep3"))
peaks

## ----seqlevels----------------------------------------------------------------
withChr <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = c(1000, 20000), width = 500))
withoutChr <- GenomicRanges::GRanges(
    "1", IRanges::IRanges(start = c(1100, 20100), width = 500))

mixed <- readPeakSets(list(a = withChr, b = withoutChr),
                      seqlevelsStyle = "UCSC")
GenomeInfoDb::seqlevels(mixed[["b"]])

## ----basic--------------------------------------------------------------------
result <- buildConsensus(peaks, verbose = FALSE)
result

## ----stats--------------------------------------------------------------------
consensusStats(result)

## ----plotRescue---------------------------------------------------------------
plotRescue(result)

## ----consensus----------------------------------------------------------------
head(consensusRanges(result), 3)

## ----weightsBam, eval = FALSE-------------------------------------------------
# weights <- computeReplicateWeights(
#     peaks,
#     method = "frip",
#     bamFiles = c("rep1.bam", "rep2.bam", "rep3.bam"))

## ----weightsFrip--------------------------------------------------------------
weights <- computeReplicateWeights(peaks, method = "frip",
                                   frip = c(0.21, 0.19, 0.06),
                                   verbose = FALSE)
weights

## ----weightsIntrinsic---------------------------------------------------------
computeReplicateWeights(peaks, method = "intrinsic", verbose = FALSE)

## ----weighted-----------------------------------------------------------------
weighted <- buildConsensus(peaks, weights = weights,
                           combinationMethod = "stouffer",
                           verbose = FALSE)
consensusStats(weighted)

## ----jaccard------------------------------------------------------------------
plotJaccard(weighted)

## ----rankProductScale, error = TRUE-------------------------------------------
try({
buildConsensus(peaks, combinationMethod = "rankProduct", verbose = FALSE)
})

## ----rankProductWorking-------------------------------------------------------
ranked <- buildConsensus(peaks, combinationMethod = "rankProduct",
                         combinedThreshold = 0.05, verbose = FALSE)
length(ranked)

## ----rankProductRelaxed-------------------------------------------------------
relaxed <- buildConsensus(peaks, combinationMethod = "rankProduct",
                          combinedThreshold = 0.05,
                          adjustmentFamily = "confirmed", verbose = FALSE)
length(relaxed)

## ----minReplicates------------------------------------------------------------
c(byTotal = length(buildConsensus(peaks, minReplicates = 3,
                                  verbose = FALSE)),
  bySupport = length(buildConsensus(peaks, minSupport = 2,
                                    verbose = FALSE)))

## ----minReplicatesPercent-----------------------------------------------------
length(buildConsensus(peaks, minReplicates = "60%", verbose = FALSE))

## ----adjustmentFamily---------------------------------------------------------
tested <- buildConsensus(peaks, verbose = FALSE)
mspcStyle <- buildConsensus(peaks, adjustmentFamily = "confirmed",
                            verbose = FALSE)

c(tested = length(tested), confirmed = length(mspcStyle))

## ----calibrate----------------------------------------------------------------
set.seed(42)
calibration <- calibrateThreshold(peaks, nPermutations = 10,
                                  targetFDR = 0.05, verbose = FALSE)

calibration$threshold

## ----plotCalibration----------------------------------------------------------
plotCalibration(calibration)

## ----calibrated---------------------------------------------------------------
calibrated <- buildConsensus(peaks,
                             combinedThreshold = calibration$threshold,
                             verbose = FALSE)
length(calibrated)

## ----parallel, eval = FALSE---------------------------------------------------
# calibration <- calibrateThreshold(peaks, nPermutations = 50, BPPARAM = 8)

## ----parallelSeed, eval = FALSE-----------------------------------------------
# calibration <- calibrateThreshold(
#     peaks, nPermutations = 50,
#     BPPARAM = BiocParallel::MulticoreParam(workers = 8, RNGseed = 42))

## ----bed3---------------------------------------------------------------------
minimalFile <- system.file("extdata", "rep1_minimal.bed",
                           package = "consensusRegions")

minimal <- readPeakSets(rep(minimalFile, 3),
                        sampleNames = c("a", "b", "c"))

minimalResult <- buildConsensus(minimal, verbose = FALSE)
analysisParameters(minimalResult)$presenceOnly

## ----merge--------------------------------------------------------------------
seeded <- buildConsensus(peaks, mergeMethod = "iterative",
                         verbose = FALSE)

summary(GenomicRanges::width(consensusRanges(result)))
summary(GenomicRanges::width(consensusRanges(seeded)))

## ----recentre-----------------------------------------------------------------
fixed <- recentrePeaks(result, width = 400, summitSource = "best")
head(fixed, 3)

## ----se-----------------------------------------------------------------------
se <- asSummarizedExperiment(result)
se

## ----blacklist, eval = FALSE--------------------------------------------------
# result <- buildConsensus(peaks, excludeRegions = blacklistGRanges)

## ----wrapper------------------------------------------------------------------
result <- runConsensus(peakFiles,
                       sampleNames = c("rep1", "rep2", "rep3"),
                       minReplicates = 2,
                       verbose = FALSE)
length(result)

## ----wrapperOptions, eval = FALSE---------------------------------------------
# result <- runConsensus(peakFiles,
#                        sampleNames = c("rep1", "rep2", "rep3"),
#                        weightMethod = "frip",
#                        bamFiles = c("rep1.bam", "rep2.bam", "rep3.bam"),
#                        calibrate = TRUE,
#                        nPermutations = 50,
#                        combinationMethod = "stouffer",
#                        minReplicates = "75%",
#                        recentre = TRUE, width = 400,
#                        excludeRegions = blacklist,
#                        outputFile = "consensus.bed",
#                        BPPARAM = 8)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

