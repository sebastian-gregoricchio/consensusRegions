#!/usr/bin/env Rscript

## Generate the example peak files shipped in inst/extdata.
##
## Three replicates over two invented chromosomes. Most positions are
## shared, a portion is deliberately left too weak for any single
## replicate to confirm so that the rescue step has something to do, and
## each replicate carries its own noise peaks so that the discard step
## does too. A BED3 copy of the first replicate is written as well, to
## exercise the path taken when the input carries no statistic at all.
##
## Run from this directory with:
##   Rscript makeExampleData.R ../extdata

set.seed(20240517)

nShared <- 220L
nWeak <- 60L
nNoise <- 90L
nChained <- 12L
chainLength <- 5L
blockSize <- 5000L
chromosomeLengths <- c(chr1 = 12000000L, chr2 = 9000000L)

## Peaks are placed one per 5 kb block, which keeps the planted positions
## clearly separated without needing a rejection loop. Blocks are drawn
## from a single shuffled pool so that no position is ever reused, and
## each chromosome contributes in proportion to its length.
blockPool <- do.call(rbind, lapply(names(chromosomeLengths), function(chr) {
    data.frame(chromosome = chr,
               block = seq_len(chromosomeLengths[[chr]] %/% blockSize) - 1L,
               stringsAsFactors = FALSE)
}))
blockPool <- blockPool[sample.int(nrow(blockPool)), ]
blocksTaken <- 0L

takePositions <- function(n, maxOffset = blockSize - 1500L) {
    selected <- blockPool[blocksTaken + seq_len(n), ]
    blocksTaken <<- blocksTaken + n

    ## a random offset inside the block, leaving room for the widest peak
    ## that will be placed there
    data.frame(chromosome = selected$chromosome,
               start = selected$block * blockSize +
                   sample.int(maxOffset, n, replace = TRUE),
               stringsAsFactors = FALSE)
}

## Independent peak calls never agree on a boundary to the base pair, so
## the shared positions are shifted a little in every replicate.
jitterPeaks <- function(start, width, amount = 90L) {
    data.frame(
        start = pmax(1L, start + sample(-amount:amount, length(start),
                                        replace = TRUE)),
        width = pmax(150L, width + sample(-40:40, length(width),
                                          replace = TRUE)))
}

writeNarrowPeak <- function(peaks, path) {
    ## the ENCODE ten-column layout: the eighth column carries -log10 p,
    ## which is what readPeakSets() picks up
    summit <- round(peaks$width / 2) +
        sample(-30:30, nrow(peaks), replace = TRUE)

    table <- data.frame(
        chromosome = peaks$chromosome,
        start = peaks$start,
        end = peaks$start + peaks$width,
        name = paste0("peak_", seq_len(nrow(peaks))),
        score = pmin(1000L, as.integer(peaks$negLog10P * 12)),
        strand = ".",
        signalValue = sprintf("%.5f", peaks$negLog10P / 3),
        pValue = sprintf("%.5f", peaks$negLog10P),
        qValue = sprintf("%.5f", pmax(0, peaks$negLog10P - 1)),
        peak = summit,
        stringsAsFactors = FALSE)

    utils::write.table(table, file = path, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = FALSE)
}

writeBed3 <- function(peaks, path) {
    table <- data.frame(chromosome = peaks$chromosome,
                        start = peaks$start,
                        end = peaks$start + peaks$width,
                        stringsAsFactors = FALSE)

    utils::write.table(table, file = path, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = FALSE)
}

## positions every replicate draws from
sharedPositions <- takePositions(nShared)
weakPositions <- takePositions(nWeak)

## A handful of loci where each replicate lays down a short run of
## partially overlapping peaks: the first and last of a run do not touch,
## so transitive merging chains them into one long region while the
## seeded rule does not. Without these the two merge methods would give
## identical output and there would be nothing to demonstrate.
chainPositions <- takePositions(nChained, maxOffset = 1000L)

## the third replicate stands in for a shallower library: the same peaks
## are there, but weaker, so it leans on the other two to confirm them
depthFactors <- c(1.00, 0.95, 0.55)
replicates <- lapply(depthFactors, function(depth) {
    sharedWidth <- sample(280:900, nShared, replace = TRUE)
    sharedShift <- jitterPeaks(sharedPositions$start, sharedWidth)
    shared <- data.frame(
        chromosome = sharedPositions$chromosome,
        start = sharedShift$start,
        width = sharedShift$width,
        negLog10P = pmax(4.2, stats::rnorm(nShared, 11, 3) * depth),
        stringsAsFactors = FALSE)

    ## these clear the weak threshold and nothing more, in any replicate
    weakWidth <- sample(250:600, nWeak, replace = TRUE)
    weakShift <- jitterPeaks(weakPositions$start, weakWidth)
    weak <- data.frame(
        chromosome = weakPositions$chromosome,
        start = weakShift$start,
        width = weakShift$width,
        negLog10P = pmax(4.05, stats::rnorm(nWeak, 5.2, 0.9) * depth),
        stringsAsFactors = FALSE)

    ## replicate-specific noise, which nothing else should support
    noisePositions <- takePositions(nNoise)
    noise <- data.frame(
        chromosome = noisePositions$chromosome,
        start = noisePositions$start,
        width = sample(200:450, nNoise, replace = TRUE),
        negLog10P = pmax(4.02, stats::rnorm(nNoise, 4.9, 0.6)),
        stringsAsFactors = FALSE)

    ## each step advances by less than one peak width, so neighbours
    ## overlap but the ends of a run do not
    chainWidth <- 600L
    chainStep <- 400L
    chained <- do.call(rbind, lapply(seq_len(chainLength), function(step) {
        data.frame(
            chromosome = chainPositions$chromosome,
            start = chainPositions$start + (step - 1L) * chainStep,
            width = chainWidth,
            negLog10P = pmax(4.2, stats::rnorm(nChained, 12, 2) * depth),
            stringsAsFactors = FALSE)
    }))

    peaks <- rbind(shared, weak, noise, chained)
    peaks[order(peaks$chromosome, peaks$start), ]
})

target <- commandArgs(trailingOnly = TRUE)
target <- if (length(target) > 0) target[1] else "../extdata"
dir.create(target, showWarnings = FALSE, recursive = TRUE)

for (i in seq_along(replicates)) {
    writeNarrowPeak(replicates[[i]],
                    file.path(target, paste0("rep", i, ".narrowPeak")))
}

## the same coordinates with every statistic stripped away
writeBed3(replicates[[1]], file.path(target, "rep1_minimal.bed"))

message("wrote ", length(replicates), " narrowPeak files to ", target)
