# Accessors for ThresholdCalibration objects

Retrieve the components of a \[ThresholdCalibration-class\] object, the
result of \[calibrateThreshold()\].

## Usage

``` r
calibrationThreshold(object)

calibrationThresholdNegLog10(object)

calibrationCurve(object)

calibrationTargetFDR(object)

calibrationAchievedFDR(object)

calibrationSeed(object)

nullStatistics(object)

observedStatistics(object)

# S4 method for class 'ThresholdCalibration'
calibrationThreshold(object)

# S4 method for class 'ThresholdCalibration'
calibrationThresholdNegLog10(object)

# S4 method for class 'ThresholdCalibration'
calibrationCurve(object)

# S4 method for class 'ThresholdCalibration'
calibrationTargetFDR(object)

# S4 method for class 'ThresholdCalibration'
calibrationAchievedFDR(object)

# S4 method for class 'ThresholdCalibration'
calibrationSeed(object)

# S4 method for class 'ThresholdCalibration'
nullStatistics(object)

# S4 method for class 'ThresholdCalibration'
observedStatistics(object)

# S4 method for class 'ThresholdCalibration'
x$name
```

## Arguments

- object:

  A \[ThresholdCalibration-class\] object.

- x:

  A \[ThresholdCalibration-class\] object.

- name:

  Slot name, for the \`\$\` shorthand.

## Value

\`calibrationThreshold()\` returns the combined p-value threshold, and
\`calibrationThresholdNegLog10()\` the same value on the -log10 scale.
\`calibrationCurve()\` returns the FDR curve as a \`data.frame\`.
\`calibrationTargetFDR()\` and \`calibrationAchievedFDR()\` return the
FDR that was asked for and the one that was reached.
\`calibrationSeed()\` returns the seed of the permutations, or \`NULL\`.
\`nullStatistics()\` and \`observedStatistics()\` return the statistics
the curve was built from.

\`\$\` reaches the slots by name, which the result of this function
carried as a plain list up to version 0.99.2.

## Author

Sebastian Gregoricchio

## Examples

``` r
peakFiles <- system.file("extdata",
    c(
        "rep1.narrowPeak", "rep2.narrowPeak",
        "rep3.narrowPeak"
    ),
    package = "consensusRegions"
)
peaks <- readPeakSets(peakFiles,
    sampleNames = c("r1", "r2", "r3"),
    verbose = FALSE
)
calibration <- calibrateThreshold(peaks,
    nPermutations = 5, seed = 42,
    verbose = FALSE
)

calibrationThreshold(calibration)
#> [1] 4.242908e-11
calibrationSeed(calibration)
#> [1] 42
```
