# ThresholdCalibration class

The result of \[calibrateThreshold()\]. It carries the threshold itself,
the FDR curve it was read off, the observed and null statistics behind
that curve, and the seed the permutations were drawn with, so that a
calibration can be inspected and repeated rather than merely trusted.

Until version 0.99.2 this was a plain list, and \`\$\` still reaches the
slots so that code written against it keeps working.

## Slots

- `threshold`:

  The combined p-value threshold.

- `thresholdNegLog10`:

  The same threshold on the -log10 scale.

- `achievedFDR`:

  The empirical FDR reached at that threshold.

- `targetFDR`:

  The FDR that was asked for.

- `fdrCurve`:

  A \`data.frame\` with one row per candidate cut, holding \`cut\`,
  \`nObserved\`, \`expectedNull\` and \`fdr\`.

- `observed`:

  The observed combined statistics.

- `null`:

  The combined statistics of the shuffled peak sets.

- `nPermutations`:

  How many shuffles the null came from.

- `seed`:

  The seed the permutations were drawn with, or \`NULL\`.

- `seedSource`:

  \`"seed"\` when the seed came from the argument, \`"BPPARAM"\` when it
  came from the backend, \`"unseeded"\` when there was none.

## Author

Sebastian Gregoricchio
