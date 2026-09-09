# ConsensusRegions class

Container for the result of \[buildConsensus()\]. It keeps the peaks as
they were submitted, the per-peak verdict of the combination step, the
consensus regions themselves, and the settings used to get there, so
that an analysis can be traced back from its output.

## Slots

- `peaks`:

  A \`GRangesList\`, one element per replicate, holding the peaks after
  classification. Each element carries the metadata columns
  \`negLog10P\`, \`class\`, \`nSupport\`, \`supportWeight\`,
  \`combinedNegLog10P\`, \`combinedNegLog10Padj\` and \`status\`.

- `consensus`:

  A \`GRanges\` of consensus regions.

- `weights`:

  A named numeric vector of replicate weights.

- `parameters`:

  A list with the call arguments.

- `stats`:

  A \`data.frame\` summarising each replicate.

- `calibration`:

  A list holding the output of \[calibrateThreshold()\] when one was
  supplied, empty otherwise.

## Author

Sebastian Gregoricchio
