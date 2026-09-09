# Combine the evidence carried by overlapping peaks

Pools a set of p-values into a single one. Four schemes are offered and
they behave differently enough that the choice matters.

Fisher's method is what MSPC uses. It treats every replicate as equally
informative and is the right default when they really are comparable.
Its weakness shows up as the number of replicates grows, because the
combined statistic keeps accumulating and a region detected faintly in
many samples starts to look highly significant.

Stouffer's Z accepts weights, so a shallow replicate can be told to
count for less instead of being dropped or trusted blindly. This is the
default here.

Lancaster's method also takes weights, folding them into the degrees of
freedom of a chi-squared sum. It is closer to Fisher in spirit and tends
to react more strongly to a single very small p-value.

The rank product ignores the p-values and works on within-replicate
ranks, which makes it the natural option when the input carries only a
score. It is also the least sensitive to a replicate whose p-values are
poorly calibrated.

## Usage

``` r
combineEvidence(
  negLog10P,
  weights = NULL,
  method = c("stouffer", "fisher", "lancaster", "rankProduct"),
  rho = NULL
)
```

## Arguments

- negLog10P:

  Numeric vector of -log10 transformed p-values.

- weights:

  Numeric vector of replicate weights, recycled to the length of
  \`negLog10P\`. Ignored by \`"fisher"\` and \`"rankProduct"\`. Defaults
  to equal weights. Default: `NULL`.

- method:

  One of \`"stouffer"\`, \`"fisher"\`, \`"lancaster"\` or
  \`"rankProduct"\`. Default: `"stouffer"`.

- rho:

  Numeric vector of relative within-replicate ranks in \`(0, 1\]\`,
  required by \`"rankProduct"\` and ignored otherwise. Default: `NULL`.

## Value

A single numeric value, the combined significance on the -log10 scale.

## Details

Everything is computed in log space. Peak callers routinely report
p-values far below the smallest representable double, and exponentiating
them first would collapse the whole set to zero.

## References

Lancaster H.O. (1961). The combination of probabilities: an application
of orthonormal functions. \*Australian Journal of Statistics\* 3:20-33.

Breitling R., Armengaud P., Amtmann A., Herzyk P. (2004). Rank products:
a simple, yet powerful, new method to detect differentially regulated
genes. \*FEBS Letters\* 573(1-3):83-92.

## Author

Sebastian Gregoricchio

## Examples

``` r
## three replicates, none convincing alone
combineEvidence(c(3, 3.2, 2.8), method = "fisher")
#> [1] 6.62626

## the same peaks, but the third replicate is known to be poor
combineEvidence(c(3, 3.2, 2.8), weights = c(1.3, 1.3, 0.4),
                method = "stouffer")
#> [1] 6.521142
```
