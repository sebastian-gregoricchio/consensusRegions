# consensusRegions — `devel-tools`

Build and check machinery for
[consensusRegions](https://github.com/sebastian-gregoricchio/consensusRegions).
The package source lives on `main`. This branch carries only the things
that check it, so nothing here ends up in the tarball.

[![R-CMD-check-bioc](https://github.com/sebastian-gregoricchio/consensusRegions/actions/workflows/R-CMD-check-bioc.yaml/badge.svg?branch=devel-tools)](https://github.com/sebastian-gregoricchio/consensusRegions/actions/workflows/R-CMD-check-bioc.yaml?query=branch%3Adevel-tools)
[![pkgdown](https://github.com/sebastian-gregoricchio/consensusRegions/actions/workflows/pkgdown.yaml/badge.svg?branch=devel-tools)](https://github.com/sebastian-gregoricchio/consensusRegions/actions/workflows/pkgdown.yaml?query=branch%3Adevel-tools)
[![Codecov](https://codecov.io/gh/sebastian-gregoricchio/consensusRegions/branch/devel-tools/graph/badge.svg)](https://app.codecov.io/gh/sebastian-gregoricchio/consensusRegions)

## Workflows

| File | What it does |
| --- | --- |
| `R-CMD-check-bioc.yaml` | `R CMD check` on Linux, macOS and Windows. `BiocCheck` and `BiocCheckGitClone` run on the R-devel job only, since that is the branch a submission is reviewed against. |
| `pkgdown.yaml` | Rebuilds the reference site and pushes it to GitHub Pages. |
| `test-coverage.yaml` | Runs `covr` and uploads to Codecov. |

Every push here starts all three. The full check takes around forty minutes,
most of it spent building Bioconductor dependencies on Windows and macOS.

## Working on the package

Code, documentation and vignettes go on `main`. If you edit them here they
will not be checked, and the two branches drift apart quietly.

Before opening a pull request against `main`:

    R CMD build .
    R CMD check --no-manual consensusRegions_*.tar.gz
    BiocCheck::BiocCheck("consensusRegions_*.tar.gz")

## Badges

Runs only ever happen on this branch, so a badge without `?branch=devel-tools`
looks at `main`, finds nothing, and renders as a failure. The README on `main`
has to spell the branch out.

## Issues

Bug reports and suggestions go in the
[issues tab](https://github.com/sebastian-gregoricchio/consensusRegions/issues),
not here.
