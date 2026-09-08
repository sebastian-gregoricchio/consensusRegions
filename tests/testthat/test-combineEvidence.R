test_that("Fisher matches the closed form computed directly", {
    pValues <- c(0.01, 0.02, 0.03)
    expected <- stats::pchisq(-2 * sum(log(pValues)), df = 6,
                              lower.tail = FALSE)

    combined <- combineEvidence(-log10(pValues), method = "fisher")

    expect_equal(10^(-combined), expected, tolerance = 1e-10)
})


test_that("Stouffer matches the closed form computed directly", {
    pValues <- c(0.01, 0.02, 0.03)
    zScores <- stats::qnorm(pValues, lower.tail = FALSE)
    expected <- stats::pnorm(sum(zScores) / sqrt(3), lower.tail = FALSE)

    combined <- combineEvidence(-log10(pValues), method = "stouffer")

    expect_equal(10^(-combined), expected, tolerance = 1e-10)
})


test_that("combining several p-values beats any of them alone", {
    negLog10P <- c(3, 3.2, 2.8)

    for (method in c("fisher", "stouffer", "lancaster")) {
        combined <- combineEvidence(negLog10P, method = method)
        expect_gt(combined, max(negLog10P))
    }
})


test_that("weights shift the result in the expected direction", {
    negLog10P <- c(8, 8, 2)

    trusting <- combineEvidence(negLog10P, method = "stouffer")
    sceptical <- combineEvidence(negLog10P,
                                 weights = c(1.3, 1.3, 0.4),
                                 method = "stouffer")

    ## down-weighting the weak replicate should raise the combined
    ## significance, since the two strong ones now carry more of it
    expect_gt(sceptical, trusting)
})


test_that("extreme p-values do not underflow", {
    ## p = 1e-320 is below the smallest representable double once
    ## exponentiated, so this only works in log space
    combined <- combineEvidence(c(320, 310, 300), method = "fisher")

    expect_true(is.finite(combined))
    expect_gt(combined, 900)
})


test_that("the rank product needs relative ranks", {
    expect_error(combineEvidence(c(3, 4), method = "rankProduct"),
                 "rho")

    combined <- combineEvidence(c(3, 4), method = "rankProduct",
                                rho = c(0.001, 0.002))
    expect_true(is.finite(combined))

    expect_error(
        combineEvidence(c(3, 4), method = "rankProduct", rho = c(0, 0.5)),
        "relative ranks")
})


test_that("an empty input is refused rather than silently returned", {
    expect_error(combineEvidence(numeric(0)), "nothing to combine")
})


test_that("the log-scale BH correction agrees with p.adjust", {
    set.seed(11)
    pValues <- runif(200, min = 1e-12, max = 0.5)

    fromLog <- consensusRegions:::.bhAdjustLog(-log10(pValues))
    fromBase <- stats::p.adjust(pValues, method = "BH")

    expect_equal(10^(-fromLog), fromBase, tolerance = 1e-10)
})
