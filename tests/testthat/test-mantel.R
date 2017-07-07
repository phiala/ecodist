context("mantel")

test_that("mantel r is the correlation", {
          set.seed(888)
          x <- runif(11175)
          y <- runif(11175)
          expect_equal(mantel(y ~ x, nperm=0, nboot=0)$r, cor(x, y))
})
