context("distance")

test_that("Euclidean distance is correct", {
          set.seed(888)
          x <- matrix(runif(50), ncol=5)
          expect_equal(distance(x, "euclidean"), dist(x))
})

test_that("Bray-Curtis distance is correct", {
          set.seed(888)
          x <- matrix(runif(50), ncol=5)
          expect_equal(distance(x, "bray"), bcdist(x))
})


test_that("Mahalanobis icov is correct", {
          set.seed(888)
          x <- matrix(runif(50), ncol=5)
          x.md <- full(distance(x, "mahal"))
          sub.md <- full(distance(x[1:5, ], icov=cov(x)))
          expect_equal(x.md[1:5, 1:5], sub.md)
})

