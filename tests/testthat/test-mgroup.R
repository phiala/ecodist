context("mgroup")

test_that("Mantel r is correct", {
          set.seed(888)
          x <- runif(11175)
          groups.char <- sample(letters[1:5], size=length(x), replace=TRUE)

          x.d <- dist(x)

          groups.factor <- factor(groups.char)
          groups.numeric <- as.numeric(groups.factor)

          groups.d <- dist(groups.numeric)
          groups.d[groups.d > 0] <- 1

          expect_equal(mantel(x.d ~ groups.d, nperm=0, nboot=0)$r, mgroup(x.d, groups.char, nperm=0))
})

