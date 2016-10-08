context("phen")

test_that("phen basics", {
  phen <- athaliana_phen()
  
  expect_true(nrow(phen) == 199)
  expect_true(ncol(phen) == 109)
})
