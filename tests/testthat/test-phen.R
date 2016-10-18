context("phen")

test_that("phen basics", {
  phen <- athaliana_phen()
  
  expect_true(nrow(phen) == 199)
  expect_true(ncol(phen) == 107 + 3)
})

test_that("order rows", {
  phen <- athaliana_phen(traits = NULL, rows_order = "snp")
  
  expect_true(all(phen$id == athaliana_ids_snp()))
})
