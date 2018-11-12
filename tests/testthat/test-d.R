context("test-d")

test_that("d produces expected value", {
  expect_equal(
    round(d(
      mlogit::mFormula(subtype ~ 1 | x1 + x2 + x3),
      "subtype",
      3,
      subtype_data
    ), 3),
    0.643
  )
})


test_that("d prints message when formula mis-specified", {
  expect_error(
    d(
      formula(subtype ~ x1 + x2 + x3),
      "subtype",
      3,
      subtype_data
    ),
    "The formula argument must be of class mFormula. Please correctly specify the model formula and try again.",
    fixed = TRUE
  )
})
