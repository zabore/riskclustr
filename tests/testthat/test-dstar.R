context("test-dstar")

test_that("dstar produces expected value", {
  expect_equal(
    round(dstar(
      mlogit::mFormula(subtype ~ 1 | x1 + x2 + x3),
      "subtype",
      4,
      subtype_data[subtype_data$subtype > 0, ]
    ), 3),
    0.402
  )
})


test_that("dstar prints message when formula mis-specified", {
  expect_error(
    dstar(
      formula(subtype ~ x1 + x2 + x3),
      "subtype",
      4,
      subtype_data[subtype_data$subtype > 0, ]
    ),
    "The formula argument must be of class mFormula. Please correctly specify the model formula and try again.",
    fixed = TRUE
  )
})


test_that("dstar prints message when label variable is not numeric", {
  expect_error(
    dstar(
      mlogit::mFormula(subtype ~ 1 | x1 + x2 + x3),
      "subtype_name",
      4,
      subtype_data[subtype_data$subtype > 0, ]
    ),
    "The argument to label must be numeric or integer. Arguments of type character and factor are not supported, please see the documentation.",
    fixed = TRUE
  )
})


test_that("dstar prints message when subtype variable doesn't start with 1", {
  expect_error(
    dstar(
      mlogit::mFormula(subtype ~ 1 | x1 + x2 + x3),
      "subtype",
      4,
      subtype_data[subtype_data$subtype > 1, ]
    ),
    "The argument to label should start with 1. Cases should be labeled with subtypes 1 through M, the total number of subtypes.",
    fixed = TRUE
  )
})


test_that("dstar prints message when there are fewer than 2 subtype levels", {
  expect_error(
    dstar(
      mlogit::mFormula(subtype ~ 1 | x1 + x2 + x3),
      "subtype",
      1,
      subtype_data[subtype_data$subtype > 0, ]
    ),
    "The argument to M, the total number of subtypes, must be a numeric value >=2.",
    fixed = TRUE
  )
})


test_that("dstar prints message when M is not equal to the number of subtypes in the data", {
  expect_error(
    dstar(
      mlogit::mFormula(subtype ~ 1 | x1 + x2 + x3),
      "subtype",
      2,
      subtype_data[subtype_data$subtype > 0, ]
    ),
    "M is not equal to the number of levels in the variable supplied to label. Please make sure M reflects the number of subtypes in the data.",
    fixed = TRUE
  )
})
