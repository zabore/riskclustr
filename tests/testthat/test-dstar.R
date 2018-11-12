context("test-dstar")

test_that("dstar produces expected value", {
  expect_equal(
    round(dstar(
      mlogit::mFormula(subtype ~ 1 | x1 + x2 + x3),
      "subtype",
      3,
      subtype_data[subtype_data$subtype > 0, ]
    ), 3),
    0.559
  )
})


test_that("dstar prints message when class variable is not numeric", {
  expect_error(
    dstar(
      mlogit::mFormula(subtype ~ 1 | x1 + x2 + x3),
      "subtype_name",
      3,
      subtype_data[subtype_data$subtype > 0, ]
    ),
    "The argument to label must be numeric or integer. Arguments of type character and factor are not supported, please see the documentation.",
    fixed = TRUE
  )
})
