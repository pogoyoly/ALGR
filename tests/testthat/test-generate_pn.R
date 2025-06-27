test_that("percentage function works", {
  expect_equal(terra::freq(generate_pn(100, 100, 1, 2, 3, 0.01, TRUE, "land_percentage", percentage = 50))[, "count"][1], 5000)

})
