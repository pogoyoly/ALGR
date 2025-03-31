set.seed(123)

r<-terra::rast(matrix(1, nrow=50, ncol=50))
map<-ALGR::establish_pac(potential_space= r,
                         cell_size=1,
                         includsion_value = 1,
                         mean_field_size = 30,
                         sd_field_size = 25,
                         distribution = "norm",
                         mean_shape_index = 1,
                         sd_shape_index = 0.3,
                         percent = 75)


crops_matrix <- data.frame(crop = c("Wheat", "Corn"),
                           percentage = c(0.5, 0.5),
                           index = c(1,2))  # Desired percentages

outcome<-ALGR::distrubute_crops(map,crops_matrix)


raster_crop<-ALGR::return_by_crop(outcome, method = 2)

test_that("Test that output works", {
  expect_equal(class(raster_crop) == 'SpatRaster', TRUE)
  expect_equal(terra::freq(raster_crop)[, "count"][2], 947)
  expect_equal(terra::freq(raster_crop)[, "count"][3], 947)

})
