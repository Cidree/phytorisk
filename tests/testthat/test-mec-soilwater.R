library(sf)
library(terra)

# Shared fixtures
dem_sr <- rast(system.file("spatial/dem_light.tiff", package = "phytorisk"))
poi_sf <- st_read(system.file("spatial/poi.geojson", package = "phytorisk"), quiet = TRUE)

# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------

test_that("mec_soilwater() errors when dem is not a SpatRaster", {
  expect_error(mec_soilwater(matrix(1:9, 3, 3), poi_sf))
})

test_that("mec_soilwater() errors when poi is not an sf object", {
  expect_error(mec_soilwater(dem_sr, data.frame(x = 1, y = 2)))
})

test_that("mec_soilwater() errors when poi has more than one row", {
  poi_multi <- rbind(poi_sf, poi_sf)
  expect_error(mec_soilwater(dem_sr, poi_multi))
})

test_that("mec_soilwater() errors when poi geometry is not POINT", {
  poi_poly <- st_buffer(poi_sf, dist = 100)
  expect_error(mec_soilwater(dem_sr, poi_poly))
})

test_that("mec_soilwater() errors when dem and poi have different CRS", {
  poi_reproj <- st_transform(poi_sf, crs = 4326)
  expect_error(mec_soilwater(dem_sr, poi_reproj))
})

test_that("mec_soilwater() errors when th is not numeric", {
  expect_error(mec_soilwater(dem_sr, poi_sf, th = "100"))
})

test_that("mec_soilwater() errors when th is negative", {
  expect_error(mec_soilwater(dem_sr, poi_sf, th = -1))
})

test_that("mec_soilwater() errors when th has length > 1", {
  expect_error(mec_soilwater(dem_sr, poi_sf, th = c(100, 200)))
})

test_that("mec_soilwater() errors when quiet is not logical", {
  expect_error(mec_soilwater(dem_sr, poi_sf, quiet = "TRUE"))
})

test_that("mec_soilwater() errors when poi does not intersect dem", {
  poi_outside <- st_sf(
    geometry = st_sfc(
      st_point(c(xmin(dem_sr) - 1e6, ymin(dem_sr) - 1e6)),
      crs = st_crs(dem_sr)
    )
  )
  expect_error(mec_soilwater(dem_sr, poi_outside))
})

# -------------------------------------------------------------------------
# assert_*() exhaustive type checks
# -------------------------------------------------------------------------

## assert_spatraster(dem) --------------------------------------------------

test_that("assert_spatraster: dem as numeric vector errors", {
  expect_error(mec_soilwater(1:9, poi_sf))
})

test_that("assert_spatraster: dem as character errors", {
  expect_error(mec_soilwater("elevation", poi_sf))
})

test_that("assert_spatraster: dem as data.frame errors", {
  expect_error(mec_soilwater(data.frame(x = 1:3, y = 1:3, z = 1:3), poi_sf))
})

test_that("assert_spatraster: dem as list errors", {
  expect_error(mec_soilwater(list(values = 1:9), poi_sf))
})

test_that("assert_spatraster: dem as NULL errors", {
  expect_error(mec_soilwater(NULL, poi_sf))
})

test_that("assert_spatraster: dem as sf object errors", {
  expect_error(mec_soilwater(poi_sf, poi_sf))
})

## assert_sf_point(poi) ----------------------------------------------------

test_that("assert_sf_point: poi as numeric vector errors", {
  expect_error(mec_soilwater(dem_sr, c(1, 2)))
})

test_that("assert_sf_point: poi as character errors", {
  expect_error(mec_soilwater(dem_sr, "poi"))
})

test_that("assert_sf_point: poi as list errors", {
  expect_error(mec_soilwater(dem_sr, list(x = 1, y = 2)))
})

test_that("assert_sf_point: poi as NULL errors", {
  expect_error(mec_soilwater(dem_sr, NULL))
})

test_that("assert_sf_point: poi as SpatRaster errors", {
  expect_error(mec_soilwater(dem_sr, dem_sr))
})

test_that("assert_sf_point: poi as sf LINESTRING errors", {
  poi_line <- st_sf(geometry = st_sfc(
    st_linestring(matrix(c(0, 0, 1, 1), ncol = 2)),
    crs = st_crs(dem_sr)
  ))
  expect_error(mec_soilwater(dem_sr, poi_line))
})

test_that("assert_sf_point: poi as sf MULTIPOINT errors", {
  poi_mpt <- st_sf(geometry = st_sfc(
    st_multipoint(matrix(c(0, 0, 1, 1), ncol = 2)),
    crs = st_crs(dem_sr)
  ))
  expect_error(mec_soilwater(dem_sr, poi_mpt))
})

## assert_positive_numeric(th) ---------------------------------------------

test_that("assert_positive_numeric: th as logical errors", {
  expect_error(mec_soilwater(dem_sr, poi_sf, th = TRUE))
})

test_that("assert_positive_numeric: th as NULL errors", {
  expect_error(mec_soilwater(dem_sr, poi_sf, th = NULL))
})

test_that("assert_positive_numeric: th as NA errors", {
  expect_error(mec_soilwater(dem_sr, poi_sf, th = NA_real_))
})

## assert_logic(quiet) -----------------------------------------------------

test_that("assert_logic: quiet as numeric errors", {
  expect_error(mec_soilwater(dem_sr, poi_sf, quiet = 1))
})

test_that("assert_logic: quiet as NULL errors", {
  expect_error(mec_soilwater(dem_sr, poi_sf, quiet = NULL))
})

test_that("assert_logic: quiet as NA errors", {
  expect_error(mec_soilwater(dem_sr, poi_sf, quiet = NA))
})

# -------------------------------------------------------------------------
# Output structure
# -------------------------------------------------------------------------

test_that("mec_soilwater() returns a SpatRaster with two named layers", {
  skip_if_not_installed("flowdem")
  result <- mec_soilwater(dem_sr, poi_sf, quiet = TRUE)
  expect_s4_class(result, "SpatRaster")
  expect_equal(nlyr(result), 2L)
  expect_equal(names(result), c("mec_soilwater", "streams"))
})

test_that("mec_soilwater() output CRS matches dem", {
  skip_if_not_installed("flowdem")
  result <- mec_soilwater(dem_sr, poi_sf, quiet = TRUE)
  expect_true(same.crs(result, dem_sr))
})

test_that("mec_soilwater() output extent matches dem", {
  skip_if_not_installed("flowdem")
  result <- mec_soilwater(dem_sr, poi_sf, quiet = TRUE)
  expect_equal(
    as.vector(ext(result)), 
    as.vector(ext(dem_sr))
  )
})

test_that("mec_soilwater layer values are binary (0/1)", {
  skip_if_not_installed("flowdem")
  result <- mec_soilwater(dem_sr, poi_sf, quiet = TRUE)
  vals <- values(result[["mec_soilwater"]], na.rm = TRUE)
  expect_true(all(vals %in% c(0, 1)))
})

test_that("streams layer values are binary (0/1)", {
  skip_if_not_installed("flowdem")
  result <- mec_soilwater(dem_sr, poi_sf, quiet = TRUE)
  vals <- values(result[["streams"]], na.rm = TRUE)
  expect_true(all(vals %in% c(0, 1)))
})

# -------------------------------------------------------------------------
# Behavioural
# -------------------------------------------------------------------------

test_that("lower th produces more stream pixels than higher th", {
  skip_if_not_installed("flowdem")
  result_low  <- mec_soilwater(dem_sr, poi_sf, th = 10,    quiet = TRUE)
  result_high <- mec_soilwater(dem_sr, poi_sf, th = 10000, quiet = TRUE)
  n_low  <- sum(values(result_low[["streams"]],  na.rm = TRUE))
  n_high <- sum(values(result_high[["streams"]], na.rm = TRUE))
  expect_gt(n_low, n_high)
})

test_that("mec_soilwater() runs without error when quiet = TRUE", {
  skip_if_not_installed("flowdem")
  expect_no_error(mec_soilwater(dem_sr, poi_sf, quiet = TRUE))
})
