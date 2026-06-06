library(sf)
library(terra)

# Shared fixtures
dem_sr <- rast(system.file("spatial/dem_light.tiff", package = "phytorisk"))
poi_sf <- st_read(system.file("spatial/poi.geojson", package = "phytorisk"), quiet = TRUE)

# Pre-compute mec_soilwater result once; NULL if flowdem is not installed
soilwater_sr <- if (requireNamespace("flowdem", quietly = TRUE)) {
  mec_soilwater(dem_sr, poi_sf, quiet = TRUE)
} else {
  NULL
}

# Structural placeholder used in validation tests where the error fires before
# the layer-names check (steps 1-8). Same extent/CRS as dem_sr, correct names.
valid_soilwater_stub <- c(dem_sr, dem_sr)
names(valid_soilwater_stub) <- c("mec_soilwater", "streams")

# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------

test_that("mec_surfacewater() errors when dem is not a SpatRaster", {
  expect_error(mec_surfacewater(matrix(1:9, 3, 3), valid_soilwater_stub, poi_sf))
})

test_that("mec_surfacewater() errors when mec_soilwater is not a SpatRaster", {
  expect_error(mec_surfacewater(dem_sr, list(mec_soilwater = 1, streams = 1), poi_sf))
})

test_that("mec_surfacewater() errors when poi is not a single-row sf POINT", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, data.frame(x = 1, y = 2)))
})

test_that("mec_surfacewater() errors when dem and mec_soilwater have different extents", {
  small_sw <- terra::crop(
    dem_sr,
    terra::ext(xmin(dem_sr), xmax(dem_sr), ymin(dem_sr), (ymin(dem_sr) + ymax(dem_sr)) / 2)
  )
  expect_error(mec_surfacewater(dem_sr, small_sw, poi_sf))
})

test_that("mec_surfacewater() errors when dem and poi have different CRS", {
  poi_reproj <- st_transform(poi_sf, crs = 4326)
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, poi_reproj))
})

test_that("mec_surfacewater() errors when buffer is not numeric", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, poi_sf, buffer = "50"))
})

test_that("mec_surfacewater() errors when buffer is negative", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, poi_sf, buffer = -1))
})

test_that("mec_surfacewater() errors when buffer has length > 1", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, poi_sf, buffer = c(50, 100)))
})

test_that("mec_surfacewater() errors when quiet is not logical", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, poi_sf, quiet = "TRUE"))
})

test_that("mec_surfacewater() errors when mec_soilwater has wrong layer names", {
  wrong_names_sr <- c(dem_sr, dem_sr)
  names(wrong_names_sr) <- c("flow", "wetlands")
  expect_error(mec_surfacewater(dem_sr, wrong_names_sr, poi_sf))
})

# -------------------------------------------------------------------------
# assert_*() exhaustive type checks
# -------------------------------------------------------------------------

## assert_spatraster(dem) --------------------------------------------------

test_that("assert_spatraster: dem as numeric vector errors", {
  expect_error(mec_surfacewater(1:9, valid_soilwater_stub, poi_sf))
})

test_that("assert_spatraster: dem as character errors", {
  expect_error(mec_surfacewater("elevation", valid_soilwater_stub, poi_sf))
})

test_that("assert_spatraster: dem as data.frame errors", {
  expect_error(mec_surfacewater(data.frame(x = 1:3, y = 1:3), valid_soilwater_stub, poi_sf))
})

test_that("assert_spatraster: dem as NULL errors", {
  expect_error(mec_surfacewater(NULL, valid_soilwater_stub, poi_sf))
})

test_that("assert_spatraster: dem as sf object errors", {
  expect_error(mec_surfacewater(poi_sf, valid_soilwater_stub, poi_sf))
})

## assert_spatraster(mec_soilwater) ----------------------------------------

test_that("assert_spatraster: mec_soilwater as numeric vector errors", {
  expect_error(mec_surfacewater(dem_sr, 1:9, poi_sf))
})

test_that("assert_spatraster: mec_soilwater as character errors", {
  expect_error(mec_surfacewater(dem_sr, "soilwater", poi_sf))
})

test_that("assert_spatraster: mec_soilwater as list errors", {
  expect_error(mec_surfacewater(dem_sr, list(mec_soilwater = dem_sr, streams = dem_sr), poi_sf))
})

test_that("assert_spatraster: mec_soilwater as NULL errors", {
  expect_error(mec_surfacewater(dem_sr, NULL, poi_sf))
})

test_that("assert_spatraster: mec_soilwater as sf object errors", {
  expect_error(mec_surfacewater(dem_sr, poi_sf, poi_sf))
})

## assert_sf_point(poi) ----------------------------------------------------

test_that("assert_sf_point: poi as numeric vector errors", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, c(1, 2)))
})

test_that("assert_sf_point: poi as character errors", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, "poi"))
})

test_that("assert_sf_point: poi as NULL errors", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, NULL))
})

test_that("assert_sf_point: poi as SpatRaster errors", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, dem_sr))
})

test_that("assert_sf_point: poi as sf POLYGON errors", {
  poi_poly <- st_buffer(poi_sf, dist = 100)
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, poi_poly))
})

## assert_positive_numeric(buffer) -----------------------------------------

test_that("assert_positive_numeric: buffer as character errors", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, poi_sf, buffer = "50"))
})

test_that("assert_positive_numeric: buffer as logical errors", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, poi_sf, buffer = TRUE))
})

test_that("assert_positive_numeric: buffer as NULL errors", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, poi_sf, buffer = NULL))
})

test_that("assert_positive_numeric: buffer as NA errors", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, poi_sf, buffer = NA_real_))
})

## assert_logic(quiet) -----------------------------------------------------

test_that("assert_logic: quiet as character errors", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, poi_sf, quiet = "FALSE"))
})

test_that("assert_logic: quiet as numeric errors", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, poi_sf, quiet = 1))
})

test_that("assert_logic: quiet as NULL errors", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, poi_sf, quiet = NULL))
})

test_that("assert_logic: quiet as NA errors", {
  expect_error(mec_surfacewater(dem_sr, valid_soilwater_stub, poi_sf, quiet = NA))
})

# -------------------------------------------------------------------------
# Output structure
# -------------------------------------------------------------------------

test_that("mec_surfacewater() returns a named list with two elements", {
  skip_if(is.null(soilwater_sr), "flowdem not available")
  result <- mec_surfacewater(dem_sr, soilwater_sr, poi_sf, quiet = TRUE)
  expect_type(result, "list")
  expect_named(result, c("mec_surfacewater", "surface_water"))
})

test_that("mec_surfacewater$mec_surfacewater is a SpatRaster with one named layer", {
  skip_if(is.null(soilwater_sr), "flowdem not available")
  result <- mec_surfacewater(dem_sr, soilwater_sr, poi_sf, quiet = TRUE)
  expect_s4_class(result$mec_surfacewater, "SpatRaster")
  expect_equal(nlyr(result$mec_surfacewater), 1L)
  expect_equal(names(result$mec_surfacewater), "mec_surfacewater")
})

test_that("mec_surfacewater$surface_water is an sf object", {
  skip_if(is.null(soilwater_sr), "flowdem not available")
  result <- mec_surfacewater(dem_sr, soilwater_sr, poi_sf, quiet = TRUE)
  expect_s3_class(result$surface_water, "sf")
})

test_that("mec_surfacewater$mec_surfacewater CRS matches dem", {
  skip_if(is.null(soilwater_sr), "flowdem not available")
  result <- mec_surfacewater(dem_sr, soilwater_sr, poi_sf, quiet = TRUE)
  expect_true(same.crs(result$mec_surfacewater, dem_sr))
})

test_that("mec_surfacewater$mec_surfacewater extent matches dem", {
  skip_if(is.null(soilwater_sr), "flowdem not available")
  result <- mec_surfacewater(dem_sr, soilwater_sr, poi_sf, quiet = TRUE)
  expect_equal(
    as.vector(ext(result$mec_surfacewater)),
    as.vector(ext(dem_sr))
  )
})

test_that("mec_surfacewater layer values are binary (0/1)", {
  skip_if(is.null(soilwater_sr), "flowdem not available")
  result <- mec_surfacewater(dem_sr, soilwater_sr, poi_sf, quiet = TRUE)
  vals <- values(result$mec_surfacewater, na.rm = TRUE)
  expect_true(all(vals %in% c(0, 1)))
})

# -------------------------------------------------------------------------
# Behavioural
# -------------------------------------------------------------------------

test_that("mec_surfacewater() runs without error when quiet = TRUE", {
  skip_if(is.null(soilwater_sr), "flowdem not available")
  expect_no_error(mec_surfacewater(dem_sr, soilwater_sr, poi_sf, quiet = TRUE))
})
