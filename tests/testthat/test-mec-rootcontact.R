library(sf)
library(terra)

# Shared fixtures
trees_sr <- rast(system.file("spatial/trees_light.tiff", package = "phytorisk"))
aoi_sf   <- st_read(system.file("spatial/tejera.geojson",  package = "phytorisk"), quiet = TRUE)
poi_sf   <- st_read(system.file("spatial/poi.geojson",     package = "phytorisk"), quiet = TRUE)

# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------

test_that("mec_rootcontact() errors when treecover is not a SpatRaster", {
  expect_error(mec_rootcontact(matrix(1:9, 3, 3), aoi_sf, poi_sf))
})

test_that("mec_rootcontact() errors when aoi is not an sf object", {
  expect_error(mec_rootcontact(trees_sr, data.frame(x = 1, y = 2), poi_sf))
})

test_that("mec_rootcontact() errors when aoi geometry is not POLYGON or MULTIPOLYGON", {
  expect_error(mec_rootcontact(trees_sr, poi_sf, poi_sf))
})

test_that("mec_rootcontact() errors when poi is not an sf object", {
  expect_error(mec_rootcontact(trees_sr, aoi_sf, data.frame(x = 1, y = 2)))
})

test_that("mec_rootcontact() errors when poi has more than one row", {
  poi_multi <- rbind(poi_sf, poi_sf)
  expect_error(mec_rootcontact(trees_sr, aoi_sf, poi_multi))
})

test_that("mec_rootcontact() errors when quiet is not logical", {
  expect_error(mec_rootcontact(trees_sr, aoi_sf, poi_sf, quiet = "TRUE"))
})

test_that("mec_rootcontact() errors when treecover and aoi have different CRS", {
  aoi_reproj <- st_transform(aoi_sf, crs = 4326)
  expect_error(mec_rootcontact(trees_sr, aoi_reproj, poi_sf))
})

test_that("mec_rootcontact() errors when treecover and poi have different CRS", {
  poi_reproj <- st_transform(poi_sf, crs = 4326)
  expect_error(mec_rootcontact(trees_sr, aoi_sf, poi_reproj))
})

test_that("mec_rootcontact() errors when poi does not intersect treecover", {
  poi_outside <- st_sf(
    geometry = st_sfc(
      st_point(c(xmin(trees_sr) - 1e6, ymin(trees_sr) - 1e6)),
      crs = st_crs(trees_sr)
    )
  )
  expect_error(mec_rootcontact(trees_sr, aoi_sf, poi_outside))
})

# -------------------------------------------------------------------------
# assert_*() exhaustive type checks
# -------------------------------------------------------------------------

## assert_spatraster(treecover) --------------------------------------------

test_that("assert_spatraster: treecover as numeric vector errors", {
  expect_error(mec_rootcontact(1:9, aoi_sf, poi_sf))
})

test_that("assert_spatraster: treecover as character errors", {
  expect_error(mec_rootcontact("treecover", aoi_sf, poi_sf))
})

test_that("assert_spatraster: treecover as data.frame errors", {
  expect_error(mec_rootcontact(data.frame(x = 1:3, y = 1:3, z = 1:3), aoi_sf, poi_sf))
})

test_that("assert_spatraster: treecover as list errors", {
  expect_error(mec_rootcontact(list(values = 1:9), aoi_sf, poi_sf))
})

test_that("assert_spatraster: treecover as NULL errors", {
  expect_error(mec_rootcontact(NULL, aoi_sf, poi_sf))
})

test_that("assert_spatraster: treecover as sf object errors", {
  expect_error(mec_rootcontact(aoi_sf, aoi_sf, poi_sf))
})

## assert_sf(aoi, geometry = c("POLYGON", "MULTIPOLYGON")) ----------------

test_that("assert_sf: aoi as numeric vector errors", {
  expect_error(mec_rootcontact(trees_sr, 1:4, poi_sf))
})

test_that("assert_sf: aoi as character errors", {
  expect_error(mec_rootcontact(trees_sr, "aoi", poi_sf))
})

test_that("assert_sf: aoi as SpatRaster errors", {
  expect_error(mec_rootcontact(trees_sr, trees_sr, poi_sf))
})

test_that("assert_sf: aoi as NULL errors", {
  expect_error(mec_rootcontact(trees_sr, NULL, poi_sf))
})

test_that("assert_sf: aoi as sf POINT errors", {
  expect_error(mec_rootcontact(trees_sr, poi_sf, poi_sf))
})

test_that("assert_sf: aoi as sf LINESTRING errors", {
  aoi_line <- st_sf(geometry = st_sfc(
    st_linestring(matrix(c(0, 0, 1, 1), ncol = 2)),
    crs = st_crs(trees_sr)
  ))
  expect_error(mec_rootcontact(trees_sr, aoi_line, poi_sf))
})

## assert_sf_point(poi) ----------------------------------------------------

test_that("assert_sf_point: poi as numeric vector errors", {
  expect_error(mec_rootcontact(trees_sr, aoi_sf, c(1, 2)))
})

test_that("assert_sf_point: poi as character errors", {
  expect_error(mec_rootcontact(trees_sr, aoi_sf, "poi"))
})

test_that("assert_sf_point: poi as NULL errors", {
  expect_error(mec_rootcontact(trees_sr, aoi_sf, NULL))
})

test_that("assert_sf_point: poi as SpatRaster errors", {
  expect_error(mec_rootcontact(trees_sr, aoi_sf, trees_sr))
})

test_that("assert_sf_point: poi as sf POLYGON errors", {
  poi_poly <- st_buffer(poi_sf, dist = 100)
  expect_error(mec_rootcontact(trees_sr, aoi_sf, poi_poly))
})

test_that("assert_sf_point: poi as sf MULTIPOINT errors", {
  poi_mpt <- st_sf(geometry = st_sfc(
    st_multipoint(matrix(c(0, 0, 1, 1), ncol = 2)),
    crs = st_crs(trees_sr)
  ))
  expect_error(mec_rootcontact(trees_sr, aoi_sf, poi_mpt))
})

## assert_logic(quiet) -----------------------------------------------------

test_that("assert_logic: quiet as character errors", {
  expect_error(mec_rootcontact(trees_sr, aoi_sf, poi_sf, quiet = "FALSE"))
})

test_that("assert_logic: quiet as numeric errors", {
  expect_error(mec_rootcontact(trees_sr, aoi_sf, poi_sf, quiet = 0))
})

test_that("assert_logic: quiet as NULL errors", {
  expect_error(mec_rootcontact(trees_sr, aoi_sf, poi_sf, quiet = NULL))
})

test_that("assert_logic: quiet as NA errors", {
  expect_error(mec_rootcontact(trees_sr, aoi_sf, poi_sf, quiet = NA))
})

# -------------------------------------------------------------------------
# Output structure
# -------------------------------------------------------------------------

test_that("mec_rootcontact() returns a SpatRaster with one named layer", {
  result <- mec_rootcontact(trees_sr, aoi_sf, poi_sf, quiet = TRUE)
  expect_s4_class(result, "SpatRaster")
  expect_equal(nlyr(result), 1L)
  expect_equal(names(result), "mec_rootcontact")
})

test_that("mec_rootcontact() output CRS matches treecover", {
  result <- mec_rootcontact(trees_sr, aoi_sf, poi_sf, quiet = TRUE)
  expect_true(same.crs(result, trees_sr))
})

test_that("mec_rootcontact layer values are binary (0/1)", {
  result <- mec_rootcontact(trees_sr, aoi_sf, poi_sf, quiet = TRUE)
  vals <- values(result, na.rm = TRUE)
  expect_true(all(vals %in% c(0, 1)))
})

# -------------------------------------------------------------------------
# Behavioural
# -------------------------------------------------------------------------

test_that("mec_rootcontact() output is cropped to aoi (fewer cells than treecover)", {
  result <- mec_rootcontact(trees_sr, aoi_sf, poi_sf, quiet = TRUE)
  expect_lte(ncell(result), ncell(trees_sr))
})

test_that("mec_rootcontact() runs without error when quiet = TRUE", {
  expect_no_error(mec_rootcontact(trees_sr, aoi_sf, poi_sf, quiet = TRUE))
})
