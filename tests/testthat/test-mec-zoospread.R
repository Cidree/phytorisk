library(sf)
library(terra)

# Shared fixtures
aoi_sf <- st_read(system.file("spatial/tejera.geojson", package = "phytorisk"), quiet = TRUE)
poi_sf <- st_read(system.file("spatial/poi.geojson",    package = "phytorisk"), quiet = TRUE)
dem_sr <- rast(system.file("spatial/dem_light.tiff",    package = "phytorisk"))

# Structural mock of mec_surfacewater() output — no flowdem dependency.
# dem_sr shares CRS with poi_sf/aoi_sf; poi_sf is a valid 1-row sf for
# surface_water (st_coordinates()[1, ] yields finite food coordinates).
mock_raster <- dem_sr
names(mock_raster) <- "mec_surfacewater"
mock_surface_lst <- list(
  mec_surfacewater = mock_raster,
  surface_water    = poi_sf
)

# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------

test_that("mec_zoospread() errors when mec_surfacewater is not a list", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, dem_sr))
})

test_that("mec_zoospread() errors when mec_surfacewater$mec_surfacewater is not a SpatRaster", {
  bad_lst <- list(mec_surfacewater = 1:4, surface_water = poi_sf)
  expect_error(mec_zoospread(aoi_sf, poi_sf, bad_lst))
})

test_that("mec_zoospread() errors when mec_surfacewater$surface_water is not an sf object", {
  bad_lst <- list(mec_surfacewater = dem_sr, surface_water = 1:4)
  expect_error(mec_zoospread(aoi_sf, poi_sf, bad_lst))
})

test_that("mec_zoospread() errors when aoi is not an sf object", {
  expect_error(mec_zoospread(data.frame(x = 1, y = 2), poi_sf, mock_surface_lst))
})

test_that("mec_zoospread() errors when aoi geometry is not POLYGON or MULTIPOLYGON", {
  expect_error(mec_zoospread(poi_sf, poi_sf, mock_surface_lst))
})

test_that("mec_zoospread() errors when poi is not a single-row sf POINT", {
  expect_error(mec_zoospread(aoi_sf, data.frame(x = 1, y = 2), mock_surface_lst))
})

test_that("mec_zoospread() errors when poi has more than one row", {
  poi_multi <- rbind(poi_sf, poi_sf)
  expect_error(mec_zoospread(aoi_sf, poi_multi, mock_surface_lst))
})

test_that("mec_zoospread() errors when poi and aoi have different CRS", {
  aoi_reproj <- st_transform(aoi_sf, crs = 4326)
  expect_error(mec_zoospread(aoi_reproj, poi_sf, mock_surface_lst))
})

test_that("mec_zoospread() errors when mec_surfacewater raster and poi have different CRS", {
  wrong_crs_raster <- dem_sr
  terra::crs(wrong_crs_raster) <- "EPSG:4326"
  names(wrong_crs_raster) <- "mec_surfacewater"
  wrong_crs_lst <- list(mec_surfacewater = wrong_crs_raster, surface_water = poi_sf)
  expect_error(mec_zoospread(aoi_sf, poi_sf, wrong_crs_lst))
})

test_that("mec_zoospread() errors when n_animals is not an integer scalar", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, n_animals = 1.5))
})

test_that("mec_zoospread() errors when n_steps is not an integer scalar", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, n_steps = 2.5))
})

test_that("mec_zoospread() errors when pixel_size is negative", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, pixel_size = -1))
})

test_that("mec_zoospread() errors when n_iter is not an integer scalar", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, n_iter = 1.5))
})

test_that("mec_zoospread() errors when dist is negative", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, dist = -1))
})

test_that("mec_zoospread() errors when quiet is not logical", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, quiet = "TRUE"))
})

test_that("mec_zoospread() errors when poi does not intersect mec_surfacewater", {
  poi_outside <- st_sf(
    geometry = st_sfc(
      st_point(c(xmin(dem_sr) - 1e6, ymin(dem_sr) - 1e6)),
      crs = st_crs(dem_sr)
    )
  )
  expect_error(mec_zoospread(aoi_sf, poi_outside, mock_surface_lst))
})

# -------------------------------------------------------------------------
# assert_*() exhaustive type checks
# -------------------------------------------------------------------------

## mec_surfacewater structure ----------------------------------------------

test_that("mec_surfacewater structure: character (not list) errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, "result"))
})

test_that("mec_surfacewater structure: numeric (not list) errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, 42))
})

test_that("mec_surfacewater structure: NULL (not list) errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, NULL))
})

test_that("mec_surfacewater structure: list with NULL $mec_surfacewater errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, list(mec_surfacewater = NULL, surface_water = poi_sf)))
})

test_that("mec_surfacewater structure: list with data.frame $surface_water errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, list(mec_surfacewater = dem_sr, surface_water = data.frame(x = 1))))
})

test_that("mec_surfacewater structure: list with unrecognised keys errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, list(raster = dem_sr, vector = poi_sf)))
})

## assert_sf(aoi, geometry = c("POLYGON", "MULTIPOLYGON")) ----------------

test_that("assert_sf: aoi as numeric errors", {
  expect_error(mec_zoospread(1:4, poi_sf, mock_surface_lst))
})

test_that("assert_sf: aoi as character errors", {
  expect_error(mec_zoospread("aoi", poi_sf, mock_surface_lst))
})

test_that("assert_sf: aoi as SpatRaster errors", {
  expect_error(mec_zoospread(dem_sr, poi_sf, mock_surface_lst))
})

test_that("assert_sf: aoi as NULL errors", {
  expect_error(mec_zoospread(NULL, poi_sf, mock_surface_lst))
})

test_that("assert_sf: aoi as sf POINT errors", {
  expect_error(mec_zoospread(poi_sf, poi_sf, mock_surface_lst))
})

test_that("assert_sf: aoi as sf LINESTRING errors", {
  aoi_line <- st_sf(geometry = st_sfc(
    st_linestring(matrix(c(0, 0, 1, 1), ncol = 2)),
    crs = st_crs(aoi_sf)
  ))
  expect_error(mec_zoospread(aoi_line, poi_sf, mock_surface_lst))
})

## assert_sf_point(poi) ----------------------------------------------------

test_that("assert_sf_point: poi as numeric errors", {
  expect_error(mec_zoospread(aoi_sf, c(1, 2), mock_surface_lst))
})

test_that("assert_sf_point: poi as character errors", {
  expect_error(mec_zoospread(aoi_sf, "poi", mock_surface_lst))
})

test_that("assert_sf_point: poi as NULL errors", {
  expect_error(mec_zoospread(aoi_sf, NULL, mock_surface_lst))
})

test_that("assert_sf_point: poi as SpatRaster errors", {
  expect_error(mec_zoospread(aoi_sf, dem_sr, mock_surface_lst))
})

test_that("assert_sf_point: poi as sf POLYGON errors", {
  poi_poly <- st_buffer(poi_sf, dist = 100)
  expect_error(mec_zoospread(aoi_sf, poi_poly, mock_surface_lst))
})

test_that("assert_sf_point: poi as sf MULTIPOINT errors", {
  poi_mpt <- st_sf(geometry = st_sfc(
    st_multipoint(matrix(c(0, 0, 1, 1), ncol = 2)),
    crs = st_crs(poi_sf)
  ))
  expect_error(mec_zoospread(aoi_sf, poi_mpt, mock_surface_lst))
})

## assert_integer_scalar(n_animals) ----------------------------------------

test_that("assert_integer_scalar: n_animals as character errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, n_animals = "5"))
})

test_that("assert_integer_scalar: n_animals as non-integer float errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, n_animals = 1.5))
})

test_that("assert_integer_scalar: n_animals as logical errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, n_animals = TRUE))
})

test_that("assert_integer_scalar: n_animals as NULL errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, n_animals = NULL))
})

## assert_integer_scalar(n_steps) ------------------------------------------

test_that("assert_integer_scalar: n_steps as character errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, n_steps = "10"))
})

test_that("assert_integer_scalar: n_steps as non-integer float errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, n_steps = 2.5))
})

test_that("assert_integer_scalar: n_steps as NULL errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, n_steps = NULL))
})

## assert_positive_numeric(pixel_size) -------------------------------------

test_that("assert_positive_numeric: pixel_size as character errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, pixel_size = "1"))
})

test_that("assert_positive_numeric: pixel_size as logical errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, pixel_size = TRUE))
})

test_that("assert_positive_numeric: pixel_size as NULL errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, pixel_size = NULL))
})

test_that("assert_positive_numeric: pixel_size as NA errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, pixel_size = NA_real_))
})

## assert_integer_scalar(n_iter) -------------------------------------------

test_that("assert_integer_scalar: n_iter as character errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, n_iter = "10"))
})

test_that("assert_integer_scalar: n_iter as non-integer float errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, n_iter = 1.5))
})

test_that("assert_integer_scalar: n_iter as NULL errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, n_iter = NULL))
})

## assert_positive_numeric(dist) -------------------------------------------

test_that("assert_positive_numeric: dist as character errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, dist = "5"))
})

test_that("assert_positive_numeric: dist as logical errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, dist = TRUE))
})

test_that("assert_positive_numeric: dist as NULL errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, dist = NULL))
})

test_that("assert_positive_numeric: dist as NA errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, dist = NA_real_))
})

## assert_logic(quiet) -----------------------------------------------------

test_that("assert_logic: quiet as character errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, quiet = "FALSE"))
})

test_that("assert_logic: quiet as numeric errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, quiet = 0))
})

test_that("assert_logic: quiet as NULL errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, quiet = NULL))
})

test_that("assert_logic: quiet as NA errors", {
  expect_error(mec_zoospread(aoi_sf, poi_sf, mock_surface_lst, quiet = NA))
})

# -------------------------------------------------------------------------
# Output structure  (n_animals=2, n_steps=3, n_iter=2 keeps tests fast)
# dist = Inf guarantees all trajectories pass the distance filter
# -------------------------------------------------------------------------

test_that("mec_zoospread() returns a SpatRaster with one named layer", {
  result <- mec_zoospread(
    aoi_sf, poi_sf, mock_surface_lst,
    n_animals = 2, n_steps = 3, n_iter = 2,
    dist = Inf, quiet = TRUE
  )
  expect_s4_class(result, "SpatRaster")
  expect_equal(nlyr(result), 1L)
  expect_equal(names(result), "mec_zoospread")
})

test_that("mec_zoospread() output CRS matches mec_surfacewater raster", {
  result <- mec_zoospread(
    aoi_sf, poi_sf, mock_surface_lst,
    n_animals = 2, n_steps = 3, n_iter = 2,
    dist = Inf, quiet = TRUE
  )
  expect_true(same.crs(result, mock_raster))
})

test_that("mec_zoospread() output extent matches mec_surfacewater raster", {
  result <- mec_zoospread(
    aoi_sf, poi_sf, mock_surface_lst,
    n_animals = 2, n_steps = 3, n_iter = 2,
    dist = Inf, quiet = TRUE
  )
  expect_equal(as.vector(ext(result)), as.vector(ext(mock_raster)))
})

test_that("mec_zoospread layer values are binary (0/1)", {
  result <- mec_zoospread(
    aoi_sf, poi_sf, mock_surface_lst,
    n_animals = 2, n_steps = 3, n_iter = 2,
    dist = Inf, quiet = TRUE
  )
  vals <- values(result, na.rm = TRUE)
  expect_true(all(vals %in% c(0, 1)))
})

# -------------------------------------------------------------------------
# Behavioural
# -------------------------------------------------------------------------

# dist = 0 → filter `trajectory$dist < 0` is always FALSE → NULL returned
test_that("mec_zoospread() returns NULL when no trajectories are within dist", {
  result <- mec_zoospread(
    aoi_sf, poi_sf, mock_surface_lst,
    n_animals = 2, n_steps = 3, n_iter = 2,
    dist = 0, quiet = TRUE
  )
  expect_null(result)
})

test_that("mec_zoospread() runs without error when quiet = TRUE", {
  expect_no_error(
    mec_zoospread(
      aoi_sf, poi_sf, mock_surface_lst,
      n_animals = 2, n_steps = 3, n_iter = 2,
      dist = Inf, quiet = TRUE
    )
  )
})
