library(sf)
library(terra)

# ── Shared fixtures ───────────────────────────────────────────────────────────

dem_sr   <- rast(system.file("spatial/dem_light.tiff",   package = "phytorisk"))
trees_sr <- rast(system.file("spatial/trees_light.tiff", package = "phytorisk"))
poi_sf   <- st_read(system.file("spatial/poi.geojson",    package = "phytorisk"), quiet = TRUE)
aoi_sf   <- st_read(system.file("spatial/tejera.geojson", package = "phytorisk"), quiet = TRUE)

# ── Synthetic mechanism outputs (no flowdem required) ────────────────────────
# A 5x5 raster used to build fake mechanism results without running real pipelines
r_base <- rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5, crs = "EPSG:25829")
set.seed(42)
values(r_base) <- runif(25)

# mec_soilwater output: SpatRaster with layers "mec_soilwater" + "streams"
soilwater_sr        <- c(r_base, r_base)
names(soilwater_sr) <- c("mec_soilwater", "streams")

# mec_rootcontact output: single-layer SpatRaster
rootcontact_sr        <- r_base
names(rootcontact_sr) <- "mec_rootcontact"

# mec_surfacewater output: list(mec_surfacewater = SpatRaster, surface_water = sf)
sw_rast          <- r_base
names(sw_rast)   <- "mec_surfacewater"
surfacewater_lst <- list(
  mec_surfacewater = sw_rast,
  surface_water    = st_sf(geometry = st_sfc(st_point(c(2.5, 2.5)), crs = "EPSG:25829"))
)

# mec_zoospread output: single-layer SpatRaster
zoospread_sr        <- r_base
names(zoospread_sr) <- "mec_zoospread"

# ── phytorisk_ensemble() ──────────────────────────────────────────────────────

# ─── Input validation ─────────────────────────────────────────────────────────

test_that("phytorisk_ensemble() errors when mec_surfacewater is not a list", {
  expect_error(
    phytorisk_ensemble(soilwater_sr, rootcontact_sr, sw_rast),
    regexp = "mec_surfacewater"
  )
})

test_that("phytorisk_ensemble() errors when mec_surfacewater list lacks required elements", {
  bad_lst <- list(foo = sw_rast, bar = 1L)
  expect_error(
    phytorisk_ensemble(soilwater_sr, rootcontact_sr, bad_lst),
    regexp = "mec_surfacewater"
  )
})

test_that("phytorisk_ensemble() errors when mec_surfacewater$mec_surfacewater is not a SpatRaster", {
  bad_lst <- list(
    mec_surfacewater = matrix(1:25, 5, 5),
    surface_water    = st_sf(geometry = st_sfc(st_point(c(2.5, 2.5)), crs = "EPSG:25829"))
  )
  expect_error(
    phytorisk_ensemble(soilwater_sr, rootcontact_sr, bad_lst),
    regexp = "mec_surfacewater"
  )
})

test_that("phytorisk_ensemble() errors when mec_surfacewater$surface_water is not sf", {
  bad_lst <- list(
    mec_surfacewater = sw_rast,
    surface_water    = data.frame(x = 1, y = 2)
  )
  expect_error(
    phytorisk_ensemble(soilwater_sr, rootcontact_sr, bad_lst),
    regexp = "mec_surfacewater"
  )
})

test_that("phytorisk_ensemble() errors when mec_soilwater is not a SpatRaster", {
  expect_error(
    phytorisk_ensemble(matrix(1:25, 5, 5), rootcontact_sr, surfacewater_lst),
    regexp = "mec_soilwater"
  )
})

test_that("phytorisk_ensemble() errors when mec_rootcontact is not a SpatRaster", {
  expect_error(
    phytorisk_ensemble(soilwater_sr, matrix(1:25, 5, 5), surfacewater_lst),
    regexp = "mec_rootcontact"
  )
})

test_that("phytorisk_ensemble() errors when mec_zoospread is provided but not a SpatRaster", {
  expect_error(
    phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst,
                       mec_zoospread = matrix(1:25, 5, 5)),
    regexp = "mec_zoospread"
  )
})

test_that("phytorisk_ensemble() errors when weights is a non-'equal' string", {
  expect_error(
    phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst, weights = "balanced"),
    regexp = "weights"
  )
})

test_that("phytorisk_ensemble() errors when weights has wrong length for 3 mechanisms", {
  expect_error(
    phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst, weights = c(0.5, 0.5)),
    regexp = "weights"
  )
})

test_that("phytorisk_ensemble() errors when weights has wrong length for 4 mechanisms", {
  expect_error(
    phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst,
                       mec_zoospread = zoospread_sr, weights = c(1/3, 1/3, 1/3)),
    regexp = "weights"
  )
})

test_that("phytorisk_ensemble() errors when weights do not sum to 1", {
  expect_error(
    phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst,
                       weights = c(0.4, 0.4, 0.4)),
    regexp = "weights"
  )
})

# ─── assert_*() exhaustive type checks ───────────────────────────────────────

## assert_spatraster(mec_soilwater) -------------------------------------------

test_that("assert_spatraster: mec_soilwater as numeric vector errors", {
  expect_error(phytorisk_ensemble(1:25, rootcontact_sr, surfacewater_lst))
})

test_that("assert_spatraster: mec_soilwater as list errors", {
  expect_error(phytorisk_ensemble(list(values = 1:25), rootcontact_sr, surfacewater_lst))
})

test_that("assert_spatraster: mec_soilwater as NULL errors", {
  expect_error(phytorisk_ensemble(NULL, rootcontact_sr, surfacewater_lst))
})

## assert_spatraster(mec_rootcontact) -----------------------------------------

test_that("assert_spatraster: mec_rootcontact as character errors", {
  expect_error(phytorisk_ensemble(soilwater_sr, "rootcontact", surfacewater_lst))
})

test_that("assert_spatraster: mec_rootcontact as data.frame errors", {
  expect_error(phytorisk_ensemble(soilwater_sr, data.frame(x = 1:5), surfacewater_lst))
})

test_that("assert_spatraster: mec_rootcontact as NULL errors", {
  expect_error(phytorisk_ensemble(soilwater_sr, NULL, surfacewater_lst))
})

## assert_spatraster(mec_zoospread) -------------------------------------------

test_that("assert_spatraster: mec_zoospread as numeric vector errors", {
  expect_error(
    phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst, mec_zoospread = 1:25)
  )
})

test_that("assert_spatraster: mec_zoospread as character errors", {
  expect_error(
    phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst, mec_zoospread = "zoo")
  )
})

# ─── Output structure ────────────────────────────────────────────────────────

test_that("phytorisk_ensemble() returns a SpatRaster", {
  result <- phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst)
  expect_s4_class(result, "SpatRaster")
})

test_that("phytorisk_ensemble() returns exactly one layer named 'ensemble_risk'", {
  result <- phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst)
  expect_equal(nlyr(result), 1L)
  expect_equal(names(result), "ensemble_risk")
})

test_that("phytorisk_ensemble() with zoospread still returns one layer named 'ensemble_risk'", {
  result <- phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst,
                               mec_zoospread = zoospread_sr)
  expect_equal(nlyr(result), 1L)
  expect_equal(names(result), "ensemble_risk")
})

test_that("phytorisk_ensemble() output has same extent as inputs", {
  result <- phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst)
  expect_equal(as.vector(ext(result)), as.vector(ext(r_base)))
})

# ─── Behavioural ─────────────────────────────────────────────────────────────

test_that("equal weights produce a simple mean across three mechanism layers", {
  result   <- phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst)
  layers   <- c(soilwater_sr$mec_soilwater, rootcontact_sr, sw_rast)
  expected <- sum(layers * rep(1 / 3, 3))
  expect_equal(as.numeric(values(result)), as.numeric(values(expected)))
})

test_that("equal weights produce a simple mean across four mechanism layers", {
  result   <- phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst,
                                 mec_zoospread = zoospread_sr)
  layers   <- c(soilwater_sr$mec_soilwater, rootcontact_sr, sw_rast, zoospread_sr)
  expected <- sum(layers * rep(1 / 4, 4))
  expect_equal(as.numeric(values(result)), as.numeric(values(expected)))
})

test_that("custom weights are applied correctly for three mechanisms", {
  w        <- c(0.5, 0.3, 0.2)
  result   <- phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst, weights = w)
  layers   <- c(soilwater_sr$mec_soilwater, rootcontact_sr, sw_rast)
  expected <- sum(layers * w)
  expect_equal(as.numeric(values(result)), as.numeric(values(expected)))
})

test_that("custom weights are applied correctly for four mechanisms", {
  w        <- c(0.4, 0.3, 0.2, 0.1)
  result   <- phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst,
                                 mec_zoospread = zoospread_sr, weights = w)
  layers   <- c(soilwater_sr$mec_soilwater, rootcontact_sr, sw_rast, zoospread_sr)
  expected <- sum(layers * w)
  expect_equal(as.numeric(values(result)), as.numeric(values(expected)))
})

test_that("mec_zoospread = NULL is accepted without error", {
  expect_no_error(
    phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst, mec_zoospread = NULL)
  )
})

test_that("weights = 'equal' is accepted without error", {
  expect_no_error(
    phytorisk_ensemble(soilwater_sr, rootcontact_sr, surfacewater_lst, weights = "equal")
  )
})

# ── phytorisk_ensemble_raw() ──────────────────────────────────────────────────

# ─── Input validation (no flowdem required) ───────────────────────────────────

test_that("phytorisk_ensemble_raw() errors when include_zoospread is not logical", {
  expect_error(
    phytorisk_ensemble_raw(aoi_sf, poi_sf, dem_sr, trees_sr, include_zoospread = "TRUE"),
    regexp = "include_zoospread"
  )
})

test_that("phytorisk_ensemble_raw() errors when include_zoospread is numeric", {
  expect_error(
    phytorisk_ensemble_raw(aoi_sf, poi_sf, dem_sr, trees_sr, include_zoospread = 1L),
    regexp = "include_zoospread"
  )
})

test_that("phytorisk_ensemble_raw() errors when include_zoospread is NULL", {
  expect_error(
    phytorisk_ensemble_raw(aoi_sf, poi_sf, dem_sr, trees_sr, include_zoospread = NULL),
    regexp = "include_zoospread"
  )
})

test_that("phytorisk_ensemble_raw() errors when append_mec is not logical", {
  expect_error(
    phytorisk_ensemble_raw(aoi_sf, poi_sf, dem_sr, trees_sr, append_mec = "TRUE"),
    regexp = "append_mec"
  )
})

test_that("phytorisk_ensemble_raw() errors when append_mec is numeric", {
  expect_error(
    phytorisk_ensemble_raw(aoi_sf, poi_sf, dem_sr, trees_sr, append_mec = 1L),
    regexp = "append_mec"
  )
})

test_that("phytorisk_ensemble_raw() errors when append_mec is NULL", {
  expect_error(
    phytorisk_ensemble_raw(aoi_sf, poi_sf, dem_sr, trees_sr, append_mec = NULL),
    regexp = "append_mec"
  )
})

# ─── Output structure (requires flowdem) ─────────────────────────────────────

test_that("phytorisk_ensemble_raw() returns a single-layer SpatRaster named 'ensemble_risk'", {
  skip_if_not_installed("flowdem")
  result <- phytorisk_ensemble_raw(aoi_sf, poi_sf, dem_sr, trees_sr, quiet = TRUE)
  expect_s4_class(result, "SpatRaster")
  expect_equal(nlyr(result), 1L)
  expect_equal(names(result), "ensemble_risk")
})

test_that("phytorisk_ensemble_raw() with append_mec = TRUE returns more than one layer", {
  skip_if_not_installed("flowdem")
  result <- phytorisk_ensemble_raw(aoi_sf, poi_sf, dem_sr, trees_sr,
                                   append_mec = TRUE, quiet = TRUE)
  expect_s4_class(result, "SpatRaster")
  expect_gt(nlyr(result), 1L)
})

test_that("phytorisk_ensemble_raw() with append_mec = TRUE includes 'ensemble_risk' layer", {
  skip_if_not_installed("flowdem")
  result <- phytorisk_ensemble_raw(aoi_sf, poi_sf, dem_sr, trees_sr,
                                   append_mec = TRUE, quiet = TRUE)
  expect_true("ensemble_risk" %in% names(result))
})

test_that("phytorisk_ensemble_raw() runs without error when quiet = TRUE", {
  skip_if_not_installed("flowdem")
  expect_no_error(
    phytorisk_ensemble_raw(aoi_sf, poi_sf, dem_sr, trees_sr, quiet = TRUE)
  )
})
