

#' Inoculum spread
#'
#' Simulates the spread of the inoculum in surface water flow
#'
#' @param dem a single-band \code{SpatRaster} with a digital elevation model
#' @param mec_soilwater the result of \link{mec_soilwater}
#' @param poi a single-point \code{sf} object denoting the point of interest
#' to run the simulations
#' @param buffer a buffer in meters to extent the spread in every direction
#' @param quiet if \code{TRUE}, suppress any message or progress bar
#'
#' @returns A \code{SpatRaster}
#' @export
#'
#' @details
#' To do...
#'
#'
#' @examples
#' ## load packages
#' library(sf)
#' library(terra)
#'
#' ## load data
#' poi_sf <- st_read(system.file("spatial/poi.geojson", package = "phytorisk"))
#' dem_sr <- rast(system.file("spatial/dem_light.tiff", package = "phytorisk"))
#' trees_sr <- rast(system.file("spatial/trees_light.tiff", package = "phytorisk"))
#'
#' ## simulate mechanism
#' mec_soilwater_sr <- mec_soilwater(dem_sr, poi_sf)
#' mec_surface_sr <- mec_surfaceflow(dem_sr, mec_soilwater_sr, poi_sf)
mec_surfaceflow <- function(dem, mec_soilwater, poi, buffer = 50, quiet = FALSE) {

  ## 0. Check for errors ----------------------
  if (!terra::same.crs(dem, mec_soilwater) | !terra::same.crs(mec_soilwater, poi)) cli::cli_abort("CRS of inputs is not the same")
  if (all(names(mec_soilwater) != c("mec_soilwater", "streams"))) cli::cli_abort("`mec_soilwater` must be the result of the {.topic phytorisk::mec_soilwater}")

  ## 1. Natural drainage network -----------------

  ## calculate the flow accumulation using a simple approximation
  if (!quiet) cli::cli_progress_step("Calculating natural drainage network...", "Natural drainage network calculated", "Natural drainage network failed")
  twi_sr <- twi(dem)

  ## union of drainage network and moist areas
  surface_water_sr <- mec_soilwater$streams + twi_sr
  surface_water_sr[surface_water_sr == 2] <- 1

  ## convert to polygons
  surface_water_vect <- terra::as.polygons(surface_water_sr)
  surface_water_vect <- surface_water_vect[surface_water_vect$streams == 1, ]

  ## 2. Surface water close to foci ------------

  ## apply limits mask to the streams raster
  ## condition: the soil humidity spot must get to a surface stream or flooding area
  if (!quiet) cli::cli_progress_step("Identifying surface water close to foci...", "Surface water close to foci identified", "Surface water close to foci failed")
  streams_masked_sr <- terra::mask(mec_soilwater$streams, mec_soilwater$mec_soilwater, maskvalues = 0)
  values <- values(streams_masked_sr)

  ## get the coordinates of pixels with value equal to 1
  coords_df <- terra::as.data.frame(streams_masked_sr, xy = TRUE)
  coords_df <- coords_df[coords_df$streams == 1, ]

  ## convert to sf
  coords_sf <- sf::st_as_sf(
    x      = coords_df,
    coords = c("x", "y"),
    crs    = terra::crs(streams_masked_sr)
  )

  ## calculate the distances from the coordinate of interest to each pixel with value 1
  distances_vec <- as.numeric(sf::st_distance(poi, coords_sf))

  ## get the coordinates of the closest pixel
  final_cell_sf <- coords_sf[which.min(distances_vec), ]
  final_cell_coord <- sf::st_coordinates(final_cell_sf)

  ## turn the coordinates into matrix indexes
  init_cell <- terra::cellFromXY(surface_water_sr, final_cell_coord)
  init_i    <- terra::rowFromCell(surface_water_sr, init_cell)
  init_j    <- terra::colFromCell(surface_water_sr, init_cell)

  ## find the pixels connected to the POI
  if (!quiet) cli::cli_progress_step("Finding connected pixels...", "Finished", "Finding connected pixels failed")
  connected_lst <- find_connected(surface_water_sr, init_i, init_j)

  ## create a new binary raster with the connected pixels marked as 1
  risk_mec_sr <- surface_water_sr
  risk_mec_sr[] <- 0
  for (pixel in connected_lst) {
    risk_mec_sr[pixel[1], pixel[2]] <- 1
  }

  ## identify the pixels with value 1 in the risk raster
  pixels_vec <- which(terra::values(risk_mec_sr, mat = FALSE) == 1)

  ## get the coordinates of pixels with value 1
  pixel_coords_mat <- terra::xyFromCell(risk_mec_sr, pixels_vec)

  ## get the heights of the DEM
  heights_vec <- terra::extract(dem, pixel_coords_mat)[, 1]

  ## find the coordinate of the pixel with value 1 and lowest height
  min_height_coord <- pixel_coords_mat[which.min(heights_vec), ]

  ## get raster extent
  risk_mec_ext <- terra::ext(risk_mec_sr)

  ## calculate the limits of the quadrats with final_cell_coord as central point
  quad1_ext <- terra::ext(risk_mec_ext[1], final_cell_coord[1], risk_mec_ext[3], final_cell_coord[2])
  quad2_ext <- terra::ext(final_cell_coord[1], risk_mec_ext[2], risk_mec_ext[3], final_cell_coord[2])
  quad3_ext <- terra::ext(risk_mec_ext[1], final_cell_coord[1], final_cell_coord[2], risk_mec_ext[4])
  quad4_ext <- terra::ext(final_cell_coord[1], risk_mec_ext[2], final_cell_coord[2], risk_mec_ext[4])

  ## determine in which quadrat is the coordinate of minimum height
  if (min_height_coord[1] <= final_cell_coord[1] && min_height_coord[2] <= final_cell_coord[2]) {
    quadrat <- quad1_ext
  } else if (min_height_coord[1] > final_cell_coord[1] && min_height_coord[2] <= final_cell_coord[2]) {
    quadrat <- quad2_ext
  } else if (min_height_coord[1] <= final_cell_coord[1] && min_height_coord[2] > final_cell_coord[2]) {
    quadrat <- quad3_ext
  } else {
    quadrat <- quad4_ext
  }

  ## extend the coordinates in every direction according to a buffer
  quadrat_bbox <- sf::st_bbox(quadrat) + c(-buffer, -buffer, buffer, buffer)

  ## crop the raster with the quadrat of the minimum's height coordinate
  risk_mec_sr <- terra::crop(risk_mec_sr, quadrat_bbox)

  ## return result
  names(risk_mec_sr) <- "mec_surfaceflow"
  return(risk_mec_sr)

}
