

#' Water inoculum movement
#'
#' Simulates inoculum transport via water in the soil matrix
#'
#' @param dem a single-band \code{SpatRaster} with a digital elevation model
#' @param poi a single-point \code{sf} object denoting the point of interest
#' to run the simulations
#' @param th threshold of flow accumulation to delineate streams
#' @param quiet if \code{TRUE}, suppress any message or progress bar
#'
#' @returns A \code{SpatRaster}
#' @export
#'
#' @details
#'
#' This function models inoculum movement in soil, considering flow direction
#' and altitude under wet and intermediate moisture conditions. It identifies
#' pixels where flow direction matches or is adjacent to the foci, then processes
#' connected pixels based on altitude differences. The result is a binary raster
#' showing the spatial pattern of inoculum dispersal, constrained by topography
#' and moisture.
#'
#' @references
#'
#' Ristaino, J., Gumpertz, M., 2000. New Frontiers in the Study of Dispersal and Spatial Analysis of Epidemics Caused by Species in the Genus Phytophthora. Annu. Rev. Phytopathol. 38, 541–576. \doi{10.1146/annurev.phyto.38.1.541}
#'
#' Vannini, A., Natili, G., Anselmi, N., Montaghi, A., Vettraino, A.M., 2010. Distribution and gradient analysis of Ink disease in chestnut forests. For. Pathol. 40, 73–86. \doi{10.1111/j.1439-0329.2009.00609.x}
#'
#' Vannini, A., Natili, G., Thomidis, T., Belli, C., Morales-Rodriguez, C., 2021. Anthropogenic and landscape features are associated with ink disease impact in Central Italy. For. Pathol. 51, e12722. \doi{10.1111/efp.12722}
#'
#'
#' @examples
#' ## load packages
#' library(sf)
#' library(terra)
#'
#' ## load data
#' dem_sr <- rast(system.file("spatial/dem_light.tiff", package = "phytorisk"))
#' poi_sf <- st_read(system.file("spatial/poi.geojson", package = "phytorisk"))
#'
#' ## simulate mechanism
#' mec_soilwater_sr <- mec_soilwater(dem_sr, poi_sf)
mec_soilwater <- function(dem, poi, th = 100, quiet = FALSE) {

  ## 0. Check for errors ----------------------
  if (!inherits(dem, "SpatRaster")) cli::cli_abort("`dem` must be a SpatRaster")
  if (!inherits(poi, "sf") | nrow(poi) != 1 | sf::st_geometry_type(poi) != "POINT")
    cli::cli_abort("`poi` must be an sf object with one POINT")

  if (!terra::same.crs(dem, poi)) cli::cli_abort("CRS of inputs is not the same")

  if (!requireNamespace("flowdem")) {
    cli::cli_abort("Package `flowdem` not available. Please install it from {.url https://github.com/KennethTM/flowdem}")
  }


  ## 1. Flow direction ------------------------

  ## Determine the flow direction, flow accumulation and drainage network
  ## https://github.com/KennethTM/flowdem/blob/main/README.md

  ## ensure dem has correct name
  names(dem) <- "dem"

  ## fill DEM and leave surfaces flat
  if (!quiet) cli::cli_progress_step("Filling DEM...", "DEM filled", "Filling DEM failed")
  dem_fill_sr <- flowdem::fill(dem, epsilon = FALSE)

  ## fill DEM and apply a gradient on flat surfaces to ensure flow
  ## NOTE: when writing DEMs filled with epsilon to file,
  ## set the datatype to flt8s
  dem_fill_eps_sr <- flowdem::fill(dem, epsilon = TRUE)

  ## fill DEM and delineate coastal drainage basins simultaneously
  if (!quiet) cli::cli_progress_step("Filling basins...", "Basins filled", "Filling basins failed")
  dem_fill_basins_sr <- flowdem::fill_basins(dem)

  ## breach DEM: resolve depression by "carving" through obstacles
  if (!quiet) cli::cli_progress_step("Removing depressions...", "Depressions removed", "Removing depressions failed")
  dem_breach_sr <- flowdem::breach(dem)

  ## use fill with epsilon on breached DEM to resolve flats and ensure drainage
  if (!quiet) cli::cli_progress_step("Filling depressions...", "Depressions filled", "Filling depressions failed")
  dem_breach_fill_eps_sr <- flowdem::fill(dem_breach_sr, epsilon = TRUE)

  ## get flow directions using the filled DEM
  if (!quiet) cli::cli_progress_step("Getting flow directions...", "Flow directions", "Flow directions failed")
  dem_dir_sr <- flowdem::dirs(dem_breach_fill_eps_sr, mode = "d8")

  ## get flow accumulation
  if (!quiet) cli::cli_progress_step("Calculating flow accumulation...", "Flow accumulation calculated", "Flow accumulation failed")
  dem_acc_sr <- flowdem::accum(dem_dir_sr, mode = "d8")

  ## flow accumulation can be used for stream delineation,
  ## e.g using a threshold of 100 contributing cells (100*100*100 = 1 km2)
  if (!quiet) cli::cli_progress_step("Delineating streams...", "Streams delineated", "Delineating streams failed")
  dem_streams_sr <- terra::ifel(dem_acc_sr > th, 1, 0)

  # 3.2. Wet front ------------------------

  ## determine the wet front from the Pc outbreak focus
  if (!quiet) cli::cli_progress_step("Determining the wet front", "Wet front determined", "Wet front determination failed")

  ## extract the flow direction in the POI
  flowdir_df <- terra::extract(dem_dir_sr, poi)
  flowdir_poi <- flowdir_df$dem

  ## create a binary raster where the values match with the flow direction of the POI
  flowdir_poi_sr <- terra::ifel(
    (dem_dir_sr == (flowdir_poi - 1)) |
      (dem_dir_sr == flowdir_poi) |
      (dem_dir_sr == (flowdir_poi + 1)),
    1, 0
  )

  flowdir_poi_sr <- terra::ifel(is.na(dem_dir_sr), NA, flowdir_poi_sr)

  ## initialization parameters
  init_cell <- terra::cellFromXY(flowdir_poi_sr, sf::st_coordinates(poi))
  init_i    <- terra::rowFromCell(flowdir_poi_sr, init_cell)
  init_j    <- terra::colFromCell(flowdir_poi_sr, init_cell)

  ## isolate the pixels that fulfil the condition
  risk_mec_sr <- isolate_pixels(flowdir_poi_sr, dem, init_i, init_j)
  terra::crs(risk_mec_sr) <- terra::crs(dem)

  ## append streams for output
  risk_mec_sr <- c(risk_mec_sr, dem_streams_sr)
  names(risk_mec_sr) <- c("mec_soilwater", "streams")

  ## return results
  return(risk_mec_sr)
}
