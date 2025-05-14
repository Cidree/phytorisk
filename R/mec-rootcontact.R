

#' Root-to-root contact
#'
#' Simulates diffusive inoculum spread through root contact
#'
#' @param treecover a single-band \code{SpatRaster} where 1 represents host
#' trees, and 0 represents background area
#' @param aoi a \code{sf} polygon representing the area of interest. Used to mask
#' the tree cover
#' @param poi a single-point \code{sf} object denoting the point of interest
#' to run the simulations
#' @param quiet if \code{TRUE}, suppress any message or progress bar
#'
#' @returns A \code{SpatRaster}
#' @export
#'
#' @details
#'
#' This function models inoculum movement due to root-to-root contact based on
#' tree continuity. It calculates vegetation continuity by counting connected
#' pixels within a 3×3 window. If continuity exists, a new raster is created with
#' connected pixels marked as 1. The function then identifies all pixels connected
#' to the foci within the study area. If no direct pixel-to-pixel contact between
#' trees exists, the infection risk is considered zero.
#'
#' @references
#'
#' Cardillo, E., Abad, E., Meyer, S., 2021. Iberian oak decline caused by Phytophthora cinnamomi: A spatiotemporal analysis incorporating the effect of host heterogeneities at landscape scale. For. Pathol. 51, e12667. \doi{10.1111/efp.12667}
#'
#' Cardillo, E., Acedo, A., Abad, E., 2018. Topographic effects on dispersal patterns of Phytophthora cinnamomi at a stand scale in a Spanish heathland. PloS One 13, e0195060.
#'
#' @examples
#' \donttest{
#' ## load packages
#' library(sf)
#' library(terra)
#'
#' ## load data
#' study_area_sf <- st_read(system.file("spatial/tejera.geojson", package = "phytorisk"))
#' poi_sf <- st_read(system.file("spatial/poi.geojson", package = "phytorisk"))
#' trees_sr <- rast(system.file("spatial/trees_light.tiff", package = "phytorisk"))
#'
#' ## simulate mechanism
#' mec_rootcontact_sr <- mec_rootcontact(trees_sr, study_area_sf, poi_sf)
#' }
mec_rootcontact <- function(treecover, aoi, poi, quiet = FALSE) {

  ## 0. Check for errors ----------------------
  if (!terra::same.crs(treecover, aoi) | !terra::same.crs(aoi, poi)) cli::cli_abort("CRS of inputs is not the same")

  ## 1. Tree cover ------------------

  ## mask tree cover to the study area
  if (!quiet) cli::cli_progress_step("Preparing tree data...", "Tree data prepared", "Tree data preparation failed")
  trees_sr <- terra::crop(treecover, aoi, mask = TRUE)

  ## 2. Spatial continuity -----------

  ## apply a buffer matrix, assuming that the roots extend beyond the crown projection
  ## TODO -> result slightly different (Antonio uses NA values)
  edges_sr <- terra::focal(trees_sr, w = matrix(1, 3, 3), fun = function(x) sum(x) - x[5])
  edges_sr[edges_sr > 0] <- 1

  ## convert coordinates to matrix indexes
  poi_cell_id <- terra::cellFromXY(edges_sr, sf::st_coordinates(poi))
  poi_index   <- terra::rowColFromCell(edges_sr, poi_cell_id)

  ## find pixels connected to the POI
  if (!quiet) cli::cli_progress_step("Finding root-to-root contact...", "Finished", "Finding root-to-root contact failed")
  conneted_list <- find_connected(edges_sr, poi_index[1], poi_index[2])

  ## create a new binary raster with the connected pixels marked as 1
  risk_mec_sr <- edges_sr
  risk_mec_sr[] <- 0
  for (pixel in conneted_list) {
    risk_mec_sr[pixel[1], pixel[2]] <- 1
  }

  ## return results
  names(risk_mec_sr) <- "mec_rootcontact"
  return(risk_mec_sr)


}
