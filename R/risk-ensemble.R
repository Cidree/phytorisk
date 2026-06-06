

#' Phytophthora ensemble risk
#'
#' Calculates Phytophthora cinnamomi's ensemble risk from mechanisms
#'
#' @param mec_soilwater The result of \link{mec_soilwater}
#' @param mec_rootcontact The result of \link{mec_rootcontact}
#' @param mec_surfacewater The result of \link{mec_surfacewater}
#' @param mec_zoospread Optional module. The result of \link{mec_zoospread}
#' @param weights Weights of the ensemble model. The default uses the same
#' weights for each model. The argument accepts a numeric vector with the
#' corresponding weights. See *Details*
#'
#' @returns A \code{SpatRaster}
#' @export
#'
#' @details
#'
#' This function calculates the ensemble risk as the weighted sum of 
#' the individual risk models as follows:
#'
#' \eqn{Risk = w_1 \cdot Mec_{soilwater} + w_2 \cdot Mec_{rootcontact} + w_3 \cdot Mec_{mec_surfacewater} + w_4  \cdot Mec_{zoospread}}.
#'
#' Being \eqn{w_i = 0.25} when \link{mec_zoospread} is not null, and 
#' \eqn{w_i = \frac{1}{3}} when \link{mec_zoospread} is null.
#' Additionally, the user can specify a numeric vector of length 
#' 3 or 4 specifying the desired weights based on expert criteria.
#'
#' @examples
#' \donttest{
#' ## load packages
#' library(phytorisk)
#' library(sf)
#' library(terra)
#' 
#' ## load data
#' poi_sf <- st_read(
#'   system.file("spatial/poi.geojson", package = "phytorisk"),
#'   quiet = TRUE
#' )
#' dem_sr <- rast(system.file("spatial/dem_light.tiff", package = "phytorisk"))
#' trees_sr <- rast(system.file("spatial/trees_light.tiff", package = "phytorisk"))
#' aoi_sf <- st_read(
#'   system.file("spatial/tejera.geojson", package = "phytorisk"),
#'   quiet = TRUE
#' )
#' 
#' ## first, calculate the individual mechanisms
#' mec_soilwater_sr <- mec_soilwater(dem_sr, poi_sf)
#' mec_surface_lst   <- mec_surfacewater(dem_sr, mec_soilwater_sr, poi_sf)
#' mec_rootcontact_sr <- mec_rootcontact(trees_sr, aoi_sf, poi_sf)
#' 
#' ## calculate ensemble risk using equal weights
#' risk_equal_sr <- phytorisk_ensemble(
#'  mec_soilwater = mec_soilwater_sr,
#'  mec_rootcontact = mec_rootcontact_sr,
#'  mec_surfacewater = mec_surface_lst
#')
#' 
#' ## assign more weight to root-to-root contact
#' risk_weighted_sr <- phytorisk_ensemble(
#'  mec_soilwater = mec_soilwater_sr,
#'  mec_rootcontact = mec_rootcontact_sr,
#'  mec_surfacewater = mec_surface_lst,
#'  weights = c(.3, .4, .3)
#')
#' }
phytorisk_ensemble <- function(
  mec_soilwater,
  mec_rootcontact,
  mec_surfacewater,
  mec_zoospread = NULL,
  weights = "equal"
) {

  # 0. Validate inputs
  ## Abort if it's not a list
  if (typeof(mec_surfacewater) != "list")
    cli::cli_abort("{.arg mec_surfacewater} argument must be the output of {.fun mec_surfacewater}")

  ## Abort if it does not contain the expected inputs
  if (!inherits(mec_surfacewater$mec_surfacewater, "SpatRaster") || !inherits(mec_surfacewater$surface_water, "sf")) 
    cli::cli_abort("{.arg mec_surfacewater} argument must be the output of {.fun mec_surfacewater}")
  
  assert_spatraster(mec_soilwater, "mec_soilwater")
  assert_spatraster(mec_rootcontact, "mec_rootcontact")
  if (!is.null(mec_zoospread)) assert_spatraster(mec_zoospread, "mec_zoospread")

  layers <- c(mec_soilwater$mec_soilwater, mec_rootcontact, mec_surfacewater$mec_surfacewater, mec_zoospread)
  n <- terra::nlyr(layers)


  # 1. Validate weights
  if (identical(weights, "equal")) {
    weights <- rep(1 / n, n)
  } else {
    if (!is.numeric(weights))  cli::cli_abort("{.arg weights} must be numeric or {.val equal}.")
    if (length(weights) != n)  cli::cli_abort("{.arg weights} must be a numeric vector of length {n}.")
    if (sum(weights) != 1)     cli::cli_abort("{.arg weights} must sum to 1.")
  }


  # 2. Weighted sum
  risk_rast <- sum(layers * weights)
  names(risk_rast) <- "ensemble_risk"
  return(risk_rast)
}






#' Raw Phytophthora ensemble risk
#'
#' Calculates Phytophthora cinnamomi's ensemble risk from raw data
#'
#' @template aoi
#' @template poi
#' @template dem
#' @template treecover
#' @param weights weights of the ensemble model. The default uses the same
#' weights for each model. The argument accepts a numeric vector with the
#' corresponding weights. See [phytorisk_ensemble]
#' @template th
#' @template buffer
#' @param include_zoospread logical. Whether to include the optional module of [mec_zoospread]
#' @param append_mec Logical. Whether to append to the results of each
#' individual module to the output SpatRaster
#' @param ... arguments passed to [mec_zoospread]
#' @template quiet
#'
#' @returns A \code{SpatRaster}
#' @export
#'
#' @details
#'
#' This function is a one-step convenience wrapper around the four dispersal
#' mechanism functions and [phytorisk_ensemble]. It runs the following pipeline
#' internally:
#'
#' \enumerate{
#'   \item \link{mec_soilwater} — models *Pc* spread through soil water pathways
#'     using flow direction and accumulation derived from \code{dem}. The
#'     \code{th} argument controls the flow-accumulation threshold used to
#'     delineate streams.
#'   \item \link{mec_rootcontact} — models root-to-root transmission across a
#'     3×3 spatial window over the \code{treecover} raster within \code{aoi}.
#'   \item \link{mec_surfacewater} — models surface water spread using the
#'     Topographic Wetness Index derived from \code{dem}. It receives the output
#'     of \link{mec_soilwater} directly; \code{buffer} extends the detected
#'     water bodies before computing risk.
#'   \item \link{mec_zoospread} (optional) — simulates animal-mediated dispersal
#'     trajectories within \code{aoi}. Only executed when
#'     \code{include_zoospread = TRUE}; additional arguments via \code{...} are
#'     forwarded to this function.
#' }
#'
#' The four mechanism outputs are then combined by [phytorisk_ensemble] into a
#' single ensemble risk surface. See [phytorisk_ensemble] for details on how
#' \code{weights} are applied.
#'
#' When \code{append_mec = TRUE} the individual mechanism rasters are
#' concatenated with the ensemble layer in the returned \code{SpatRaster},
#' allowing inspection of each component alongside the final risk surface.
#'
#' @examples
#' \donttest{
#' ## load packages
#' library(phytorisk)
#' library(sf)
#' library(terra)
#' 
#' ## load data
#' poi_sf <- st_read(
#'   system.file("spatial/poi.geojson", package = "phytorisk"),
#'   quiet = TRUE
#' )
#' dem_sr <- rast(system.file("spatial/dem_light.tiff", package = "phytorisk"))
#' trees_sr <- rast(system.file("spatial/trees_light.tiff", package = "phytorisk"))
#' aoi_sf <- st_read(
#'   system.file("spatial/tejera.geojson", package = "phytorisk"),
#'   quiet = TRUE
#' )
#' 
#' ## calculate ensemble risk, returning individual mechanisms
#' risk_equal_sr <- phytorisk_ensemble_raw(
#'   aoi = aoi_sf,
#'   poi = poi_sf,
#'   dem = dem_sr,
#'   treecover = trees_sr,
#'   append_mec = TRUE
#' )
#' 
#' ## visualize results
#' plot(risk_equal_sr)
#' }
phytorisk_ensemble_raw <- function(
  aoi,
  poi,
  dem,
  treecover,
  weights = "equal",
  th = 100,
  buffer = 50,
  include_zoospread = FALSE,
  append_mec = FALSE,
  ...,
  quiet = FALSE) {
  
  # 0. Validate speciific inputs (rest are handled downstream)
  assert_logic(include_zoospread, "include_zoospread")
  assert_logic(append_mec, "append_mec")

  # 1. Calculate individual components

  if (!quiet) cli::cli_h1("Mec Ii - Spread in soil water")
  mec_soilwater <- mec_soilwater(dem, poi, th = th, quiet = quiet)

  if (!quiet) cli::cli_h1("Mec Iii - Root-to-root contact")
  mec_rootcontact <- mec_rootcontact(treecover, aoi, poi, quiet = quiet)

  if (!quiet) cli::cli_h1("Mec II - Spread in surface water")
  mec_surfacewater <- mec_surfacewater(dem, mec_soilwater, poi, buffer = buffer, quiet = quiet)

  if (include_zoospread) {
    if (!quiet) cli::cli_h1("Mec III - Dispersion by animals")
    mec_zoospread   <- mec_zoospread(aoi, poi, mec_surfacewater, ..., quiet = quiet)
  } else {
    mec_zoospread <- NULL
  }

  # 2. Calculate ensemble risk
  final_risk <- phytorisk_ensemble(
    mec_soilwater    = mec_soilwater, 
    mec_rootcontact  = mec_rootcontact, 
    mec_surfacewater = mec_surfacewater, 
    mec_zoospread    = mec_zoospread, 
    weights          = weights
  )

  # 3. Return results based on append_mec
  if (append_mec) {
    c(final_risk, mec_soilwater, mec_rootcontact, mec_zoospread, mec_zoospread)
  } else {
    final_risk
  }

}
