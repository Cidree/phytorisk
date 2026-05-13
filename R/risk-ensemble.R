

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
#' ## TODO
phytorisk_ensemble <- function(
  mec_soilwater,
  mec_rootcontact,
  mec_surfacewater,
  mec_zoospread = NULL,
  weights = "equal"
) {

  # 0. Validate inputs
  assert_spatraster(mec_soilwater, "mec_soilwater")
  assert_spatraster(mec_rootcontact, "mec_rootcontact")
  assert_spatraster(mec_surfacewater, "mec_surfacewater")
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
  return(sum(layers * weights))
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
#' @param ... arguments passed to [mec_zoospread]
#' @template quiet
#'
#' @returns A \code{SpatRaster}
#' @export
#'
#' @details
#'
#' #TODO
#'
phytorisk_ensemble_raw <- function(
  aoi,
  poi,
  dem,
  treecover,
  weights = "equal",
  th = 100,
  buffer = 50,
  include_zoospread = FALSE,
  ...,
  quiet = FALSE) {

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

  final_risk <- phytorisk_ensemble(
    mec_soilwater    = mec_soilwater, 
    mec_rootcontact  = mec_rootcontact, 
    mec_surfacewater = mec_surfacewater, 
    mec_zoospread    = mec_zoospread, 
    weights          = weights
  )

  return(final_risk)

}
