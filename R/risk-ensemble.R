

#' Phytophthora ensemble risk
#'
#' Calculates Phytophthora cinnamomi's ensemble risk from mechanisms
#'
#' @param mec_soilwater the result of \link{mec_soilwater}
#' @param mec_rootcontact the result of \link{mec_rootcontact}
#' @param mec_surfacewater the result of \link{mec_surfacewater}
#' @param mec_zoospread optional module. The result of \link{mec_zoospread}
#' @param weights weights of the ensemble model. The default uses the same
#' weights for each model. The argument accepts a numeric vector with the
#' corresponding weights. See *Details*
#'
#' @returns A \code{SpatRaster}
#' @export
#'
#' @details
#'
#' This function calculates the ensemble risk as the weighted sum of the individual
#' risk models as follows:
#'
#' \eqn{Risk = w_1 \cdot Mec_{soilwater} + w_2 \cdot Mec_{rootcontact} + w_3 \cdot Mec_{mec_surfacewater} + w_4  \cdot Mec_{zoospread}}.
#'
#' Being \eqn{w_i = 0.25} when \link{mec_zoospread} is not null, and \eqn{w_i = \frac{1}{3}} when \link{mec_zoospread} is null.
#' Additionally, the user can specify a numeric vector of length 3 or 4 specifying the desired weights based
#' on expert criteria.
#'
#' @examples
#' ## TODO
phytorisk_ensemble <- function(mec_soilwater,
                               mec_rootcontact,
                               mec_surfacewater,
                               mec_zoospread = NULL,
                               weights       = "equal") {

  # 0. Check errors

  # 1. Get weights
  if (weights == "equal") {
    if (is.null(mec_zoospread))  weights <- rep(1 / 3, 3)
    if (!is.null(mec_zoospread)) weights <- rep(1 / 4, 4)
  }

  ## weights must sum 1
  if (sum(weights) != 1) cli::cli_abort("'weights' must sum 1")

  # 2. Apply function depending if mec_zoospread is NULL
  if (is.null(mec_zoospread)) {

    ## check if weights are correct
    if (length(weights) != 3) cli::cli_abort("Wrong length of 'weights'. It must have 4 values")

    ## calculate risk
    risk_sr <- sum(c(mec_soilwater$mec_soilwater, mec_rootcontact, mec_surfacewater$mec_surfacewater) * weights)

  } else {

    ## check if weights are correct
    if (length(weights) != 4) cli::cli_abort("Wrong length of 'weights'. It must have 4 values")

    ## calculate risk
    risk_sr <- sum(c(mec_soilwater$mec_soilwater, mec_rootcontact, mec_surfacewater$mec_surfacewater, mec_zoospread) * weights)
  }

  return(risk_sr)


}






#' Raw Phytophthora ensemble risk
#'
#' Calculates Phytophthora cinnamomi's ensemble risk from raw data
#'
#' @param aoi a \code{sf} polygon representing the area of interest. Used to mask
#' the tree cover
#' @param poi a single-point \code{sf} object denoting the point of interest
#' to run the simulations
#' @param dem a single-band \code{SpatRaster} with a digital elevation model
#' @param treecover a single-band \code{SpatRaster} where 1 represents host
#' trees, and 0 represents background area
#' @param weights weights of the ensemble model. The default uses the same
#' weights for each model. The argument accepts a numeric vector with the
#' corresponding weights. See \link{phytorisk_ensemble}
#' @param th threshold of flow accumulation to delineate streams
#' @param buffer a buffer in meters to extend the spread in every direction
#' @param include_zoospread logical. Whether to include the optional module of \link{mec_zoospread}
#' @param ... arguments passed to \link{mec_zoospread}
#' @param quiet if \code{TRUE}, suppress any message or progress bar
#'
#' @returns A \code{SpatRaster}
#' @export
#'
#' @details
#'
#' TO DO
#'
phytorisk_ensemble_raw <- function(aoi,
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
  mec_surfacewater   <- mec_surfacewater(dem, mec_soilwater, poi, buffer = buffer, quiet = quiet)

  if (include_zoospread) {
    if (!quiet) cli::cli_h1("Mec III - Dispersion by animals")
    mec_zoospread   <- mec_zoospread(aoi, poi, mec_surfacewater, ..., quiet = quiet)
  } else {
    mec_zoospread <- NULL
  }

  final_risk <- phytorisk_ensemble(mec_soilwater, mec_rootcontact, mec_surfacewater, mec_zoospread, weights = weights)

  return(final_risk)

}
