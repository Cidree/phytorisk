
#' Dispersion by animals
#'
#' Optional Module. Simulates dispersion by domestic and wild animal movement
#'
#' @param aoi a \code{sf} polygon representing the area of interest
#' @param poi a single-point \code{sf} object denoting the point of interest
#' to run the simulations
#' @param mec_surface the result of \link{mec_surfacewater}
#' @param n_animals number of animals
#' @param n_steps number of steps of each animal for the simulations
#' @param pixel_size size of the movement in pixels
#' @param n_iter number of random iterations for each animal
#' @param dist filter trajectories less than the specified meters away of the
#' POI (susceptible to inoculum)
#' @param quiet if \code{TRUE}, suppress any message or progress bar
#'
#' @returns A \code{SpatRaster}
#' @importFrom stats aggregate
#' @export
#'
#' @details
#' This function models animal movement as a series of straight-line steps
#' towards resource sources, while avoiding exclusion areas and staying within
#' defined boundaries. It randomly assigns initial animal positions and calculates
#' movement direction towards resources. The direction is normalized and randomized
#' using a rotation matrix. Animal trajectories are stored in vector format, and
#' those crossing the foci are identified for further analysis. Simulation
#' parameters, such as the number of animals and steps, are set at the start.
#'
#' @references
#'
#' Kliejunas, J.T., Ko, W.H., 1976. Dispersal of Phytophthora cinnamomi on the island of Hawaii. Phytopathology 66, 457–460.
#'
#' Li, A.Y., Williams, N., Fenwick, S.G., Hardy, G.E.St.J., Adams, P.J., 2014b. Potential for dissemination of Phytophthora cinnamomi by feral pigs via ingestion of infected plant material. Biol. Invasions 16, 765–774. \doi{10.1007/s10530-013-0535-7}
#'
#' Cardillo, E., Acedo, A., Abad, E., 2018. Topographic effects on dispersal patterns of Phytophthora cinnamomi at a stand scale in a Spanish heathland. PloS One 13, e0195060.
#'
#' @examples
#' ## load packages
#' # TODO
mec_zoospread <- function(aoi,
                          poi,
                          mec_surface,
                          n_animals = 5,
                          n_steps = 100,
                          pixel_size = 1,
                          n_iter = 10,
                          dist = 5,
                          quiet = FALSE) {

  ## 1. Generate animal trajectories --------------

  trajectory_lst <- list()
  if (!quiet) cli::cli_h2("Starting animal movement simulation")
  if (!quiet) cli::cli_progress_bar("Simulating trajectories", total = n_iter)
  for (iteration in 1:n_iter) {

    ## user feedback
    if (!quiet) cli::cli_progress_update()

    ## initialize the animal's position randomly within the study area
    withr::with_seed(iteration, {
      random_points <- sf::st_sample(aoi, size = n_animals)
      animals <- lapply(1:n_animals, function(i) {
        coords <- sf::st_coordinates(random_points[i])
        Animal(coords[1], coords[2], sf::st_crs(aoi))
      })

      ## randomly select a source of food for each animal
      food_coords_list <- lapply(1:n_animals, function(i) {
        food_coords <- sf::st_coordinates(mec_surface$surface_water[sample(nrow(mec_surface$surface_water), 1), ])[1, ]
        return(food_coords)
      })
    })

    ## store animal positions in each step
    animal_positions_df <- data.frame(
      step      = rep(1:n_steps, each = n_animals),
      x         = numeric(n_animals * n_steps),
      y         = numeric(n_animals * n_steps),
      animal_id = rep(1:n_animals, each = n_steps)
    )

    ## simulate animals movement
    for (step in 1:n_steps) {
      for (i in 1:n_animals) {
        ## move animals towards food source in pixels
        animals[[i]] <- move_towards_food(animals[[i]], food_coords_list[[i]], pixel_size)

        ## stay within the study area
        animals[[i]] <- stay_within_area(animals[[i]], aoi)

        ## store animal's positions in the data frame
        index <- (step - 1) * n_animals + i
        animal_positions_df[index, ] <- c(step, sf::st_coordinates(animals[[i]])[1], sf::st_coordinates(animals[[i]])[2], i)
      }
    }

    ## convert positions to points
    animal_positions_sf <- sf::st_as_sf(animal_positions_df, coords = c("x", "y"), crs = sf::st_crs(aoi))

    ## group by animal and create trajectory lines
    trajectory_union_sf <- aggregate(
      animal_positions_sf["geometry"],
      by = list(animal_positions_sf$animal_id),
      FUN = function(x) sf::st_union(x)
    )

    ## rename grouping column
    names(trajectory_union_sf)[1] <- "animal_id"

    ## cast geometries to LINESTRING
    trajectory_lines_sf <- sf::st_cast(trajectory_union_sf, "LINESTRING")
    trajectory_lst[[iteration]] <- trajectory_lines_sf

  }

  ## completed simulation
  if (!quiet) cli::cli_progress_done()
  if (!quiet) cli::cli_alert_success("Simulation completed. {n_iter} trajectories generated.")


  ## 2. Prepare data -----------------------
  if (!quiet) cli::cli_progress_step("Preparing results", "Success", "Failed")

  ## bind in a single file
  trajectory_sf <- do.call(rbind, trajectory_lst)

  ## calculate distances between trajectories and POI
  trajectory_sf$dist <- as.numeric(sf::st_distance(trajectory_sf, poi))

  ## filter trajectories less than X meters away of the POI (susceptible to inoculum)
  trajectory_sf <- trajectory_sf[trajectory_sf$dist < dist, ]

  ## if there are no trajectories, exit function
  if (nrow(trajectory_sf) == 0) {
    cli::cli_alert_warning("There are no trajectories within {dist} meters of the inoculum")
    return(NULL)
  }

  ## create a raster template for the study area
  risk_pc_sr <- terra::rast(mec_surface$mec_surfaceflow, vals = 0)

  ## rasterize trajectories
  for (i in 1:nrow(trajectory_sf)) {
    new_line_sr <- terra::rasterize(trajectory_sf[i,], risk_pc_sr, field = 1, update = TRUE, fun = 'sum')
    risk_pc_sr <- max(risk_pc_sr, new_line_sr, na.rm = TRUE)
  }

  ## stablish the pixel values no intercepted at 0
  risk_pc_sr[risk_pc_sr[] > 1] <- 1

  ## return
  names(risk_pc_sr) <- "mec_zoospread"
  return(risk_pc_sr)


}


