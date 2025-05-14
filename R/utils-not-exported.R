
#'  (Internal) Extracts cell neighbours
#'
#' @param i row id
#' @param j column id
#' @param nrows total number of rows
#' @param ncols total number of columns
#'
#' @return a list of pixel neighbours
#' @keywords internal
extract_neighbours <- function(i, j, nrows, ncols) {

  neighbours <- list(
    N = c(i-1, j), S = c(i+1, j), E = c(i, j+1), O = c(i, j-1),
    NE = c(i-1, j+1), NO = c(i-1, j-1), SE = c(i+1, j+1), SO = c(i+1, j-1)
  )

  valid_neighbours <- lapply(neighbours, function(coord) {
    if (coord[1] > 0 && coord[1] <= nrows && coord[2] > 0 && coord[2] <= ncols) return(coord)
    else return(NULL)
  })

  return(valid_neighbours[!sapply(valid_neighbours, is.null)])

}





#'  (Internal) Isolate pixels
#'
#' @param flowdir a \code{SpatRaser} with the flow directions
#' @param dem a \code{SpatRaster} with a digital elevation model
#' @param init_i ith row to start
#' @param init_j jth column to stat
#'
#' @return a SpatRaster
#' @keywords internal
isolate_pixels <- function(flowdir, dem, init_i, init_j) {

  ## create a binary matrix to store the cells that fullfill the condition
  pixels_condition <- matrix(0, terra::nrow(flowdir), terra::ncol(flowdir))

  ## process the cells from the begining
  cells_to_visit <- list(c(init_i, init_j))
  visited <- matrix(FALSE, terra::nrow(flowdir), terra::ncol(flowdir))
  visited[init_i, init_j] <- TRUE

  while (length(cells_to_visit) > 0) {
    ## extract the first cell
    current_cell <- cells_to_visit[[1]]
    cells_to_visit <- cells_to_visit[-1]

    ## get the current row and col
    i <- current_cell[1]
    j <- current_cell[2]

    ## get the current's cell height
    current_height <- dem[i, j]

    ## get the current's cell neighbours
    neighbours <- extract_neighbours(i, j, nrow(flowdir), ncol(flowdir))

    for (neighbour in neighbours) {
      neighbour_i <- neighbour[1]
      neighbour_j <- neighbour[2]

      if (!visited[neighbour_i, neighbour_j]) {
        ## get the height of the neighbour
        height_neighbour <- dem[neighbour_i, neighbour_j]

        ## verify the condition of flow direction in the neighbour
        flow_cell <- terra::cellFromRowCol(flowdir, neighbour_i, neighbour_j)
        if (terra::values(flowdir)[flow_cell]) {
          ## if the neighbour's height is lower and fullfills the flow's condition
          if (height_neighbour < current_height) {
            pixels_condition[neighbour_i, neighbour_j] <- 1
            ## add the neighbour cell to the list of cells to visit
            cells_to_visit <- append(cells_to_visit, list(c(neighbour_i, neighbour_j)))
            visited[neighbour_i, neighbour_j] <- TRUE
          } else {
            pixels_condition[neighbour_i, neighbour_j] <- 0
          }
        }
      }
    }
  }

  return(terra::rast(pixels_condition, extent = flowdir))
}






#'  (Internal) Verify if a pixel is connected
#'
#' @param mat a matrix
#' @param x a row id
#' @param y a column id
#'
#' @return logical
#' @keywords internal
is_connected <- function(mat, x, y) {
  return(!is.na(mat[x, y]) && mat[x, y] == 1)
}





#'  (Internal) Finds connected pixels
#'
#' @param r a SpatRaster
#' @param init_i ith row to start
#' @param init_j jth row to start
#'
#' @return a list
#' @keywords internal
find_connected <- function(r, init_i, init_j) {

  ## convert raster to matrix
  mat <- as.matrix(r, wide = TRUE)

  visited_mat <- matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))
  connected_lst <- list()
  pile_lst <- list(c(init_i, init_j))

  while (length(pile_lst) > 0) {
    current <- pile_lst[[length(pile_lst)]]
    pile_lst <- pile_lst[-length(pile_lst)]
    x <- current[1]
    y <- current[2]

    if (x < 1 || y < 1 || x > nrow(mat) || y > ncol(mat)) {
      next
    }
    if (visited_mat[x, y]) {
      next
    }
    if (!is_connected(mat, x, y)) {
      next
    }

    visited_mat[x, y] <- TRUE
    connected_lst <- append(connected_lst, list(c(x, y)))

    ## add adjacent positions (including diagonals)
    adjacent_lst <- list(
      c(x + 1, y),     # down
      c(x - 1, y),     # up
      c(x, y + 1),     # right
      c(x, y - 1),     # left
      c(x + 1, y + 1), # down-right
      c(x + 1, y - 1), # down-left
      c(x - 1, y + 1), # up-right
      c(x - 1, y - 1)  # up-left
    )

    for (adjacent in adjacent_lst) {
      pile_lst <- append(pile_lst, list(adjacent))
    }
  }

  return(connected_lst)
}



#'  (Internal) Calculates Topographic Water Index
#'
#' @param dem a SpatRaster
#'
#' @return a SpatRaster
#' @keywords internal
twi <- function(dem) {
  ## get terrain slope
  slope <- terra::terrain(dem, v = "slope")

  ## get TWI
  twi <- log(1 / slope)

  ## use percentile 95 (based on bibliography)
  twi_p95 <- stats::quantile(twi, probs = 0.95, na.rm = TRUE)

  ## create a binary raster
  terra::as.int(twi > twi_p95)

}



#'  (Internal) Creates an amial
#'
#' @param x coordinate X
#' @param y coordinate Y
#' @param crs coordinate reference system
#'
#' @return a new animal
#' @keywords internal
## utils-not-exported
Animal <- function(x, y, crs) {
  animal <- sf::st_sfc(sf::st_point(c(x, y)))
  sf::st_crs(animal) <- crs
  return(animal)
}





#'  (Internal) Move animals towards food
#'
#'  Moves animals towards a food source specified in pixels based on randomness
#'
#' @param animal location of the animal
#' @param food_coords coordinates of food
#' @param pixel_size size of the pixel
#' @param exclusion_area_sf area of exclusion
#'
#' @importFrom stats runif
#'
#' @return a new animal
#' @keywords internal
move_towards_food <- function(animal, food_coords, pixel_size, exclusion_area_sf) {

  ## maximum number of iterations to avoid exclusion area
  max_attempts <- 10
  attempts <- 0

  repeat {
    attempts <- attempts + 1

    ## calculate direction towards food source
    direction_vector <- c(food_coords[1] - sf::st_coordinates(animal)[1], food_coords[2] - sf::st_coordinates(animal)[2])

    ## normalize direction vector
    direction_vector <- direction_vector / sqrt(sum(direction_vector^2))

    ## aggregate a random component to the movement
    random_angle <- runif(1, -pi / 4, pi / 4)  # change angle randomly within a range
    rotation_matrix <- matrix(c(cos(random_angle), -sin(random_angle), sin(random_angle), cos(random_angle)), nrow = 2)
    random_direction <- rotation_matrix %*% direction_vector

    ## move towards the food source with the given pixel size, and random component
    new_coords <- c(
      sf::st_coordinates(animal)[1] + random_direction[1] * pixel_size,
      sf::st_coordinates(animal)[2] + random_direction[2] * pixel_size
    )

    ## create a new point with the proposed coordinates
    new_animal <- sf::st_sfc(sf::st_point(new_coords))
    sf::st_crs(new_animal) <- sf::st_crs(animal)  # maintain original CRS

    ## verify if the new point is within the exclusion area
    if (length(sf::st_intersects(new_animal, exclusion_area_sf)) == 0 || attempts >= max_attempts) {
      ## exit the loop if the animal is not within the exclusion area or if we reached the maximum number of attempts
      break
    }
  }

  return(new_animal)
}





#'  (Internal) Do not leave study area
#'
#' @param animal location of the animal
#' @param aoi area of interest (sf)
#'
#' @return animal
#' @keywords internal
stay_within_area <- function(animal, aoi) {

  ## verify if CRS if the same
  if (sf::st_crs(animal) != sf::st_crs(aoi)) cli::cli_abort("CRS is not the same")

  if (!sf::st_within(animal, aoi, sparse = FALSE)) {

    ## if the animal is out of the study area, move it inside
    nearest_point <- sf::st_nearest_points(animal, aoi)

    ## coords of the closest point to the study area
    nearest_coords <- sf::st_coordinates(nearest_point)[2, ]
    direction_vector <- c(nearest_coords[1] - sf::st_coordinates(animal)[1],
                          nearest_coords[2] - sf::st_coordinates(animal)[2])
    norm_vector <- sqrt(sum(direction_vector^2))
    new_coords <- c(sf::st_coordinates(animal)[1] + direction_vector[1] / norm_vector,
                    sf::st_coordinates(animal)[2] + direction_vector[2] / norm_vector)
    animal <- sf::st_sfc(sf::st_point(new_coords))
    sf::st_crs(animal) <- sf::st_crs(aoi)
  }

  return(animal)
}
