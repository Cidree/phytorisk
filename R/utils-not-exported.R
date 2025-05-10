
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
