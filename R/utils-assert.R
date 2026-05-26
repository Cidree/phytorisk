
assert_positive_numeric <- function(arg, ref) { # nocov start

    if (!is.numeric(arg) || length(arg) != 1 || arg < 0) {
        cli::cli_abort(
            "{.arg {ref}} must be a single positive numeric value (>= 0).",
            .frame = parent.frame()
        )
    }
} # nocov end


assert_logic <- function(arg, ref = "quiet") { # nocov start

    if (!is.logical(arg)) {
        cli::cli_abort(
            "{.arg {ref}} must be either TRUE or FALSE.",
            .frame = parent.frame()
            )
        }
} # nocov end

assert_integer_scalar <- function(arg, ref) { # nocov start

    if (!is.numeric(arg) || length(arg) != 1 || arg != as.integer(arg)) {
        cli::cli_abort(
            "{.arg {ref}} must be a single integer value.",
            .frame = parent.frame()
        )
    }
} # nocov end


assert_spatraster <- function(arg, ref) { # nocov start
    if (!inherits(arg, "SpatRaster")) {
        cli::cli_abort(
            "{.arg {ref}} must be a {.cls SpatRaster} object.",
            .frame = parent.frame()
        )
    }
} # nocov end

assert_sf <- function(arg, ref, geometry = NULL) { # nocov start
    if (!inherits(arg, "sf")) {
        cli::cli_abort(
            "{.arg {ref}} must be a {.cls sf} object.",
            .frame = parent.frame()
        )
    }

    if (!is.null(geometry)) {
        geom_types <- as.character(sf::st_geometry_type(arg, by_geometry = FALSE))
        if (!geom_types %in% geometry) {
            cli::cli_abort(
                "{.arg {ref}} must have {.val {geometry}} geometry, not {.val {geom_types}}.",
                .frame = parent.frame()
            )
        }
    }
} # nocov end

assert_sf_point <- function(arg, ref) { # nocov start
    if (!inherits(arg, "sf") || nrow(arg) != 1 || sf::st_geometry_type(arg) != "POINT") {
        cli::cli_abort(
            "{.arg {ref}} must be a single-row {.cls sf} object with a {.cls POINT} geometry.",
            .frame = parent.frame()
        )
    }
} # nocov end


assert_same_crs <- function(arg1, ref1, arg2, ref2) { # nocov start
    if (!terra::same.crs(arg1, arg2)) {
        cli::cli_abort(
            "{.arg {ref1}} and {.arg {ref2}} must have the same CRS.",
            .frame = parent.frame()
        )
    }
} # nocov end


assert_same_ext <- function(arg1, ref1, arg2, ref2) { # nocov start
    if (!terra::compareGeom(arg1, arg2, stopOnError = FALSE)) {
        cli::cli_abort(
            "{.arg {ref1}} and {.arg {ref2}} must have the same extent and resolution.",
            .frame = parent.frame()
        )
    }
} # nocov end
