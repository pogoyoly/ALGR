# establish_tes
#' Establish fields by a tesselation algorithm
#'
#' @param potential_space a raster including a potential space category for field placement
#' @param includsion_value inclusion value for the potential space raster
#' @param mean_field_size mean field size counted by number of cells
#' @param sd_field_size sd field size counted by number of cells
#' @param mean_shape_index mean shape index calculated by a relation between width/length of placement
#' @param sd_shape_index sd shape index calculated by a relation between width/length of placement
#'
#' @return ALGR object
#' @export
#'
#' @examples
#'
#' # Create a raster with some potential space
#'test<-ALGR::generate_pn(50,50,1,2,3,0.01,TRUE, "land_percentage", percetange = 75)
#'terra::plot(test)
#'
#'# Run the simulation
#'result <- establish_tes(
#'  potential_space = test,
#'  includsion_value = 1,
#'  mean_field_size = 25,
#'  sd_field_size = 25,
#'  mean_shape_index = 2,
#'  sd_shape_index = 2)
#'
#'return_by_field(result, method = 1)
#'
establish_tes <- function(potential_space, includsion_value,
                            mean_field_size, sd_field_size,
                            mean_shape_index, sd_shape_index) {


  value_to_count <- includsion_value
  freq_table <- terra::freq(potential_space)
  count <- freq_table[freq_table$value == value_to_count, "count"]
  n_fields <- count / mean_field_size
  # Identify potential cells
  potential_cells <- which(values(potential_space) == includsion_value)
  field_raster <- terra::rast(potential_space)
  values(field_raster) <- NA
  resolution <- terra::res(field_raster)

  field_list <- list()
  field_id <- 1

  # Track used cells
  used_cells <- rep(FALSE, ncell(field_raster))

  #try each field 10 times

  # Step 1: Generate candidate field templates
  field_templates <- vector("list", n_fields)
  for (i in seq_len(n_fields)) {
    area <- max(4, round(abs(rnorm(1, mean_field_size, sd_field_size))))
    shape <- abs(rnorm(1, mean_shape_index, sd_shape_index))
    width <- max(1, round(sqrt(area * shape)))
    height <- max(1, round(area / width))
    field_templates[[i]] <- list(area = width * height, width = width, height = height)
  }

  # Sort largest fields first
  ordered_indices <- order(sapply(field_templates, function(x) x$area), decreasing = TRUE)
  field_templates <- field_templates[ordered_indices]
  placed_templates <- rep(FALSE, length(field_templates))

  # Get coordinates for all potential cells in row-major order
  candidate_rc <- terra::rowColFromCell(field_raster, potential_cells)

  for (k in seq_along(potential_cells)) {
    row <- candidate_rc[k, 1]
    col <- candidate_rc[k, 2]

    for (i in which(!placed_templates)) {
      template <- field_templates[[i]]
      dims <- list(
        list(height = template$height, width = template$width),
        list(height = template$width, width = template$height)  # Try both orientations
      )

      for (dim in dims) {
        rows <- row:(row + dim$height - 1)
        cols <- col:(col + dim$width - 1)

        if (max(rows) > nrow(field_raster) || max(cols) > ncol(field_raster)) next

        # Get all cell indices in the proposed rectangle
        cell_indices <- as.vector(outer(rows, cols, function(r, c) cellFromRowCol(field_raster, r, c)))

        if (all(cell_indices %in% potential_cells) && all(!used_cells[cell_indices])) {
          # Place field
          used_cells[cell_indices] <- TRUE
          values(field_raster)[cell_indices] <- field_id

          coords <- terra::xyFromCell(field_raster, cell_indices)
          ymin <- terra::ymin(field_raster)
          ymax <- terra::ymax(field_raster)
          coords[,2] <- ymax - (coords[,2] - ymin)
          loc <- list( x = coords[,2], y = coords[,1])



          new_field <- new("Field",
                           number = field_id,
                           location = loc,
                           farmer = NA_real_,
                           crop = NA_real_)
          field_list[[field_id]] <- new_field
          field_id <- field_id + 1

          placed_templates[i] <- TRUE
          break
        }
      }

      if (placed_templates[i]) break  # Move to next potential cell
    }

    if (all(placed_templates)) break  # All templates placed
  }


  return(list(
    map = field_raster,
    field_list = field_list
  ))
}



# Define the Field class
setClass("Field", slots = list(
  number = "numeric",
  location = "list",
  farmer = "numeric",
  crop = "numeric"
))


