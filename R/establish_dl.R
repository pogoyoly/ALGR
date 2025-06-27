# Function to generate dead leaves texture
#' Establish fields by a dead leaves algorithm algorithm
#'
#' @param potential_space a raster including a potential space category for field placement
#' @param cell_size cell size for output
#' @param includsion_value inclusion value for the potential space raster
#' @param mean_field_size mean field size counted by number of cells
#' @param sd_field_size sd field size counted by number of cells
#' @param distribution distribution type for field size
#' @param mean_shape_index mean shape index calculated by a relation between width/length of placement
#' @param sd_shape_index sd shape index calculated by a relation between width/length of placement
#' @param percent percent of the potential space to be filled with fields
#'
#' @return ALGR object
#' @export
#' @importFrom stats na.omit rlnorm rnorm runif setNames
#' @importFrom methods new
#'
#' @examples
#' r<-terra::rast(matrix(1, nrow=100, ncol=100))
#' output <- establish_dl(potential_space = r,
#'                                                 cell_size = 1,
#'                                                 includsion_value = 1,
#'                                                 mean_field_size = 50,
#'                                                 sd_field_size = 25,
#'                                                 distribution = "norm",
#'                                                 mean_shape_index = .5,
#'                                                 sd_shape_index = .1,
#'                                                 percent = 0.95)
#' return_by_field(output)


establish_dl <- function(potential_space,
                                     cell_size,
                                     includsion_value,
                                     mean_field_size,
                                     sd_field_size,
                                     distribution = "norm",
                                     mean_shape_index,
                                     sd_shape_index,
                                     percent
) {
  # Initialize empty canvas

  checkmate::assertClass(potential_space, "SpatRaster")
  checkmate::assert_numeric(cell_size)
  checkmate::assert_numeric(includsion_value)
  checkmate::assert_numeric(mean_field_size)
  checkmate::assert_numeric(sd_field_size)
  checkmate::assert_numeric(mean_shape_index)
  checkmate::assert_numeric(sd_shape_index)
  checkmate::assert_numeric(percent)



  #setup matrix to be filled
  canvas <- matrix(0, nrow = nrow(potential_space), ncol = ncol(potential_space))
  potential_space_matrix <- terra::as.matrix(potential_space, wide=TRUE)
  potential_space_matrix[potential_space_matrix != includsion_value] <- 0
  potential_space_matrix[potential_space_matrix == includsion_value] <- 1

  #creat empty field list
  field_list <- list()

  #calculate size of potential space and set realized/to be filled
  num_potential_patches <- length(potential_space_matrix[potential_space_matrix == 1])
  realized_patches <- 0
  to_be_filled <- num_potential_patches

  #set tracker for field numbers
  i <- 1

  # Generate patches
  while (realized_patches < to_be_filled) {

    #set field size and shape index
    if (distribution == "norm") {
      field_size <- max(1, round(rnorm(1, mean = mean_field_size, sd = sd_field_size)))
    }
    if (distribution == "lnorm") {
      mu_N <- log(mean_field_size^2 / sqrt(sd_field_size^2 + mean_field_size^2))
      sigma_N <- sqrt(log(1 + (sd_field_size^2 / mean_field_size^2)))

      field_size <- max(1, round(rlnorm(1, meanlog = mu_N, sdlog = sigma_N)))

    }
    shape_index <- rnorm(1, mean = mean_shape_index, sd = sd_shape_index)
    if (shape_index <= 1) {
      shape_index <- 1
    }
    if (field_size <= 0) {
      field_size <- 1
    }

    ratio <- (shape_index - 1) / 4
    axis_exs <- c("horiz", "vert")
    axis_ex <- sample(axis_exs,1)
    if (axis_ex == "horiz") {
      field_col_size <- round(sqrt(field_size) * (1 + ratio))  # More columns as shape_index increases
      field_row_size <- round(field_size / field_col_size)  # Adjust rows to maintain total_cells

    }
    if (axis_ex == "vert") {
      field_row_size <- round(sqrt(field_size) * (1 + ratio))  # More columns as shape_index increases
      field_col_size <- round(field_size / field_row_size)  # Adjust rows to maintain total_cells

    }


    #calcualte how much is realized
    realized_patches <- length(canvas[canvas > 0])




    if (field_row_size > nrow(potential_space) | field_col_size > ncol(potential_space)) {
      next
    }

    #choose start location
    patch_size <- field_size
    x <- sample(1:(ncol(potential_space) - field_row_size + 1), 1)  # Ensure patch stays within canvas bounds
    y <- sample(1:(nrow(potential_space) - field_col_size + 1), 1) # Ensure patch stays within canvas bounds

    # Add patch to canvas
    # Get the number of rows and columns in the canvas
    n_rows <- nrow(canvas)
    n_cols <- ncol(canvas)

    # Calculate the bounds for x and y
    x_max <- x + field_row_size - 1
    y_max <- y + field_col_size - 1

    # Ensure that the indices stay within bounds
    if (x_max < n_cols && y_max < n_rows) {
      # If within bounds, update the canvas as intended

      canvas[y:y_max, x:x_max] <- i
    } else {
      # If not fully within bounds, calculate the largest allowed region
      valid_x_max <- min(x_max, n_cols)
      valid_y_max <- min(y_max, n_rows)




      # Update the canvas within the valid region
      if(y <= nrow(canvas) && y >= 0 && x <= ncol(canvas) && x >= 0){
        canvas[y:valid_y_max, x:valid_x_max] <- i
      }
    }
    #update field number
    i <- i + 1

    #remove all sections outside of potential space
    canvas <- matrixcalc::hadamard.prod(canvas, potential_space_matrix)


  }

  #final removel of all sections outside of potential space
  output <- matrixcalc::hadamard.prod(canvas, potential_space_matrix)



  #convert to raster format to be used with landscape metrics tools
  dead_leves_rast <- terra::rast(output)


  #remove small patches till you reach the percentage that needs to be filled
  patched_raster <- landscapemetrics::get_patches(dead_leves_rast)[[1]]
  patchedf <- landscapemetrics::lsm_p_area(patched_raster, directions = 8)
  # patchedf<- patchedf %>% arrange(value)
  patchedf <- patchedf[order(patchedf$value), ]

  unique_list <- terra::unique(dead_leves_rast)
  realized_patches <- length(dead_leves_rast[dead_leves_rast > 0])
  to_be_filled <- realized_patches * percent

  #remove the smallest patch
  while (realized_patches > to_be_filled) {
    val <- patchedf[1,1]
    val <- as.numeric(val)
    replacement_value <- as.numeric(unique_list[val,1])
    dead_leves_rast[dead_leves_rast == replacement_value] <- 0
    realized_patches <- length(dead_leves_rast[dead_leves_rast > 0])
    patchedf <- patchedf[-1,]
  }

  #reindex
  #sapply(unique_list, function(unique_list) dead_leves_rast[dead_leves_rast == unique_list] <- unique_list)
  for (j in 1:length(unique_list)) {
    dead_leves_rast[dead_leves_rast == unique_list[j,1]] <- j

  }


  #get patched raster for saving field locations in field list
  patched_raster <- landscapemetrics::get_patches(dead_leves_rast)[[1]]

  for (i in 2:length(patched_raster)) {

    mati <- terra::as.matrix(patched_raster[[i]], wide = TRUE)
    dimnames(mati) <- list(x = 1:nrow(mati), y = 1:ncol(mati))
    #mydf <- reshape2::melt(mati,level = 1, varnames = c("x", "y"))  # Ensure proper names
    #print(mydf)
    #names(mydf) <- c("x", "y", "Z")
    #newdata <- mydf[!is.na(mydf$Z),]  # Filter out NAs
    newdata <- expand.grid(x = 1:nrow(mati), y = 1:ncol(mati))
    newdata$Z <- as.vector(mati)
    newdata <- newdata[!is.na(newdata$Z),]

    field_obj <- new("Field", number = (i - 1), location = list(newdata$x,newdata$y), farmer = 1)
    field_list <- c(field_list,field_obj)

    # Give names to each row and column as well as names of each dimension of the matrix itself.

  }


  #set extent by input raster
  terra::ext(dead_leves_rast) <- c(0, cell_size*ncol(dead_leves_rast), 0, cell_size*nrow(dead_leves_rast))

  #finalize result into single object
  result <- list(map = dead_leves_rast, field_list = field_list)

  return(result)
}
