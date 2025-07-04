

#' Establish fields by a place and conquer algorithm
#'
#' @param potential_space a raster including a potential space category for field placement
#' @param cell_size cell size for output
#' @param includsion_value inclusion value for the potential space raster
#' @param additional_lim for use in case of inclusion of a raod raster
#' @param mean_field_size mean field size counted by number of cells
#' @param sd_field_size sd field size counted by number of cells
#' @param distribution either norm or lnorm
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
#' r<-terra::rast(matrix(1, nrow=50, ncol=50))
#' output<-establish_pac(potential_space= r,
#'                          cell_size=1,
#'                          includsion_value = 1,
#'                          mean_field_size = 200,
#'                          sd_field_size = 100,
#'                          distribution = "norm",
#'                          mean_shape_index = 3,
#'                          sd_shape_index = 0.3,
#'                          percent = 70)
#'
#'
#' return_by_field(output)


establish_pac<-function(potential_space,
                                     cell_size,
                                     includsion_value,
                                     additional_lim = NULL,
                                     mean_field_size,
                                     sd_field_size,
                                     distribution = "norm",
                                     mean_shape_index,
                                     sd_shape_index,
                                     percent
){

  checkmate::assertClass(potential_space, "SpatRaster")
  checkmate::assert_numeric(cell_size)
  checkmate::assert_numeric(includsion_value)
  checkmate::assert_numeric(mean_field_size)
  checkmate::assert_numeric(sd_field_size)
  checkmate::assert_numeric(mean_shape_index)
  checkmate::assert_numeric(sd_shape_index)
  checkmate::assert_numeric(percent)

  if (mean_shape_index < 1 || mean_shape_index > 5) {
    rlang::abort("mean_shape_index must be between 1 and 5.")
  }


  # Compute raster dimensions
  raster_nrow <- nrow(potential_space)
  raster_ncol <- ncol(potential_space)
  raster_total_cells <- raster_nrow * raster_ncol

  # Estimate maximum expected field size (conservative ~99.7% of normal range)
  expected_max_field_size <- mean_field_size + 3 * sd_field_size

  # Check if expected field is larger than available space
  if (expected_max_field_size > raster_total_cells) {
    rlang::abort(paste0(
      "Mean field size (plus 3*SD) is too large for the raster.\n",
      "Expected max field size: ", expected_max_field_size,
      " vs raster cells: ", raster_total_cells, ".\n",
      "Reduce mean_field_size or adjust raster extent."
    ))
  }

  if(is.null(additional_lim) == FALSE){
    #make sure all rasters are the same size
    additional_lim <- terra::rast(additional_lim)
    template <- potential_space
    additional_lim <- terra::resample(additional_lim, template, method = "bilinear") # Use "near" for categorical data
  }


  #set land raster on which to save the fields
  potential_space[is.na(potential_space[])] <- 0
  potential_space<-as.matrix(potential_space, wide=TRUE)
  land = matrix(0, nrow(potential_space), ncol(potential_space))
  land_raster<-terra::rast(land)

  #setup place holders
  achieved_percent <- 0
  field_num <- 1
  field_list <- list()


  #add any additional limitations to a seperate raster
  if(is.null(additional_lim) == FALSE){

    #denote the unaccecable land
    additional_lim[which(is.na(additional_lim[]) == FALSE)] <- 255
    additional_lim[is.na(additional_lim[])] <- 0


    additional_lim<-as.matrix(additional_lim, wide=TRUE)

    #add the roats to the potential space limiter
    potential_space <- potential_space + additional_lim

  }

  field_num <- 1

  #main placement loop
  while(achieved_percent < percent){


    field_placed <- FALSE
    while(field_placed == FALSE){

      #set field size and shape index
      if (distribution == "norm") {
        field_size <- max(1, round(rnorm(1, mean = mean_field_size, sd = sd_field_size)))
      } else if (distribution == "lnorm") {
        mu_N <- log(mean_field_size^2 / sqrt(sd_field_size^2 + mean_field_size^2))
        sigma_N <- sqrt(log(1 + (sd_field_size^2 / mean_field_size^2)))
        field_size <- max(1, round(rlnorm(1, meanlog = mu_N, sdlog = sigma_N)))
      } else if (distribution == "unif") {
        min_val <- max(1, mean_field_size - sd_field_size)
        max_val <- mean_field_size + sd_field_size
        field_size <- max(1, round(runif(1, min = min_val, max = max_val)))
      } else {
        rlang::abort(paste0("Unsupported distribution type: '", distribution,
                            "'. Choose one of 'norm', 'lnorm', or 'unif'."))
      }


      shape_index <- rnorm(1, mean=mean_shape_index, sd=sd_shape_index)
      if(shape_index <= 1){
        shape_index<-1
      }
      if(field_size <= 0){
        field_size<-1
      }

      ratio <- (shape_index - 1) / 4
      axis_exs <- c("horiz", "vert")
      axis_ex<-sample(axis_exs,1)
      if(axis_ex == "horiz"){
        field_col_size <- round(sqrt(field_size) * (1 + ratio))  # More columns as shape_index increases
        field_row_size <- round(field_size / field_col_size)  # Adjust rows to maintain total_cells

      }
      if(axis_ex == "vert"){
        field_row_size <- round(sqrt(field_size) * (1 + ratio))  # More columns as shape_index increases
        field_col_size <- round(field_size / field_row_size)  # Adjust rows to maintain total_cells

      }



      #choose random start location and check if it can be used
      check<-999
      check2 <- 1
      while(check != includsion_value || check2 != 0){
        start_col <- sample(1:ncol(potential_space), 1)
        start_row <- sample(1:nrow(potential_space), 1)
        if(is.na(potential_space[start_row,start_col]) == FALSE){
          check<-potential_space[start_row,start_col]
          check2<-land[start_row,start_col]
        }
      }



      #choose direction for starting placment
      directions <- c("one", "two")
      direction<-sample(directions,1)
      switch(direction,
             "one" = {
               dir<- 1
             },
             "two" = {
               dir<- - 1
             },
      )


      #define max size
      max_size<- field_row_size * field_col_size
      changes <- 0
      cur_col <- start_col
      cur_row <- start_row
      placed_cells <- 0
      #start expansion according to direction

      while(placed_cells < max_size && changes != 2){
        suppressWarnings({

          #create block vectors
          vec1<- seq(start_row, start_row + field_row_size, by=1)
          vec2 <- rep(cur_col, field_col_size)

          #function to check if there is overlap
          coverlap<-function(a,b){
            if(a > nrow(potential_space) || b > ncol(potential_space) || a < 1 || b < 1 ||
               is.na(potential_space[a,b]) == TRUE ||  land[a,b] != 0 ||  potential_space[a,b] != includsion_value){
              return("overlap")
            }
          }

          #count how many overlaps
          co_vec <- mapply(coverlap, vec1, cur_col, SIMPLIFY = TRUE)
          num_co <- sum(co_vec == "overlap")

          #if no overlaps place
          if(num_co == 0){
            for(i in 1:length(vec1)){
              land[vec1[i], cur_col] <- field_num

            }

            placed_cells <- placed_cells + length(vec1)
            cur_col <- cur_col + dir

          }
          if(num_co > 0){
            #if there are overlaps we try to move the block left and right to see if we can minimize

            #fist we create a matrix of the possible shift
            shift_vec<-seq( - field_row_size, field_row_size, by = 1)
            shift_func <- function(vec, shift) {
              vec<- vec + shift
            }
            shifted_vectors <- mapply(function(shift) shift_func(vec1, shift), shift_vec, SIMPLIFY = FALSE)
            result_matrix <- do.call(rbind, shifted_vectors)


            #count how many overlaps for each possible shifted vector
            coverlap_test<-function(a){

              if(a > nrow(potential_space) || cur_col > ncol(potential_space) || a < 1 || cur_col < 1 ||
                 is.na(potential_space[a,cur_col]) == TRUE || land[a,cur_col] != 0 ||  potential_space[a,cur_col] != includsion_value){
                return(1)
              }
              else{
                return(0)
              }
            }

            co_vec <- rowSums(apply(FUN = coverlap_test, MARGIN = c(1,2), result_matrix))

            if(sum(co_vec)== 0){
              dir <- dir * -1
              cur_col <- start_col + dir
              changes <- changes + 1
            }
            else{
              #choose the one with the lowest value of co_vec without contributing to the shift
              find_threshold_with_sapply <- function(n_vector, m_vector) {
                # Ensure that n_vector and m_vector are sorted in the same order based on the absolute value of n_vector
                sorted_indices <- order(abs(n_vector))
                n_vector <- n_vector[sorted_indices]
                m_vector <- m_vector[sorted_indices]

                # Find the index where the minimum m value is found
                min_m_index <- which.min(m_vector)

                # If there are multiple minimum values, find the last occurrence
                if(length(min_m_index) > 1) {
                  min_m_index <- max(min_m_index)
                }

                # Return the n value at the threshold point
                return(n_vector[min_m_index])
              }
              #This is the shift that is now choosen
              threshold_n <- find_threshold_with_sapply(shift_vec, co_vec)

              #shift by the choosen
              right_vec <- shift_func(vec1,threshold_n)
              for(i in 1:length(right_vec)){
                if(right_vec[i] > 0 && right_vec[i] < nrow(land) && cur_col > 0 && cur_col < ncol(land) && potential_space[right_vec[i], cur_col] == includsion_value){
                  land[right_vec[i], cur_col] <- field_num
                }

              }

              placed_cells <- placed_cells +length(vec1)

              cur_col <- cur_col + dir

              co_vec <- mapply(coverlap, right_vec, cur_col, SIMPLIFY = TRUE)
              num_co <- sum(co_vec == "overlap")


              if(cur_col >= ncol(land) || num_co == length(right_vec)){
                dir <- dir * -1
                cur_col <- start_col + dir
                changes <- changes + 1

              }




            }




            #otherwise place field with current shift

          }

          #now add to cur_col and start again


        })

      }


      #check block one direction

      #count how many you cannot place

      #expand block to other direction

      #add to land raster with field number indexed

      #repeat till max size




      field_placed <- TRUE
      field_num <- field_num + 1


    }

    #calculate how much of the space has been filled
    calculate_part_potential <- function(potential_space, inclusion_value) {
      # Ensure potential_space is a matrix and inclusion_value is a vector or scalar
      part_potential <- length(which(potential_space %in% inclusion_value))
      return(part_potential)
    }
    part_potential <- calculate_part_potential(potential_space,includsion_value )
    part_filled <- length(land[land > 0])
    achieved_percent <- (part_filled / part_potential) * 100

  }

  land_raster <- terra::rast(land)
  patched_raster <- landscapemetrics::get_patches(land_raster)[[1]]


  for(i in 2:length(patched_raster)){
    mati <- terra::as.matrix(patched_raster[[i]], wide=TRUE)
    dimnames(mati) <- list(x = 1:nrow(mati), y = 1:ncol(mati))
    newdata <- expand.grid(x = 1:nrow(mati), y = 1:ncol(mati))
    newdata$Z <- as.vector(mati)
    newdata <- newdata[!is.na(newdata$Z),]

    field_obj <- new("Field", number = (i-1), location = list(newdata$x,newdata$y), farmer = 1)
    field_list<-c(field_list,field_obj)

    # Give names to each row and column as well as names of each dimension of the matrix itself.

  }



  #distribute fields between farmers

  #make sure the output extent is defined by cell size and has origin in 00
  terra::ext(land_raster)<-c(0, cell_size*ncol(land_raster), 0, cell_size*nrow(land_raster))

  #unify resutls in list
  result<-list(map = land_raster, field_list = field_list)

  return(result)
}


setClass("Field", slots=list(number="numeric",location="list",farmer="numeric", crop = "numeric"))
